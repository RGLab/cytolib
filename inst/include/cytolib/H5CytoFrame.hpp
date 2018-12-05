/*
 * H5CytoFrame.hpp
 *
 *  Created on: Sep 18, 2017
 *      Author: wjiang2
 */

#ifndef INST_INCLUDE_CYTOLIB_H5CYTOFRAME_HPP_
#define INST_INCLUDE_CYTOLIB_H5CYTOFRAME_HPP_
#include <cytolib/MemCytoFrame.hpp>


namespace cytolib
{
/**
 * The class represents the H5 version of cytoFrame
 * It doesn't store and own the event data in memory.
 * Instead, data is read from H5 file on demand, which is more memory efficient.
 */
class H5CytoFrame:public CytoFrame{
protected:
	string filename_;
	/*EDIT: We now do not maintain these handlers, instead treat each IO as atomic operations
	 * Because it is not easy to accomplish the resource sharing among multiple H5CytoFrame objects solely depending on H5's mechanisms.
	 * e.g. a second openFile call with H5F_ACC_RDWR will overwrite the previous H5F_ACC_RDONLY, thus cause the unexpected data tampering
	 * these H5 handlers remain open during the life cycle of H5CytoFrame
	 * for faster accessing the data
	 */
	H5File file;
	DataSet dataset;
	DataSpace dataspace;
	hsize_t dims[2];              // dataset dimensions
	//flags indicating if cached meta data needs to be flushed to h5
	bool is_dirty_params;
	bool is_dirty_keys;
	bool is_dirty_pdata;
	EVENT_DATA_VEC read_data(uvec col_idx) const
	{
		unsigned nrow = n_rows();
		unsigned ncol = col_idx.size();
		/*
		 * Define the memory dataspace.
		 */
		hsize_t dimsm[] = {ncol, nrow};
		DataSpace memspace(2,dimsm);
		hsize_t      offset_mem[2];
		hsize_t      count_mem[2];
		EVENT_DATA_VEC data(nrow, ncol);
		//reach one col at a time
		for(unsigned i = 0; i < ncol; i++)
		{
			//select slab for h5 data space
			unsigned idx = col_idx[i];
			hsize_t      offset[] = {idx, 0};   // hyperslab offset in the file
			hsize_t      count[] = {1, nrow};    // size of the hyperslab in the file
			dataspace.selectHyperslab( H5S_SELECT_SET, count, offset );

			//select slab for mem
			offset_mem[0] = i;
			offset_mem[1] = 0;
			count_mem[0]  = 1;
			count_mem[1]  = nrow;
			memspace.selectHyperslab( H5S_SELECT_SET, count_mem, offset_mem );

			dataset.read(data.memptr(), h5_datatype_data(DataTypeLocation::MEM) ,memspace, dataspace);
		}

		return data;
	}


public:

	~H5CytoFrame(){
		flush_meta();

	};
	/*
	 * for simplicity, we don't want to handle the object that has all the h5 handler closed
	 * because it will require lots of validity checks before each disk IO operations
	 */
//	void close_h5(){
//		dataspace.close();
//		dataset.close();
//		file.close();
//	}
	void flush_meta(){
		//flush the cached meta data from CytoFrame into h5
		if(is_dirty_params)
			flush_params();
		if(is_dirty_keys)
			flush_keys();
		if(is_dirty_pdata)
			flush_pheno_data();
	}
	void flush_params()
	{
		CompType param_type = get_h5_datatype_params(DataTypeLocation::MEM);
		DataSet ds = file.openDataSet("params");
		ds.write(params.data(), param_type );
		is_dirty_params = false;
	}

	void flush_keys()
	{
		CompType key_type = get_h5_datatype_keys();
		DataSet ds = file.openDataSet("keywords");

		//convert to vector
		vector<KEY_WORDS_SIMPLE> keyVec;
		for(const auto & e : keys_)
		{
			keyVec.push_back(KEY_WORDS_SIMPLE(e.first, e.second));
		}


		ds.write(&keyVec[0], key_type );
		is_dirty_keys = false;
	}
	void flush_pheno_data()
	{
		CompType key_type = get_h5_datatype_keys();
		DataSet ds = file.openDataSet("pdata");

		//convert to vector

		vector<KEY_WORDS_SIMPLE> keyVec;
		for(std::pair<std::string, string> e : pheno_data_)
		{
			keyVec.push_back(KEY_WORDS_SIMPLE(e.first, e.second));
		}


		ds.write(&keyVec[0], key_type );
		is_dirty_pdata = false;
	}

	H5CytoFrame(const H5CytoFrame & frm):CytoFrame(frm)
	{
		filename_ = frm.filename_;
		file = frm.file;//safe to copy due to refcount during copy constructor provided by h5
		dataset = frm.dataset;//safe to copy due to refcount during copy constructor provided by h5
		dataspace = frm.dataspace;//safe to copy due to explicit copy through its assignment operator provided by h5
		is_dirty_params = frm.is_dirty_params;
		is_dirty_keys = frm.is_dirty_keys;
		is_dirty_pdata = frm.is_dirty_pdata;
		memcpy(dims, frm.dims, sizeof(dims));

	}
	H5CytoFrame(H5CytoFrame && frm):CytoFrame(frm)
	{
//		swap(pheno_data_, frm.pheno_data_);
//		swap(keys_, frm.keys_);
//		swap(params, frm.params);
//		swap(channel_vs_idx, frm.channel_vs_idx);
//		swap(marker_vs_idx, frm.marker_vs_idx);
		swap(filename_, frm.filename_);
		swap(file, frm.file);
		swap(dataset, frm.dataset);
		swap(dataspace, frm.dataspace);
		swap(dims, frm.dims);

		swap(is_dirty_params, frm.is_dirty_params);
		swap(is_dirty_keys, frm.is_dirty_keys);
		swap(is_dirty_pdata, frm.is_dirty_pdata);
	}
	H5CytoFrame & operator=(const H5CytoFrame & frm)
	{
		CytoFrame::operator=(frm);
		filename_ = frm.filename_;
		file = frm.file;
		dataset = frm.dataset;
		dataspace = frm.dataspace;
		is_dirty_params = frm.is_dirty_params;
		is_dirty_keys = frm.is_dirty_keys;
		is_dirty_pdata = frm.is_dirty_pdata;
		memcpy(dims, frm.dims, sizeof(dims));
		return *this;
	}
	H5CytoFrame & operator=(H5CytoFrame && frm)
	{
		CytoFrame::operator=(frm);
		swap(filename_, frm.filename_);
		swap(file, frm.file);
		swap(dataset, frm.dataset);
		swap(dataspace, frm.dataspace);
		swap(dims, frm.dims);

		swap(is_dirty_params, frm.is_dirty_params);
		swap(is_dirty_keys, frm.is_dirty_keys);
		swap(is_dirty_pdata, frm.is_dirty_pdata);
		return *this;
	}
	unsigned n_rows() const{
				return dims[1];
		}
	void set_channel(const string & oldname, const string &newname, bool is_update_keywords = true)
	{
		CytoFrame::set_channel(oldname, newname, is_update_keywords);
		is_dirty_params = true;
		if(is_update_keywords)
			is_dirty_keys = true;
	}
	/**
	 * constructor from FCS
	 * @param fcs_filename
	 * @param h5_filename
	 */
	H5CytoFrame(const string & fcs_filename, FCS_READ_PARAM & config, const string & h5_filename, unsigned int flags = H5F_ACC_RDWR):filename_(h5_filename), is_dirty_params(false), is_dirty_keys(false), is_dirty_pdata(false)
	{
		MemCytoFrame fr(fcs_filename, config);
		fr.read_fcs();
		fr.write_h5(h5_filename);
		*this = H5CytoFrame(h5_filename, flags);
	}
	/**
	 * constructor from the H5
	 * @param _filename H5 file path
	 */
	H5CytoFrame(const string & h5_filename, unsigned int flags = H5F_ACC_RDONLY):CytoFrame(),filename_(h5_filename), is_dirty_params(false), is_dirty_keys(false), is_dirty_pdata(false)
	{
		file.openFile(filename_, flags);
		load_meta();
		//open dataset for event data

		dataset = file.openDataSet(DATASET_NAME);
		dataspace = dataset.getSpace();
		dataspace.getSimpleExtentDims(dims);

	}
	/**
	 * abandon the changes to the meta data in cache by reloading them from disk
	 */
	void load_meta(){

		DataSet ds_param = file.openDataSet("params");
	//	DataType param_type = ds_param.getDataType();

		hsize_t dim_param[1];
		DataSpace dsp_param = ds_param.getSpace();
		dsp_param.getSimpleExtentDims(dim_param);
		int nParam = dim_param[0];

		StrType str_type(0, H5T_VARIABLE);	//define variable-length string data type

		FloatType datatype = h5_datatype_data(DataTypeLocation::MEM);

		/*
		 * read params as array of compound type
		 */


		hsize_t dim_pne[] = {2};
		ArrayType pne(datatype, 1, dim_pne);

		/*
		 * have to redefine cytoParam to use char * since H5 can't directly read into string member of the compound type
		 */
		struct cytoParam1{
			char * channel;
			char * marker;
			EVENT_DATA_TYPE min, max, PnG;
			EVENT_DATA_TYPE PnE[2];
			int PnB;

		};
		vector<cytoParam1> pvec(nParam);

		CompType param_type(sizeof(cytoParam1));
		param_type.insertMember("channel", HOFFSET(cytoParam1, channel), str_type);
		param_type.insertMember("marker", HOFFSET(cytoParam1, marker), str_type);
		param_type.insertMember("min", HOFFSET(cytoParam1, min), datatype);
		param_type.insertMember("max", HOFFSET(cytoParam1, max), datatype);
		param_type.insertMember("PnG", HOFFSET(cytoParam1, PnG), datatype);
		param_type.insertMember("PnE", HOFFSET(cytoParam1, PnE), pne);
		param_type.insertMember("PnB", HOFFSET(cytoParam1, PnB), PredType::NATIVE_INT8);
		ds_param.read(pvec.data(),param_type);
		//cp back to param
		params.resize(nParam);
		for(auto i = 0; i < nParam; i++)
		{
			params[i].channel = pvec[i].channel;
			delete [] pvec[i].channel;//reclaim vlen char
			params[i].marker = pvec[i].marker;
			delete [] pvec[i].marker;
			params[i].min = pvec[i].min;
			params[i].max = pvec[i].max;
			params[i].PnG = pvec[i].PnG;
			params[i].PnE[0] = pvec[i].PnE[0];
			params[i].PnE[1] = pvec[i].PnE[1];
			params[i].PnB = pvec[i].PnB;
		}
		build_hash();
		is_dirty_params = false;
		/*
		 * read keywords
		 */
		struct key_t{
				char * key;
				char * value;
	//			key_t(const char * k, const char * v):key(k),value(v){};
			};
		CompType key_type(sizeof(key_t));
		key_type.insertMember("key", HOFFSET(key_t, key), str_type);
		key_type.insertMember("value", HOFFSET(key_t, value), str_type);


		DataSet ds_key = file.openDataSet( "keywords");
		DataSpace dsp_key = ds_key.getSpace();
		hsize_t dim_key[1];
		dsp_key.getSimpleExtentDims(dim_key);
		int nKey = dim_key[0];

		vector<key_t> keyVec(nKey);
		ds_key.read(keyVec.data(), key_type);
//		keys.resize(nKey);
		for(auto i = 0; i < nKey; i++)
		{
			keys_[keyVec[i].key] = keyVec[i].value;
			delete [] keyVec[i].key;
			delete [] keyVec[i].value;

//			keys[i].first = keyVec[i].key;
//			delete [] keyVec[i].key;
//			keys[i].second = keyVec[i].value;
//			delete [] keyVec[i].value;
		}
		is_dirty_keys = false;
		/*
		 *
		 * read pdata
		 */
		DataSet ds_pd = file.openDataSet( "pdata");
		DataSpace dsp_pd = ds_pd.getSpace();
		hsize_t dim_pd[1];
		dsp_pd.getSimpleExtentDims(dim_pd);
		int nPd = dim_pd[0];

		keyVec.resize(nPd);
		ds_pd.read(keyVec.data(), key_type);
		for(auto i = 0; i < nPd; i++)
		{
			pheno_data_[keyVec[i].key] = keyVec[i].value;
			delete [] keyVec[i].key;
			delete [] keyVec[i].value;
		}
		is_dirty_pdata = false;

	}

	string get_h5_file_path() const{
		return filename_;
	}

	void convertToPb(pb::CytoFrame & fr_pb, const string & h5_filename, H5Option h5_opt) const
	{
		fr_pb.set_is_h5(true);

		if(!fs::equivalent(fs::path(filename_).parent_path(), fs::path(h5_filename).parent_path()))
		{
			switch(h5_opt)
			{
			case H5Option::copy:
				fs::copy(filename_, h5_filename);
				break;
			case H5Option::move:
				fs::rename(filename_, h5_filename);
				break;
			case H5Option::link:
				fs::create_hard_link(filename_, h5_filename);
				break;
			case H5Option::symlink:
				fs::create_symlink(filename_, h5_filename);
				break;
			case H5Option::skip:
				break;
			}
		}
	}

	/**
	 * Read data from disk.
	 * The caller will directly receive the data vector without the copy overhead thanks to the move semantics supported for vector container in c++11
	 * @return
	 */


	EVENT_DATA_VEC get_data() const{
		unsigned n = n_cols();
		uvec col_idx(n);
		for(unsigned i = 0; i < n; i++)
			col_idx[i] = i;
		return read_data(col_idx);
	}
	/**
	 * Partial IO
	 * @param col_idx
	 * @return
	 */
	EVENT_DATA_VEC get_data(uvec col_idx) const
	{
		return read_data(col_idx);
	}

	CytoFramePtr copy(const string & h5_filename = "") const
	{
		string new_filename = h5_filename;
		if(new_filename == "")
		{
			new_filename = generate_unique_filename(fs::temp_directory_path(), "", ".h5");
			fs::remove(new_filename);
		}
		fs::copy_file(filename_, new_filename);
		CytoFramePtr ptr(new H5CytoFrame(new_filename, H5F_ACC_RDWR));
		//copy cached meta
		ptr->set_params(get_params());
		ptr->set_keywords(get_keywords());
		ptr->set_pheno_data(get_pheno_data());
		return ptr;
	}
	CytoFramePtr copy_realized(uvec row_idx, uvec col_idx, const string & h5_filename = "") const
	{

		string new_filename = h5_filename;
		if(new_filename == "")
		{
			new_filename = generate_unique_filename(fs::temp_directory_path(), "", ".h5");
			fs::remove(new_filename);
		}		//if view is empty, then simply invoke copy
		if(row_idx.size() == 0 && col_idx.size() == 0)
			return copy(new_filename);
		//otherwise, realize it to memory and write back to new file
		MemCytoFrame fr(*this);
		fr.copy_realized(row_idx, col_idx)->write_h5(new_filename);//this flushes the meta data as well
		return CytoFramePtr(new H5CytoFrame(new_filename, H5F_ACC_RDWR));
	}

	/**
	 * copy setter
	 * @param _data
	 */
	void set_data(const EVENT_DATA_VEC & _data)
	{

		dataset.write(_data.mem, h5_datatype_data(DataTypeLocation::MEM));

	}

	void set_data(EVENT_DATA_VEC && _data)
	{
//		EVENT_DATA_VEC data = _data;
		set_data(_data);
	}
};

};



#endif /* INST_INCLUDE_CYTOLIB_H5CYTOFRAME_HPP_ */
