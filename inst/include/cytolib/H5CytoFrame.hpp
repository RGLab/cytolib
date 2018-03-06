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
	/*
	 * these H5 handlers remain open during the life cycle of H5CytoFrame
	 * for faster accessing the data
	 */
	H5File file;
	DataSet dataset;
	DataSpace dataspace;
	hsize_t dims[2];              // dataset dimensions
	arma::uvec row_idx_;
	bool is_row_indexed;

public:
	~H5CytoFrame(){};
	H5CytoFrame():is_row_indexed(false){};
	/**
	 * constructor from FCS
	 * @param fcs_filename
	 * @param h5_filename
	 */
	H5CytoFrame(const string & fcs_filename, FCS_READ_PARAM & config, const string & h5_filename):filename_(h5_filename)
	{
		MemCytoFrame fr(fcs_filename, config);
		fr.read_fcs();
		fr.write_h5(h5_filename);
		*this = H5CytoFrame(h5_filename);
	}
	/**
	 * constructor from the H5
	 * @param _filename H5 file path
	 */
	H5CytoFrame(const string & h5_filename, unsigned int flags = H5F_ACC_RDONLY):filename_(h5_filename),is_row_indexed(false)
	{
		file.openFile(filename_, flags);

		DataSet ds_param = file.openDataSet("params");
	//	DataType param_type = ds_param.getDataType();

		hsize_t dim_param[1];
		DataSpace dsp_param = ds_param.getSpace();
		dsp_param.getSimpleExtentDims(dim_param);
		int nParam = dim_param[0];

		StrType str_type(0, H5T_VARIABLE);	//define variable-length string data type

		FloatType datatype( PredType::NATIVE_FLOAT );
		datatype.setOrder(is_host_big_endian()?H5T_ORDER_BE:H5T_ORDER_LE );

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
			unsigned original_col_idx;

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
		param_type.insertMember("original_col_idx", HOFFSET(cytoParam1, original_col_idx), PredType::NATIVE_INT8);
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
			params[i].original_col_idx = pvec[i].original_col_idx;
		}
		build_hash();

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
		/*
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
		//open dataset for event data

		dataset = file.openDataSet(DATASET_NAME);
		dataspace = dataset.getSpace();
		dataspace.getSimpleExtentDims(dims);

	}

	const string & get_h5_file_path(){
		return filename_;
	}
	unsigned n_rows() const{
		//read nEvents
		if(is_row_indexed)
			return row_idx_.size();
		else
			return dims[1];
	}

	void convertToPb(pb::CytoFrame & fr_pb, const string & h5_filename, H5Option h5_opt)
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
		unsigned nrow = n_rows();
		unsigned ncol = n_cols();
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
			unsigned idx = params[i].original_col_idx;
			hsize_t      offset[] = {idx, 0};   // hyperslab offset in the file
			hsize_t      count[] = {1, nrow};    // size of the hyperslab in the file
			dataspace.selectHyperslab( H5S_SELECT_SET, count, offset );

			//select slab for mem
			offset_mem[0] = i;
			offset_mem[1] = 0;
			count_mem[0]  = 1;
			count_mem[1]  = nrow;
			memspace.selectHyperslab( H5S_SELECT_SET, count_mem, offset_mem );

			dataset.read(data.memptr(), PredType::NATIVE_FLOAT ,memspace, dataspace);
		}

		if(is_row_indexed)
			data = data(row_idx_);
		return data;
	}

	/*
	 * This overloaded function exists so that data subsetting can be operated directly on abstract CytoFrame object
	 */
	EVENT_DATA_VEC get_data(vector<string> colnames, ColType col_type) const
	{
		return cols(colnames, col_type).get_data();
	}

	void cols_(vector<string> colnames, ColType col_type)
	{

		uvec col_idx = get_col_idx(colnames, col_type);

		//update params
		CytoFrame::cols_(col_idx);
	}


	void rows_(vector<unsigned> row_idx)
	{
		if(is_row_indexed && row_idx.size()!=row_idx_.size())
			throw(domain_error("The size of the new row index is not the same as the total number of events!"));
		row_idx_ = arma::conv_to<uvec>::from(row_idx);
	}

	H5CytoFrame rows(vector<unsigned> row_idx) const
	{
		H5CytoFrame res(*this);
		res.rows_(row_idx);
		return res;
	}

	H5CytoFrame cols(vector<string> colnames, ColType col_type) const
	{
		H5CytoFrame res(*this);
		res.cols_(colnames, col_type);
		return res;
	}

	/**
	 * Doesn't create the H5 copy so the new H5CytoFrame shares the same H5 data with the original object.
	 * @return
	 */
	CytoFrame * copy() const
	{
		return new H5CytoFrame(*this);
	}
	/**
	 * copy setter
	 * @param _data
	 */
	void set_data(const EVENT_DATA_VEC & _data)
	{

		dataset.write(_data.mem, PredType::NATIVE_FLOAT );

	}

	void set_data(EVENT_DATA_VEC && _data)
	{
//		EVENT_DATA_VEC data = _data;
		set_data(_data);
	}
};

};



#endif /* INST_INCLUDE_CYTOLIB_H5CYTOFRAME_HPP_ */
