#include <cytolib/H5CytoFrame.hpp>
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#include <cytolib/global.hpp>


namespace cytolib
{
	EVENT_DATA_VEC H5CytoFrame::read_data(uvec col_idx) const
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


	/*
	 * for simplicity, we don't want to handle the object that has all the h5 handler closed
	 * because it will require lots of validity checks before each disk IO operations
	 */
//	void close_h5(){
//		dataspace.close();
//		dataset.close();
//		file.close();
//	}
	void H5CytoFrame::flush_meta(){
		//flush the cached meta data from CytoFrame into h5
		if(is_dirty_params)
			flush_params();
		if(is_dirty_keys)
			flush_keys();
		if(is_dirty_pdata)
			flush_pheno_data();
	}
	void H5CytoFrame::flush_params()
	{
		CompType param_type = get_h5_datatype_params(DataTypeLocation::MEM);
		DataSet ds = file.openDataSet("params");
		hsize_t size[1] = {params.size()};
		ds.extend(size);
		auto params_char = params_c_str();

		ds.write(&params_char[0], param_type );
		ds.flush(H5F_SCOPE_LOCAL);
		is_dirty_params = false;
	}

	void H5CytoFrame::flush_keys()
	{
		CompType key_type = get_h5_datatype_keys();
		DataSet ds = file.openDataSet("keywords");
		auto keyVec = to_kw_vec<KEY_WORDS>(keys_);

		hsize_t size[1] = {keyVec.size()};
		ds.extend(size);
		ds.write(&keyVec[0], key_type );
		ds.flush(H5F_SCOPE_LOCAL);

		is_dirty_keys = false;
	}
	void H5CytoFrame::flush_pheno_data()
	{
		CompType key_type = get_h5_datatype_keys();
		DataSet ds = file.openDataSet("pdata");

		auto keyVec = to_kw_vec<PDATA>(pheno_data_);
		hsize_t size[1] = {keyVec.size()};
		ds.extend(size);

		ds.write(&keyVec[0], key_type );
		ds.flush(H5F_SCOPE_LOCAL);

		is_dirty_pdata = false;
	}

	H5CytoFrame::H5CytoFrame(const H5CytoFrame & frm):CytoFrame(frm)
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
	H5CytoFrame::H5CytoFrame(H5CytoFrame && frm):CytoFrame(frm)
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
	H5CytoFrame & H5CytoFrame::operator=(const H5CytoFrame & frm)
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
	H5CytoFrame & H5CytoFrame::operator=(H5CytoFrame && frm)
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
	/**
	 * constructor from FCS
	 * @param fcs_filename
	 * @param h5_filename
	 */
	H5CytoFrame::H5CytoFrame(const string & fcs_filename, FCS_READ_PARAM & config
			, const string & h5_filename, bool readonly):filename_(h5_filename), is_dirty_params(false), is_dirty_keys(false), is_dirty_pdata(false)
	{
		MemCytoFrame fr(fcs_filename, config);
		fr.read_fcs();
		fr.write_h5(h5_filename);
		*this = H5CytoFrame(h5_filename, readonly);
	}
	/**
	 * constructor from the H5
	 * @param _filename H5 file path
	 */
	H5CytoFrame::H5CytoFrame(const string & h5_filename, bool readonly):CytoFrame(readonly),filename_(h5_filename), is_dirty_params(false), is_dirty_keys(false), is_dirty_pdata(false)
	{


		file.openFile(filename_, H5F_ACC_RDWR);//always use the same flag and keep lock at cf level to avoid h5 open error caused conflicting h5 flags among cf objects that points to the same h5
		load_meta();


		//open dataset for event data

		dataset = file.openDataSet(DATASET_NAME);
		dataspace = dataset.getSpace();
		dataspace.getSimpleExtentDims(dims);

	}
	/**
	 * abandon the changes to the meta data in cache by reloading them from disk
	 */
	void H5CytoFrame::load_meta(){

		DataSet ds_param = file.openDataSet("params");
	//	DataType param_type = ds_param.getDataType();

		hsize_t dim_param[1];
		DataSpace dsp_param = ds_param.getSpace();
		dsp_param.getSimpleExtentDims(dim_param);
		int nParam = dim_param[0];

		StrType str_type(H5::PredType::C_S1, H5T_VARIABLE);	//define variable-length string data type

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
		vector<cytoParam> params(nParam);
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
		CytoFrame::set_params(params, true);
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

	void H5CytoFrame::convertToPb(pb::CytoFrame & fr_pb, const string & h5_filename, H5Option h5_opt) const
	{
		fr_pb.set_is_h5(true);
		if(h5_opt != H5Option::skip)
		{
			auto dest = fs::path(h5_filename).parent_path();
			if(!fs::exists(dest))
				throw(logic_error(dest.string() + "doesn't exist!"));

			if(!fs::equivalent(fs::path(filename_).parent_path(), dest))
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
				default:
					throw(logic_error("invalid h5_opt!"));
				}
			}
		}
	}

	/**
	 * Read data from disk.
	 * The caller will directly receive the data vector without the copy overhead thanks to the move semantics supported for vector container in c++11
	 * @return
	 */


	EVENT_DATA_VEC H5CytoFrame::get_data() const{
		unsigned n = n_cols();
		uvec col_idx(n);
		for(unsigned i = 0; i < n; i++)
			col_idx[i] = i;
		return read_data(col_idx);
	}

	CytoFramePtr H5CytoFrame::copy(const string & h5_filename) const
	{
		string new_filename = h5_filename;
		if(new_filename == "")
		{
			new_filename = generate_unique_filename(fs::temp_directory_path(), "", ".h5");
			fs::remove(new_filename);
		}
		fs::copy_file(filename_, new_filename);
		CytoFramePtr ptr(new H5CytoFrame(new_filename, false));
		//copy cached meta
		ptr->set_params(get_params());
		ptr->set_keywords(get_keywords());
		ptr->set_pheno_data(get_pheno_data());
		return ptr;
	}
	CytoFramePtr H5CytoFrame::copy_realized(uvec row_idx, uvec col_idx, const string & h5_filename) const
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
		return CytoFramePtr(new H5CytoFrame(new_filename, false));
	}

	/**
	 * copy setter
	 * @param _data
	 */
	void H5CytoFrame::set_data(const EVENT_DATA_VEC & _data)
	{
		check_write_permission();
		hsize_t dims_data[2] = {_data.n_cols, _data.n_rows};
		dataset.extend(dims_data);
		//refresh data space and dims
		dataspace = dataset.getSpace();
		dataspace.getSimpleExtentDims(dims);

		dataset.write(_data.mem, h5_datatype_data(DataTypeLocation::MEM));
		dataset.flush(H5F_SCOPE_LOCAL);

	}



};

