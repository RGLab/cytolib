// Copyright 2019 Fred Hutchinson Cancer Research Center
// See the included LICENSE file for details on the licence that is granted to the user of this software.
#include <cytolib/H5CytoFrame.hpp>

namespace cytolib
{
	EVENT_DATA_VEC H5CytoFrame::read_data(uvec col_idx) const
	{
		H5File file(filename_, h5_flags(), FileCreatPropList::DEFAULT, access_plist_);
		auto dataset = file.openDataSet(DATASET_NAME);
		auto dataspace = dataset.getSpace();

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
		check_write_permission();
		H5File file(filename_, h5_flags(), FileCreatPropList::DEFAULT, access_plist_);

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
		check_write_permission();
		H5File file(filename_, h5_flags(), FileCreatPropList::DEFAULT, access_plist_);
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
		check_write_permission();
		H5File file(filename_, h5_flags(), FileCreatPropList::DEFAULT, access_plist_);
		CompType key_type = get_h5_datatype_keys();
		DataSet ds = file.openDataSet("pdata");

		auto keyVec = to_kw_vec<PDATA>(pheno_data_);
		hsize_t size[1] = {keyVec.size()};
		ds.extend(size);

		ds.write(&keyVec[0], key_type );
		ds.flush(H5F_SCOPE_LOCAL);

		is_dirty_pdata = false;
	}



	/**
	 * abandon the changes to the meta data in cache by reloading them from disk
	 */
	void H5CytoFrame::load_meta(){
		H5File file(filename_, h5_flags(), FileCreatPropList::DEFAULT, access_plist_);
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
		CytoFrame::set_params(params);
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





	/**
	 * copy setter
	 * @param _data
	 */
	void H5CytoFrame::set_data(const EVENT_DATA_VEC & _data)
	{
		H5File file(filename_, h5_flags(), FileCreatPropList::DEFAULT, access_plist_);
		check_write_permission();
		hsize_t dims_data[2] = {_data.n_cols, _data.n_rows};
		auto dataset = file.openDataSet(DATASET_NAME);

		dataset.extend(dims_data);
		//refresh data space and dims
		auto dataspace = dataset.getSpace();
		dataspace.getSimpleExtentDims(dims);

		dataset.write(_data.mem, h5_datatype_data(DataTypeLocation::MEM));
		dataset.flush(H5F_SCOPE_LOCAL);

	}



};

