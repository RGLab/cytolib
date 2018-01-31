/*
 * H5CytoFrame.hpp
 *
 *  Created on: Sep 18, 2017
 *      Author: wjiang2
 */

#ifndef INST_INCLUDE_CYTOLIB_H5CYTOFRAME_HPP_
#define INST_INCLUDE_CYTOLIB_H5CYTOFRAME_HPP_
#include <cytolib/CytoFrame.hpp>


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
public:
	~H5CytoFrame(){};
	/**
	 * constructor from FCS
	 * @param fcs_filename
	 * @param h5_filename
	 */
	H5CytoFrame(const string & fcs_filename, FCS_READ_PARAM & config,  bool is_read_data, const string & h5_filename):filename_(h5_filename)
	{
		MemCytoFrame fr(fcs_filename, config);
		fr.read_fcs();
		fr.writeH5(h5_filename);
		*this = H5CytoFrame(h5_filename);
	}
	/**
	 * constructor from the H5
	 * @param _filename H5 file path
	 */
	H5CytoFrame(const string & h5_filename, unsigned int flags = H5F_ACC_RDONLY):filename_(h5_filename)
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
		buildHash();

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
			keys[keyVec[i].key] = keyVec[i].value;
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
			pd[keyVec[i].key] = keyVec[i].value;
			delete [] keyVec[i].key;
			delete [] keyVec[i].value;
		}
		//open dataset for event data

		dataset = file.openDataSet(DATASET_NAME);
		dataspace = dataset.getSpace();
		dataspace.getSimpleExtentDims(dims);

	}

	unsigned nRow() const{
		//read nEvents
		return dims[1];
	}
	/**
	 * Read data from disk.
	 * The caller will directly receive the data vector without the copy overhead thanks to the move semantics supported for vector container in c++11
	 * @return
	 */
	EVENT_DATA_VEC get_data() const{
		EVENT_DATA_VEC data(nCol() * nRow());
		dataset.read(data.data(), PredType::NATIVE_FLOAT);

		return data;
	}
	EVENT_DATA_VEC get_data(const string & colname, ColType type) const{
		int idx = getColId(colname, type);
		if(idx<0)
			throw(domain_error("colname not found: " + colname));
		EVENT_DATA_VEC data(nRow());
	//	hsize_t dim[1];
		return data;

	}
	/**
	 * copy setter
	 * @param _data
	 */
	void setData(const EVENT_DATA_VEC & _data)
	{

		dataset.write(_data.data(), PredType::NATIVE_FLOAT );

	}

	void setData(EVENT_DATA_VEC && _data)
	{
		EVENT_DATA_VEC data = _data;
		setData(data);
	}
};

};



#endif /* INST_INCLUDE_CYTOLIB_H5CYTOFRAME_HPP_ */
