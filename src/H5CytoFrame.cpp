/*
 * H5CytoFrame.cpp
 *
 *  Created on: Sep 18, 2017
 *      Author: wjiang2
 */

#include <cytolib/H5CytoFrame.hpp>
void H5CytoFrame::compensate(const compensation &){

}

H5CytoFrame::H5CytoFrame(const string & _filename):filename(_filename)
{
	H5File file( filename, H5F_ACC_RDONLY);

	DataSet ds_param = file.openDataSet("params");
	DataType param_type = ds_param.getDataType();

	hsize_t dim_param[1];
	DataSpace dsp_param = ds_param.getSpace();
	dsp_param.getSimpleExtentDims(dim_param);
	CytoFrame::params.resize(dim_param[0]);
//	StrType str_type(0, H5T_VARIABLE);	//define variable-length string data type
//
//	FloatType datatype( PredType::NATIVE_FLOAT );
//	datatype.setOrder(is_host_big_endian()?H5T_ORDER_BE:H5T_ORDER_LE );

	/*
	 * read params as array of compound type
	 */

//
//	hsize_t dim_pne[] = {2};
//	ArrayType pne(datatype, 1, dim_pne);
//
//	CompType param_type(sizeof(cytoParam));
//	param_type.insertMember("channel", HOFFSET(cytoParam, channel), str_type);
//	param_type.insertMember("marker", HOFFSET(cytoParam, marker), str_type);
//	param_type.insertMember("min", HOFFSET(cytoParam, min), datatype);
//	param_type.insertMember("max", HOFFSET(cytoParam, max), datatype);
//	param_type.insertMember("PnG", HOFFSET(cytoParam, PnG), datatype);
//	param_type.insertMember("PnE", HOFFSET(cytoParam, PnE), pne);
//	param_type.insertMember("PnB", HOFFSET(cytoParam, PnB), PredType::NATIVE_INT8);

	ds_param.read(&CytoFrame::params[0],param_type);
	/*
	 * read keywords
	 */
	//convert to vector

//	struct key_t{
//		string key, value;
//		key_t(const string & k, const string & v):key(k),value(v){};
//	};
//	vector<key_t> keyVec;
//
//	CompType key_type(sizeof(key_t));
//	key_type.insertMember("key", HOFFSET(key_t, key), str_type);
//	key_type.insertMember("value", HOFFSET(key_t, value), str_type);
//
//
//	DataSet ds_key = file.createDataSet( "keywords", key_type, dsp_key);
//	ds_key.write(&keyVec[0], key_type );
//	for(std::pair<std::string, string> e : keys)
//	{
//			keyVec.push_back(key_t(e.first, e.second));
//		}
	hsize_t dimsf[2];              // dataset dimensions
	DataSet dataset = file.openDataSet(DATASET_NAME);
	DataSpace dataspace = dataset.getSpace();
	dataspace.getSimpleExtentDims(dimsf);
	nEvents = dimsf[1];
}

EVENT_DATA_TYPE * H5CytoFrame::getData(){
	EVENT_DATA_TYPE * data = new EVENT_DATA_TYPE[nCol() * nRow()];
	return data;
}
EVENT_DATA_TYPE * H5CytoFrame::getData(const string & colname, ColType type){
	int idx = getColId(colname, type);
	EVENT_DATA_TYPE * data = new EVENT_DATA_TYPE[nRow()];
	return data;
//	return data.get() + idx * nRow();
}
