/*
 * CytoFrame.hpp
 *
 *  Created on: Sep 18, 2017
 *      Author: wjiang2
 */

#ifndef INST_INCLUDE_CYTOLIB_CYTOFRAME_HPP_
#define INST_INCLUDE_CYTOLIB_CYTOFRAME_HPP_

#include "readFCSHeader.hpp"
#include "compensation.hpp"

#include <c++/H5Cpp.h>

namespace cytolib
{
enum class ColType {channel, marker, unknown};
enum class RangeType {instrument, data};
enum class FrameType {FCS, H5};

typedef unordered_map<string, string> PDATA;

using namespace H5;

const H5std_string  DATASET_NAME( "data");


/**
 * The class representing a single FCS file
 */
class CytoFrame{
protected:
	PDATA pd;
	KEY_WORDS keys;//keyword pairs parsed from FCS Text section
	vector<cytoParam> params;// parameters coerced from keywords and computed from data for quick query
	unordered_map<string, int> channel_vs_idx;//hash map for query by channel
	unordered_map<string, int> marker_vs_idx;//hash map for query by marker
public:
	virtual ~CytoFrame(){};

	compensation get_compensation(const string & key = "SPILL")
	{
		compensation comp;

		if(keys.find(key)!=keys.end())
		{
			string val = keys[key];
			vector<string> valVec;
			boost::split(valVec, val, boost::is_any_of(","));
			int n = boost::lexical_cast<int>(valVec[0]);
			if(n > 0)
			{
				comp.spillOver.resize(n*n);
				comp.marker.resize(n);
				for(int i = 0; i < n; i++)//param name
				{
					comp.marker[i] = valVec[i+1];
				}
				for(unsigned i = n + 1; i < valVec.size(); i++)//param name
				{
					comp.spillOver[i-n-1] = boost::lexical_cast<double>(valVec[i]);
				}

			}

		}
		return comp;
	}


	virtual void compensate(const compensation & comp)
	{

		int nMarker = comp.marker.size();
		EVENT_DATA_VEC dat = getData();
		arma::mat A(dat.data(), nRow(), nCol(), false, true);//point to the original data
//		A.rows(1,3).print("\ndata");
		mat B = comp.get_spillover_mat();
//		B.print("comp");
		B = inv(B);
//		B.print("comp inverse");
		uvec indices(nMarker);
		for(int i = 0; i < nMarker; i++)
		{
			int id = getColId(comp.marker[i], ColType::channel);
			if(id < 0)
				throw(domain_error(params[i].channel + " not found in compensation parameters!"));

			indices[i] = id;
		}
//		uvec rind={1,2};
//		A.submat(rind,indices).print("data ");
		A.cols(indices) = A.cols(indices) * B;
//		A.submat(rind,indices).print("data comp");
		setData(dat);
	}
	/**
	 * getter from cytoParam vector
	 * @return
	 */
	const vector<cytoParam> & getParams () const
	{
		return params;
	}
//	virtual void writeFCS(const string & filename);
	/**
	 * save the CytoFrame as HDF5 format
	 *
	 * @param filename the path of the output H5 file
	 */
	virtual void writeH5(const string & filename)
	{
		H5File file( filename, H5F_ACC_TRUNC );
		StrType str_type(0, H5T_VARIABLE);	//define variable-length string data type

		FloatType datatype( PredType::NATIVE_FLOAT );
		datatype.setOrder(is_host_big_endian()?H5T_ORDER_BE:H5T_ORDER_LE );

		/*
		 * write params as array of compound type
		 */


		hsize_t dim_pne[] = {2};
		ArrayType pne(datatype, 1, dim_pne);

		CompType param_type(sizeof(cytoParam));
		param_type.insertMember("channel", HOFFSET(cytoParam, channel), str_type);
		param_type.insertMember("marker", HOFFSET(cytoParam, marker), str_type);
		param_type.insertMember("min", HOFFSET(cytoParam, min), datatype);
		param_type.insertMember("max", HOFFSET(cytoParam, max), datatype);
		param_type.insertMember("PnG", HOFFSET(cytoParam, PnG), datatype);
		param_type.insertMember("PnE", HOFFSET(cytoParam, PnE), pne);
		param_type.insertMember("PnB", HOFFSET(cytoParam, PnB), PredType::NATIVE_INT8);

		hsize_t dim_param[] = {nCol()};
		DataSpace dsp_param(1, dim_param);
	//	vector<const char *> cvec;
	//	for(auto & c : params)
	//	{
	//		cvec.push_back(c.channel.c_str());
	////		cout << c.channel << ":" <<cvec.back() << endl;
	//	}
	//
	//	for(auto c : cvec)
	//		cout << c << endl;
		DataSet ds_param = file.createDataSet( "params", param_type, dsp_param);
		ds_param.write(params.data(), param_type );

		/*
		 * write keywords
		 */
		//convert to vector

		struct key_t{
			string key, value;
			key_t(const string & k, const string & v):key(k),value(v){};
		};
		vector<key_t> keyVec;
		for(std::pair<std::string, string> e : keys)
		{
			keyVec.push_back(key_t(e.first, e.second));
		}


		CompType key_type(sizeof(key_t));
		key_type.insertMember("key", HOFFSET(key_t, key), str_type);
		key_type.insertMember("value", HOFFSET(key_t, value), str_type);

		hsize_t dim_key[] = {keyVec.size()};
		DataSpace dsp_key(1, dim_key);

		DataSet ds_key = file.createDataSet( "keywords", key_type, dsp_key);
		ds_key.write(&keyVec[0], key_type );

		/*
		 * write pdata
		 */
		keyVec.clear();//reuse keyVec and key_type for pd
		for(std::pair<std::string, string> e : pd)
		{
			keyVec.push_back(key_t(e.first, e.second));
		}


		hsize_t dim_pd[] = {keyVec.size()};
		DataSpace dsp_pd(1, dim_pd);

		DataSet ds_pd = file.createDataSet( "pdata", key_type, dsp_pd);
		ds_pd.write(&keyVec[0], key_type );

		 /*
		* store events data as fixed
		* size dataset.
		*/
		unsigned nEvents = nRow();
		hsize_t dimsf[2] = {nCol(), nEvents};              // dataset dimensions
		DSetCreatPropList plist;
		hsize_t	chunk_dims[2] = {1, nEvents};
		plist.setChunk(2, chunk_dims);
	//	plist.setFilter()
		DataSpace dataspace( 2, dimsf);
		DataSet dataset = file.createDataSet( DATASET_NAME, datatype, dataspace, plist);
		/*
		* Write the data to the dataset using default memory space, file
		* space, and transfer properties.
		*/
		dataset.write(&getData()[0], PredType::NATIVE_FLOAT );
	}

	/**
	 * get the data of entire event matrix
	 * @return
	 */
	virtual EVENT_DATA_VEC getData() const=0;
	/**
	 * get the data for the single channel
	 *
	 * @param colname the channel for marker name
	 * @param type enum class indicates the type of colname, can be either ColType::channel or ColType::marker or ColType::unknown
	 * when ColType::unknown, both types will be tried for the column match.
	 * @return
	 */
	virtual EVENT_DATA_VEC getData(const string & colname, ColType type) const=0;
	virtual void setData(const EVENT_DATA_VEC &)=0;
	virtual void setData(EVENT_DATA_VEC &&)=0;
	/**
	 * extract all the keyword pairs
	 *
	 * @return a vector of pairs of strings
	 */
	 virtual const KW_PAIR & getKeywords() const{
		return keys.getPairs();
	}
	/**
	 * extract the value of the single keyword by keyword name
	 *
	 * @param key keyword name
	 * @return keyword value as a string
	 */
	virtual string getKeyword(const string & key) const
	{
		string res="";
		auto it = keys.find(key);
		if(it!=keys.end())
			res = it->second;
		return res;
	}

	/**
	 * set the value of the single keyword
	 * @param key keyword name
	 * @param value keyword value
	 */
	virtual void setKeyword(const string & key, const string & value)
	{
		keys[key] = value;
	}

	/**
	 * get the number of columns(or parameters)
	 *
	 * @return
	 */
	virtual unsigned nCol() const
	{
		return params.size();
	}

	/**
	 * get the number of rows(or events)
	 * @return
	 */
	virtual unsigned nRow() const=0;
	/**
	 * check if the hash map for channel and marker has been built
	 * @return
	 */
	virtual bool isHashed() const
	{
		return channel_vs_idx.size()==nCol();
	}

	/**
	 * build the hash map for channel and marker for the faster query
	 *
	 */
	virtual void buildHash()
	{
		for(unsigned i = 0; i < nCol(); i++)
		{
			channel_vs_idx[params[i].channel] = i;
			marker_vs_idx[params[i].marker] = i;
		}
	}

	/**
	 * get all the channel names
	 * @return
	 */
	virtual vector<string> getChannels() const
	{
		vector<string> res(nCol());
		for(unsigned i = 0; i < nCol(); i++)
			res[i] = params[i].channel;
		return res;
	}

	virtual void updateChannels(const CHANNEL_MAP & chnl_map)
	{
		for(auto & it : chnl_map)
		{
			setChannel(it.first, it.second);
		}

	}
	/**
	 * get all the marker names
	 * @return
	 */
	virtual vector<string> getMarkers() const
	{
		vector<string> res(nCol());
			for(unsigned i = 0; i < nCol(); i++)
				res[i] = params[i].marker;
		return res;
	}

	/**
	 * get the numeric index for the given column
	 * @param colname column name
	 * @param type the type of column
	 * @return
	 */
	virtual int getColId(const string & colname, ColType type) const
	{
		if(!isHashed())
			throw(domain_error("please call buildHash() first to build the hash map for column index!"));

		switch(type)
		{
		case ColType::channel:
			{
				auto it1 = channel_vs_idx.find(colname);
				if(it1==channel_vs_idx.end())
					return -1;
				else
					return it1->second;

			}
		case ColType::marker:
			{
				auto it2 = marker_vs_idx.find(colname);
				if(it2==marker_vs_idx.end())
					return -1;
				else
					return it2->second;
			}
		case ColType::unknown:
			{
				auto it1 = channel_vs_idx.find(colname);
				auto it2 = marker_vs_idx.find(colname);
				if(it1==channel_vs_idx.end()&&it2==marker_vs_idx.end())
					return -1;
				else if(it1!=channel_vs_idx.end()&&it2!=marker_vs_idx.end())
					throw(domain_error("ambiguous colname without colType: " + colname ));
				else if(it1!=channel_vs_idx.end())
					return it1->second;
				else
					return it2->second;
			}
		default:
			throw(domain_error("invalid col type"));
		}

	}

	virtual void setChannel(const string & oldname, const string &newname)
	{
		int id = getColId(oldname, ColType::channel);
		if(id<0)
			throw(domain_error("colname not found: " + oldname));
		if(g_loglevel>=GATING_HIERARCHY_LEVEL)
			PRINT(oldname + "-->"  + newname + "\n");
		params[id].channel=newname;
		channel_vs_idx.erase(oldname);
		channel_vs_idx[newname] = id;
	}

	virtual void setMarker(const string & oldname, const string & newname)
	{
		int id = getColId(oldname, ColType::marker);
		if(id<0)
			throw(domain_error("colname not found: " + oldname));
		params[id].marker=newname;
		marker_vs_idx.erase(oldname);
		marker_vs_idx[newname] = id;
	}

	/**
	 * the range of a specific column
	 * @param colname
	 * @param ctype the type of column
	 * @param rtype either RangeType::data or RangeType::instrument
	 * @return
	 */
	virtual pair<EVENT_DATA_TYPE, EVENT_DATA_TYPE> getRange(const string & colname, ColType ctype, RangeType rtype) const
	{

		switch(rtype)
		{
		case RangeType::data:
			{

				EVENT_DATA_VEC vec = getData(colname, ctype);
				EVENT_DATA_TYPE * data = &vec[0];
				auto res = minmax_element(data, data + nRow());
				return make_pair(*res.first, *res.second);
			}
		case RangeType::instrument:
		{
			int idx = getColId(colname, ctype);
			if(idx<0)
				throw(domain_error("colname not found: " + colname));
			return make_pair(params[idx].min, params[idx].max);
		}
		default:
			throw(domain_error("invalid range type"));
		}
	}

	const PDATA & getPData() const {return pd;}
	string getPData(const string & name) const {
		auto it = pd.find(name);
		if(it==pd.end())
			return "";
		else
			return it->second;

	}
	void setPData(const string & name, const string & value){
		pd[name] = value;
	}
	void delPData(const string & name){pd.erase(name);}
};
};


#endif /* INST_INCLUDE_CYTOLIB_CYTOFRAME_HPP_ */
