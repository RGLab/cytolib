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

#include <H5Cpp.h>

namespace cytolib
{
enum class ColType {channel, marker, unknown};
enum class RangeType {instrument, data};
enum class FrameType {FCS, H5};
enum class H5Option {copy, move, skip, link, symlink};

typedef unordered_map<string, string> PDATA;

using namespace H5;

const H5std_string  DATASET_NAME( "data");



class CytoFrame;
typedef shared_ptr<CytoFrame> CytoFramePtr;

/**
 * The class representing a single FCS file
 */
class CytoFrame{
protected:
	PDATA pheno_data_;
	KEY_WORDS keys_;//keyword pairs parsed from FCS Text section
	vector<cytoParam> params;// parameters coerced from keywords and computed from data for quick query
	unordered_map<string, int> channel_vs_idx;//hash map for query by channel
	unordered_map<string, int> marker_vs_idx;//hash map for query by marker


	CytoFrame (){};
	virtual bool is_hashed() const
	{
		return channel_vs_idx.size()==n_cols();
	}

	/**
	 * build the hash map for channel and marker for the faster query
	 *
	 */
	virtual void build_hash()
	{
		channel_vs_idx.clear();
		marker_vs_idx.clear();
		for(unsigned i = 0; i < n_cols(); i++)
		{
			channel_vs_idx[params[i].channel] = i;
			marker_vs_idx[params[i].marker] = i;
		}
	}
public:
	virtual ~CytoFrame(){};

	CytoFrame(const CytoFrame & frm)
	{
//		cout << "copy CytoFrame member" << endl;
		pheno_data_ = frm.pheno_data_;
		keys_ = frm.keys_;
		params = frm.params;
		channel_vs_idx = frm.channel_vs_idx;
		marker_vs_idx = frm.marker_vs_idx;
	}

	CytoFrame & operator=(const CytoFrame & frm)
	{
		pheno_data_ = frm.pheno_data_;
		keys_ = frm.keys_;
		params = frm.params;
		channel_vs_idx = frm.channel_vs_idx;
		marker_vs_idx = frm.marker_vs_idx;
		return *this;

	}

	CytoFrame & operator=(CytoFrame && frm)
	{
		swap(pheno_data_, frm.pheno_data_);
		swap(keys_, frm.keys_);
		swap(params, frm.params);
		swap(channel_vs_idx, frm.channel_vs_idx);
		swap(marker_vs_idx, frm.marker_vs_idx);
		return *this;
	}


	CytoFrame(CytoFrame && frm)
	{
		swap(pheno_data_, frm.pheno_data_);
		swap(keys_, frm.keys_);
		swap(params, frm.params);
		swap(channel_vs_idx, frm.channel_vs_idx);
		swap(marker_vs_idx, frm.marker_vs_idx);
	}

	compensation get_compensation(const string & key = "SPILL")
	{
		compensation comp;

		if(keys_.find(key)!=keys_.end())
		{
			string val = keys_[key];
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

	virtual void convertToPb(pb::CytoFrame & fr_pb, const string & h5_filename, H5Option h5_opt) const = 0;


	virtual void compensate(const compensation & comp)
	{

		int nMarker = comp.marker.size();
		EVENT_DATA_VEC dat = get_data();
		arma::mat A(dat.memptr(), n_rows(), n_cols(), false, true);//point to the original data
//		A.rows(1,3).print("\ndata");
		mat B = comp.get_spillover_mat();
//		B.print("comp");
		B = inv(B);
//		B.print("comp inverse");
		uvec indices(nMarker);
		for(int i = 0; i < nMarker; i++)
		{
			int id = get_col_idx(comp.marker[i], ColType::channel);
			if(id < 0)
				throw(domain_error("compensation parameter '" + comp.marker[i] + " not found in cytoframe parameters!"));

			indices[i] = id;
		}
//		uvec rind={1,2};
//		A.submat(rind,indices).print("data ");
		A.cols(indices) = A.cols(indices) * B;
//		A.submat(rind,indices).print("data comp");
		set_data(dat);
	}
	/**
	 * getter from cytoParam vector
	 * @return
	 */
	const vector<cytoParam> & get_params() const
	{
		return params;
	}

	void set_params(const vector<cytoParam> & _params)
	{
		params = _params;
	}
//	virtual void writeFCS(const string & filename);
	/**
	 * save the CytoFrame as HDF5 format
	 *
	 * @param filename the path of the output H5 file
	 */
	virtual void write_h5(const string & filename) const
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


		hsize_t dim_param[] = {n_cols()};
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
		for(const auto & e : keys_)
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
		for(std::pair<std::string, string> e : pheno_data_)
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
		unsigned nEvents = n_rows();
		hsize_t dimsf[2] = {n_cols(), nEvents};              // dataset dimensions
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
		dataset.write(&get_data()[0], PredType::NATIVE_FLOAT );
	}

	/**
	 * get the data of entire event matrix
	 * @return
	 */
	virtual EVENT_DATA_VEC get_data() const=0;
	virtual EVENT_DATA_VEC get_data(uvec col_idx) const=0;
	virtual EVENT_DATA_VEC get_data(vector<string>cols, ColType col_type) const
	{
		return get_data(get_col_idx(cols, col_type));
	}

	virtual void set_data(const EVENT_DATA_VEC &)=0;
	virtual void set_data(EVENT_DATA_VEC &&)=0;
	/**
	 * extract all the keyword pairs
	 *
	 * @return a vector of pairs of strings
	 */
	 virtual const KEY_WORDS & get_keywords() const{
		return keys_;
	}
	/**
	 * extract the value of the single keyword by keyword name
	 *
	 * @param key keyword name
	 * @return keyword value as a string
	 */
	virtual string get_keyword(const string & key) const
	{
		string res="";
		auto it = keys_.find(key);
		if(it!=keys_.end())
			res = it->second;
		return res;
	}

	/**
	 * set the value of the single keyword
	 * @param key keyword name
	 * @param value keyword value
	 */
	virtual void set_keyword(const string & key, const string & value)
	{
		keys_[key] = value;
	}

	/**
	 * get the number of columns(or parameters)
	 *
	 * @return
	 */
	virtual unsigned n_cols() const
	{
		return params.size();
	}

	/**
	 * get the number of rows(or events)
	 * @return
	 */
	virtual unsigned n_rows() const=0;

	virtual void subset_parameters(uvec col_idx)
	{
			unsigned n = col_idx.size();
			vector<cytoParam> params_new(n);
			for(unsigned i = 0; i < n; i++)
			{
					params_new[i] = params[col_idx[i]];
			}
			params = params_new;
			build_hash();//update idx table

//			set_data(get_data(col_idx));
			//update keywords PnX
//          for(unsigned i = 0; i < n; i++)
//          {
//                  if(it.first)
//          }
	}
	/**
	 * check if the hash map for channel and marker has been built
	 * @return
	 */

	/**
	 * get all the channel names
	 * @return
	 */
	virtual vector<string> get_channels() const
	{
		vector<string> res(n_cols());
		for(unsigned i = 0; i < n_cols(); i++)
			res[i] = params[i].channel;
		return res;
	}

	virtual void update_channels(const CHANNEL_MAP & chnl_map)
	{
		for(auto & it : chnl_map)
		{
			set_channel(it.first, it.second);
		}

	}
	/**
	 * get all the marker names
	 * @return
	 */
	virtual vector<string> get_markers() const
	{
		vector<string> res(n_cols());
			for(unsigned i = 0; i < n_cols(); i++)
				res[i] = params[i].marker;
		return res;
	}

	/**
	 * get the numeric index for the given column
	 * @param colname column name
	 * @param type the type of column
	 * @return
	 */
	virtual int get_col_idx(const string & colname, ColType type) const
	{
		if(!is_hashed())
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
	uvec get_col_idx(vector<string> colnames, ColType col_type) const
	{

		unsigned n = colnames.size();
		uvec col_idx(n);
		for(unsigned i = 0; i < n; i++)
		{
			int idx = get_col_idx(colnames[i], col_type);
			if(idx<0)
				throw(domain_error("colname not found: " + colnames[i]));
			col_idx[i] = idx;

		}
		return col_idx;
	}
	virtual void set_channel(const string & oldname, const string &newname, bool is_update_keywords = true)
	{
		int id = get_col_idx(oldname, ColType::channel);
		if(id<0)
			throw(domain_error("colname not found: " + oldname));
		if(g_loglevel>=GATING_HIERARCHY_LEVEL)
			PRINT(oldname + "-->"  + newname + "\n");
		params[id].channel=newname;
		channel_vs_idx.erase(oldname);
		channel_vs_idx[newname] = id;

		//update keywords(linear time, not sure how to improve it other than optionally skip it
		if(is_update_keywords)
		{
			for(auto & it : keys_)
				if(it.second == oldname)
					it.second = newname;
		}

	}

	virtual void set_marker(const string & oldname, const string & newname)
	{
		int id = get_col_idx(oldname, ColType::marker);
		if(id<0)
			throw(domain_error("colname not found: " + oldname));
		params[id].marker=newname;
		marker_vs_idx.erase(oldname);
		marker_vs_idx[newname] = id;
	}

	/**
	 * Update the instrument range (typically after data transformation)
	 * @param colname
	 * @param ctype
	 * @param new_range
	 */
	virtual void set_range(const string & colname, ColType ctype, pair<EVENT_DATA_TYPE, EVENT_DATA_TYPE> new_range){
		int idx = get_col_idx(colname, ctype);
		if(idx<0)
			throw(domain_error("colname not found: " + colname));
		params[idx].min = new_range.first;
		params[idx].max = new_range.second;
	}
	/**
	 * the range of a specific column
	 * @param colname
	 * @param ctype the type of column
	 * @param rtype either RangeType::data or RangeType::instrument
	 * @return
	 */
	virtual pair<EVENT_DATA_TYPE, EVENT_DATA_TYPE> get_range(const string & colname, ColType ctype, RangeType rtype) const
	{

		switch(rtype)
		{
		case RangeType::data:
			{

				EVENT_DATA_VEC vec = get_data({colname}, ctype);

				return make_pair(vec.min(), vec.max());
			}
		case RangeType::instrument:
		{
			int idx = get_col_idx(colname, ctype);
			if(idx<0)
				throw(domain_error("colname not found: " + colname));
			return make_pair(params[idx].min, params[idx].max);
		}
		default:
			throw(domain_error("invalid range type"));
		}
	}
	/**
	 * Compute the time step from keyword either "$TIMESTEP" or "$BTIM", "$TIMESTEP" is preferred when it is present
	 * This is used to convert time channel to the meaningful units later on during data transforming
	 * @param
	 * @return
	 */
	EVENT_DATA_TYPE get_time_step(const string time_channel)
	{

	  //check if $TIMESTEP is available
		EVENT_DATA_TYPE ts;
		auto it_time = keys_.find("$TIMESTEP");
		if(it_time != keys_.end())
				ts = boost::lexical_cast<EVENT_DATA_TYPE>(it_time->second);
		else
		{
		  auto it_btime = keys_.find("$BTIM");
		  auto it_etime = keys_.find("$ETIM");
		  if(it_btime == keys_.end() || it_etime == keys_.end())
			  ts = 1;
		  else
		  {
			  TM_ext btime = parse_time_with_fractional_seconds(it_btime->second);
			  TM_ext etime = parse_time_with_fractional_seconds(it_etime->second);

			  ts = difftime(mktime(&btime.time),mktime(&btime.time));
			  ts = ts + etime.fractional_secs/100 - btime.fractional_secs/100;

			  const auto time_range = get_range(time_channel, ColType::channel, RangeType::data);
			  ts /= (time_range.second - time_range.first);
//	      as.numeric(time.total)/diff(unit.range)

	    }
	  }
		return ts;
	}



	virtual CytoFramePtr copy(const string & h5_filename = "") const=0;
	virtual string get_h5_file_path() const=0;
	virtual CytoFramePtr copy(uvec row_idx, uvec col_idx, const string & h5_filename = "") const=0;

	const PDATA & get_pheno_data() const {return pheno_data_;}
	string get_pheno_data(const string & name) const {
		auto it = pheno_data_.find(name);
		if(it==pheno_data_.end())
			return "";
		else
			return it->second;

	}
	void set_pheno_data(const string & name, const string & value){
		pheno_data_[name] = value;
	}
	void set_pheno_data(const PDATA & _pd)
	{
		pheno_data_ = _pd;
	}
	void del_pheno_data(const string & name){pheno_data_.erase(name);}
};
};


#endif /* INST_INCLUDE_CYTOLIB_CYTOFRAME_HPP_ */
