#include <cytolib/CytoFrame.hpp>
#include <boost/lexical_cast.hpp>
#include <cytolib/global.hpp>


namespace cytolib
{

	/**
	 * build the hash map for channel and marker for the faster query
	 *
	 */
	void CytoFrame::build_hash()
	{
		channel_vs_idx.clear();
		marker_vs_idx.clear();
		for(unsigned i = 0; i < n_cols(); i++)
		{
			channel_vs_idx[params[i].channel] = i;
			marker_vs_idx[params[i].marker] = i;
		}
	}
	void CytoFrame::check_write_permission(){
		if(readonly_)
			throw(domain_error("Can't modify the read-only CytoFrame object!"));

	}
//	void close_h5() =0;
	CytoFrame::CytoFrame(const CytoFrame & frm)
	{
//		cout << "copy CytoFrame member" << endl;
		pheno_data_ = frm.pheno_data_;
		keys_ = frm.keys_;
		params = frm.params;
		channel_vs_idx = frm.channel_vs_idx;
		marker_vs_idx = frm.marker_vs_idx;
		readonly_ = frm.readonly_;
	}

	CytoFrame & CytoFrame::operator=(const CytoFrame & frm)
	{
		pheno_data_ = frm.pheno_data_;
		keys_ = frm.keys_;
		params = frm.params;
		channel_vs_idx = frm.channel_vs_idx;
		marker_vs_idx = frm.marker_vs_idx;
		readonly_ = frm.readonly_;
		return *this;

	}

	CytoFrame & CytoFrame::operator=(CytoFrame && frm)
	{
		swap(pheno_data_, frm.pheno_data_);
		swap(keys_, frm.keys_);
		swap(params, frm.params);
		swap(channel_vs_idx, frm.channel_vs_idx);
		swap(marker_vs_idx, frm.marker_vs_idx);
		swap(readonly_, frm.readonly_);
		return *this;
	}


	CytoFrame::CytoFrame(CytoFrame && frm)
	{
		swap(pheno_data_, frm.pheno_data_);
		swap(keys_, frm.keys_);
		swap(params, frm.params);
		swap(channel_vs_idx, frm.channel_vs_idx);
		swap(readonly_, frm.readonly_);
	}

	compensation CytoFrame::get_compensation(const string & key)
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



	void CytoFrame::compensate(const compensation & comp)
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

	void CytoFrame::scale_time_channel(string time_channel){
		auto idx = get_col_idx(time_channel, ColType::channel);
		if(idx >= 0){
			//special treatment for time channel
			EVENT_DATA_TYPE timestep = get_time_step(time_channel);
			if(g_loglevel>=GATING_HIERARCHY_LEVEL)
				PRINT("multiplying "+time_channel+" by :"+ to_string(timestep) + "\n");
			//TODO:partial IO
			EVENT_DATA_VEC data = get_data();;//get_data(idx);
			EVENT_DATA_TYPE * x = data.colptr(idx);
			int nEvents = n_rows();
			for(int i = 0; i < nEvents; i++)
				x[i] = x[i] * timestep;
			set_data(data);
			//TODO:update instrument range as well
			auto param_range = get_range(time_channel, ColType::channel, RangeType::data);
			param_range.first = param_range.first * timestep;
			param_range.second = param_range.second * timestep;
			set_range(time_channel, ColType::channel, param_range);
		}

	}


	void CytoFrame::set_params(const vector<cytoParam> & _params, bool force)
	{
		if(!force)
			check_write_permission();
		params = _params;
		build_hash();//update idx table

	}
//	void writeFCS(const string & filename);

	FloatType CytoFrame::h5_datatype_data(DataTypeLocation storage_type) const
	{
		if(storage_type)
		{
			FloatType datatype( PredType::NATIVE_FLOAT );
			datatype.setOrder(is_host_big_endian()?H5T_ORDER_BE:H5T_ORDER_LE );
			return datatype;
		}
		else
			return FloatType(PredType::NATIVE_DOUBLE);
	}
	CompType CytoFrame::get_h5_datatype_params(DataTypeLocation storage_type) const
	{
		StrType str_type(H5::PredType::C_S1, H5T_VARIABLE);	//define variable-length string data type

		/*
		 * define params as array of compound type
		 */

		FloatType datatype= h5_datatype_data(storage_type);
		hsize_t dim_pne[] = {2};
		ArrayType pne(datatype, 1, dim_pne);

		CompType param_type(sizeof(cytoParam_cstr));
		param_type.insertMember("channel", HOFFSET(cytoParam_cstr, channel), str_type);
		param_type.insertMember("marker", HOFFSET(cytoParam_cstr, marker), str_type);
		param_type.insertMember("min", HOFFSET(cytoParam_cstr, min), datatype);
		param_type.insertMember("max", HOFFSET(cytoParam_cstr, max), datatype);
		param_type.insertMember("PnG", HOFFSET(cytoParam_cstr, PnG), datatype);
		param_type.insertMember("PnE", HOFFSET(cytoParam_cstr, PnE), pne);
		param_type.insertMember("PnB", HOFFSET(cytoParam_cstr, PnB), PredType::NATIVE_INT8);

		return param_type;

	}
	CompType CytoFrame::get_h5_datatype_keys() const
	{
		StrType str_type(H5::PredType::C_S1, H5T_VARIABLE);	//define variable-length string data type

		CompType key_type(sizeof(KEY_WORDS_SIMPLE));
		key_type.insertMember("key", HOFFSET(KEY_WORDS_SIMPLE, key), str_type);
		key_type.insertMember("value", HOFFSET(KEY_WORDS_SIMPLE, value), str_type);
		return key_type;

	}
	void CytoFrame::write_h5_params(H5File file) const
	{
		hsize_t dim_param[] = {n_cols()};
		hsize_t dim_max[] = {H5S_UNLIMITED};
		DataSpace dsp_param(1, dim_param, dim_max);
		DSetCreatPropList plist;
		plist.setChunk(1, dim_param);
		DataSet ds = file.createDataSet( "params", get_h5_datatype_params(DataTypeLocation::H5), dsp_param, plist);
		auto params_char = params_c_str();
		ds.write(&params_char[0], get_h5_datatype_params(DataTypeLocation::MEM));
	}
	/**
	 * Convert string to cstr in params for writing to h5
	 * @return
	 */
	vector<cytoParam_cstr> CytoFrame::params_c_str() const{
		auto nParams = params.size();
		vector<cytoParam_cstr> res(nParams);
		for(unsigned i = 0; i < nParams; i++)
		{
			res[i].channel = params[i].channel.c_str();
			res[i].marker = params[i].marker.c_str();
			res[i].min = params[i].min;
			res[i].max = params[i].max;
			res[i].PnG = params[i].PnG;
			res[i].PnE[0] = params[i].PnE[0];
			res[i].PnE[1] = params[i].PnE[1];
			res[i].PnB = params[i].PnB;
		}
		return res;
	}

	void CytoFrame::write_h5_keys(H5File file) const
	{
		CompType key_type = get_h5_datatype_keys();
		hsize_t dim_key[] = {keys_.size()};
		hsize_t dim_max[] = {H5S_UNLIMITED};
		DataSpace dsp_key(1, dim_key, dim_max);
		DSetCreatPropList plist;
		plist.setChunk(1, dim_key);
		DataSet ds = file.createDataSet( "keywords", key_type, dsp_key, plist);

		auto keyVec = to_kw_vec<KEY_WORDS>(keys_);
		ds.write(&keyVec[0], key_type );

	}
	void CytoFrame::write_h5_pheno_data(H5File file) const
	{
		CompType key_type = get_h5_datatype_keys();
		hsize_t nSize = pheno_data_.size();
		if(nSize==0)
			throw runtime_error("CytoFrame requires non-empty pdata to write to h5!");
		hsize_t dim_pd[] = {nSize};
		hsize_t dim_max[] = {H5S_UNLIMITED};
		DataSpace dsp_pd(1, dim_pd, dim_max);
		DSetCreatPropList plist;
		plist.setChunk(1, dim_pd);

		DataSet ds = file.createDataSet( "pdata", key_type, dsp_pd, plist);

		auto keyVec = to_kw_vec<PDATA>(pheno_data_);
		ds.write(&keyVec[0], key_type );

	}
	/**
	 * save the CytoFrame as HDF5 format
	 *
	 * @param filename the path of the output H5 file
	 */
	void CytoFrame::write_h5(const string & filename) const
	{
		H5File file( filename, H5F_ACC_TRUNC );

		write_h5_params(file);

		write_h5_keys(file);

		write_h5_pheno_data(file);


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
		hsize_t dim_max[] = {H5S_UNLIMITED, H5S_UNLIMITED};

		DataSpace dataspace( 2, dimsf, dim_max);
		DataSet dataset = file.createDataSet( DATASET_NAME, h5_datatype_data(DataTypeLocation::H5), dataspace, plist);
		/*
		* Write the data to the dataset using default memory space, file
		* space, and transfer properties.
		*/
		EVENT_DATA_VEC dat = get_data();
		dataset.write(dat.mem, h5_datatype_data(DataTypeLocation::MEM));
	}

	/**
	 * extract the value of the single keyword by keyword name
	 *
	 * @param key keyword name
	 * @return keyword value as a string
	 */
	string CytoFrame::get_keyword(const string & key) const
	{
		string res="";
		auto it = keys_.find(key);
		if(it!=keys_.end())
			res = it->second;
		return res;
	}

	void CytoFrame::subset_parameters(uvec col_idx)
	{
			unsigned n = col_idx.size();
			vector<cytoParam> params_new(n);
			for(unsigned i = 0; i < n; i++)
			{
					params_new[i] = params[col_idx[i]];
			}
			set_params(params_new);

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
	vector<string> CytoFrame::get_channels() const
	{
		vector<string> res(n_cols());
		for(unsigned i = 0; i < n_cols(); i++)
			res[i] = params[i].channel;
		return res;
	}

	void CytoFrame::set_channels(const CHANNEL_MAP & chnl_map)
	{
		check_write_permission();
		for(auto & it : chnl_map)
		{
			try//catch the unmatched col error so that it can proceed the rest
			{
				set_channel(it.first, it.second);
			}catch (const domain_error & e){
				string msg = e.what();
				if(msg.find("colname not found") == string::npos)
					throw(e);
			}

		}

	}
	/**
	 * get all the marker names
	 * @return
	 */
	vector<string> CytoFrame::get_markers() const
	{
		vector<string> res(n_cols());
			for(unsigned i = 0; i < n_cols(); i++)
				res[i] = params[i].marker;
		return res;
	}
	/**
	 * Look up marker by channel info
	 * @param channel
	 * @return
	 */
	string CytoFrame::get_marker(const string & channel)
	{
		int idx = get_col_idx(channel, ColType::channel);
		if(idx<0)
			throw(domain_error("colname not found: " + channel));
		return params[idx].marker;
	}


	/**
	 * get the numeric index for the given column
	 * @param colname column name
	 * @param type the type of column
	 * @return
	 */
	int CytoFrame::get_col_idx(const string & colname, ColType type) const
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
	uvec CytoFrame::get_col_idx(vector<string> colnames, ColType col_type) const
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
	void CytoFrame::set_channel(const string & oldname, const string &newname, bool is_update_keywords)
	{
		check_write_permission();
		int id = get_col_idx(oldname, ColType::channel);
		if(id<0)
			throw(domain_error("colname not found: " + oldname));
		if(oldname!=newname)
		{

			if(g_loglevel>=GATING_HIERARCHY_LEVEL)
				PRINT(oldname + "-->"  + newname + "\n");
			int id1 = get_col_idx(newname, ColType::channel);
			if(id1>=0&&id1!=id)
				throw(domain_error("colname already exists: " + newname));
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
	}

	void CytoFrame::set_marker(const string & channelname, const string & markername)
	{
		check_write_permission();
		int id = get_col_idx(channelname, ColType::channel);
		if(id<0)
			throw(domain_error("channel not found: " + channelname));
                string oldmarkername = params[id].marker;
		if(oldmarkername!=markername)
		{
			int id1 = get_col_idx(markername, ColType::marker);
			if(id1>=0&&id1!=id)
				throw(domain_error("marker already exists: " + markername));
			params[id].marker=markername;
			marker_vs_idx.erase(oldmarkername);
			marker_vs_idx[markername] = id;
		}
	}

	/**
	 * Update the instrument range (typically after data transformation)
	 * @param colname
	 * @param ctype
	 * @param new_range
	 */
	void CytoFrame::set_range(const string & colname, ColType ctype, pair<EVENT_DATA_TYPE, EVENT_DATA_TYPE> new_range
			, bool is_update_keywords){
		check_write_permission();
		int idx = get_col_idx(colname, ctype);
		if(idx<0)
			throw(domain_error("colname not found: " + colname));
		params[idx].min = new_range.first;
		params[idx].max = new_range.second;

		if(is_update_keywords)
		{
			string pid = to_string(idx+1);
			set_keyword("flowCore_$P" + pid + "Rmin", boost::lexical_cast<string>(new_range.first));
			set_keyword("flowCore_$P" + pid + "Rmax", boost::lexical_cast<string>(new_range.second));
		}
	}
	/**
	 * the range of a specific column
	 * @param colname
	 * @param ctype the type of column
	 * @param rtype either RangeType::data or RangeType::instrument
	 * @return
	 */
	pair<EVENT_DATA_TYPE, EVENT_DATA_TYPE> CytoFrame::get_range(const string & colname, ColType ctype, RangeType rtype) const
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
	EVENT_DATA_TYPE CytoFrame::get_time_step(const string time_channel) const
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


	string CytoFrame::get_pheno_data(const string & name) const {
		auto it = pheno_data_.find(name);
		if(it==pheno_data_.end())
			return "";
		else
			return it->second;

	}

};


