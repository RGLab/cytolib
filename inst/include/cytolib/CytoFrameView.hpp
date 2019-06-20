/*
 * CytoFrameView.hpp
 *
 *  Created on: Apr 3, 2018
 *      Author: wjiang2
 */

#ifndef INST_INCLUDE_CYTOLIB_CYTOFRAMEVIEW_HPP_
#define INST_INCLUDE_CYTOLIB_CYTOFRAMEVIEW_HPP_
#include "CytoFrame.hpp"
namespace cytolib
{
class CytoFrameView{
	CytoFramePtr ptr_;
	arma::uvec row_idx_;
	arma::uvec col_idx_;
public:
	CytoFrameView(){};
	CytoFrameView(CytoFramePtr ptr):ptr_(ptr){};
	CytoFramePtr get_cytoframe_ptr() const{
		if(ptr_)
			return ptr_;
		else
			throw(domain_error("Empty CytoFrameView!"));
	}
	bool is_row_indexed() const{return row_idx_.size()>0;};
	bool is_col_indexed() const{return col_idx_.size()>0;};

	/*forwarded apis
	 *
	 */
	void scale_time_channel(string time_channel = "time"){
		get_cytoframe_ptr()->scale_time_channel(time_channel);
			}
	void set_readonly(bool flag){
		get_cytoframe_ptr()->set_readonly(flag);
	}
	virtual void flush_meta(){
		get_cytoframe_ptr()->flush_meta();
	};
	virtual void load_meta(){
		get_cytoframe_ptr()->load_meta();
	};
	string get_h5_file_path() const{
			return get_cytoframe_ptr()->get_h5_file_path();
		}
	vector<string> get_channels() const{
		vector<string> orig = get_cytoframe_ptr()->get_channels();
		unsigned n = col_idx_.size();
		if(n == 0)
			return orig;
		else
		{
			vector<string> res(n);
			for(unsigned i = 0; i < n; i++)
				res[i] = orig[col_idx_[i]];
			return res;
		}
	}
	vector<string> get_markers() const{
		vector<string> orig = get_cytoframe_ptr()->get_markers();
		unsigned n = col_idx_.size();
		if(n == 0)
			return orig;
		else
		{
			vector<string> res(n);
			for(unsigned i = 0; i < n; i++)
				res[i] = orig[col_idx_[i]];
			return res;
		}
	}
	void set_channels(const CHANNEL_MAP & chnl_map){get_cytoframe_ptr()->set_channels(chnl_map);}
	void convertToPb(pb::CytoFrame & fr_pb, const string & h5_filename, H5Option h5_opt) const{get_cytoframe_ptr()->convertToPb(fr_pb, h5_filename, h5_opt);};
	void set_channel(const string & oldname, const string &newname)
	{
		get_cytoframe_ptr()->set_channel(oldname, newname);
	}
	string get_marker(const string & channel)
	{
		return get_cytoframe_ptr()->get_marker(channel);
	}

	void set_marker(const string & channelname, const string & markername)
	{
		get_cytoframe_ptr()->set_marker(channelname, markername);
	}
	void compensate(const compensation & comp)
	{
			get_cytoframe_ptr()->compensate(comp);
	}
	compensation get_compensation(const string & key = "SPILL")
	{
		return	get_cytoframe_ptr()->get_compensation(key);
	}
	void write_h5(const string & filename) const
	{
		get_cytoframe_ptr()->write_h5(filename);
	}

	KEY_WORDS get_keywords() const{
		KEY_WORDS res = get_cytoframe_ptr()->get_keywords();
		//TODO: we currently do this filtering in R
//		if(is_col_indexed())//filter out the PnX based on the col_idx
//		{
//
//			    pat <- paste(to.del, collapse = "|")
//			    pat <- paste0("^\\$P(", pat , ")[A-Z]$")
//			    sel <- grep(pat, names(kw))
//		}
		return res;
	}
	/**
	 * extract the value of the single keyword by keyword name
	 *
	 * @param key keyword name
	 * @return keyword value as a string
	 */
	string get_keyword(const string & key) const
	{
		return get_cytoframe_ptr()->get_keyword(key);

	}

	/**
	 * set the value of the single keyword
	 * @param key keyword name
	 * @param value keyword value
	 */
	void set_keyword(const string & key, const string & value)
	{
		get_cytoframe_ptr()->set_keyword(key, value);
	}
	void set_keywords(const KEY_WORDS & keys){
		get_cytoframe_ptr()->set_keywords(keys);
	}

	void set_range(const string & colname, ColType ctype, pair<EVENT_DATA_TYPE, EVENT_DATA_TYPE> new_range){
		get_cytoframe_ptr()->set_range(colname, ctype, new_range);
	}
	/**
	 * the range of a specific column
	 * @param colname
	 * @param ctype the type of column
	 * @param rtype either RangeType::data or RangeType::instrument
	 * @return
	 */
	pair<EVENT_DATA_TYPE, EVENT_DATA_TYPE> get_range(const string & colname, ColType ctype, RangeType rtype) const
	{

		return get_cytoframe_ptr()->get_range(colname, ctype	, rtype);
	}

	const PDATA & get_pheno_data() const {return get_cytoframe_ptr()->get_pheno_data();}
	void set_pheno_data(const string & name, const string & value) {get_cytoframe_ptr()->set_pheno_data(name, value);}
	void set_pheno_data(const PDATA & _pd) {get_cytoframe_ptr()->set_pheno_data(_pd);}
	const vector<cytoParam> & get_params() const
	{
		return get_cytoframe_ptr()->get_params();
	}

	void set_params(const vector<cytoParam> & _params)
	{
		get_cytoframe_ptr()->set_params(_params);
	}

	/*subsetting*/

	void cols_(vector<string> colnames, ColType col_type)
	{

		uvec col_idx = get_cytoframe_ptr()->get_col_idx(colnames, col_type);//idx from original frame
		if(is_col_indexed())//convert abs idex to relative indx
		{
			for(unsigned i = 0; i < col_idx.size(); i++)
			{

				auto it = std::find(col_idx_.begin(), col_idx_.end(), col_idx[i]);
				if(it == col_idx_.end())
					throw(domain_error(colnames[i] + " not present in the current cytoframeView!"));

				col_idx[i] = it - col_idx_.begin();
			}

		}
		//update params
		CytoFrameView::cols_(col_idx);

	}

	/**
	 *
	 * @param col_idx column index relative to view
	 */
	void cols_(uvec col_idx)
	{
		unsigned max_idx = col_idx.max();
		unsigned min_idx = col_idx.min();
		if(max_idx >= n_cols() || min_idx < 0)
			throw(domain_error("The size of the new col index is not within the original mat size!"));
		if(is_col_indexed())//covert relative idx to abs idx
		{
//			cout << "indexing " << endl;
			for(auto & i : col_idx)
			{
//				cout << "relative idx: " << i << " abs: " << col_idx_[i] << endl;
				i = col_idx_[i];
			}


		}
		col_idx_ = col_idx;

	}
	void cols_(vector<unsigned> col_idx)
	{
		cols_(arma::conv_to<uvec>::from(col_idx));
	}

	void rows_(vector<unsigned> row_idx)
	{
		rows_(arma::conv_to<uvec>::from(row_idx));
	}

	void rows_(uvec row_idx)
	{
		unsigned max_idx = row_idx.max();
		unsigned min_idx = row_idx.min();
		if(max_idx >= n_rows() || min_idx < 0)
			throw(domain_error("The size of the new row index is not within the original mat size!"));
		if(is_row_indexed())
		{
			for(auto & i : row_idx)
				i = row_idx_[i];

		}
		row_idx_ = row_idx;

	}
	/**
	 * Corresponding to the original $Pn FCS TEXT
	 * @return
	 */
	vector<unsigned> get_original_col_ids() const
	{
		unsigned n = n_cols();
		vector<unsigned> res(n);
		for(unsigned i = 0; i < n; i++)
			if(is_col_indexed())
				res[i] = col_idx_[i];
			else
				res[i] = i;
		return res;

	}
	unsigned n_cols() const
	{
		if(is_col_indexed())
			return col_idx_.size();
		else
			return get_cytoframe_ptr()->n_cols();
	}
	/**
	 * get the number of rows(or events)
	 * @return
	 */
	unsigned n_rows() const
	{
		//read nEvents
		if(is_row_indexed())
			return row_idx_.size();
		else
			return get_cytoframe_ptr()->n_rows();
	}
	/**
	 * clear the row and column index
	 */
	void reset_view(){
		row_idx_.clear();
		col_idx_.clear();

	}
	/**
	 * Realize the delayed subsetting (i.e. cols() and rows()) operations to the underlying data
	 * and clear the view
	 */
	CytoFrameView copy_realized(const string & h5_filename = "") const
	{
		return get_cytoframe_ptr()->copy_realized(row_idx_, col_idx_, h5_filename);;
	}
	void set_data(const EVENT_DATA_VEC & data_in){
		//fetch the original view of data
		EVENT_DATA_VEC data_orig = get_cytoframe_ptr()->get_data();
		//update it
		if(is_col_indexed()&&is_row_indexed())
			data_orig.submat(row_idx_, col_idx_) = data_in;
		else if(is_row_indexed())
			data_orig.rows(row_idx_) = data_in;
		else if(is_col_indexed())
			data_orig.cols(col_idx_) = data_in;
		else
			if(data_orig.n_cols!=data_in.n_cols||data_orig.n_rows!=data_in.n_rows)
				throw(domain_error("The size of theinput data is different from the cytoframeview!"));
			else
				data_orig = data_in;


		//write back to ptr_
		get_cytoframe_ptr()->set_data(data_orig);
	}
	EVENT_DATA_VEC get_data() const
	{
		EVENT_DATA_VEC data;

		if(is_col_indexed())
			data = get_cytoframe_ptr()->get_data(col_idx_);
		else
			data = get_cytoframe_ptr()->get_data();

		if(is_row_indexed())
			data = data.rows(row_idx_);

		return data;
	}

	CytoFrameView copy(const string & h5_filename = "") const
	{
		CytoFrameView cv(*this);
		cv.ptr_ = get_cytoframe_ptr()->copy(h5_filename);
		return cv;

	}
};
}



#endif /* INST_INCLUDE_CYTOLIB_CYTOFRAMEVIEW_HPP_ */
