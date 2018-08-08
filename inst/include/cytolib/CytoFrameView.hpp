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
	CytoFramePtr get_cytoframe_ptr(){return ptr_;}
	bool is_row_indexed() const{return row_idx_.size()>0;};
	bool is_col_indexed() const{return col_idx_.size()>0;};

	/*forwarded apis
	 *
	 */
	vector<string> get_channels() const{
		vector<string> orig = ptr_->get_channels();
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
		vector<string> orig = ptr_->get_markers();
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
	void update_channels(const CHANNEL_MAP & chnl_map){ptr_->update_channels(chnl_map);}
	void convertToPb(pb::CytoFrame & fr_pb, const string & h5_filename, H5Option h5_opt) const{ptr_->convertToPb(fr_pb, h5_filename, h5_opt);};
	void set_channel(const string & oldname, const string &newname)
	{
		ptr_->set_channel(oldname, newname);
	}

	void set_marker(const string & oldname, const string & newname)
	{
		ptr_->set_marker(oldname, newname);
	}
	void compensate(const compensation & comp)
	{
			ptr_->compensate(comp);
	}
	compensation get_compensation(const string & key = "SPILL")
	{
		return	ptr_->get_compensation(key);
	}
	void write_h5(const string & filename) const
	{
		ptr_->write_h5(filename);
	}
	const KEY_WORDS & get_keywords() const{
		return ptr_->get_keywords();
	}
	/**
	 * extract the value of the single keyword by keyword name
	 *
	 * @param key keyword name
	 * @return keyword value as a string
	 */
	string get_keyword(const string & key) const
	{
		return ptr_->get_keyword(key);

	}

	/**
	 * set the value of the single keyword
	 * @param key keyword name
	 * @param value keyword value
	 */
	void set_keyword(const string & key, const string & value)
	{
		ptr_->set_keyword(key, value);
	}
	void set_range(const string & colname, ColType ctype, pair<EVENT_DATA_TYPE, EVENT_DATA_TYPE> new_range){
		ptr_->set_range(colname, ctype, new_range);
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

		return ptr_->get_range(colname, ctype	, rtype);
	}

	const PDATA & get_pheno_data() const {return ptr_->get_pheno_data();}
	void set_pheno_data(const string & name, const string & value) {ptr_->set_pheno_data(name, value);}
	/*subsetting*/

	void cols_(vector<string> colnames, ColType col_type)
	{

		uvec col_idx = ptr_->get_col_idx(colnames, col_type);

		//update params
		CytoFrameView::cols_(col_idx);

	}

	void cols_(uvec col_idx)
	{
		unsigned max_idx = col_idx.max();
		unsigned min_idx = col_idx.min();
		if(max_idx >= n_cols() || min_idx < 0)
			throw(domain_error("The size of the new row index is not within the original mat size!"));
		if(is_col_indexed())
		{
			for(auto & i : col_idx)
				i = col_idx_[i];

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

	unsigned n_cols() const
	{
		if(is_col_indexed())
			return col_idx_.size();
		else
			return ptr_->n_cols();
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
			return ptr_->n_rows();
	}
	/**
	 * clear the row and column index
	 */
	void reset_view(){
		row_idx_.clear();
		col_idx_.clear();

	}
	/**
	 * Realize the delayed subsetting (i.e. cols() and rows()) operations
	 * and clear the view
	 */
	void flush_view()
	{

		ptr_->set_data(get_data());

		reset_view();
	}

	EVENT_DATA_VEC get_data() const
	{
		EVENT_DATA_VEC data;

		if(is_col_indexed())
			data = ptr_->get_data(col_idx_);
		else
			data = ptr_->get_data();

		if(is_row_indexed())
			data = data.rows(row_idx_);

		return data;
	}

	CytoFramePtr copy(const string & h5_filename = "") const
	{
		CytoFramePtr ptr = ptr_->copy(row_idx_, col_idx_, h5_filename);
//		//TODO:optimize into one step for h5
//		if(is_row_indexed())
//			ptr->rows_(row_idx_);
//		if(is_col_indexed())
//			ptr->cols_(col_idx_);
		return ptr;

	}
};
}



#endif /* INST_INCLUDE_CYTOLIB_CYTOFRAMEVIEW_HPP_ */
