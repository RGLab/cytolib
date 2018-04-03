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


	/*subsetting*/

	void cols_(vector<string> colnames, ColType col_type)
	{

		uvec col_idx = ptr_->get_col_idx(colnames, col_type);

		//update params
		CytoFrameView::cols_(col_idx);

	}

	void cols_(uvec col_idx)
	{
		if(col_idx.max() >= n_cols() || col_idx.min() < 0)
			throw(domain_error("The size of the new row index is not within the original mat size!"));
		if(is_col_indexed())
		{
			for(auto & i : col_idx)
				i = col_idx_[i];

		}
		col_idx_ = col_idx;

	}

	void rows_(vector<unsigned> row_idx)
	{
		rows_(arma::conv_to<uvec>::from(row_idx));
	}

	void rows_(uvec row_idx)
	{
		if(row_idx.max() >= n_rows() || row_idx.min() < 0)
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
		return ptr_->copy(h5_filename);

	}
};
}



#endif /* INST_INCLUDE_CYTOLIB_CYTOFRAMEVIEW_HPP_ */
