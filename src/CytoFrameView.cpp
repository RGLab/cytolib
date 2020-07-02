// Copyright 2019 Fred Hutchinson Cancer Research Center
// See the included LICENSE file for details on the licence that is granted to the user of this software.
#include <cytolib/CytoFrameView.hpp>
#include <cytolib/global.hpp>
namespace cytolib
{
	CytoFramePtr CytoFrameView::get_cytoframe_ptr() const{
		if(ptr_)
			return ptr_;
		else
			throw(domain_error("Empty CytoFrameView!"));
	}
	vector<string> CytoFrameView::get_channels() const{
		vector<string> orig = get_cytoframe_ptr()->get_channels();
		unsigned n = col_idx_.size();
		if(!is_col_indexed_)
			return orig;
		else if(n == 0)
			return vector<string>();
		else
		{
			vector<string> res(n);
			for(unsigned i = 0; i < n; i++)
				res[i] = orig[col_idx_[i]];
			return res;
		}
	}
	vector<string> CytoFrameView::get_markers() const{
		vector<string> orig = get_cytoframe_ptr()->get_markers();
		unsigned n = col_idx_.size();
		if(!is_col_indexed_)
			return orig;
		else if(n == 0)
			return vector<string>();
		else
		{
			vector<string> res(n);
			for(unsigned i = 0; i < n; i++)
				res[i] = orig[col_idx_[i]];
			return res;
		}
	}
	/*subsetting*/

	void CytoFrameView::cols_(vector<string> colnames, ColType col_type)
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
	void CytoFrameView::cols_(uvec col_idx)
	{
		if(col_idx.is_empty()){
			col_idx_.reset();
		}else{
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
		is_col_indexed_ = true;
	}
	void CytoFrameView::cols_(vector<unsigned> col_idx)
	{
		cols_(arma::conv_to<uvec>::from(col_idx));
	}

	void CytoFrameView::rows_(vector<unsigned> row_idx)
	{
		rows_(arma::conv_to<uvec>::from(row_idx));
	}

	void CytoFrameView::rows_(uvec row_idx)
	{
		if(row_idx.is_empty()){
			row_idx_.reset();
		}else{
			unsigned max_idx = row_idx.max();
			unsigned min_idx = row_idx.min();
			if(max_idx >= n_rows() || min_idx < 0)
				throw(domain_error("The size of the new row index ("
						+ to_string(min_idx) + "," + to_string(min_idx)
						+ ") is not within the original mat size (0, " + to_string(n_rows()) + ")"
						)
					);
			if(is_row_indexed())
			{
				for(auto & i : row_idx)
					i = row_idx_[i];
				
			}
			row_idx_ = row_idx;
		}
		is_row_indexed_ = true;

	}
	/**
	 * Corresponding to the original $Pn FCS TEXT
	 * @return
	 */
	vector<unsigned> CytoFrameView::get_original_col_ids() const
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
	unsigned CytoFrameView::n_cols() const
	{
		if(is_col_indexed_)
			return col_idx_.size();
		else
			return get_cytoframe_ptr()->n_cols();
	}
	/**
	 * get the number of rows(or events)
	 * @return
	 */
	unsigned CytoFrameView::n_rows() const
	{
		//read nEvents
		if(is_row_indexed_)
			return row_idx_.size();
		else
			return get_cytoframe_ptr()->n_rows();
	}

	void CytoFrameView::set_data(const EVENT_DATA_VEC & data_in){
		if(is_empty()){
			// Setting empty to empty is an allowed no-op, but not setting empty to non-empty
			if(!data_in.is_empty()){
				throw(domain_error("Cannot assign non-empty input data to empty CytoFrameView!"));	
			}
		}else{
			//fetch the original view of data
			EVENT_DATA_VEC data_orig = get_cytoframe_ptr()->get_data();
			//update it
			if(is_col_indexed_&&is_row_indexed_)
				data_orig.submat(row_idx_, col_idx_) = data_in;
			else if(is_row_indexed_)
				data_orig.rows(row_idx_) = data_in;
			else if(is_col_indexed_)
				data_orig.cols(col_idx_) = data_in;
			else
				if(data_orig.n_cols!=data_in.n_cols||data_orig.n_rows!=data_in.n_rows)
					throw(domain_error("The size of the input data is different from the cytoframeview!"));
				else
					data_orig = data_in;
				
				
			//write back to ptr_
			get_cytoframe_ptr()->set_data(data_orig);
		}
	}
	EVENT_DATA_VEC CytoFrameView::get_data() const
	{
		EVENT_DATA_VEC data;
		if(is_empty()){
			data = EVENT_DATA_VEC(n_rows(), n_cols());
		}else{
			auto ptr = get_cytoframe_ptr();
			if(is_col_indexed()&&is_row_indexed())
				data = ptr->get_data(row_idx_, col_idx_);
			else if(is_col_indexed())
			{
				data = ptr->get_data(col_idx_, true);
			}else if(is_row_indexed())
				data = ptr->get_data(row_idx_, false);
			else
				data = ptr->get_data();
		}
		return data;
	}

	CytoFrameView CytoFrameView::copy(const string & h5_filename) const
	{
		CytoFrameView cv(*this);
		cv.ptr_ = get_cytoframe_ptr()->copy(h5_filename);
		return cv;

	}
};

