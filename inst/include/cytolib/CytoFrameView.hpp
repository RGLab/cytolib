/* Copyright 2019 Fred Hutchinson Cancer Research Center
 * See the included LICENSE file for details
 * on the licence that is granted to the user of this software.
 * CytoFrameView.hpp
 *
 *  Created on: Apr 3, 2018
 *      Author: wjiang2
 */

#ifndef INST_INCLUDE_CYTOLIB_CYTOFRAMEVIEW_HPP_
#define INST_INCLUDE_CYTOLIB_CYTOFRAMEVIEW_HPP_
#include "MemCytoFrame.hpp"
#include "H5CytoFrame.hpp"
#include "TileCytoFrame.hpp"
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

namespace cytolib
{
class CytoFrameView{
	CytoFramePtr ptr_;
	arma::uvec row_idx_;
	arma::uvec col_idx_;
	bool is_row_indexed_ = false;
	bool is_col_indexed_ = false;
public:
	CytoFrameView(){};
	CytoFrameView(CytoFramePtr ptr):ptr_(ptr){};
	CytoFramePtr get_cytoframe_ptr() const;

	bool is_row_indexed() const{return is_row_indexed_;};
	bool is_col_indexed() const{return is_col_indexed_;};
	bool is_empty() const{return (is_row_indexed_ && row_idx_.is_empty()) || (is_col_indexed_ && col_idx_.is_empty());};

	/*forwarded apis
	 *
	 */
	void scale_time_channel(string time_channel = "time"){
		get_cytoframe_ptr()->scale_time_channel(time_channel);
			}
	void set_readonly(bool flag){
		get_cytoframe_ptr()->set_readonly(flag);
	}
	bool get_readonly(){
			return get_cytoframe_ptr()->get_readonly();
		}
	virtual void flush_meta(){
		get_cytoframe_ptr()->flush_meta();
	};
	virtual void load_meta(){
		get_cytoframe_ptr()->load_meta();
	};
	string get_uri() const{
			return get_cytoframe_ptr()->get_uri();
		}
	FileFormat get_backend_type()  const{
		return get_cytoframe_ptr()->get_backend_type();
	};
	vector<string> get_channels() const;
	vector<string> get_markers() const;
	void set_channels(const CHANNEL_MAP & chnl_map){get_cytoframe_ptr()->set_channels(chnl_map);}
	void convertToPb(pb::CytoFrame & fr_pb
			, const string & cf_filename
			, CytoFileOption h5_opt
			, const CytoCtx & ctx = CytoCtx()) const{
		if(is_row_indexed_ || is_col_indexed_)
		{
			if(h5_opt == CytoFileOption::copy||h5_opt == CytoFileOption::move)
			{
				//realize view
				auto cfv = copy_realized(cf_filename, true);
				//trigger archive logic on the new cfv (which will skip overwriting itself)
				cfv.convertToPb(fr_pb, cf_filename, h5_opt, ctx);
				auto oldh5 = get_uri();
				if(h5_opt == CytoFileOption::move&&oldh5!="")
				{
					if(!fs::equivalent(fs::path(oldh5), fs::path(cf_filename)))
						fs::remove_all(oldh5);

				}
			}
			else
				throw(domain_error("Only 'copy' or 'move' option is supported for the indexed CytoFrameView object!"));
		}
		else
			get_cytoframe_ptr()->convertToPb(fr_pb, cf_filename, h5_opt, ctx);

	};
	void set_channel(const string & oldname, const string &newname)
	{
		get_cytoframe_ptr()->set_channel(oldname, newname);
	}
	int set_channels(const vector<string> & channels)
	{
		auto n1 = n_cols();
		auto n2 = channels.size();
		if(n2!=n1)
			throw(domain_error("The size of the input of 'set_channels' (" + to_string(n2) + ") is different from the original one (" + to_string(n1) + ")"));
		//update the original vector of channels
		auto idx = get_original_col_ids();
		auto old = get_cytoframe_ptr()->get_channels();
		for(unsigned i = 0; i < n1; i++)
		{
			old[idx[i]] = channels[i];
		}
		//update the cf
		return get_cytoframe_ptr()->set_channels(old);
	}

	// Realize the view to the underlying cytoframe (if necessary) and then append extra columns.
	// Realizing the view is necessary for unambiguous Pn indices for the added columns.
	void append_columns(const vector<string> & new_colnames, const EVENT_DATA_VEC & new_cols){
		if(is_row_indexed() || is_col_indexed()){
			// Realize to the original file and reset the view
			ptr_ = realize(ptr_, get_uri(), true);
			reset_view();
		}
		ptr_->append_columns(new_colnames, new_cols);
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
	compensation get_compensation(const string & key = "$SPILLOVER")
	{
		return	get_cytoframe_ptr()->get_compensation(key);
	}
	void write_to_disk(const string & filename, FileFormat format = FileFormat::TILE
			, const CytoCtx ctx = CytoCtx()) const
	{
		//create a mem-based cfv to avoid extra disk write IO from realization call
		CytoFrameView cv(*this);
		cv.ptr_ = CytoFramePtr(new MemCytoFrame(*(get_cytoframe_ptr())));
		//TODO:it would less overhead if we could have in-place realize method without creating the copy
		auto cv1 = cv.copy_realized();
		auto ptr = cv1.get_cytoframe_ptr();

		ptr->write_to_disk(filename, format, ctx);

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

	void rename_keyword(const string & old_key, const string & new_key){
		get_cytoframe_ptr()->rename_keyword(old_key, new_key);
	}
	void remove_keyword(const string & key){
		get_cytoframe_ptr()->remove_keyword(key);
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

	void cols_(vector<string> colnames, ColType col_type);

	/**
	 *
	 * @param col_idx column index relative to view
	 */
	void cols_(uvec col_idx);
	void cols_(vector<unsigned> col_idx);

	void rows_(vector<unsigned> row_idx);

	void rows_(uvec row_idx);
	/**
	 * Corresponding to the original $Pn FCS TEXT
	 * @return
	 */
	vector<unsigned> get_original_col_ids() const;
	unsigned n_cols() const;
	/**
	 * get the number of rows(or events)
	 * @return
	 */
	unsigned n_rows() const;
	/**
	 * clear the row and column index
	 */
	void reset_view(){
		row_idx_.clear();
		col_idx_.clear();
		is_row_indexed_ = is_col_indexed_ = false;
	}
	/**
	 * Realize the delayed subsetting (i.e. cols() and rows()) operations to the underlying data
	 * and clear the view
	 */
	CytoFramePtr realize(CytoFramePtr ptr, const string & cf_filename = "", bool overwrite = false) const
	{
		if(is_row_indexed_ && is_col_indexed_){
			return ptr->copy(row_idx_, col_idx_, cf_filename, overwrite);
		}else if(is_row_indexed_){
			return ptr->copy(row_idx_, true, cf_filename, overwrite);
		}else if(is_col_indexed_){
			return ptr->copy(col_idx_, false, cf_filename, overwrite);
		}else{
			return ptr->copy(cf_filename, overwrite);
		}
	}
	CytoFrameView copy_realized(const string & cf_filename = "", bool overwrite = false) const
	{
		return CytoFrameView(realize(get_cytoframe_ptr(), cf_filename, overwrite));
	}
	/*
	 * this API exists to bypass potential writing the on-disk cf
	 */
	shared_ptr<MemCytoFrame> get_realized_memcytoframe() const{

		shared_ptr<MemCytoFrame> ptr(new MemCytoFrame(*get_cytoframe_ptr()));
		auto res = realize(ptr);
		return dynamic_pointer_cast<MemCytoFrame>(res);
	}
	void set_data(const EVENT_DATA_VEC & data_in);
	EVENT_DATA_VEC get_data() const;

	CytoFrameView copy(const string & cf_filename = "") const;
};



}



#endif /* INST_INCLUDE_CYTOLIB_CYTOFRAMEVIEW_HPP_ */
