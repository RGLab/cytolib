/* Copyright 2019 Fred Hutchinson Cancer Research Center
 * See the included LICENSE file for details on the license that is granted to the
 * user of this software.
 * H5CytoFrame.hpp
 *
 *  Created on: Sep 18, 2017
 *      Author: wjiang2
 */

#ifndef INST_INCLUDE_CYTOLIB_H5CYTOFRAME_HPP_
#define INST_INCLUDE_CYTOLIB_H5CYTOFRAME_HPP_
#include <cytolib/MemCytoFrame.hpp>
#include <cytolib/global.hpp>
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

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
	hsize_t dims[2];              // dataset dimensions
	bool readonly_;//whether allow the public API to modify it, can't rely on h5 flag mechanism since
//					its behavior is uncerntain for multiple opennings
	//flags indicating if cached meta data needs to be flushed to h5
	bool is_dirty_params;
	bool is_dirty_keys;
	bool is_dirty_pdata;
	FileAccPropList access_plist_;//used to custom fapl, especially for s3 backend
	EVENT_DATA_VEC read_data(uvec col_idx) const;
	int h5_flags() const{
		if(get_readonly())
			return H5F_ACC_RDONLY;
		else
			return H5F_ACC_RDWR;
	};
public:
	void flush_meta();
	void flush_params();

	void flush_keys();
	void flush_pheno_data();
	void set_readonly(bool flag){
		readonly_ = flag;
	}
	bool get_readonly() const{
		return readonly_ ;
	}
	FileFormat get_backend_type() const{
			return FileFormat::H5;
		};
	H5CytoFrame(const H5CytoFrame & frm):CytoFrame(frm)
	{
		filename_ = frm.filename_;
		is_dirty_params = frm.is_dirty_params;
		is_dirty_keys = frm.is_dirty_keys;
		is_dirty_pdata = frm.is_dirty_pdata;
		readonly_ = frm.readonly_;
		access_plist_ = frm.access_plist_;
		memcpy(dims, frm.dims, sizeof(dims));

	}
	H5CytoFrame(H5CytoFrame && frm):CytoFrame(frm)
	{
//		swap(pheno_data_, frm.pheno_data_);
//		swap(keys_, frm.keys_);
//		swap(params, frm.params);
//		swap(channel_vs_idx, frm.channel_vs_idx);
//		swap(marker_vs_idx, frm.marker_vs_idx);
		swap(filename_, frm.filename_);
		swap(dims, frm.dims);
		swap(access_plist_, frm.access_plist_);

		swap(readonly_, frm.readonly_);
		swap(is_dirty_params, frm.is_dirty_params);
		swap(is_dirty_keys, frm.is_dirty_keys);
		swap(is_dirty_pdata, frm.is_dirty_pdata);
	}
	H5CytoFrame & operator=(const H5CytoFrame & frm)
	{
		CytoFrame::operator=(frm);
		filename_ = frm.filename_;
		is_dirty_params = frm.is_dirty_params;
		is_dirty_keys = frm.is_dirty_keys;
		is_dirty_pdata = frm.is_dirty_pdata;
		readonly_ = frm.readonly_;
		access_plist_ = frm.access_plist_;
		memcpy(dims, frm.dims, sizeof(dims));
		return *this;
	}
	H5CytoFrame & operator=(H5CytoFrame && frm)
	{
		CytoFrame::operator=(frm);
		swap(filename_, frm.filename_);
		swap(dims, frm.dims);
		swap(is_dirty_params, frm.is_dirty_params);
		swap(is_dirty_keys, frm.is_dirty_keys);
		swap(is_dirty_pdata, frm.is_dirty_pdata);
		swap(readonly_, frm.readonly_);
		swap(access_plist_, frm.access_plist_);
		return *this;
	}

	unsigned n_rows() const{
				return dims[1];
		}

	void set_params(const vector<cytoParam> & _params)
	{
		CytoFrame::set_params(_params);
		is_dirty_params = true;

	}
	void set_keywords(const KEY_WORDS & keys){
		CytoFrame::set_keywords(keys);
		is_dirty_keys = true;

	}
	void set_keyword(const string & key, const string & value){
		CytoFrame::set_keyword(key, value);
		is_dirty_keys = true;

	}
	void rename_keyword(const string & old_key, const string & new_key){
		CytoFrame::rename_keyword(old_key, new_key);
		is_dirty_keys = true;
	}
	void remove_keyword(const string & key){
		CytoFrame::remove_keyword(key);
		is_dirty_keys = true;
	}
	void set_channel(const string & oldname, const string &newname, bool is_update_keywords = true)
	{
		CytoFrame::set_channel(oldname, newname, is_update_keywords);
		is_dirty_params = true;
		if(is_update_keywords)
			is_dirty_keys = true;
	}
	int set_channels(const vector<string> & channels)
	{
		int res = CytoFrame::set_channels(channels);
		is_dirty_params = true;
		is_dirty_keys = true;
		return res;
	}

	void append_data_columns(const EVENT_DATA_VEC & new_cols)
	{
		EVENT_DATA_VEC this_data = get_data();
		this_data.insert_cols(this_data.n_cols, new_cols);
		set_data(this_data);
	}

	void set_marker(const string & channelname, const string & markername)
	{
		CytoFrame::set_marker(channelname, markername);
		is_dirty_params = true;
	}
	void set_range(const string & colname, ColType ctype, pair<EVENT_DATA_TYPE, EVENT_DATA_TYPE> new_range, bool is_update_keywords = true){
		CytoFrame::set_range(colname, ctype, new_range, is_update_keywords);
		is_dirty_params = true;
		if(is_update_keywords)
			is_dirty_keys = true;

	}
	void set_pheno_data(const string & name, const string & value){
		CytoFrame::set_pheno_data(name, value);
		is_dirty_pdata = true;
	}
	void set_pheno_data(const PDATA & _pd)
	{
		CytoFrame::set_pheno_data(_pd);
		is_dirty_pdata = true;
	}
	void del_pheno_data(const string & name){
		CytoFrame::del_pheno_data(name);
		is_dirty_pdata = true;
	}
	/**
	 * constructor from FCS
	 * @param fcs_filename
	 * @param h5_filename
	 */
	H5CytoFrame(const string & fcs_filename, FCS_READ_PARAM & config, const string & h5_filename
			, bool readonly = false):filename_(h5_filename), is_dirty_params(false), is_dirty_keys(false), is_dirty_pdata(false)
	{
		MemCytoFrame fr(fcs_filename, config);
		fr.read_fcs();
		fr.write_h5(h5_filename);
		*this = H5CytoFrame(h5_filename, readonly);
	}
	/**
	 * constructor from the H5
	 * @param _filename H5 file path
	 */
	H5CytoFrame(const string & h5_filename, bool readonly = true, bool init = true):CytoFrame(),filename_(h5_filename), readonly_(readonly), is_dirty_params(false), is_dirty_keys(false), is_dirty_pdata(false)
	{
		access_plist_ = FileAccPropList::DEFAULT;
		if(init)//optionally delay load for the s3 derived cytoframe which needs to reset fapl before load
			init_load();
	}
	void init_load(){
		//always use the same flag and keep lock at cf level to avoid h5 open error caused conflicting h5 flags among cf objects that points to the same h5
		H5File file(filename_, h5_flags(), FileCreatPropList::DEFAULT, access_plist_);
		load_meta();


		//open dataset for event data

		auto dataset = file.openDataSet(DATASET_NAME);
		auto dataspace = dataset.getSpace();
		dataspace.getSimpleExtentDims(dims);

	}
	/**
	 * abandon the changes to the meta data in cache by reloading them from disk
	 */
	void load_meta();;;

	string get_uri() const{
		return filename_;
	}
	void check_write_permission() const{
		if(readonly_)
			throw(domain_error("Can't write to the read-only H5CytoFrame object!"));

	}

	void convertToPb(pb::CytoFrame & fr_pb
			, const string & h5_filename
			, CytoFileOption h5_opt
			, const CytoCtx & ctx = CytoCtx()) const
	{
			fr_pb.set_is_h5(true);
			if(h5_opt != CytoFileOption::skip)
			{
				auto h5path = fs::path(h5_filename);
				auto dest = h5path.parent_path();
				if(!fs::exists(dest))
					throw(logic_error(dest.string() + "doesn't exist!"));

				if(!fs::equivalent(fs::path(filename_).parent_path(), dest))
				{
					switch(h5_opt)
					{
					case CytoFileOption::copy:
						{
							if(fs::exists(h5path))
								fs::remove(h5path);
							fs::copy(filename_, h5_filename);
							break;
						}
					case CytoFileOption::move:
						{
							if(fs::exists(h5path))
								fs::remove(h5path);
							fs::rename(filename_, h5_filename);
							break;
						}
					case CytoFileOption::link:
						{
							throw(logic_error("'link' option for H5CytoFrame is no longer supported!"));
							fs::create_hard_link(filename_, h5_filename);
							break;
						}
					case CytoFileOption::symlink:
						{
							if(fs::exists(h5path))
								fs::remove(h5path);
							fs::create_symlink(filename_, h5_filename);
							break;
						}
					default:
						throw(logic_error("invalid h5_opt!"));
					}
				}
			}
		}
	/**
	 * Read data from disk.
	 * The caller will directly receive the data vector without the copy overhead thanks to the move semantics supported for vector container in c++11
	 * @return
	 */


	EVENT_DATA_VEC get_data() const
	{
		unsigned n = n_cols();
		uvec col_idx(n);
		for(unsigned i = 0; i < n; i++)
			col_idx[i] = i;
		return read_data(col_idx);
	}
	/**
	 * Partial IO
	 * @param col_idx
	 * @return
	 */
	EVENT_DATA_VEC get_data(uvec idx, bool is_col) const
	{
		if(is_col)
			return read_data(idx);
		else
			return get_data().rows(idx);
	}
	EVENT_DATA_VEC get_data(uvec row_idx, uvec col_idx) const
	{
		return read_data(col_idx).rows(row_idx);
	}
	/*
	 * protect the h5 from being overwritten accidentally
	 * which will make the original cf object invalid
	 */
	void copy_overwrite_check(const string & dest, bool overwrite) const
	{
		// Check if trying to write to same file
		if(fs::equivalent(fs::path(filename_), fs::path(dest))){
			// First check if that is even allowed by this CytoFrame
			check_write_permission();
			// Then make sure it has been explicitly approved by the caller
			if(!overwrite)
				throw(domain_error("Copying H5CytoFrame to itself is not supported! "+ dest));
		}
	}

	CytoFramePtr copy(const string & h5_filename = "", bool overwrite = false) const
	{
		copy_overwrite_check(h5_filename, overwrite);
		string new_filename = h5_filename;
		if(new_filename == "")
		{
			new_filename = generate_unique_filename(fs::temp_directory_path().string(), "", ".h5");
			fs::remove(new_filename);
		}
		fs::copy_file(filename_, new_filename);
		CytoFramePtr ptr(new H5CytoFrame(new_filename, false));
		//copy cached meta
		ptr->set_params(get_params());
		ptr->set_keywords(get_keywords());
		ptr->set_pheno_data(get_pheno_data());
		return ptr;
	}
	CytoFramePtr copy(uvec row_idx, uvec col_idx, const string & h5_filename = "", bool overwrite = false) const
	{
		copy_overwrite_check(h5_filename, overwrite);

		string new_filename = h5_filename;
		if(new_filename == "")
		{
			new_filename = generate_unique_filename(fs::temp_directory_path().string(), "", ".h5");
			fs::remove(new_filename);
		}
		MemCytoFrame fr(*this);
		fr.copy(row_idx, col_idx)->write_h5(new_filename);//this flushes the meta data as well
		return CytoFramePtr(new H5CytoFrame(new_filename, false));
	}

	CytoFramePtr copy(uvec idx, bool is_row_indexed, const string & h5_filename = "", bool overwrite = false) const
	{
		copy_overwrite_check(h5_filename, overwrite);

		string new_filename = h5_filename;
		if(new_filename == "")
		{
			new_filename = generate_unique_filename(fs::temp_directory_path().string(), "", ".h5");
			fs::remove(new_filename);
		}
		MemCytoFrame fr(*this);
		fr.copy(idx, is_row_indexed)->write_h5(new_filename);//this flushes the meta data as well
		return CytoFramePtr(new H5CytoFrame(new_filename, false));
	}

	/**
	 * copy setter
	 * @param _data
	 */
	void set_data(const EVENT_DATA_VEC & _data);

	void set_data(EVENT_DATA_VEC && _data)
	{
		check_write_permission();
		set_data(_data);
	}
};

};



#endif /* INST_INCLUDE_CYTOLIB_H5CYTOFRAME_HPP_ */
