/*
 * TileCytoFrame.hpp
 *
 *  Created on: Apr 21, 2020
 *      Author: wjiang2
 */

#ifndef INST_INCLUDE_CYTOLIB_TILECYTOFRAME_HPP_
#define INST_INCLUDE_CYTOLIB_TILECYTOFRAME_HPP_
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
class TileCytoFrame:public CytoFrame{
protected:
	string uri_;
	hsize_t dims[2];              // dataset dimensions
	bool readonly_;//whether allow the public API to modify it, can't rely on h5 flag mechanism since its behavior is uncerntain for multiple opennings
	//flags indicating if cached meta data needs to be flushed to h5
	bool is_dirty_params;
	bool is_dirty_keys;
	bool is_dirty_pdata;
	FileAccPropList access_plist_;//used to custom fapl, especially for s3 backend
//	EVENT_DATA_VEC read_data(uvec col_idx) const{return EVENT_DATA_VEC();};
	tiledb::Context ctx_;

public:
	unsigned int default_flags = H5F_ACC_RDWR;
	void flush_meta(){};
	void flush_params(){};

	void flush_keys(){};
	void flush_pheno_data(){};
	void set_readonly(bool flag){
		readonly_ = flag;
	}
	bool get_readonly(){
		return readonly_ ;
	}
	TileCytoFrame(const TileCytoFrame & frm):CytoFrame(frm)
	{
		uri_ = frm.uri_;
		is_dirty_params = frm.is_dirty_params;
		is_dirty_keys = frm.is_dirty_keys;
		is_dirty_pdata = frm.is_dirty_pdata;
		readonly_ = frm.readonly_;
		access_plist_ = frm.access_plist_;
		memcpy(dims, frm.dims, sizeof(dims));

	}
	TileCytoFrame(TileCytoFrame && frm):CytoFrame(frm)
	{
//		swap(pheno_data_, frm.pheno_data_);
//		swap(keys_, frm.keys_);
//		swap(params, frm.params);
//		swap(channel_vs_idx, frm.channel_vs_idx);
//		swap(marker_vs_idx, frm.marker_vs_idx);
		swap(uri_, frm.uri_);
		swap(dims, frm.dims);
		swap(access_plist_, frm.access_plist_);

		swap(readonly_, frm.readonly_);
		swap(is_dirty_params, frm.is_dirty_params);
		swap(is_dirty_keys, frm.is_dirty_keys);
		swap(is_dirty_pdata, frm.is_dirty_pdata);
	}
	TileCytoFrame & operator=(const TileCytoFrame & frm)
	{
		CytoFrame::operator=(frm);
		uri_ = frm.uri_;
		is_dirty_params = frm.is_dirty_params;
		is_dirty_keys = frm.is_dirty_keys;
		is_dirty_pdata = frm.is_dirty_pdata;
		readonly_ = frm.readonly_;
		access_plist_ = frm.access_plist_;
		memcpy(dims, frm.dims, sizeof(dims));
		return *this;
	}
	TileCytoFrame & operator=(TileCytoFrame && frm)
	{
		CytoFrame::operator=(frm);
		swap(uri_, frm.uri_);
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
	TileCytoFrame(const string & fcs_filename, FCS_READ_PARAM & config, const string & uri
			, bool readonly = false, const S3Cred & cred = S3Cred()):uri_(uri), is_dirty_params(false), is_dirty_keys(false), is_dirty_pdata(false)
	{
		MemCytoFrame fr(fcs_filename, config);
		fr.read_fcs();
		tiledb::Config cfg;
		cfg["vfs.s3.aws_access_key_id"] = cred.access_key_id_;
		cfg["vfs.s3.aws_secret_access_key"] = cred.access_key_;
		cfg["vfs.s3.region"] = cred.region_;
		ctx_ = tiledb::Context(cfg);
		fr.write_tile(uri, ctx_);
		*this = TileCytoFrame(uri, readonly, true, cred);
	}
	/**
	 * constructor from the H5
	 * @param _filename H5 file path
	 */
	TileCytoFrame(const string & uri, bool readonly = true, bool init = true
			, const S3Cred & cred = S3Cred()):CytoFrame(),uri_(uri), readonly_(readonly), is_dirty_params(false), is_dirty_keys(false), is_dirty_pdata(false)
	{
		access_plist_ = FileAccPropList::DEFAULT;
		tiledb::Config cfg;
		cfg["vfs.s3.aws_access_key_id"] = cred.access_key_id_;
		cfg["vfs.s3.aws_secret_access_key"] = cred.access_key_;
		cfg["vfs.s3.region"] = cred.region_;
		ctx_ = tiledb::Context(cfg);
		if(init)//optionally delay load for the s3 derived cytoframe which needs to reset fapl before load
			init_load();
	}
	void init_load(){
		load_meta();

		fs::path arraypath(uri_);
		auto mat_uri = (arraypath / "mat").string();

		//open dataset for event data
		tiledb::Array array(ctx_, mat_uri, TILEDB_READ);
		tiledb::ArraySchema schema(ctx_, mat_uri);
		auto dm = schema.domain();
		auto nch = dm.dimension("channel").domain<int>().second;
		auto nevent = dm.dimension("cell").domain<int>().second;
		dims[0] = nevent;
		dims[1] = nch;
//		const std::vector<int> subarray = {1, nch};
//
//		tiledb::Query query(ctx, array);
//		query.set_subarray(subarray);
//		query.set_layout(TILEDB_GLOBAL_ORDER);
//		vector<float> min_vec(nch);
//		query.set_buffer("min", min_vec);
//		vector<float> max_vec(nch);
//		query.set_buffer("max", max_vec);
//		query.submit();
//		query.finalize();

		uint64_t npd = array.metadata_num();
		uint32_t v_num;
		tiledb_datatype_t v_type;
		for (uint64_t i = 0; i < npd; ++i) {
			string pn,pv;

			const void* v;
			array.get_metadata_from_index(i, &pn, &v_type, &v_num, &v);
			if(v)
			{
				pv = string(static_cast<const char *>(v));
				pv.resize(v_num);
			}
			pheno_data_[pn] = pv;
		}

//		auto dataset = file.openDataSet(DATASET_NAME);
//		auto dataspace = dataset.getSpace();
//		dataspace.getSimpleExtentDims(dims);

	}
	/**
	 * abandon the changes to the meta data in cache by reloading them from disk
	 */
	void load_meta()
	{
		fs::path arraypath(uri_);

		auto params_uri = (arraypath / "params").string();
		load_params(params_uri);

		auto kw_uri = (arraypath / "keywords").string();
		load_kw(kw_uri);

	}
	void load_kw(const string & uri)
	{
		tiledb::Array array(ctx_, uri, TILEDB_READ);
		uint64_t nkw = array.metadata_num();

		keys_.resize(nkw);

		uint32_t v_num;
		tiledb_datatype_t v_type;
		for (uint64_t i = 0; i < nkw; ++i) {
			const void* v;
			array.get_metadata_from_index(i, &keys_[i].first, &v_type, &v_num, &v);
			if(v)
			{
				keys_[i].second = string(static_cast<const char *>(v));
				keys_[i].second.resize(v_num);
			}

		}

	}
	void load_params(const string & uri)
	{
		tiledb::Array array(ctx_, uri, TILEDB_READ);
//		tiledb::ArraySchema schema(ctx, uri);
//		auto dm = schema.domain();
//		auto nch = dm.dimension("params").domain<int>().second;
		int nch = array.metadata_num();
		const std::vector<int> subarray = {1, nch};

		tiledb::Query query(ctx_, array);
		query.set_subarray(subarray);
		query.set_layout(TILEDB_GLOBAL_ORDER);
		vector<float> min_vec(nch);
		query.set_buffer("min", min_vec);
		vector<float> max_vec(nch);
		query.set_buffer("max", max_vec);
		query.submit();
		query.finalize();

		params.resize(nch);

		uint32_t v_num;
		tiledb_datatype_t v_type;
		for (int i = 0; i < nch; ++i) {
			const void* v;
			array.get_metadata_from_index(i, &params[i].channel, &v_type, &v_num, &v);
			if(v)
			{
				params[i].marker = string(static_cast<const char *>(v));
				params[i].marker.resize(v_num);
			}

			params[i].min = min_vec[i];
			params[i].max = max_vec[i];

		}


	}
	string get_uri() const{
		return uri_;
	}
	void check_write_permission(){
		if(readonly_)
			throw(domain_error("Can't write to the read-only TileCytoFrame object!"));

	}

	void convertToPb(pb::CytoFrame & fr_pb, const string & h5_filename, H5Option h5_opt) const
	{
//			fr_pb.set_is_h5(true);
//			if(h5_opt != H5Option::skip)
//			{
//				auto h5path = fs::path(h5_filename);
//				auto dest = h5path.parent_path();
//				if(!fs::exists(dest))
//					throw(logic_error(dest.string() + "doesn't exist!"));
//
//				if(!fs::equivalent(fs::path(filename_).parent_path(), dest))
//				{
//					switch(h5_opt)
//					{
//					case H5Option::copy:
//						{
//							if(fs::exists(h5path))
//								fs::remove(h5path);
//							fs::copy(filename_, h5_filename);
//							break;
//						}
//					case H5Option::move:
//						{
//							if(fs::exists(h5path))
//								fs::remove(h5path);
//							fs::rename(filename_, h5_filename);
//							break;
//						}
//					case H5Option::link:
//						{
//							throw(logic_error("'link' option for TileCytoFrame is no longer supported!"));
//							fs::create_hard_link(filename_, h5_filename);
//							break;
//						}
//					case H5Option::symlink:
//						{
//							if(fs::exists(h5path))
//								fs::remove(h5path);
//							fs::create_symlink(filename_, h5_filename);
//							break;
//						}
//					default:
//						throw(logic_error("invalid h5_opt!"));
//					}
//				}
//			}
		}
	/**
	 * Read data from disk.
	 * The caller will directly receive the data vector without the copy overhead thanks to the move semantics supported for vector container in c++11
	 * @return
	 */


	EVENT_DATA_VEC get_data() const
	{
		fs::path arraypath(uri_);
		auto mat_uri = (arraypath / "mat").string();

		//open dataset for event data
		tiledb::Array array(ctx_, mat_uri, TILEDB_READ);

		int ncol = dims[1];
		int nrow = dims[0];

		const std::vector<int> subarray = {1, nrow, 1, ncol};

		tiledb::Query query(ctx_, array);
		query.set_subarray(subarray);
		query.set_layout(TILEDB_GLOBAL_ORDER);
		EVENT_DATA_VEC data(nrow, ncol);

		query.set_buffer<double>("mat", data.memptr(), nrow * ncol);
		query.submit();
		query.finalize();


		return data;
	}
	/**
	 * Partial IO
	 * @param col_idx
	 * @return
	 */
	EVENT_DATA_VEC get_data(uvec idx, bool is_col) const
	{
		fs::path arraypath(uri_);
		auto mat_uri = (arraypath / "mat").string();

		//open dataset for event data
		tiledb::Array array(ctx_, mat_uri, TILEDB_READ);
		tiledb::Query query(ctx_, array);
		query.set_layout(TILEDB_COL_MAJOR);
		int ncol,nrow, dim_idx;
		if(is_col)
		{
			ncol = idx.size();
			nrow = dims[0];
			dim_idx = 1;
			query.add_range<int>(0, 1, nrow);//select all rows
		}
		else
		{
			nrow = idx.size();
			ncol = dims[1];
			dim_idx  = 0;
			query.add_range<int>(1, 1, ncol);//select all cols
		}

		//tiledb idx starting from 1
		for(int i : idx)
			query.add_range<int>(dim_idx, i+1, i+1);


		EVENT_DATA_VEC data(nrow, ncol);

		query.set_buffer<double>("mat", data.memptr(), nrow * ncol);
		query.submit();
		query.finalize();


		return data;

	}
	EVENT_DATA_VEC get_data(uvec row_idx, uvec col_idx) const
	{
		fs::path arraypath(uri_);
		auto mat_uri = (arraypath / "mat").string();

		//open dataset for event data
		tiledb::Array array(ctx_, mat_uri, TILEDB_READ);
		auto ncol = col_idx.size();
		auto nrow = row_idx.size();


		tiledb::Query query(ctx_, array);
		query.set_layout(TILEDB_COL_MAJOR);
		//tiledb idx starting from 1
		for(int i : row_idx)
			query.add_range<int>(0, i+1, i+1);
		for(int i : col_idx)
			query.add_range<int>(1, i+1, i+1);

		EVENT_DATA_VEC data(nrow, ncol);

		query.set_buffer<double>("mat", data.memptr(), nrow * ncol);
		query.submit();
		query.finalize();


		return data;
	}
	/*
	 * protect the h5 from being overwritten accidentally
	 * which will make the original cf object invalid
	 */
	void copy_overwrite_check(const string & dest) const
	{
		if(fs::equivalent(fs::path(uri_), fs::path(dest)))
			throw(domain_error("Copying TileCytoFrame to itself is not supported! "+ dest));
	}

	CytoFramePtr copy(const string & h5_filename = "", bool overwrite = false) const
	{
		if(!overwrite)
			copy_overwrite_check(h5_filename);
		string new_filename = h5_filename;
		if(new_filename == "")
		{
			new_filename = generate_unique_filename(fs::temp_directory_path().string(), "", ".h5");
			fs::remove(new_filename);
		}
		fs::copy_file(uri_, new_filename);
		CytoFramePtr ptr(new TileCytoFrame(new_filename, false));
		//copy cached meta
		ptr->set_params(get_params());
		ptr->set_keywords(get_keywords());
		ptr->set_pheno_data(get_pheno_data());
		return ptr;
	}
	CytoFramePtr copy(uvec row_idx, uvec col_idx, const string & h5_filename = "", bool overwrite = false) const
	{

		if(!overwrite)
			copy_overwrite_check(h5_filename);
		string new_filename = h5_filename;
		if(new_filename == "")
		{
			new_filename = generate_unique_filename(fs::temp_directory_path().string(), "", ".h5");
			fs::remove(new_filename);
		}
		MemCytoFrame fr(*this);
		fr.copy(row_idx, col_idx)->write_h5(new_filename);//this flushes the meta data as well
		return CytoFramePtr(new TileCytoFrame(new_filename, false));
	}

	CytoFramePtr copy(uvec idx, bool is_row_indexed, const string & h5_filename = "", bool overwrite = false) const
	{

		if(!overwrite)
			copy_overwrite_check(h5_filename);
		string new_filename = h5_filename;
		if(new_filename == "")
		{
			new_filename = generate_unique_filename(fs::temp_directory_path().string(), "", ".h5");
			fs::remove(new_filename);
		}
		MemCytoFrame fr(*this);
		fr.copy(idx, is_row_indexed)->write_h5(new_filename);//this flushes the meta data as well
		return CytoFramePtr(new TileCytoFrame(new_filename, false));
	}

	/**
	 * copy setter
	 * @param _data
	 */
	void set_data(const EVENT_DATA_VEC & _data)
	{

	}

	void set_data(EVENT_DATA_VEC && _data)
	{
		check_write_permission();
		set_data(_data);
	}
};

};





#endif /* INST_INCLUDE_CYTOLIB_TILECYTOFRAME_HPP_ */
