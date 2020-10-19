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

#ifdef HAVE_TILEDB
#include <tiledb/tiledb>

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
//	EVENT_DATA_VEC read_data(uvec col_idx) const{return EVENT_DATA_VEC();};
	/*
	 * storing context object itself doesn't make  mat_array_ptr_ safely copiable
	 * because array store Context reference, once copied, still point to the original
	 * ctx ref so we want to ensure it lives some cycle with mat_array
	 * TODO:ideally treat them as atomic unit
	 */
	CytoCtx ctx_;

	shared_ptr<tiledb::Array> mat_array_ptr_;
	tiledb::Array & get_mat_array_ref() const {
		if(mat_array_ptr_)
			return *mat_array_ptr_;
		else
			throw(domain_error("Empty mat_array_ptr_!"));
	}
public:

	void flush_meta(){
		//flush the cached meta data from CytoFrame into h5
		if(is_dirty_params)
			flush_params();
		if(is_dirty_keys)
			flush_keys();
		if(is_dirty_pdata)
			flush_pheno_data();

	};
	void flush_params(){
		check_write_permission();
		write_tile_params(uri_, ctx_);
		is_dirty_params = false;

	};

	void flush_keys(){
		check_write_permission();
		write_tile_kw(uri_, ctx_);
		is_dirty_keys = false;

	};
	void flush_pheno_data(){
		check_write_permission();
		write_tile_pd(uri_, ctx_);
		is_dirty_pdata = false;

	};
	void set_readonly(bool flag){
		readonly_ = flag;
	}
	bool get_readonly() const{
		return readonly_ ;
	}


	unsigned n_rows() const{
				return dims[0];
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
	 * @param uri
	 */
	TileCytoFrame(const string & fcs_filename, FCS_READ_PARAM & config, const string & uri
			, bool readonly = false, CytoCtx ctx = CytoCtx()):uri_(uri)
	, is_dirty_params(false), is_dirty_keys(false), is_dirty_pdata(false), ctx_(ctx)
	{
		MemCytoFrame fr(fcs_filename, config);
		fr.read_fcs();
		fr.write_tile(uri, ctx_);
		*this = TileCytoFrame(uri, readonly, true, ctx_);
	}
	/**
	 * constructor from the H5
	 * @param _filename H5 file path
	 */
	TileCytoFrame(const string & uri, bool readonly = true
			, bool init = true, CytoCtx ctx = CytoCtx()):CytoFrame(),uri_(uri)
	, readonly_(readonly), is_dirty_params(false), is_dirty_keys(false), is_dirty_pdata(false), ctx_(ctx)
	{

		if(init)//optionally delay load for the s3 derived cytoframe which needs to reset fapl before load
			init_load();
	}
	void init_load(){
		load_meta();

		fs::path arraypath(uri_);
		auto mat_uri = (arraypath / "mat").string();
		auto ctxptr_ = static_pointer_cast<tiledb::Context>(ctx_.get_ctxptr());
		//open dataset for event data
		mat_array_ptr_ = shared_ptr<tiledb::Array>(new tiledb::Array(*ctxptr_, mat_uri, TILEDB_READ));

//		tiledb::ArraySchema schema(ctx_, mat_uri);
//		auto dm = schema.domain();
		auto nch = mat_array_ptr_->non_empty_domain<int>("channel").second;//dm.dimension("channel").domain<int>().second;
		auto nevent = mat_array_ptr_->non_empty_domain<int>("cell").second;//dm.dimension("cell").domain<int>().second;
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

		uint64_t npd = mat_array_ptr_->metadata_num();
		uint32_t v_num;
		tiledb_datatype_t v_type;
		for (uint64_t i = 0; i < npd; ++i) {
			string pn,pv;

			const void* v;
			mat_array_ptr_->get_metadata_from_index(i, &pn, &v_type, &v_num, &v);
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

		load_params();
		load_kw();
		load_pd();
	}
	void load_kw()
	{
		auto ctxptr_ = static_pointer_cast<tiledb::Context>(ctx_.get_ctxptr());

		tiledb::Array array(*ctxptr_, (fs::path(uri_) / "keywords").string(), TILEDB_READ);
		uint64_t nkw = array.metadata_num();

		keys_.resize(nkw);

		uint32_t v_num;
		tiledb_datatype_t v_type;
		for (uint64_t i = 0; i < nkw; ++i) {
			const void* v;
			array.get_metadata_from_index(i, &(keys_[i].first), &v_type, &v_num, &v);
			if(v)
			{
				keys_[i].second = string(static_cast<const char *>(v));
				keys_[i].second.resize(v_num);
			}
			else
				keys_[i].second = "";


		}

	}
	void load_pd()
	{
		auto ctxptr_ = static_pointer_cast<tiledb::Context>(ctx_.get_ctxptr());
		tiledb::Array array(*ctxptr_, (fs::path(uri_) / "pdata").string(), TILEDB_READ);
		uint64_t nkw = array.metadata_num();

		uint32_t v_num;
		tiledb_datatype_t v_type;
		for (uint64_t i = 0; i < nkw; ++i) {
			const void* v;
			string key, value;
			array.get_metadata_from_index(i, &key, &v_type, &v_num, &v);
			if(v)
			{
				value = string(static_cast<const char *>(v));
				value.resize(v_num);
			}
			pheno_data_[key] = value;
		}

	}
	void load_params()
	{
		auto ctxptr_ = static_pointer_cast<tiledb::Context>(ctx_.get_ctxptr());
		tiledb::Array array(*ctxptr_, (fs::path(uri_) / "params").string(), TILEDB_READ);
		int nch = array.metadata_num();
		if(nch > 0)
		{
			const std::vector<int> subarray = {1, nch};

			tiledb::Query query(*ctxptr_, array);
			query.set_subarray(subarray);
			query.set_layout(TILEDB_GLOBAL_ORDER);
			vector<float> min_vec(nch);
			query.set_buffer("min", min_vec);
			vector<float> max_vec(nch);
			query.set_buffer("max", max_vec);
	//		auto max_buf = array.max_buffer_elements(subarray);
	//		auto nbuf_size = max_buf["channel"].second;
	//		vector<char> ch_vec(nbuf_size);
	//		if(max_buf["channel"].first!=nch)
	//			throw(domain_error("channel attribute length is not consistent with its marker meta data size!"));
	//		vector<uint64_t> off(nch);
	//		query.set_buffer("channel", off, ch_vec);

			query.submit();
			query.finalize();

			params.resize(nch);
			//get channel in order
			tiledb::Array array1(*ctxptr_, (fs::path(uri_) / "channel_idx").string(), TILEDB_READ);
			uint32_t v_num;
			tiledb_datatype_t v_type;
			for (int i = 0; i < nch; ++i) {
					const void* v;
					string channel;
					array1.get_metadata_from_index(i, &channel, &v_type, &v_num, &v);
					int idx = *(static_cast<const int *>(v));
					params[idx].channel = channel;

				}


			for (int i = 0; i < nch; ++i) {
	//			auto start = off[i];
	//			auto nsize = (i< (nch - 1)?off[i+1]:nbuf_size) - start;
	//			params[i].channel = string(&ch_vec[start], nsize);
				const void* v;
				array.get_metadata(params[i].channel, &v_type, &v_num, &v);
				if(v)
				{
					params[i].marker = string(static_cast<const char *>(v));
					params[i].marker.resize(v_num);
				}

				params[i].min = min_vec[i];
				params[i].max = max_vec[i];

			}

			build_hash();
		}
	}
	string get_uri() const{
		return uri_;
	}
	void check_write_permission() const{
		if(readonly_)
			throw(domain_error("Can't write to the read-only TileCytoFrame object!"));

	}

	void convertToPb(pb::CytoFrame & fr_pb
			, const string & uri
			, CytoFileOption file_opt
			, const CytoCtx & ctx = CytoCtx()) const
	{
			fr_pb.set_is_h5(true);
			if(file_opt != CytoFileOption::skip)
			{
				auto ctxptr_ = static_pointer_cast<tiledb::Context>(ctx_.get_ctxptr());

				tiledb::VFS vfs(*ctxptr_);
				auto filepath = fs::path(uri);
				auto dest = filepath.parent_path();
				auto src = fs::path(uri_).parent_path();
				bool is_eq;
				if(is_remote_path(uri)||is_remote_path(uri_))
					is_eq = uri == uri_;
				else
					is_eq = fs::equivalent(src, dest);
				if(!is_eq)
				{
					if(vfs.is_dir(uri))
						vfs.remove_dir(uri);

					switch(file_opt)
					{
					case CytoFileOption::copy:
						{
							write_tile(uri, ctx_);//todo:cp dir for more efficient operations
							break;
						}
					case CytoFileOption::move:
						{
							vfs.move_dir(uri_, uri);
							break;
						}
					case CytoFileOption::link:
						{
							throw(logic_error("'link' option for TileCytoFrame is no longer supported!"));
							break;
						}
					case CytoFileOption::symlink:
						{
							if(is_remote_path(uri)||is_remote_path(uri_))
								throw(logic_error("'symlink' option for remote TileCytoFrame is not supported!"));
							fs::create_symlink(uri_, uri);
							break;
						}
					default:
						throw(logic_error("invalid file_opt!"));
					}
				}
			}
		}
	FileFormat get_backend_type()  const{
				return FileFormat::TILE;
			};
	/**
	 * Read data from disk.
	 * The caller will directly receive the data vector without the copy overhead thanks to the move semantics supported for vector container in c++11
	 * @return
	 */


	EVENT_DATA_VEC get_data() const
	{
		auto & array = get_mat_array_ref();
		if(!array.is_open()||array.query_type() != TILEDB_READ)
		{
			array.open(TILEDB_READ);
		}
		int ncol = dims[1];
		int nrow = dims[0];

		if(ncol*nrow>0)
		{

			const std::vector<int> subarray = {1, nrow, 1, ncol};
			auto ctxptr_ = static_pointer_cast<tiledb::Context>(ctx_.get_ctxptr());
			tiledb::Query query(*ctxptr_, *mat_array_ptr_);
			query.set_subarray(subarray);
			query.set_layout(TILEDB_GLOBAL_ORDER);

			arma::Mat<float> buf(nrow, ncol);

			query.set_buffer("mat", buf.memptr(), nrow*ncol);
			query.submit();
			query.finalize();
	//		double runtime = (gettime() - start);
	//		cout << "get all: " << runtime << endl;

			return arma::conv_to<arma::mat>::from(buf);
		}
		else
			return EVENT_DATA_VEC(nrow, ncol);
	}

	EVENT_DATA_VEC read_cols(uvec cidx) const
	{
		auto & array = get_mat_array_ref();
		if(!array.is_open()||array.query_type() != TILEDB_READ)
		{
			array.open(TILEDB_READ);
		}
		auto ctxptr_ = static_pointer_cast<tiledb::Context>(ctx_.get_ctxptr());
		tiledb::Query query(*ctxptr_, *mat_array_ptr_);
		query.set_layout(TILEDB_COL_MAJOR);
		int ncol,nrow, dim_idx;
		ncol = cidx.size();
		nrow = dims[0];
		dim_idx = 1;
		query.add_range<int>(0, 1, nrow);//select all rows

		//tiledb idx starting from 1
		for(int i : cidx)
			query.add_range<int>(dim_idx, i+1, i+1);


		arma::Mat<float> buf(nrow, ncol);



		query.set_buffer("mat", buf.memptr(), nrow * ncol);
		query.submit();
		query.finalize();

		return arma::conv_to<arma::mat>::from(buf);

	}
	/**
	 * Partial IO
	 * @param col_idx
	 * @return
	 */
	EVENT_DATA_VEC get_data(uvec idx, bool is_col) const
	{
		EVENT_DATA_VEC data;

		if(is_col)
		{
			data = read_cols(idx);
		}
		else
		{
			data = get_data().rows(idx);

		}




		return data;

	}
	EVENT_DATA_VEC get_data(uvec row_idx, uvec col_idx) const
	{


		return read_cols(col_idx).rows(row_idx);

	}
	/*
	 * protect the h5 from being overwritten accidentally
	 * which will make the original cf object invalid
	 */
	void copy_overwrite_check(const string & dest, bool overwrite) const
	{
		// Check if trying to write to same file
		if(fs::equivalent(fs::path(uri_), fs::path(dest))){
			// First check if that is even allowed by this CytoFrame
			check_write_permission();
			// Then make sure it has been explicitly approved by the caller
			if(!overwrite)
				throw(domain_error("Copying TileCytoFrame to itself is not supported! "+ dest));
		}
	}

	CytoFramePtr copy(const string & uri = "", bool overwrite = false) const
	{
		copy_overwrite_check(uri, overwrite);

		string new_uri = uri;
		if(new_uri == "")
		{
			new_uri = generate_unique_filename(fs::temp_directory_path().string(), "", ".tile");
			fs::remove(new_uri);
		}
		recursive_copy(uri_, new_uri);
		CytoFramePtr ptr(new TileCytoFrame(new_uri, false, true, ctx_));
		//copy cached meta
		ptr->set_params(get_params());
		ptr->set_keywords(get_keywords());
		ptr->set_pheno_data(get_pheno_data());
		return ptr;
	}
	CytoFramePtr copy(uvec row_idx, uvec col_idx, const string & uri = "", bool overwrite = false) const
	{
		copy_overwrite_check(uri, overwrite);

		string new_uri = uri;
		if(new_uri == "")
		{
			new_uri = generate_unique_filename(fs::temp_directory_path().string(), "", ".tile");
			fs::remove(new_uri);
		}
		MemCytoFrame fr(*this);
		fr.copy(row_idx, col_idx)->write_tile(new_uri);//this flushes the meta data as well
		return CytoFramePtr(new TileCytoFrame(new_uri, false, true, ctx_));
	}

	CytoFramePtr copy(uvec idx, bool is_row_indexed, const string & uri = "", bool overwrite = false) const
	{
		copy_overwrite_check(uri, overwrite);

		string new_uri = uri;
		if(new_uri == "")
		{
			new_uri = generate_unique_filename(fs::temp_directory_path().string(), "", ".tile");
			fs::remove(new_uri);
		}
		MemCytoFrame fr(*this);
		fr.copy(idx, is_row_indexed)->write_tile(new_uri);//this flushes the meta data as well
		return CytoFramePtr(new TileCytoFrame(new_uri, false, true, ctx_));
	}

	/**
	 * copy setter
	 * @param _data
	 */
	void set_data(const EVENT_DATA_VEC & _data)
	{
		check_write_permission();

		// write_tile_data() will resize if necessary but it involves a delete and re-write
		write_tile_data(uri_, _data, ctx_);

		auto & array = get_mat_array_ref();

		if((dims[0] != _data.n_rows) || (dims[1] != _data.n_cols)){
			// If an array re-size was necessary, close() followed by open()
			// seems necessary
			array.close();
			array.open(TILEDB_READ);
			// And the dimensions need to be updated
			dims[0] = _data.n_rows;
			dims[1] = _data.n_cols;
		}else{
			// Otherwise, it is a little more efficient to use reopen()
			if(array.is_open())
				array.reopen();
		}

	}

	void set_data(EVENT_DATA_VEC && _data)
	{
		set_data(_data);
	}

	void append_data_columns(const EVENT_DATA_VEC & new_cols)
	{
		EVENT_DATA_VEC this_data = get_data();
		this_data.insert_cols(this_data.n_cols, new_cols);
		set_data(this_data);
	}
};

};




#endif
#endif /* INST_INCLUDE_CYTOLIB_TILECYTOFRAME_HPP_ */
