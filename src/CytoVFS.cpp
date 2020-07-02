/*
 * CytoVFS.cpp
 *
 *  Created on: Jun 22, 2020
 *      Author: wjiang2
 */

#include <cytolib/CytoVFS.hpp>
#ifdef HAVE_TILEDB
#include <tiledb/tiledb>
#else
#include <fstream>
#endif
namespace cytolib
{
#ifdef HAVE_TILEDB
	const bool have_tiledb = true;
/**
 *  implementation tiledb::VFS
 * @param file
 * @param buf
 */
	CytoCtx::CytoCtx():access_key_id_(""),access_key_(""),region_("us-west-1"),num_threads_(1){
		init_ctxptr();
	};

	CytoCtx::CytoCtx(const string & secret_id
						, const string & secret_key
						, const string & aws_region
						, int num_threads):access_key_id_(secret_id)
			, access_key_(secret_key), region_(aws_region), num_threads_(num_threads){
		init_ctxptr();
	};

	void CytoCtx::init_ctxptr(){
		tiledb::Config cfg;
		cfg["vfs.s3.aws_access_key_id"] = access_key_id_;
		cfg["vfs.s3.aws_secret_access_key"] =  access_key_;
		cfg["vfs.s3.region"] =  region_;
		cfg["sm.num_reader_threads"] = num_threads_;

		ctxptr_.reset(new tiledb::Context(cfg));

	}

	CytoVFS::CytoVFS(CytoCtx ctx){
		auto ctxptr = static_pointer_cast<tiledb::Context>(ctx.get_ctxptr());
		vfsptr_.reset(new tiledb::VFS(*ctxptr));
	};
	void CytoVFS::write_buf(const string & file, const string & buf){
		auto & vfs = *(static_pointer_cast<tiledb::VFS>(vfsptr_));

		tiledb::VFS::filebuf sbuf(vfs);
		sbuf.open(file, ios::out);
		ostream output(&sbuf);
		output.write(&buf[0], buf.size());

	}
	vector<char> CytoVFS::read_buf(const string & file){
		auto & vfs = *(static_pointer_cast<tiledb::VFS>(vfsptr_));
		tiledb::VFS::filebuf sbuf(vfs);
		sbuf.open(file, ios::in);
		istream input(&sbuf);
		if (!input)
			throw(invalid_argument("File not found.." + file));
		 auto length = file_size(file);
		 vector<char> buf(length);
		 input.read(buf.data(), length);
		 return buf;
		}
/**
	 * forward the apis to tiledb::VFS
	 */
	vector<string> CytoVFS::ls(string p)
	{
		auto & vfs = *(static_pointer_cast<tiledb::VFS>(vfsptr_));
		return vfs.ls(p);
	}
	bool CytoVFS::is_dir(string p) const{
		auto & vfs = *(static_pointer_cast<tiledb::VFS>(vfsptr_));
		return vfs.is_dir(p);
	}
	bool CytoVFS::is_file(string p){
		auto & vfs = *(static_pointer_cast<tiledb::VFS>(vfsptr_));
		return vfs.is_file(p);
	}
	void CytoVFS::remove_dir(string p){
		auto & vfs = *(static_pointer_cast<tiledb::VFS>(vfsptr_));
		vfs.remove_dir(p);
	}
	void CytoVFS::create_dir(string p){
		auto & vfs = *(static_pointer_cast<tiledb::VFS>(vfsptr_));
		vfs.create_dir(p);
	}
	void CytoVFS::move_dir(string p, string p1){
		auto & vfs = *(static_pointer_cast<tiledb::VFS>(vfsptr_));
		vfs.move_dir(p, p1);
	}
	int CytoVFS::file_size(string p){
		auto & vfs = *(static_pointer_cast<tiledb::VFS>(vfsptr_));
		return vfs.file_size(p);
	}

#else
	const bool have_tiledb = false;
	/**
	 *  implementation based on std::fs
	 */

	//dummy
	CytoCtx::CytoCtx():access_key_id_(""),access_key_(""),region_("us-west-1"),num_threads_(1){};
	//dummy
	CytoCtx::CytoCtx(const string & secret_id
						, const string & secret_key
						, const string & aws_region
						, int num_threads):access_key_id_(secret_id)
			, access_key_(secret_key), region_(aws_region), num_threads_(num_threads){};
	//dummy
	void CytoCtx::init_ctxptr(){}
	//dummy
	CytoVFS::CytoVFS(CytoCtx ctx){};

	/**
	 * forward the apis to std::fs
	 */
	void CytoVFS::write_buf(const string & file, const string & buf){
		ofstream output(file, ios::out | ios::binary);
		output.write(&buf[0], buf.size());

	}
	vector<char> CytoVFS::read_buf(const string & file){
		ifstream input(file, ios::in | ios::binary);
		if (!input)
			throw(invalid_argument("File not found.." + file));
		 auto length = file_size(file);
		 vector<char> buf(length);
		 input.read(buf.data(), length);
		 return buf;
	}

	vector<string> CytoVFS::ls(string p)
		{

			vector<string>res;
			for(auto e : fs::directory_iterator(p))
				res.push_back(fs::path(e).string());
			return res;
		}
	bool CytoVFS::is_dir(string p) const{return fs::is_directory(p);}
	bool CytoVFS::is_file(string p){return !fs::is_directory(p)&&fs::exists(p);}
	void CytoVFS::remove_dir(string p){fs::remove_all(p);}
	void CytoVFS::create_dir(string p){ fs::create_directory(p);}
	void CytoVFS::move_dir(string p, string p1){fs::rename(p, p1);}
	int CytoVFS::file_size(string p){return fs::file_size(p);}

#endif

}



