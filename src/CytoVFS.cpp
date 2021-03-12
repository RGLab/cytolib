/*
 * CytoVFS.cpp
 *
 *  Created on: Jun 22, 2020
 *      Author: wjiang2
 */

#include <cytolib/CytoVFS.hpp>
#include <fstream>
namespace cytolib
{
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

	vector<string> CytoVFS::ls(string p) const
		{

			vector<string>res;
			for(auto e : fs::directory_iterator(p))
				res.push_back(fs::path(e).string());
			return res;
		}
	bool CytoVFS::is_dir(string p) const{return fs::is_directory(p);}
	bool CytoVFS::is_file(string p) const{return !fs::is_directory(p)&&fs::exists(p);}
	void CytoVFS::remove_dir(string p){fs::remove_all(p);}
	void CytoVFS::create_dir(string p){ fs::create_directory(p);}
	void CytoVFS::move_dir(string p, string p1){fs::rename(p, p1);}
	int CytoVFS::file_size(string p){return fs::file_size(p);}

}



