/*Copyright 2019 Fred Hutchinson Cancer Research Center
 * See the included LICENSE file for details on the license that is granted to the
 * user of this software.
 * CytoVFS.hpp
 *
 *  Created on: Jun 22, 2020
 *      Author: wjiang2
 */

#ifndef INST_INCLUDE_CYTOLIB_CytoVFS_HPP_
#define INST_INCLUDE_CYTOLIB_CytoVFS_HPP_
#include <string>
#include <vector>
#include <memory>
#include <map>
using namespace std;
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

namespace cytolib
{
/**
 * abstract layer so that the  tiledb:Context and tiledb::VFS is separated from
 * the rest code base and they can work when tiledb is disabled at compile time
 * no longer needed since tiledb support is dropped
 */
	class CytoCtx
	{
		string access_key_id_;
		string access_key_;
		string region_;
		int num_threads_;
		shared_ptr<void> ctxptr_;
		void init_ctxptr();
	public:
			CytoCtx();

			/**
			 * API to be exposed to R
			 * @return
			 */
			map<string, string> get_config() const{
				map<string, string> res;
				res["access_key_id"] = access_key_id_;
				res["access_key"] = access_key_;
				res["region"] = region_;
				res["num_threads"] = to_string(num_threads_);
				return res;
			}
			CytoCtx(const string & secret_id
						, const string & secret_key
						, const string & aws_region
						, int num_threads = 1);
			shared_ptr<void> get_ctxptr() const{return ctxptr_;};

	};

//	typedef shared_ptr<CytoCtx> CtxPtr;

	/**
	 *
	 * dummy wrapper class around std::fs to mimic tiledb::VFS api
	 */
	class CytoVFS
	{
		shared_ptr<void> vfsptr_;

	public:
		CytoVFS(CytoCtx ctx);
		void write_buf(const string & file, const string & buf);
		vector<char> read_buf(const string & file);
		vector<string> ls (string p) const;
		bool is_dir(string p) const;
		bool is_file(string p) const;
		void remove_dir(string p);
		void create_dir(string p);
		void move_dir(string p, string p1);
		int file_size(string p);
	};


}
#endif /* INST_INCLUDE_CYTOLIB_CytoVFS_HPP_ */
