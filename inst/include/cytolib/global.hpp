/* Copyright 2019 Fred Hutchinson Cancer Research Center
 * See the included LICENSE file for details on the license that is granted to the
 * user of this software.
 * global.hpp
 *
 *  Created on: Mar 31, 2014
 *      Author: wjiang2
 */

#ifndef GLOBAL_HPP_
#define GLOBAL_HPP_

#include <iostream>
#include <string>
#include <memory>
#include <ctype.h>
#include <vector>
#include <chrono>
#include <unordered_set>
#include "datatype.hpp"
#include "CytoVFS.hpp"
using namespace std;

namespace cytolib
{
	enum class FileFormat {TILE, H5, MEM};
	inline string fmt_to_str(FileFormat fmt)
	{
		switch(fmt)
		{
		case FileFormat::H5:
			return "h5";
		case FileFormat::TILE:
			return "tile";
		default:
			return "mem";
		}

	}
	FileFormat uri_backend_type(const string & path, const CytoVFS & vfs);
	string s3_to_http(string uri);
	void check_sample_guid(const string & sample_guid);
	bool is_remote_path(const string &);

	#define GATING_SET_LEVEL 1
	#define GATING_HIERARCHY_LEVEL 2
	#define POPULATION_LEVEL 3
	#define GATE_LEVEL 4

	void PRINT(string a);
	void PRINT(const char * a);

	extern vector<string> spillover_keys;
	extern unsigned short g_loglevel;// debug print is turned off by default
	extern bool my_throw_on_error;//can be toggle off to get a partially parsed gating tree for debugging purpose

	const int bsti = 1;  // Byte swap test integer
	#define is_host_big_endian() ( (*(char*)&bsti) == 0 )

	enum class ColType {channel, marker, unknown};

	#define PRT true
	string fs_tmp_path();

	string generate_unique_filename(const string & dir, const string & prefix, const string & suffix);

	string generate_unique_dir(const string & dir, const string & prefix);
	void recursive_copy(const fs::path &src, const fs::path &dst);

	/**
	 * Generate time stamp as string
	 * @return
	 */
	string generate_timestamp();
	/**
	 * Generate uniquely identifiable id
	 * @return
	 */
	string generate_uid();
	string generate_uid_old(int len);
	string path_dir_name(const string & full_path);
	string path_base_name(const string & full_path);

	/**
	 * validity check and reordering(if needed) channels for newv
	 * @param oldv the existing data
	 * @param newv the new data to be checked
	 * @param sample_uid the sample name of the new data
	 */
	template<class T1, class T2> T2 channel_consistency_check(const T1 & oldv, const T2 & newv, const string & sample_uid)
	{
		auto res = newv;
		if(oldv.size()>0)
		{
			string msg = "Found channel inconsistency across samples. ";
			auto c1 = oldv.get_channels();
			unordered_set<string> old_ch(c1.begin(), c1.end());
			auto c2 = newv.get_channels();
			unordered_set<string> new_ch(c2.begin(), c2.end());
			for(auto ch : c1)
			{
				if(new_ch.find(ch) == new_ch.end())
					throw(domain_error(msg + "'" + ch + "' is missing from "  + sample_uid));

			}
			for(auto ch : c2)
			{
				if(old_ch.find(ch) == old_ch.end())
					throw(domain_error(msg + sample_uid + " has the channel '" + ch + "' that is not found in other samples!"));

			}
			//check if need to order newv
			for(unsigned i = 0; i < c1.size(); i++)
			{
				if(c1[i]!=c2[i])
				{
					res.cols_(c1, ColType::channel);
					break;
				}
			}
		}
		return res;
	}
	struct TM_ext
	{
		tm _time;
		EVENT_DATA_TYPE fractional_secs;
		TM_ext():fractional_secs(0){
			time_t rawtime;
			time(&rawtime);
			struct tm * timeinfo = localtime (&rawtime);//The returned value points to an internal object
			//init time member to avoid random values for day,month,year
			_time = *timeinfo;
		}
	};
	/**
	 * Parse the time string with fractional seconds
	 * std lib doesn't handle and boost::posix_time is not header-only
	 * @param s_time time string "H:M:S.ss"
	 * @return
	 */
	TM_ext parse_time_with_fractional_seconds(const string s_time);

	#ifdef _OPENMP
#define gettime() omp_get_wtime()
#else
#define gettime() clock()/(double)(CLOCKS_PER_SEC / 1000)
#endif
};

#endif /* GLOBAL_HPP_ */
