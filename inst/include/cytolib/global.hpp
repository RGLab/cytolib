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
#include "datatype.hpp"
using namespace std;

namespace cytolib
{
	#define GATING_SET_LEVEL 1
	#define GATING_HIERARCHY_LEVEL 2
	#define POPULATION_LEVEL 3
	#define GATE_LEVEL 4

	void PRINT(string a);
	void PRINT(const char * a);


	extern unsigned short g_loglevel;// debug print is turned off by default
	extern bool my_throw_on_error;//can be toggle off to get a partially parsed gating tree for debugging purpose

	const int bsti = 1;  // Byte swap test integer
	#define is_host_big_endian() ( (*(char*)&bsti) == 0 )


	#define PRT true
	string fs_tmp_path();

	string generate_unique_filename(const string & dir, const string & prefix, const string & suffix);

	string generate_unique_dir(const string & dir, const string & prefix);
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
