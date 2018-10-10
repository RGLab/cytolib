/*
 * global.hpp
 *
 *  Created on: Mar 31, 2014
 *      Author: wjiang2
 */

#ifndef GLOBAL_HPP_
#define GLOBAL_HPP_

#define GATING_SET_LEVEL 1
#define GATING_HIERARCHY_LEVEL 2
#define POPULATION_LEVEL 3
#define GATE_LEVEL 4


#ifdef ROUT
#include <R_ext/Print.h>
#endif

#include <iostream>
#include <string>
#include <memory>
#include <ctype.h>
#include <vector>

#include <armadillo>
using namespace arma;

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <chrono>
#include <H5Cpp.h>

using namespace std;
extern unsigned short g_loglevel;// debug print is turned off by default
extern bool my_throw_on_error;//can be toggle off to get a partially parsed gating tree for debugging purpose
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
using namespace H5;

namespace cytolib
{
	#define CYTOLIB_INIT() \
			bool my_throw_on_error = true;\
			unsigned short g_loglevel = 0;\

	inline void PRINT(string a){
	#ifdef ROUT
	 Rprintf(a.c_str());
	#else
	 cout << a;
	#endif

	}
	inline void PRINT(const char * a){
	#ifdef ROUT
	 Rprintf(a);
	#else
	 cout << a;
	#endif

	}



	const int bsti = 1;  // Byte swap test integer
	#define is_host_big_endian() ( (*(char*)&bsti) == 0 )

	typedef double EVENT_DATA_TYPE;
	typedef arma::Mat<EVENT_DATA_TYPE> EVENT_DATA_VEC;
	#define EVENT_DATA_TYPE_IN_MEM_H5 PredType::NATIVE_DOUBLE

	#define PRT true

	inline string generate_temp_filename(string dir = "/tmp")
	{
		char tmp[15] = "/tmp/XXXXXX.h5";
		int fid = mkstemps(tmp, 3);
		if(fid == -1)
			throw(domain_error("Can't create the unique temp file: " + string(tmp)));

		close(fid);
		return tmp;
	}

	/**
	 * Generate time stamp as string
	 * @return
	 */
	inline string generate_timestamp()
	{
		time_t rawtime;
		time(&rawtime);
		char out[13];
		strftime(out, 13, "%y%m%d%H%M%S", localtime(&rawtime));
		out[12] ='\0';
		return string(out);
	}
	/**
	 * Generate uniquely identifiable id (pseudo guid)
	 * @return
	 */
	inline string generate_uid(int len)
	{
//		int t = time(NULL);//this only returns second-wise precision, not sufficient for distinguishing multiple gatingsets that are generated within short period

		chrono::milliseconds ms = chrono::duration_cast< chrono::milliseconds >(
													chrono::system_clock::now().time_since_epoch()
													);
		int t = ms.count();
//		cout << "random seed value: " << t << endl;
		srand (t);
		char s[len+1];
		static const char alphanum[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
		s[0] = alphanum[rand() % (sizeof(alphanum) - 11) + 10];
		for (int i = 1; i < len; ++i) {
			s[i] = alphanum[rand() % (sizeof(alphanum) - 1)];
		}
		s[len] = '\0';
		return string(s);
	}
	inline string path_dir_name(const string & full_path)
	{
//		//TODO:handle windows path
//		size_t i_last_slash = full_path.find_last_of('/');
//		return full_path.substr(0,i_last_slash);
		return fs::path(full_path).parent_path().string();
	}
	inline string path_base_name(const string & full_path)
	{
//		//TODO:handle windows path
//		size_t i_last_slash = full_path.find_last_of('/');
//		return full_path.substr(i_last_slash);
		return fs::path(full_path).filename().string();
	}

	struct TM_ext
	{
		tm time;
		EVENT_DATA_TYPE fractional_secs;
	};
	/**
	 * Parse the time string with fractional seconds
	 * std lib doesn't handle and boost::posix_time is not header-only
	 * @param s_time time string "H:M:S.ss"
	 * @return
	 */
	inline TM_ext parse_time_with_fractional_seconds(const string s_time){
		TM_ext res;
		vector<string> time_vec;
		//split the H:M:S.ms by .
		boost::split(time_vec, s_time, boost::is_any_of("."));
		//using std lib to parse the first half
		strptime(time_vec[0].c_str(), "%H:%M:%S", &(res.time));
		 //parse the second half as fractional seconds
		if(time_vec.size()==2)
		{
			res.fractional_secs = boost::lexical_cast<EVENT_DATA_TYPE>(time_vec[1]);
		}
		else
			res.fractional_secs = 0;
		return res;
	}
};

#endif /* GLOBAL_HPP_ */
