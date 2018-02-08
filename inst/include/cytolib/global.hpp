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
using namespace std;
extern unsigned short g_loglevel;// debug print is turned off by default
extern bool my_throw_on_error;//can be toggle off to get a partially parsed gating tree for debugging purpose

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
	typedef vector<EVENT_DATA_TYPE> EVENT_DATA_VEC;
	typedef unique_ptr<EVENT_DATA_TYPE[] > EVENT_DATA_PTR;

	#define PRT true

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
	inline string generate_guid(int len)
	{
		char s[len];
		static const char alphanum[] =
				"0123456789"
				"ABCDEFGHIJKLMNOPQRSTUVWXYZ"
				"abcdefghijklmnopqrstuvwxyz";

		for (int i = 0; i < len; ++i) {
			s[i] = alphanum[rand() % (sizeof(alphanum) - 1)];
		}

		return string(s);
	}

};

#endif /* GLOBAL_HPP_ */
