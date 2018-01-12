/*
 * global.hpp
 *
 *  Created on: Mar 31, 2014
 *      Author: wjiang2
 */

#ifndef GLOBAL_HPP_
#define GLOBAL_HPP_

enum loglevel_t{
	NO_LOG = 0,
	GATING_SET_LEVEL = 1,
	GATING_HIERARCHY_LEVEL = 2,
	POPULATION_LEVEL = 3,
	GATE_LEVEL = 4
	};

#ifdef ROUT
#include <R_ext/Print.h>
#endif

#include <iostream>
#include <string>
#include <memory>
#include <ctype.h>
#include <vector>
using namespace std;


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

#endif /* GLOBAL_HPP_ */
