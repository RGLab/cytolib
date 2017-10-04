/*
 * global.cpp
 *
 *  Created on: Mar 31, 2014
 *      Author: wjiang2
 */

#include <cytolib/global.hpp>

unsigned short g_loglevel = 0;
bool my_throw_on_error = true;//can be toggle off to get a partially parsed gating tree for debugging purpose

void PRINT(string a){
#ifdef ROUT
 Rprintf(a.c_str());
#else
 cout << a;
#endif

}
void PRINT(const char * a){
#ifdef ROUT
 Rprintf(a);
#else
 cout << a;
#endif

}

