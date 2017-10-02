/*
 * config.hpp
 *
 *  Created on: Sep 29, 2017
 *      Author: wjiang2
 */

#ifndef INST_INCLUDE_CYTOLIB_CONFIG_HPP_
#define INST_INCLUDE_CYTOLIB_CONFIG_HPP_
#ifdef ROUT
//#include <Rcpp.h>
//using namespace Rcpp;
#define COUT cout //TODO: originally it was Rout, to get around Rcpp, we will have to use Rprintf/printf to rewrite all the console print
#endif


#ifndef ROUT
#define COUT cout
#endif





#endif /* INST_INCLUDE_CYTOLIB_CONFIG_HPP_ */
