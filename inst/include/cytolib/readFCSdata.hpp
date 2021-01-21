/* Copyright 2019 Fred Hutchinson Cancer Research Center
 * See the included LICENSE file for details on the license that is granted to the
 * user of this software.
 * readFCSdata.hpp
 *
 *  Created on: Sep 21, 2017
 *      Author: wjiang2
 */

#ifndef INST_INCLUDE_CYTOLIB_READFCSDATA_HPP_
#define INST_INCLUDE_CYTOLIB_READFCSDATA_HPP_
#include "readFCSHeader.hpp"


#ifdef _OPENMP
#include <omp.h>
#endif


namespace cytolib
{
typedef unsigned char BYTE;


/**
 * The struct stores all the parsing arguments for events data
 */
struct FCS_READ_DATA_PARAM{
	 bool scale, truncate_max_range, truncate_min_val;
	 EVENT_DATA_TYPE decades, min_limit;
	 TransformType transform;
	 int num_threads; //number of cores to be used for parallel-read of data (channel / core)
	 vector<int64_t> which_lines; //select rows to be read in
	 int seed;
	 bool isTransformed;//record the outcome after parsing
	 FCS_READ_DATA_PARAM(){
		 scale = false;
		 truncate_max_range = true;
		 truncate_min_val = false;
		 decades = 0;
		 min_limit=-111;
		 transform =  TransformType::linearize;
		 num_threads = 1;
		 isTransformed = false;
		 seed = 1;
	 }


};


};



#endif /* INST_INCLUDE_CYTOLIB_READFCSDATA_HPP_ */
