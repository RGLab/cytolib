/*
 * readFCSdata.hpp
 *
 *  Created on: Sep 21, 2017
 *      Author: wjiang2
 */

#ifndef INST_INCLUDE_CYTOLIB_READFCSDATA_HPP_
#define INST_INCLUDE_CYTOLIB_READFCSDATA_HPP_
#include "readFCSHeader.hpp"
#include <memory>
#include <ctype.h>

typedef unsigned char BYTE;
const int bsti = 1;  // Byte swap test integer

#define is_host_big_endian() ( (*(char*)&bsti) == 0 )


struct FCS_READ_DATA_PARAM{
	 bool scale, truncate_max_range, truncate_min_val;
	 EVENT_DATA_TYPE decades, min_limit;
	 TransformType transform;
	 int num_threads;
	 FCS_READ_DATA_PARAM(){
		 scale = false;
		 truncate_max_range = true;
		 truncate_min_val = false;
		 decades = 0;
		 min_limit=-111;
		 transform =  TransformType::linearize;
		 num_threads = 1;
	 }


};
EVENT_DATA_PTR readFCSdata(ifstream &in, const FCS_Header & header,KEY_WORDS & keys, vector<cytoParam> & params,int & nEvents, FCS_READ_DATA_PARAM &);



#endif /* INST_INCLUDE_CYTOLIB_READFCSDATA_HPP_ */
