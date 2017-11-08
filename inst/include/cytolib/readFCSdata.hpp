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
#ifdef _OPENMP
#include <omp.h>
#endif
typedef unsigned char BYTE;


/**
 * The struct stores all the parsing arguments for events data
 */
struct FCS_READ_DATA_PARAM{
	 bool scale, truncate_max_range, truncate_min_val;
	 EVENT_DATA_TYPE decades, min_limit;
	 TransformType transform;
	 int num_threads; //number of cores to be used for parallel-read of data (channel / core)
	 vector<int> which_lines; //select rows to be read in
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
	 }


};
/**
 * parse the data segment of FCS
 *
 * @param in (input) file stream object opened from FCS file
 * @param output (output) the container to store the parsed data
 * @param header (input) the header parsed by readHeaderAndText
 * @param keys (input&output) the keywords parsed by readHeaderAndText
 * @param params (input&output) the params parsed by readHeaderAndText
 * @param config (input) the parsing arguments for data
 */
void readFCSdata(ifstream &in, EVENT_DATA_VEC & output, const FCS_Header & header,KEY_WORDS & keys, vector<cytoParam> & params,  FCS_READ_DATA_PARAM & config);



#endif /* INST_INCLUDE_CYTOLIB_READFCSDATA_HPP_ */
