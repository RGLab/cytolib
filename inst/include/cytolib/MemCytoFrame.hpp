/*
 * MemCytoFrame.hpp
 *
 *  Created on: Sep 18, 2017
 *      Author: wjiang2
 */

#ifndef INST_INCLUDE_CYTOLIB_MEMCYTOFRAME_HPP_
#define INST_INCLUDE_CYTOLIB_MEMCYTOFRAME_HPP_
#include "CytoFrame.hpp"
#include "readFCSdata.hpp"


struct FCS_READ_PARAM{
	FCS_READ_HEADER_PARAM header;
	FCS_READ_DATA_PARAM data;
};
class MemCytoFrame: public CytoFrame{
	EVENT_DATA_PTR data;

public:
	MemCytoFrame(const string &filename, FCS_READ_PARAM &, bool onlyTxt);
	void compensate(const compensation &);
//	void transform(const transformation &);
	void save(const string & filename, FrameType type);
	EVENT_DATA_TYPE * getData();
	EVENT_DATA_TYPE * getData(const string &, ColType);


};

#ifdef _OPENMP
#define gettime() omp_get_wtime()
#else
#define gettime() clock()/(double)(CLOCKS_PER_SEC / 1000)
#endif




#endif /* INST_INCLUDE_CYTOLIB_MEMCYTOFRAME_HPP_ */
