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

/**
 * the container that stores the different FCS parse arguments
 */
struct FCS_READ_PARAM{
	FCS_READ_HEADER_PARAM header;
	FCS_READ_DATA_PARAM data;
};

/**
 * The class represents the in-memory version of CytoFrame, which stores and owns the events data
 */
class MemCytoFrame: public CytoFrame{
	EVENT_DATA_VEC data;

public:
	/**
	 * Constructor from the FCS file
	 *
	 * @param filename FCS file path
	 * @param config the parse arguments.
	 * @param onlyTxt flag indicates whether to only parse text segment (which contains the keywords)
	 */
	MemCytoFrame(const string &filename, FCS_READ_PARAM & config, bool onlyTxt);
	void compensate(const compensation &);
	int nRow();

	EVENT_DATA_VEC getData();
	EVENT_DATA_VEC getData(const string &, ColType);


};

#ifdef _OPENMP
#define gettime() omp_get_wtime()
#else
#define gettime() clock()/(double)(CLOCKS_PER_SEC / 1000)
#endif




#endif /* INST_INCLUDE_CYTOLIB_MEMCYTOFRAME_HPP_ */
