/*
 * MemCytoFrame.hpp
 *
 *  Created on: Sep 18, 2017
 *      Author: wjiang2
 */

#ifndef INST_INCLUDE_CYTOLIB_MEMCYTOFRAME_HPP_
#define INST_INCLUDE_CYTOLIB_MEMCYTOFRAME_HPP_
#include "CytoFrame.hpp"
#include <fstream>
#include <cstring>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

struct FCS_Header{
	float FCSversion;
	int textstart, textend, datastart, dataend, anastart, anaend, additional;
};

class MemCytoFrame: public CytoFrame{
	EVENT_DATA_TYPE * data;

public:
	~MemCytoFrame();
	MemCytoFrame(){data = NULL;};
	MemCytoFrame(const string &filename, bool isEmptyKeyValue, int nDataset, bool scale, double decades, double min_limit, bool truncate_max_range, bool ignoreTextOffset, bool onlyTxt);
	void compensate(const compensation &);
//	void transform(const transformation &);
	void save(const string & filename, FrameType type);
	EVENT_DATA_TYPE * getData();
	EVENT_DATA_TYPE * getData(const string &, ColType);


};






#endif /* INST_INCLUDE_CYTOLIB_MEMCYTOFRAME_HPP_ */
