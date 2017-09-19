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

struct FCS_Header{
	float FCSversion;
	int textstart, textend, datastart, dataend, anastart, anaend, additional;
};

class MemCytoFrame{
	EVENT_DATA_TYPE * data;

public:
	~MemCytoFrame();
	MemCytoFrame():data(NULL){};
	MemCytoFrame(const string &filename, bool emptyValue, int nDataset, bool scale, double decades, double min_limit, bool truncate_max_range, bool onlyTxt);
	void compensate(const compensation &);
//	void transform(const transformation &);
	void save(const string & filename);
	EVENT_DATA_TYPE * getData();
	EVENT_DATA_TYPE * getData(const string &, ColType);
	KEY_WORDS getKeywords();
	string getKeyword(const string &);
	void setKeyword(const string &, const string &);
	int nCol();
	int nRow();
	vector<string> getChannels();
	vector<string> getMarkers();
	int getColId(const string &, ColType);
	void setChannel(const string &, const string &);
	void setMarker(const string &, const string &);

};






#endif /* INST_INCLUDE_CYTOLIB_MEMCYTOFRAME_HPP_ */
