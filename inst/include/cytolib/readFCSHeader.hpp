/*
 * readFCSHeader.hpp
 *
 *  Created on: Sep 21, 2017
 *      Author: wjiang2
 */

#ifndef INST_INCLUDE_CYTOLIB_READFCSHEADER_HPP_
#define INST_INCLUDE_CYTOLIB_READFCSHEADER_HPP_

#include <fstream>
#include <cstring>
#include "global.hpp"
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <vector>
#include <numeric>
#include <unordered_map>
#include <iostream>
#include <algorithm>
using namespace std;


enum class TransformType {none, linearize, scale,linearize_with_PnG_scaling};
enum class endianType {big, small, mixed};
struct FCS_Header{
	float FCSversion;
	int textstart, textend, datastart, dataend, anastart, anaend, additional;
};

typedef unordered_map <string, string> KEY_WORDS;
typedef float EVENT_DATA_TYPE;
typedef unique_ptr<EVENT_DATA_TYPE[] > EVENT_DATA_PTR;
struct FCS_READ_HEADER_PARAM{
	bool isEmptyKeyValue, ignoreTextOffset;
	int nDataset;
	FCS_READ_HEADER_PARAM(){
		isEmptyKeyValue = false;
		ignoreTextOffset = false;
		nDataset = 1;
	};

};
/**
 * parameters parsed from keywords, but may be accessed frequently later thus saved as the copy
 * in this struct for fast and easy query
 */
struct cytoParam{
	string channel, marker;
	EVENT_DATA_TYPE min, max, PnG;
	pair<EVENT_DATA_TYPE, EVENT_DATA_TYPE> PnE;
	int PnB;

};
void readHeaderAndText(ifstream &in,FCS_Header & header, KEY_WORDS & keys, vector<cytoParam> & params, const FCS_READ_HEADER_PARAM & config);


#endif /* INST_INCLUDE_CYTOLIB_READFCSHEADER_HPP_ */
