/*
 * CytoFrame.hpp
 *
 *  Created on: Sep 18, 2017
 *      Author: wjiang2
 */

#ifndef INST_INCLUDE_CYTOLIB_CYTOFRAME_HPP_
#define INST_INCLUDE_CYTOLIB_CYTOFRAME_HPP_


#include "readFCSHeader.hpp"
#include "compensation.hpp"
#include <c++/H5Cpp.h>
enum class ColType {channel, marker, unknown};
enum class RangeType {instrument, data};
enum class FrameType {FCS, H5};

using namespace H5;

const H5std_string  DATASET_NAME( "data");

/**
 * The class representing a single FCS file
 */
class CytoFrame{
protected:
//	int nEvents;
	KEY_WORDS keys;//keyword pairs parsed from FCS Text section
	vector<cytoParam> params;// parameters coerced from keywords and computed from data for quick query
	unordered_map<string, int> channel_vs_idx;//hash map for query by channel
	unordered_map<string, int> marker_vs_idx;//hash map for query by marker
public:
	virtual ~CytoFrame(){};
	virtual void compensate(const compensation &)=0;
//	virtual void transform(const transformation &)=0;
	virtual void writeFCS(const string & filename);
	/**
	 * save the CytoFrame as HDF5 format
	 *
	 * @param filename the path of the output H5 file
	 */
	virtual void writeH5(const string & filename);
	/**
	 * get the data of entire event matrix
	 * @return
	 */
	virtual EVENT_DATA_VEC getData()=0;
	/**
	 * get the data for the single channel
	 *
	 * @param colname the channel for marker name
	 * @param type enum class indicates the type of colname, can be either ColType::channel or ColType::marker or ColType::unknown
	 * when ColType::unknown, both types will be tried for the column match.
	 * @return
	 */
	virtual EVENT_DATA_VEC getData(const string & colname, ColType type)=0;
	/**
	 * extract all the keyword pairs
	 *
	 * @return a vector of pairs of strings
	 */
	virtual vector<pair <string, string>> getKeywords();
	/**
	 * extract the value of the single keyword by keyword name
	 *
	 * @param key keyword name
	 * @return keyword value as a string
	 */
	virtual string getKeyword(const string & key);
	/**
	 * set the value of the single keyword
	 * @param key keyword name
	 * @param value keyword value
	 */
	virtual void setKeyword(const string & key, const string & value);
	/**
	 * get the number of columns(or parameters)
	 *
	 * @return
	 */
	virtual int nCol();
	/**
	 * get the number of rows(or events)
	 * @return
	 */
	virtual int nRow()=0;
	/**
	 * check if the hash map for channel and marker has been built
	 * @return
	 */
	virtual bool isHashed();
	/**
	 * build the hash map for channel and marker for the faster query
	 *
	 */
	virtual void buildHash();
	/**
	 * get all the channel names
	 * @return
	 */
	virtual vector<string> getChannels();
	/**
	 * get all the marker names
	 * @return
	 */
	virtual vector<string> getMarkers();
	/**
	 * get the numeric index for the given column
	 * @param colname column name
	 * @param type the type of column
	 * @return
	 */
	virtual int getColId(const string & colname, ColType type);
	virtual void setChannel(const string &, const string &);
	virtual void setMarker(const string &, const string &);
	/**
	 * the range of a specific column
	 * @param colname
	 * @param ctype the type of column
	 * @param rtype either RangeType::data or RangeType::instrument
	 * @return
	 */
	virtual pair<EVENT_DATA_TYPE, EVENT_DATA_TYPE> getRange(const string & colname, ColType ctype, RangeType rtype);
};



#endif /* INST_INCLUDE_CYTOLIB_CYTOFRAME_HPP_ */
