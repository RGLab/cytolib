/*
 * H5CytoFrame.hpp
 *
 *  Created on: Sep 18, 2017
 *      Author: wjiang2
 */

#ifndef INST_INCLUDE_CYTOLIB_H5CYTOFRAME_HPP_
#define INST_INCLUDE_CYTOLIB_H5CYTOFRAME_HPP_
#include <cytolib/CytoFrame.hpp>



class H5CytoFrame{
	~H5CytoFrame(){};
	void compensate(const compensation &);
//	void transform(const transformation &);
	H5CytoFrame(const string & filename, FrameType type);
	void save(const string & filename, FrameType type);
	double * getData();
	double * getData(const string &, ColType);
	KEY_WORDS getKeywords();
	string getKeyword(const string &);
	void setKeyword(const string &, const string &);
	int nCol();
	int nRow();
	vector<string> getChannels();
	vector<string> getMarkers();
	void setChannel(const string &, const string &);
	void setMarker(const string &, const string &);

};




#endif /* INST_INCLUDE_CYTOLIB_H5CYTOFRAME_HPP_ */
