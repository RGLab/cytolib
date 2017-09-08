/*
 * flowData.hpp
 *
 *  Created on: Apr 13, 2012
 *      Author: wjiang2
 */

#ifndef FLOWDATA_HPP_
#define FLOWDATA_HPP_
#include <vector>
#include <iostream>
#include <string>
#include <stdexcept>

using namespace std;

#include <boost/config.hpp>
#include <boost/graph/adj_list_serialize.hpp>
#include <boost/algorithm/string.hpp>

unsigned find_pos(vector<string> s,string pattern, bool ignore_case);
/*
 * representing one FCS data
 *
 * does not own the actual data, instead it holds the pointer to the external data array
 * so that there is zero overhead when passing data from R
 */
class flowData{

public:
	vector<string> params;
	unsigned sampleID;//it is only valid when access cdf version of flowdata, used as index for sample dimension
	double * data;
	unsigned nEvents;
	bool ignore_case; //whether channel-searching is case sensitive


	flowData();
	flowData(double* mat,vector<string>,unsigned _nEvents,unsigned _sampleID, bool _ignore_case = false);
	double * subset(string channel) const;
	/*
	 * accessors
	 */
	void setParams(vector<string> _params);
	vector<string> getParams(){return params;};
	void setEventCount(unsigned _nEvents){nEvents=_nEvents;};
	unsigned getEventsCount(){return nEvents;};
	void setSampleID(unsigned _sampleID){sampleID=_sampleID;};
	unsigned getSampleID(){return sampleID;};

	double * getData();
};




#endif /* FLOWDATA_HPP_ */
