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
#include "global.hpp"
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


	/*
	 * accessors
	 */
	vector<string> getParams(){return params;};
	void setEventCount(unsigned _nEvents){nEvents=_nEvents;};
	unsigned getEventsCount(){return nEvents;};
	//keep the redundant APIs to be consistent with cytoframe
	unsigned nRow(){return nEvents;};
	vector<string> getChannels(){return params;};

	void setSampleID(unsigned _sampleID){sampleID=_sampleID;};
	unsigned getSampleID(){return sampleID;};

	flowData():data(NULL),nEvents(0),ignore_case(false){};

	flowData(double* mat,vector<string> _params,unsigned _nEvents,unsigned _sampleID, bool _ignore_case = false){


		params=_params;

		nEvents=_nEvents;
		sampleID=_sampleID;

		data=mat;
		ignore_case = _ignore_case;
	}
	double * getData(){
		return data;
	}

	void setParams(vector<string> _params){
		if(_params.size()!=params.size())
			throw(domain_error("the number of parameters is not consistent with cdf file!"));
		params=_params;
	}
	struct InsensitiveCompare
	{
	  std::string comp;

	  InsensitiveCompare( std::string const &s ) : comp(s) {}

	  bool operator() ( std::string const &test ) const
	  {
	//	  string icomp = boost::to_lower_copy(comp);
	//	  string itest = boost::to_lower_copy(test);
	//	  return (icomp.compare(itest) == 0);
	    return (boost::iequals(comp, test));
	  }
	};
	unsigned find_pos(vector<string> s,string pattern, bool ignore_case) const{
		vector<string>::iterator it1,it2,res;
		it1=s.begin();
		it2=s.end();
		if(ignore_case)
		{
			res=find_if(it1,it2,InsensitiveCompare(pattern));
		}
		else
			res=find(it1,it2,pattern);
		if(res==it2)
			throw(domain_error(pattern.append(" not found in flowData!")));
		return (res-it1);
	}



	double * subset(string channel) const{
		unsigned paramInd=find_pos(params,channel, ignore_case);

				return data + paramInd*nEvents;
	}

};




#endif /* FLOWDATA_HPP_ */
