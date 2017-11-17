/*
 * flowData.cpp
 *
 *  Created on: Apr 13, 2012
 *      Author: wjiang2
 */

#include <cytolib/flowData.hpp>

#include <algorithm>


flowData::flowData():data(NULL),nEvents(0),ignore_case(false){};

flowData::flowData(double* mat,vector<string> _params,unsigned _nEvents,unsigned _sampleID, bool _ignore_case){


	params=_params;

	nEvents=_nEvents;
	sampleID=_sampleID;

	data=mat;
	ignore_case = _ignore_case;
}
double * flowData::getData(){
	return data;
}

void flowData::setParams(vector<string> _params){
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
unsigned find_pos(vector<string> s,string pattern, bool ignore_case){
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



double * flowData::subset(string channel) const{
	unsigned paramInd=find_pos(params,channel, ignore_case);

			return data + paramInd*nEvents;
}



