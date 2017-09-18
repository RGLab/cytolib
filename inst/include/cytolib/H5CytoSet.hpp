/*
 * H5CytoSet.hpp
 *
 *  Created on: Sep 15, 2017
 *      Author: wjiang2
 */

#ifndef H5CYTOSET_HPP_
#define H5CYTOSET_HPP_
#include <cytolib/CytoSet.hpp>


class H5CytoSet:public CytoSet{
	string filename;
public:
	~H5CytoSet(){};
	H5CytoSet(const string & _filename):filename(_filename){};
	CytoFrame getFrame(const string & sampleName)=0;
	void setFrame(string sampleName, const CytoFrame & frm)=0;
	CytoSet subset(const vector<string> & sampleNames)=0;
	vector<string> getColnames()=0;
	void setColname(const string & _old, const string & _new)=0;
	vector<string> getMarkernames()=0;
	vector<string> getColnames()=0;
	pData getPhenoData()=0;
	void setPhenoData(const pData & pd)=0;
	vector<string> getSampleNames()=0;
	void setSampleName(const string & _old, const string & _new)=0;


};


#endif /* H5CYTOSET_HPP_ */
