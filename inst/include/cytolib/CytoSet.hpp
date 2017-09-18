/*
 * cytoData.hpp
 *
 *  Created on: Sep 15, 2017
 *      Author: wjiang2
 */

#ifndef CYTOSET_HPP_
#define CYTOSET_HPP_


#include "CytoFrame.hpp"



class pData{
	vector<vector<string>> tbl;
};


/**
 * A container of the collection of CytoFrames
 */
class CytoSet{

public:
	virtual ~CytoSet()=0;
	virtual CytoFrame getFrame(const string & sampleName)=0;
	virtual void setFrame(string sampleName, const CytoFrame & frm)=0;
	virtual CytoSet subset(const vector<string> & sampleNames)=0;
	virtual vector<string> getColnames()=0;
	virtual void setColname(const string & _old, const string & _new)=0;
	virtual vector<string> getMarkernames()=0;
	virtual vector<string> getColnames()=0;
	virtual pData getPhenoData()=0;
	virtual void setPhenoData(const pData & pd)=0;
	virtual vector<string> getSampleNames()=0;
	virtual void setSampleName(const string & _old, const string & _new)=0;

};


#endif /* CYTOSET_HPP_ */
