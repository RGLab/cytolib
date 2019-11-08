/* Copyright 2019 Fred Hutchinson Cancer Research Center
 * See the included LICENSE file for details on the license that is granted to the
 * user of this software.
 * POPINDICES.hpp
 *
 *  Created on: Oct 8, 2012
 *      Author: wjiang2
 */


#ifndef POPINDICES_HPP_
#define POPINDICES_HPP_

#include "gate.hpp"

namespace cytolib
{
void packToBytes(const vector <bool> & x, vector<unsigned char> & bytes);
void unpackFromBytes(vector <bool> & x, const vector<unsigned char>& x_bytes);

/**
 * \class POPINDICES
 * \brief the event indices for the subpopulation
 *
 *  It is an abstract class that hides the details of implementation of the actual indices type.
 */
class POPINDICES{
protected:
	unsigned nEvents;
public:
	POPINDICES():nEvents(0){};
	POPINDICES(unsigned _nEvents):nEvents(_nEvents){};
	virtual ~POPINDICES(){};
	/**
	 * convert the POPINDICES to bool vector
	 *
	 */
	virtual vector<bool> getIndices()=0;
	virtual vector<unsigned> getIndices_u()=0;
	/**
	 * compute the event count from the event indices
	 */
	virtual unsigned getCount()=0;
	/**
	 * retrieve the total number of events for the original ungated events data
	 */
	unsigned getTotal(){return nEvents;}
	virtual POPINDICES * clone()=0;
	virtual void convertToPb(pb::POPINDICES & ind_pb) = 0;

};
/*
 * bool vector
 */
class BOOLINDICES:public POPINDICES{
private:
	vector <bool> x;
public:
	BOOLINDICES():POPINDICES(){};

	BOOLINDICES(vector <unsigned> _ind, unsigned _nEvent);
	BOOLINDICES(vector <bool> _ind);
	vector<bool> getIndices(){
		return x;
	}
	vector<unsigned> getIndices_u();


	unsigned getCount(){
		return count(x.begin(),x.end(),true);
	}




	POPINDICES * clone(){
		BOOLINDICES * res=new BOOLINDICES(*this);
		return res;
	}
	void convertToPb(pb::POPINDICES & ind_pb);
	BOOLINDICES(const pb::POPINDICES & ind_pb);

};

/*
 * int vector
 */
class INTINDICES:public POPINDICES{
private:
	vector <unsigned> x;
public:
	INTINDICES():POPINDICES(){};

	INTINDICES(vector <bool> _ind);

	INTINDICES(vector <unsigned> _ind, unsigned _nEvent):POPINDICES(_nEvent),x(_ind){};

	vector<bool> getIndices();

	vector<unsigned> getIndices_u(){return x;};
	unsigned getCount(){

		return x.size();
	}


	POPINDICES * clone(){

		INTINDICES * res=new INTINDICES(*this);
		return res;
	}
	void convertToPb(pb::POPINDICES & ind_pb);

	INTINDICES(const pb::POPINDICES & ind_pb);


};

/*
 * root node
 */
class ROOTINDICES:public POPINDICES{

public:
	ROOTINDICES():POPINDICES(){};
	ROOTINDICES(unsigned _nEvents):POPINDICES(_nEvents){};

	vector<bool> getIndices(){

		vector<bool> res(nEvents,true);

		return res;
	}
	vector<unsigned> getIndices_u();


	unsigned getCount(){

		return nEvents;
	}

	POPINDICES * clone(){

		ROOTINDICES * res=new ROOTINDICES(*this);
		return res;
	}
	void convertToPb(pb::POPINDICES & ind_pb){
		ind_pb.set_indtype(pb::ROOT);
		ind_pb.set_nevents(nEvents);
	}
	ROOTINDICES(const pb::POPINDICES & ind_pb){
		nEvents = ind_pb.nevents();
	}
};

};

#endif /* POPINDICES_HPP_ */
