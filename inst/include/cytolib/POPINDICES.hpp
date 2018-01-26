/*
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
inline void packToBytes(const vector <bool> & x, vector<unsigned char> & bytes){
	/*
	 * pack bits into bytes
	 */

	for(unsigned i =0 ; i < x.size(); i++) {
		unsigned byteIndex = i / 8;
		unsigned bitIndex = i % 8;
		if(x[i]) {
			// set bit
			unsigned char mask  = 1 << bitIndex;
			bytes[byteIndex] = bytes[byteIndex] | mask;
		}
	}

}
inline void unpackFromBytes(vector <bool> & x, const vector<unsigned char>& x_bytes){

	/*
	 * unpack bytes into bits
	 */
	for(unsigned i =0 ; i < x.size(); i++) {
		unsigned byteIndex = i / 8;
		unsigned bitIndex = i % 8;

		x[i] = x_bytes[byteIndex] & (1 << bitIndex);
	}

}

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

	BOOLINDICES(vector <unsigned> _ind, unsigned _nEvent){
		nEvents = _nEvent;
		x.resize(nEvents);
		for(auto i : _ind)
			x[i] = true;
	}
	BOOLINDICES(vector <bool> _ind){
		x=_ind;
		nEvents=_ind.size();

	}
	vector<bool> getIndices(){
		return x;
	}
	vector<unsigned> getIndices_u(){
		vector<unsigned> res;

		for(unsigned i = 0; i < x.size(); i++){
			if(x[i])
				res.push_back(i);
		}
		return res;
	}


	unsigned getCount(){
		return count(x.begin(),x.end(),true);
	}




	POPINDICES * clone(){
		BOOLINDICES * res=new BOOLINDICES(*this);
		return res;
	}
	void convertToPb(pb::POPINDICES & ind_pb){
		ind_pb.set_indtype(pb::BOOL);
		unsigned nBits=x.size();
		unsigned nBytes=ceil(float(nBits)/8);
		vector<unsigned char> bytes(nBytes,0);
		packToBytes(x, bytes);
		string * byte_pb = ind_pb.mutable_bind();
		for(unsigned i = 0; i < bytes.size(); i++){
			unsigned char byte = bytes[i];
			byte_pb->append(string(1, byte));
		}
		ind_pb.set_nevents(nEvents);
	}
	BOOLINDICES(const pb::POPINDICES & ind_pb){
		nEvents = ind_pb.nevents();
		//fetch byte stream from pb
		vector<unsigned char> bytes(ind_pb.bind().begin(),ind_pb.bind().end());
		//convert it to bit vector
		x.resize(nEvents,false);
		unpackFromBytes(x, bytes);
	}

};

/*
 * int vector
 */
class INTINDICES:public POPINDICES{
private:
	vector <unsigned> x;
public:
	INTINDICES():POPINDICES(){};

	INTINDICES(vector <bool> _ind){

		for(vector<bool>::iterator it=_ind.begin();it!=_ind.end();it++)
		{
			unsigned i=it-_ind.begin();
			if(*it)
				x.push_back(i);
		}
		nEvents=_ind.size();
	}

	INTINDICES(vector <unsigned> _ind, unsigned _nEvent):POPINDICES(_nEvent),x(_ind){};

	vector<bool> getIndices(){

		vector<bool> res(nEvents,false);

		for(vector<unsigned>::iterator it=x.begin();it!=x.end();it++){
			unsigned i=*it;
			res[i]=true;
		}
		return res;
	}

	vector<unsigned> getIndices_u(){return x;};
	unsigned getCount(){

		return x.size();
	}


	POPINDICES * clone(){

		INTINDICES * res=new INTINDICES(*this);
		return res;
	}
	void convertToPb(pb::POPINDICES & ind_pb){
		ind_pb.set_indtype(pb::INT);
		BOOST_FOREACH(vector<unsigned>::value_type & it, x){
			ind_pb.add_iind(it);
		}
		ind_pb.set_nevents(nEvents);
	}

	INTINDICES(const pb::POPINDICES & ind_pb){
		nEvents = ind_pb.nevents();
		unsigned nSize = ind_pb.iind_size();
		x = vector<unsigned>(nSize);
		for(unsigned i = 0; i < nSize; i++)
			x[i] = ind_pb.iind(i);
	}


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
	vector<unsigned> getIndices_u(){

		vector<unsigned> res(nEvents);
		for(unsigned i = 0; i < nEvents; i++)
			res[i]=i;
		return res;
	}


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
