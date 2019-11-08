// Copyright 2019 Fred Hutchinson Cancer Research Center
// See the included LICENSE file for details on the licence that is granted to the user of this software.
#include <cytolib/POPINDICES.hpp>
#include <boost/foreach.hpp>

namespace cytolib
{
void packToBytes(const vector <bool> & x, vector<unsigned char> & bytes){
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
void unpackFromBytes(vector <bool> & x, const vector<unsigned char>& x_bytes){

	/*
	 * unpack bytes into bits
	 */
	for(unsigned i =0 ; i < x.size(); i++) {
		unsigned byteIndex = i / 8;
		unsigned bitIndex = i % 8;

		x[i] = x_bytes[byteIndex] & (1 << bitIndex);
	}

}


	BOOLINDICES::BOOLINDICES(vector <unsigned> _ind, unsigned _nEvent){
		nEvents = _nEvent;
		x.resize(nEvents);
		for(auto i : _ind)
			x[i] = true;
	}
	BOOLINDICES::BOOLINDICES(vector <bool> _ind){
		x=_ind;
		nEvents=_ind.size();

	}
	vector<unsigned> BOOLINDICES::getIndices_u(){
		vector<unsigned> res;

		for(unsigned i = 0; i < x.size(); i++){
			if(x[i])
				res.push_back(i);
		}
		return res;
	}

	void BOOLINDICES::convertToPb(pb::POPINDICES & ind_pb){
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
	BOOLINDICES::BOOLINDICES(const pb::POPINDICES & ind_pb){
		nEvents = ind_pb.nevents();
		//fetch byte stream from pb
		vector<unsigned char> bytes(ind_pb.bind().begin(),ind_pb.bind().end());
		//convert it to bit vector
		x.resize(nEvents,false);
		unpackFromBytes(x, bytes);
	}

	INTINDICES::INTINDICES(vector <bool> _ind){

		for(vector<bool>::iterator it=_ind.begin();it!=_ind.end();it++)
		{
			unsigned i=it-_ind.begin();
			if(*it)
				x.push_back(i);
		}
		nEvents=_ind.size();
	}


	vector<bool> INTINDICES::getIndices(){

		vector<bool> res(nEvents,false);

		for(vector<unsigned>::iterator it=x.begin();it!=x.end();it++){
			unsigned i=*it;
			res[i]=true;
		}
		return res;
	}

	void INTINDICES::convertToPb(pb::POPINDICES & ind_pb){
		ind_pb.set_indtype(pb::INT);
		BOOST_FOREACH(vector<unsigned>::value_type & it, x){
			ind_pb.add_iind(it);
		}
		ind_pb.set_nevents(nEvents);
	}

	INTINDICES::INTINDICES(const pb::POPINDICES & ind_pb){
		nEvents = ind_pb.nevents();
		unsigned nSize = ind_pb.iind_size();
		x = vector<unsigned>(nSize);
		for(unsigned i = 0; i < nSize; i++)
			x[i] = ind_pb.iind(i);
	}


	vector<unsigned> ROOTINDICES::getIndices_u(){

		vector<unsigned> res(nEvents);
		for(unsigned i = 0; i < nEvents; i++)
			res[i]=i;
		return res;
	}


};

