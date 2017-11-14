/*
 * cytoData.hpp
 *
 *  Created on: Sep 15, 2017
 *      Author: wjiang2
 */

#ifndef CYTOSET_HPP_
#define CYTOSET_HPP_


#include "CytoFrame.hpp"

/**
 * A container of the collection of CytoFrames
 */
class CytoSet{
	typedef unordered_map<string, CytoFrame> FRAMES;
	FRAMES frames;
	typedef FRAMES::iterator iterator;
	typedef FRAMES::const_iterator const_iterator;
public:
	/*
	 * forwarding APIs
	 */
	 size_t size(){return frames.size();}
	 iterator end(){return frames.end();}
	 iterator begin(){return frames.begin();}
	 iterator find(const string &sampleName){
	         return frames.find(sampleName);
	 }
	 size_type erase ( const string& k ){return frames.erase(k);}
	 /**
	  * insert
	  * @param sampleName
	  * @return
	  */
	 string & operator [](const string & sampleName){
	         return frames[sampleName];
	   }

	 /**
	  * forward to the first element's getChannels
	  */
	vector<string> getChannels(){return begin()->second.getChannels();};
	/**
	 * modify the channels for each individual frame
	 * @param _old
	 * @param _new
	 */
	void setChannels(const string & _old, const string & _new){
		for(auto & p : frames)
			p.second.setChannels(_old, _new);
	};

	//* forward to the first element's getChannels
	vector<string> getMarkers(){return begin()->second.getMarkers();};

	void setMarkers(const string & _old, const string & _new){
		for(auto & p : frames)
			p.second.setMarkers(_old, _new);
	};

	int nCol(){return begin()->second.nCol();}
	/**
	 * iterate through hash map to extract sample names
	 * @return
	 */
	vector<string> getSampleNames(){
		vector<string> res;
		for(const auto & f : frames)
			res.push_back(f.first);
		return res;

	};
	/**
	 * modify the name of one sample , which involves delete/insert the existing frame
	 * @param _old
	 * @param _new
	 */
	void setSampleName(const string & _old, const string & _new){
		auto it = find(_new);
		if(it!=end())
			throw(range_error(_new + " already exists!"));
		it = find(_old);
		if(it==end())
			throw(range_error(_old + " not found!"));

		auto frm = it->second;
		erase(_old);
		frames[_new] = frm;

	};

	/**
	 * validity checks on the frame to see if its data structure is consistent with cytoset
	 *
	 * @param frm
	 * @return
	 */
	int isNotValidFrame(const CytoFrame & frm){
		//validity check the channels against the existing frms
		if(nCol() != frm.nCol())
			return -1;

		//check channel in linear time(taking advantage of the hash map)
		auto frm1 = begin()->second;
		for(const auto & c : frm.getChannels())
		{
			if(frm1.getColId(c, ColType::channel) <0 )
				return -2;
		}

		//check the pdata
		const auto & pd1 = frm1.getPData();
		const auto & pd2 = frm.getPData();
		if(pd1.size()!=pd2.size())
			return -3;

		for(const auto & p : pd2)
		{
			if(pd1.find(p->first)==pd1.end())
				return -4;
		}
		return 0;
	}

	/**
	 * insert a new frame with some validity checks
	 * @param sampleName
	 * @param frm
	 */
	void addFrame(const string & sampleName, const CytoFrame & frm){
		auto it = find(sampleName);
		if(it!=end())
			throw(range_error(sampleName + " already exists!"));
		int flag = isNotValidFrame(frm);
		if(!flag)
		{
			frames[sampleName] = frm;
		}
		else
		{
			switch(flag)
			{
			case -1:
				throw(domain_error(sampleName + " the number of channels is not consistent with existing frames!"));
				break;
			case -2:
				throw(domain_error(sampleName + "  some channels are not present in existing frames!"));
				break;
			case -3:
				throw(domain_error(sampleName + " the number of pData columns is not consistent with existing frames!"));
				break;
			case -4:
				throw(domain_error(sampleName + " some columns of pData are not present in existing frames!"));
				break;
			default:
				throw(domain_error("invalid error code!"));
			}

		}

	}
	CytoSet subset(const vector<string> & sampleNames);
	vector<PDATA> getPData();
	void setPData(const pData & pd);
};




#endif /* CYTOSET_HPP_ */
