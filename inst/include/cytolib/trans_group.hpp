/* Copyright 2019 Fred Hutchinson Cancer Research Center
 * See the included LICENSE file for details on the license that is granted to the
 * user of this software.
 * trans_group.hpp
 *
 *  Created on: Mar 9, 2018
 *      Author: wjiang2
 */

#ifndef INST_INCLUDE_CYTOLIB_TRANS_GROUP_HPP_
#define INST_INCLUDE_CYTOLIB_TRANS_GROUP_HPP_

#include "transformation.hpp"

namespace cytolib
{

typedef map<string, TransPtr, ciLessBoost> trans_map;/* we always do case-insensitive searching for transformation lookup
due to some of channel name discrepancies occured in flowJo workspaces*/
struct PARAM{
		string param;
		bool log;
		unsigned range;
		unsigned highValue;
		unsigned calibrationIndex;
		//EDIT: can't trust this info from xml
//		EVENT_DATA_TYPE timestep;//only meaningful for time channel which is used to scale time channel (only for data, not for gates since gates are already stored at scaled value)
		PARAM(){};
		void update_channels(const CHANNEL_MAP & chnl_map);
		PARAM(const pb::PARAM & param_pb);
		 void convertToPb(pb::PARAM & param_pb);
		};
typedef vector<PARAM> PARAM_VEC;

/**
 *
 * @param pVec
 * @param name  channel name (possibly prefixed by comp.prefix)
 * @param comp
 * @return
 */
PARAM_VEC::const_iterator findTransFlag(const PARAM_VEC & pVec, const string & name
		, const string & prefix, const string & suffix);


class trans_local{
private:
	trans_map tp;
public:
//	trans_map getTransMap(){return tp;};
	const trans_map getTransMap()const{return tp;};
	void setTransMap(trans_map _tp){tp=_tp;};

	TransPtr getTran(string channel)const;

	trans_map cloneTransMap();

	void addTrans(string tName, TransPtr trans){tp[tName]=trans;};
	trans_local(){};

	void convertToPb(pb::trans_local & lg_pb);

	trans_local(const pb::trans_local & lg_pb);
	//legacy archive
	trans_local(const pb::trans_local & lg_pb, map<intptr_t, TransPtr> & trans_tbl);

	void update_channels(const CHANNEL_MAP & chnl_map);

	trans_local copy() const;
};

class trans_global:public trans_local{
private:
	string groupName;
	vector<int> sampleIDs;
public:
	void setSampleIDs(vector<int> _sampleIDs){sampleIDs=_sampleIDs;}
	vector<int> getSampleIDs()const{return sampleIDs;}
	string getGroupName()const {return groupName;}
	void setGroupName(string _groupName){groupName=_groupName;};
	trans_global(){};
};

typedef vector<trans_global> trans_global_vec;

};

#endif /* INST_INCLUDE_CYTOLIB_TRANS_GROUP_HPP_ */
