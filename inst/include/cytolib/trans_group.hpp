/*
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
		void update_channels(const CHANNEL_MAP & chnl_map){
			CHANNEL_MAP::const_iterator itChnl = chnl_map.find(param);
			if(itChnl!=chnl_map.end())
				param = itChnl->second;
		};
		PARAM(const pb::PARAM & param_pb){
			param = param_pb.param();
			log = param_pb.log();
			range = param_pb.range();
			highValue = param_pb.highvalue();
			calibrationIndex = param_pb.calibrationindex();
//			timestep = param_pb.timestep();
		};
		 void convertToPb(pb::PARAM & param_pb){
			 param_pb.set_param(param);
			 param_pb.set_log(log);
			 param_pb.set_range(range);
			 param_pb.set_highvalue(highValue);
			 param_pb.set_calibrationindex(calibrationIndex);
//			 param_pb.set_timestep(timestep);
		 };
		};
typedef vector<PARAM> PARAM_VEC;

/**
 *
 * @param pVec
 * @param name  channel name (possibly prefixed by comp.prefix)
 * @param comp
 * @return
 */
inline PARAM_VEC::const_iterator findTransFlag(const PARAM_VEC & pVec, const string & name, const string & prefix, const string & suffix){
	PARAM_VEC::const_iterator it;
	for(it=pVec.begin();it!=pVec.end();it++)
	{
		//	try both the bare and prefixed chnl ma,e
		string chnl = it->param;
		string chnl_comp = prefix + chnl + suffix;//append prefix
		if(chnl.compare(name)==0||chnl_comp.compare(name)==0)
			break;
	}

	return it;
}


class trans_local{
private:
	trans_map tp;
public:
//	trans_map getTransMap(){return tp;};
	const trans_map getTransMap()const{return tp;};
	void setTransMap(trans_map _tp){tp=_tp;};

	TransPtr getTran(string channel)const{
		TransPtr res(NULL);
		if(channel!="Time"&&channel!="time")
		{
			trans_map::const_iterator it=tp.find(channel);
			if(it!=tp.end())
				res=it->second;
		}
		return res;
	}

	trans_map cloneTransMap(){

		trans_map res;
		/*
		 * clone trans map
		 */

		for(trans_map::iterator it=tp.begin();it!=tp.end();it++)
		{
			TransPtr curTran=it->second;
			if(curTran!=NULL)
			{
				if(g_loglevel>=POPULATION_LEVEL)
					PRINT("cloning transformatioin:"+curTran->getChannel()+"\n");
				res[it->first] = curTran->clone();
			}
		}
		return res;
	}

	void addTrans(string tName, TransPtr trans){tp[tName]=trans;};
	trans_local(){};

	void convertToPb(pb::trans_local & lg_pb){

		for(const auto & it : tp){
			pb::trans_pair * tp = lg_pb.add_tp();
			tp->set_name(it.first);
			pb::transformation * trans_pb = tp->mutable_trans();

			it.second->convertToPb(*trans_pb);
		}
	}

	trans_local(const pb::trans_local & lg_pb){

		for(int i = 0; i < lg_pb.tp_size(); i ++)
		{
			const pb::trans_pair & tp_pb = lg_pb.tp(i);
				const pb::transformation & trans_pb = tp_pb.trans();
				switch(trans_pb.trans_type())
				{
				case pb::PB_CALTBL:
					tp[tp_pb.name()].reset(new transformation(trans_pb));
					break;
				case pb::PB_BIEXP:
					tp[tp_pb.name()].reset(new biexpTrans(trans_pb));
					break;
				case pb::PB_FASIGNH:
					tp[tp_pb.name()].reset(new fasinhTrans(trans_pb));
					break;
				case pb::PB_FLIN:
					tp[tp_pb.name()].reset(new flinTrans(trans_pb));
					break;
				case pb::PB_LIN:
					tp[tp_pb.name()].reset(new linTrans(trans_pb));
					break;
				case pb::PB_LOG:
					tp[tp_pb.name()].reset(new logTrans(trans_pb));
					break;
				case pb::PB_LOGICLE:
					tp[tp_pb.name()].reset(new logicleTrans(trans_pb));
					break;
				default:
					throw(domain_error("unknown type of transformation archive!"));
				}
			}
	}
	//legacy archive
	trans_local(const pb::trans_local & lg_pb, map<intptr_t, TransPtr> & trans_tbl){

		for(int i = 0; i < lg_pb.tp_size(); i ++){
			const pb::trans_pair & tp_pb = lg_pb.tp(i);
			intptr_t old_address = (intptr_t)tp_pb.trans_address();
			//look up from the tbl for the new pointer
			map<intptr_t, TransPtr>::iterator it = trans_tbl.find(old_address);
			if(it!=trans_tbl.end()){
				tp[tp_pb.name()] = it->second;
			}
			else
				throw(domain_error("the current archived transformation is not found in the global table!"));

		}
	}

	void update_channels(const CHANNEL_MAP & chnl_map){

		//iterate throiugh chnl_map instead of tp since tp iterator will be invalidated when erased
		for(CHANNEL_MAP::const_iterator itChnl = chnl_map.begin(); itChnl != chnl_map.end(); itChnl++)
		{

			string oldN = itChnl->first;
			trans_map::iterator itTp = tp.find(oldN);

			if(itTp!=tp.end())
			{
				string newN = itChnl->second;
				if(g_loglevel>=GATING_SET_LEVEL)
					PRINT("update transformation: "+ oldN + "-->" + newN +"\n");

				TransPtr curTran = itTp->second;
				curTran->setChannel(newN);
				/*
				 *
				 * have to delete the old one before adding newN
				 * because tp[newN] could be just updating existing tp[oldN]
				 * instead of inserting new entry due to the fact tp trans_map is set to case insensitive key searching
				 * otherwise, tp.erase would lead to losing the entry when oldN and newN are equivalent
				 */
				tp.erase(itTp);
				tp[newN] = curTran; //add new entry

			}

		}
	}

	trans_local copy() const
	{
		trans_local res;
		for(const auto & it : tp)
		{
			auto trans = it.second;
			string chnl = it.first;
			if(g_loglevel>=GATING_HIERARCHY_LEVEL)
			{
				string type;
				if(!trans)
					throw(domain_error("Empty trans func: " + chnl));

				trans->getType(type);
				PRINT("copy transformation type = " + type +" : " + chnl+"\n");
			}
			res.tp[chnl] = trans->clone();
		}
		return res;
	}
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
