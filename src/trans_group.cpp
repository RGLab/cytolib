// Copyright 2019 Fred Hutchinson Cancer Research Center
// See the included LICENSE file for details on the licence that is granted to the user of this software.
#include <cytolib/trans_group.hpp>
#include <cytolib/global.hpp>

namespace cytolib
{

		void PARAM::update_channels(const CHANNEL_MAP & chnl_map){
			CHANNEL_MAP::const_iterator itChnl = chnl_map.find(param);
			if(itChnl!=chnl_map.end())
				param = itChnl->second;
		};
		PARAM::PARAM(const pb::PARAM & param_pb){
			param = param_pb.param();
			log = param_pb.log();
			range = param_pb.range();
			highValue = param_pb.highvalue();
			calibrationIndex = param_pb.calibrationindex();
//			timestep = param_pb.timestep();
		};
		 void PARAM::convertToPb(pb::PARAM & param_pb){
			 param_pb.set_param(param);
			 param_pb.set_log(log);
			 param_pb.set_range(range);
			 param_pb.set_highvalue(highValue);
			 param_pb.set_calibrationindex(calibrationIndex);
//			 param_pb.set_timestep(timestep);
		 };
/**
 *
 * @param pVec
 * @param name  channel name (possibly prefixed by comp.prefix)
 * @param comp
 * @return
 */
 PARAM_VEC::const_iterator findTransFlag(const PARAM_VEC & pVec, const string & name, const string & prefix, const string & suffix){
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



	TransPtr trans_local::getTran(string channel)const{
		TransPtr res(NULL);
		trans_map::const_iterator it=tp.find(channel);
		if(it!=tp.end())
			res=it->second;
		return res;
	}

	trans_map trans_local::cloneTransMap(){

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


	void trans_local::convertToPb(pb::trans_local & lg_pb){

		for(const auto & it : tp){
			pb::trans_pair * tp = lg_pb.add_tp();
			tp->set_name(it.first);
			pb::transformation * trans_pb = tp->mutable_trans();

			it.second->convertToPb(*trans_pb);
		}
	}

	trans_local::trans_local(const pb::trans_local & lg_pb){

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
				case pb::PB_SCALE:
					tp[tp_pb.name()].reset(new scaleTrans(trans_pb));
					break;
				default:
					throw(domain_error("unknown type of transformation archive!"));
				}
			}
	}
	//legacy archive
	trans_local::trans_local(const pb::trans_local & lg_pb, map<intptr_t, TransPtr> & trans_tbl){

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

	void trans_local::update_channels(const CHANNEL_MAP & chnl_map){

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

	trans_local trans_local::copy() const
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

