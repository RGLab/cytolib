// Copyright 2019 Fred Hutchinson Cancer Research Center
// See the included LICENSE file for details on the licence that is granted to the user of this software.
#include <cytolib/nodeProperties.hpp>
#include <cytolib/global.hpp>
#include <boost/foreach.hpp>

namespace cytolib
{

	/*
	 * convert pb object to internal structure
	 * @param np_pb
	 */
	nodeProperties::nodeProperties(const pb::nodeProperties & np_pb):thisGate(NULL),hidden(false){
		thisName = np_pb.thisname();
		if(g_loglevel>=POPULATION_LEVEL)
				PRINT("loading node: "+thisName+"\n");;

		hidden = np_pb.hidden();
		for(int i = 0; i < np_pb.fcstats_size(); i++){
		   const pb::POPSTATS &	stat_pb = np_pb.fcstats(i);
		   fcStats[stat_pb.stattype()] = stat_pb.statval();
		}
		for(int i = 0; i < np_pb.fjstats_size(); i++){
		   const pb::POPSTATS & stat_pb = np_pb.fjstats(i);
		   fjStats[stat_pb.stattype()] = stat_pb.statval();
		}
		if(np_pb.has_indices()){
			const pb::POPINDICES & ind_pb = np_pb.indices();
			switch(ind_pb.indtype()){
			case pb::BOOL:
				indices.reset(new BOOLINDICES(ind_pb));
				break;
			case pb::INT:
				indices.reset(new INTINDICES(ind_pb));
				break;
			case pb::ROOT:
				indices.reset(new ROOTINDICES(ind_pb));
				break;
			default:
				throw(domain_error("unknown type of event indices archive!"));
			}
		}

		/*
		 * parse gate
		 */
		if(np_pb.has_thisgate()){
			const pb::gate & gate_pb = np_pb.thisgate();
			switch(gate_pb.type())
			{
			case pb::RANGE_GATE:
				thisGate = gatePtr(new rangeGate(gate_pb));
				break;
			case pb::BOOL_GATE:
				thisGate = gatePtr(new boolGate(gate_pb));
				break;
			case pb::POLYGON_GATE:
				thisGate = gatePtr(new polygonGate(gate_pb));
				break;
			case pb::RECT_GATE:
				thisGate = gatePtr(new rectGate(gate_pb));
				break;
			case pb::ELLIPSE_GATE:
				thisGate = gatePtr(new ellipseGate(gate_pb));
				break;
			case pb::ELLIPSOID_GATE:
				thisGate = gatePtr(new ellipsoidGate(gate_pb));
				break;
			case pb::LOGICAL_GATE:
				thisGate = gatePtr(new logicalGate(gate_pb));
				break;
			case pb::CLUSTER_GATE:
				thisGate = gatePtr(new clusterGate(gate_pb));
				break;
			case pb::QUAD_GATE:
				thisGate = gatePtr(new quadGate(gate_pb));
				break;
			default:
				throw(domain_error("unknown type of gate archive!"));
			}

		}
	}


	void nodeProperties::convertToPb(pb::nodeProperties & np_pb, bool isRoot){

		np_pb.set_thisname(thisName);
		np_pb.set_hidden(hidden);
		//copy fj stats
		BOOST_FOREACH(POPSTATS::value_type & it, fjStats){
			pb::POPSTATS *fj = np_pb.add_fjstats();
			fj->set_stattype(it.first);
			fj->set_statval(it.second);
		}
		//copy fc stats
		BOOST_FOREACH(POPSTATS::value_type & it, fcStats){
			pb::POPSTATS *fc = np_pb.add_fcstats();
			fc->set_stattype(it.first);
			fc->set_statval(it.second);
		}

		bool isGated = indices != NULL;

		if(isRoot){
			if(isGated){
				pb::POPINDICES * ind_pb = np_pb.mutable_indices();
				indices->convertToPb(*ind_pb);
			}

		}
		else
		{
			//cp gate
			if(thisGate!=NULL){
				pb::gate * gate_pb = np_pb.mutable_thisgate();
				thisGate->convertToPb(*gate_pb);

				//only archive event index for gated non-bool gate to save disk)
				bool nonBool = thisGate->getType() != BOOLGATE;

				if(isGated&&nonBool){
					pb::POPINDICES * ind_pb = np_pb.mutable_indices();
					indices->convertToPb(*ind_pb);
				}

			}

		}


	}

	/* since nodeProperties contains non-copyable scope_ptr member
	 * , customized copy and assignment constructor is required
	 *
	 */
	nodeProperties::nodeProperties(const nodeProperties& np){
		thisName=np.thisName;

		thisGate=np.thisGate==NULL?NULL:np.thisGate->clone();
		if(np.indices.get()!=NULL)
			indices.reset(np.indices->clone());
		fjStats=np.fjStats;
		fcStats=np.fcStats;
		hidden=np.hidden;


	}
	nodeProperties & nodeProperties::operator=(nodeProperties np){
		std::swap(thisName, np.thisName);
		std::swap(thisGate, np.thisGate);
		if(np.indices.get()!=NULL)
			indices.reset(np.indices->clone());
		std::swap(fjStats, np.fjStats);
		std::swap(fcStats, np.fcStats);
		std::swap(hidden, np.hidden);

		return *this;

	}



	/**
	 * setter method for the private member of pop stats
	 * @param s POPSTATS
	 * @param isFlowCore flag indicates if the stats is for flowJo workspace.
	 */
	void nodeProperties::setStats(POPSTATS s,bool isFlowCore){
		if(isFlowCore)
			fcStats=s;
		else
			fjStats=s;

	}
	/**
	 * getter for the private member of gate
	 * @return the pointer to an abstract base \link<gate> object
	 */
	gatePtr nodeProperties::getGate(){
		if(thisGate==NULL)
			throw(logic_error("gate is not parsed!"));
		return(thisGate);
	}
	/**
	 * setter for the private member of population name
	 */

	void nodeProperties::setName(const char * popName){
		if(string(popName).find('/') != std::string::npos){
			throw(domain_error("pop name contains '/' character!"));
		}
		thisName=popName;
	}


	/**
	 * Retrieve the event indices (relative to the whole ungated event data) from the node
	 * @return a bool vector
	 */
	vector<bool> nodeProperties::getIndices(){
			if(!this->isGated())
				throw(domain_error("trying to get Indices for unGated node!"));
			return indices->getIndices();
			}
	vector<unsigned> nodeProperties::getIndices_u(){
			if(!this->isGated())
				throw(domain_error("trying to get Indices for unGated node!"));
			return indices->getIndices_u();
			}

	/**
	 * update the node with the new indices
	 *
	 */
	void nodeProperties::setIndices(vector<bool> _ind){
		unsigned nEvents=count(_ind.begin(),_ind.end(),true);
		unsigned nSizeInt=sizeof(unsigned)*nEvents;
		unsigned nSizeBool=_ind.size()/8;

		if(nSizeInt<nSizeBool)
			indices.reset(new INTINDICES(_ind));
		else
			indices.reset(new BOOLINDICES(_ind));

	}

	void nodeProperties::setIndices(INDICE_TYPE _ind, unsigned nTotal){
		unsigned nEvents=_ind.size();;
		unsigned nSizeInt=sizeof(unsigned)*nEvents;
		unsigned nSizeBool=nTotal/8;

		if(nSizeInt<nSizeBool)
			indices.reset(new INTINDICES(_ind, nTotal));
		else
			indices.reset(new BOOLINDICES(_ind, nTotal));

	}
	/**
	 * calculate the cell count for the current population.
	 *
	 * It may or may not be the same as the value retrieved from getStats especially when the gate indices are updated but computeStats has not been called.
	 */
	unsigned nodeProperties::getCounts(){
		if(!this->isGated())
			throw(domain_error("trying to get counts for unGated node!"));
		return indices->getCount();

	}

};

