/*
 * populationNode.hpp
 *
 *  Created on: Mar 16, 2012
 *      Author: wjiang2
 */

#ifndef NODEPROPERTIES_HPP_
#define NODEPROPERTIES_HPP_

#include "POPINDICES.hpp"

using namespace std;


/*! population stats */
typedef map<string,float> POPSTATS;
/*
 *TODO: this class should exist apart from populationTree object
 *so all its constructor and desctuctor functions should be private
 */
/* gate is polymorphic member,so has to to be pointer
 * can't be reference either, because this member should belong to nodeProperties
 * and destroyed by nodeProperties's destructor ,reference member means refer to the object
 * outside, So it is not possible to instantiate it since we may not need to parse gate if only
 * stats from flowJo are needed.
 * Also reference member's life cycle is different from its host object, which could be problematic.
 *
 */

typedef boost::scoped_ptr<POPINDICES> popIndPtr;/*! the pointer to the event indices*/
/**
 * \class nodeProperties
 * \brief The container that holds gate and population information
 *
 * It has the population name, the pointer to the base \link<gate> object, gate indices and population stats.
 */
class nodeProperties{
private:
	string thisName;
	gate * thisGate; /**< pointer to the abstract gate object */
	popIndPtr indices;/**< scoped_ptr to the POPINDICES */
	POPSTATS fjStats,fcStats;
	bool hidden;


public:

	bool isGated(){return indices.get()!=NULL;};
	int getTotal(){return indices->getTotal();};

	nodeProperties():thisGate(NULL),hidden(false){}

	/*
	 * convert pb object to internal structure
	 * @param np_pb
	 */
	nodeProperties(const pb::nodeProperties & np_pb):thisGate(NULL),hidden(false){
		thisName = np_pb.thisname();
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
				thisGate = new rangeGate(gate_pb);
				break;
			case pb::BOOL_GATE:
				thisGate = new boolGate(gate_pb);
				break;
			case pb::POLYGON_GATE:
				thisGate = new polygonGate(gate_pb);
				break;
			case pb::RECT_GATE:
				thisGate = new rectGate(gate_pb);
				break;
			case pb::ELLIPSE_GATE:
				thisGate = new ellipseGate(gate_pb);
				break;
			case pb::ELLIPSOID_GATE:
				thisGate = new ellipsoidGate(gate_pb);
				break;
			case pb::LOGICAL_GATE:
				thisGate = new logicalGate(gate_pb);
				break;
			default:
				throw(domain_error("unknown type of gate archive!"));
			}

		}
	}


	void convertToPb(pb::nodeProperties & np_pb, bool isRoot){

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
	nodeProperties(const nodeProperties& np){
		thisName=np.thisName;

		thisGate=np.thisGate==NULL?NULL:np.thisGate->clone();
		if(np.indices.get()!=NULL)
			indices.reset(np.indices->clone());
		fjStats=np.fjStats;
		fcStats=np.fcStats;
		hidden=np.hidden;


	}
	nodeProperties & operator=(nodeProperties np){
		std::swap(thisName, np.thisName);
		std::swap(thisGate, np.thisGate);
		if(np.indices.get()!=NULL)
			indices.reset(np.indices->clone());
		std::swap(fjStats, np.fjStats);
		std::swap(fcStats, np.fcStats);
		std::swap(hidden, np.hidden);

		return *this;

	}

	/*
	 * gate is dynamically created,so they are freed here in destroy method
	 */
	~nodeProperties(){

	//	PRINT("entring the destructor of nodeProperties\n");

		if(thisGate!=NULL)
		{
			if(g_loglevel>=GATE_LEVEL)
				PRINT("free gate:"+this->thisName+"\n");
			delete thisGate;
		}
	}

	/**
	 * retrieve the pop stats that was pre-calculated and stored in the node.
	 *
	 * @param isFlowCore when true, return the calculated stats; if false, returns the stats parsed from xml (only relevant for the GatingHierarchy parsed from flowJo)
	  */
	POPSTATS getStats(bool isFlowCore=false){

		return(isFlowCore?this->fcStats:this->fjStats);
	}

	/**
	 * setter method for the private member of pop stats
	 * @param s POPSTATS
	 * @param isFlowCore flag indicates if the stats is for flowJo workspace.
	 */
	void setStats(POPSTATS s,bool isFlowCore=false){
		if(isFlowCore)
			fcStats=s;
		else
			fjStats=s;

	}
	/**
	 * getter for the private member of gate
	 * @return the pointer to an abstract base \link<gate> object
	 */
	gate * getGate(){
		if(thisGate==NULL)
			throw(logic_error("gate is not parsed!"));
		return(thisGate);
	}
	/**
	 * getter for the private member of population name
	 */
	string getName(){
		return(this->thisName);
	}
	/**
	 * setter for the private member of population name
	 */

	void setName(const char * popName){
		if(string(popName).find('/') != std::string::npos){
			throw(domain_error("pop name contains '/' character!"));
		}
		thisName=popName;
	}
	void setHiddenFlag(bool _value){
		hidden=_value;
	}
	bool getHiddenFlag(){
		return (hidden);
	}

	/**
	 * setter for the private member of gate
	 */
	void setGate(gate *gate){
		thisGate=gate;
	}

	/**
	 * Retrieve the event indices (relative to the whole ungated event data) from the node
	 * @return a bool vector
	 */
	vector<bool> getIndices(){
			if(!this->isGated())
				throw(domain_error("trying to get Indices for unGated node!"));
			return indices->getIndices();
			}
	vector<unsigned> getIndices_u(){
			if(!this->isGated())
				throw(domain_error("trying to get Indices for unGated node!"));
			return indices->getIndices_u();
			}

	void setIndices(unsigned _nEvent){
			indices.reset(new ROOTINDICES(_nEvent));
	}

	/**
	 * update the node with the new indices
	 *
	 */
	void setIndices(vector<bool> _ind){
		unsigned nEvents=count(_ind.begin(),_ind.end(),true);
		unsigned nSizeInt=sizeof(unsigned)*nEvents;
		unsigned nSizeBool=_ind.size()/8;

		if(nSizeInt<nSizeBool)
			indices.reset(new INTINDICES(_ind));
		else
			indices.reset(new BOOLINDICES(_ind));

	}

	void setIndices(INDICE_TYPE _ind, unsigned nTotal){
		unsigned nEvents=_ind.size();;
		unsigned nSizeInt=sizeof(unsigned)*nEvents;
		unsigned nSizeBool=nTotal/8;

		if(nSizeInt<nSizeBool)
			indices.reset(new INTINDICES(_ind, nTotal));
		else
			indices.reset(new BOOLINDICES(_ind, nTotal));

	}
	/*
	 * potentially it is step can be done within the same loop in gating
	 * TODO:MFI can be calculated here as well
	 */
	/**
	 * update the pop stats
	 *
	 * It is important to call this function after gate indices are updated.
	 */
	void computeStats(){
			fcStats["count"]=getCounts();
	}
	/**
	 * calculate the cell count for the current population.
	 *
	 * It may or may not be the same as the value retrieved from getStats especially when the gate indices are updated but computeStats has not been called.
	 */
	unsigned getCounts(){
		if(!this->isGated())
			throw(domain_error("trying to get counts for unGated node!"));
		return indices->getCount();

	}

};

#endif /* NODEPROPERTIES_HPP_ */
