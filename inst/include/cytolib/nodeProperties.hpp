/* Copyright 2019 Fred Hutchinson Cancer Research Center
 * See the included LICENSE file for details on the license that is granted to the
 * user of this software.
 * populationNode.hpp
 *
 *  Created on: Mar 16, 2012
 *      Author: wjiang2
 */

#ifndef NODEPROPERTIES_HPP_
#define NODEPROPERTIES_HPP_

#include "POPINDICES.hpp"

using namespace std;

namespace cytolib
{

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

typedef unique_ptr<POPINDICES> popIndPtr;/*! the pointer to the event indices*/
/**
 * \class nodeProperties
 * \brief The container that holds gate and population information
 *
 * It has the population name, the pointer to the base \link<gate> object, gate indices and population stats.
 */
class nodeProperties{
private:
	string thisName;
	gatePtr thisGate; /**< pointer to the abstract gate object */
	popIndPtr indices;/**< ptr to the POPINDICES */
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
	nodeProperties(const pb::nodeProperties & np_pb);


	void convertToPb(pb::nodeProperties & np_pb, bool isRoot);
	/* since nodeProperties contains non-copyable scope_ptr member
	 * , customized copy and assignment constructor is required
	 *
	 */
	nodeProperties(const nodeProperties& np);
	nodeProperties & operator=(nodeProperties np);


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
	void setStats(POPSTATS s,bool isFlowCore=false);
	/**
	 * getter for the private member of gate
	 * @return the pointer to an abstract base \link<gate> object
	 */
	gatePtr getGate();
	/**
	 * getter for the private member of population name
	 */
	string getName(){
		return(this->thisName);
	}
	/**
	 * setter for the private member of population name
	 */

	void setName(const char * popName);
	void setHiddenFlag(bool _value){
		hidden=_value;
	}
	bool getHiddenFlag(){
		return (hidden);
	}

	/**
	 * setter for the private member of gate
	 */
	void setGate(gatePtr gate){
		thisGate=gate;
	}

	/**
	 * Retrieve the event indices (relative to the whole ungated event data) from the node
	 * @return a bool vector
	 */
	vector<bool> getIndices();
	vector<unsigned> getIndices_u();

	void setIndices(unsigned _nEvent){
			indices.reset(new ROOTINDICES(_nEvent));
	}

	/**
	 * update the node with the new indices
	 *
	 */
	void setIndices(vector<bool> _ind);

	void setIndices(INDICE_TYPE _ind, unsigned nTotal);
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
	unsigned getCounts();

};
};

#endif /* NODEPROPERTIES_HPP_ */
