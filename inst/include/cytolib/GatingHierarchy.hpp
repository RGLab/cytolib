/*Copyright 2019 Fred Hutchinson Cancer Research Center
 * See the included LICENSE file for details on the license that is granted to the
 * user of this software.
 * GatingHierarchy.hpp
 *
 *  Created on: Mar 17, 2012
 *      Author: mike
 */

#ifndef GATINGHIERARCHY_HPP_
#define GATINGHIERARCHY_HPP_
#include <iostream>
#include <string>
#include <vector>
#include <queue>
#include <unordered_set>
#include "populationTree.hpp"
#include <fstream>
#include <algorithm>
#include "MemCytoFrame.hpp"
#include "CytoFrameView.hpp"
#include "TileCytoFrame.hpp"
#include "H5CytoFrame.hpp"
using namespace std;

namespace cytolib
{
#define REGULAR 0
#define TSORT 1
#define BFS 2
/*
 * because __SIZE_TYPE__ is long long unsigned int by gcc on win64 (mingw64)
 * we cast it to unsigned int before pass it to Rcpp::wrap to avoid error
 */
typedef unsigned int NODEID;


typedef map<string,VertexID> VertexID_map;
typedef vector<VertexID> VertexID_vec;
typedef vector<pair<VertexID,VertexID>> VertexID_edge_vec; //represents edges as <start, end> pairs
typedef vector<string> StringVec;
typedef vector<EVENT_DATA_TYPE> DoubleVec;
typedef vector<bool> BoolVec;
typedef vector<NODEID> NODEID_vec;

struct phylo{
	VertexID_edge_vec edges;
	VertexID_vec leaf_nodes; // aka "tips"
	vector<string> leaf_names;
	VertexID_vec internal_nodes; // all non-leaf nodes
	vector<string> internal_names;
};


struct OurVertexPropertyWriterR {

	OurVertexPropertyWriterR(populationTree &g_) : g(g_) {}

    template <class Vertex>
    void operator() (std::ostream &out, Vertex u) {
    	nodeProperties &curNode=g[u];
    	bool isBoolGate=false;
    	bool hidden = false;
    	if(u!=0)
    	{
    		unsigned short gateType=curNode.getGate()->getType();
    		isBoolGate=(gateType==BOOLGATE);
    		hidden=curNode.getHiddenFlag();
    	}
    	out<<"[shape=record,label=\""<<curNode.getName()<<"\",isBool="<<isBoolGate<<",hidden="<<hidden<<"]";


    }

    populationTree &g;
};

class GatingHierarchy;
typedef shared_ptr<GatingHierarchy> GatingHierarchyPtr;

CytoFramePtr load_cytoframe(const string & uri, bool readonly = true
			, CytoCtx ctxptr = CytoCtx());
/**
 ** \class GatingHierarchy
 **
 ** \brief The class that holds the gating tree.
 **
 ** It stores the transformation, compensation and gating information as well as the flow data
	Each FCS file is associated with one GatingHierarchy object.
	It can also serves as a gating template when data is empty.
 */
class GatingHierarchy{
private:
	compensation comp; /*< compensation object */

	/* compensation is currently done in R due to the linear Algebra
						e[, cols] <- t(solve(t(spillover))%*%t(e[,cols]))
						we can try uBlas for this simple task, but when cid=="-1",we still need to
						do this in R since comp is extracted from FCS keyword (unless it can be optionally extracted from workspace keyword)
	 	 	 	 	  */
	populationTree tree; /**< the gating tree */

	PARAM_VEC transFlag; /*< for internal use of parse flowJo workspace */
	trans_local trans; /*< the transformation used for this particular GatingHierarchy object */
	CytoFrameView frame_;
public:
	bool is_cytoFrame_only() const{return tree.m_vertices.size()==1;};
	CytoFrameView & get_cytoframe_view_ref(){return frame_;}
	CytoFrameView get_cytoframe_view() const{return frame_;}
	void set_cytoframe_view(CytoFrameView fr){
		frame_ = fr;
	}
	/*
	 * forwarding cytoFrameView APIs
	 */
	unsigned n_rows() const{return get_cytoframe_view().n_rows();}
	int n_cols() const{return get_cytoframe_view().n_cols();}
	vector<string> get_markers() const{return get_cytoframe_view().get_markers();};
	vector<string> get_channels() const{return get_cytoframe_view().get_channels();};
	void set_marker(const string & _channel, const string & _marker){get_cytoframe_view_ref().set_marker(_channel, _marker);}
	void set_channel(const string & _old, const string & _new){get_cytoframe_view_ref().set_channel(_old, _new);}
	const PDATA & get_pheno_data() const {return get_cytoframe_view().get_pheno_data();}


	/**
	 * setter for channels (dedicated fro Rcpp API and  it won't throw on the unmatched old channel name)
	 * @param chnl_map
	 */
	void set_channels(const CHANNEL_MAP & chnl_map);
	/**
	 * compensate the data by the spillover provided by workspace
	 * or FCS TEXT keyword
	 * add prefix (e.g. Comp_ or <>) to channel name of the data
	 */
	void compensate(CytoFrame & cytoframe);
	trans_local getLocalTrans() const{return trans;}

	/**
	 * transform the data
	 * The reason we pass in MemCytoFrame is because the data member frame_ may not be finalized yet at this stage of parsing.
	 */
	void transform_data(MemCytoFrame & cytoframe);

	void calgate(MemCytoFrame & cytoframe, VertexID u, bool computeTerminalBool, INTINDICES &parentIndice);
	void extendGate(MemCytoFrame & cytoframe, float extend_val);

	/**
	 * default constructor that creates an empty gating tree
	 *
	 * examples:
	 * \code
	 * 	GatingHierarchy *curGh=new GatingHierarchy();
	 * \endcode
	 */
	GatingHierarchy(){
		addRoot();//add default root to avoid the risk of crashing on accessing the empty boost graph
		}
	GatingHierarchy(const CytoFrameView & frame_view):frame_(frame_view){addRoot();};
	GatingHierarchy(compensation _comp, PARAM_VEC _transFlag, trans_local _trans):comp(_comp), transFlag(_transFlag),trans(_trans) {addRoot();};
	/**
	 *
	 * @param gh_pb
	 * @param uri
	 * @param h5_opt
	 * @param is_skip_data whether to skip writing cytoframe data to pb. It is typically remain as default unless for debug purpose (e.g. re-writing gs that is loaded from legacy pb archive without actual data associated)
	 */
	void convertToPb(pb::GatingHierarchy & gh_pb, string uri, CytoFileOption h5_opt
			, bool is_skip_data = false
			, const CytoCtx & ctx = CytoCtx());
	GatingHierarchy(CytoCtx ctx, pb::GatingHierarchy & pb_gh, string uri, bool is_skip_data
			, bool readonly = true){
			const pb::populationTree & tree_pb =  pb_gh.tree();
			int nNodes = tree_pb.node_size();

			tree = populationTree(nNodes);
			for(int i = 0; i < nNodes; i++){
				const pb::treeNodes & node_pb = tree_pb.node(i);
				const pb::nodeProperties & np_pb = node_pb.node();

				VertexID curChildID = i;
				tree[curChildID] = nodeProperties(np_pb);

				if(node_pb.has_parent()){
					VertexID parentID = node_pb.parent();
					boost::add_edge(parentID,curChildID,tree);
				}

			}
			//restore comp
			comp = compensation(pb_gh.comp());
			//restore trans flag
			for(int i = 0; i < pb_gh.transflag_size(); i++){
				transFlag.push_back(PARAM(pb_gh.transflag(i)));
			}
			//restore trans local
			trans = trans_local(pb_gh.trans());

			//restore fr
			if(!is_skip_data)
			{
				CytoFramePtr ptr = load_cytoframe(uri, readonly, ctx);
				pb::CytoFrame fr = *pb_gh.mutable_frame();
				if(!fr.is_h5())
					ptr.reset(new MemCytoFrame(*ptr));
				frame_ = CytoFrameView(ptr);
			}
		}

	//load legacy pb
	GatingHierarchy(pb::GatingHierarchy & pb_gh, map<intptr_t, TransPtr> & trans_tbl);
	/**
	 * add empty root node to the gating tree with the name set to 'root'
	 *
	 * \return the newly added root node Id
	 *
	 * For example:
	 * \code
	 * 	GatingHierarchy *curGh=new GatingHierarchy();
	 *	curGh->addRoot();
	 * \endcode
	 *
	 */
	VertexID addRoot();
	/*
	 * this is for semi-automated pipeline to add population node associated with gate
	 * assuming gate split the parent population into two subpops, one of which is to keep
	 * depends on isNegate flag of the gate
	 */
	VertexID addGate(gatePtr g,VertexID parentID,string popName);
	/*
	 * remove the node along with associated population properities including indices and gates
	 * can't do recursive removal here since vid is subject to change during the removal precoess
	 */
	void removeNode(VertexID nodeID)
	{
		if(nodeID>0)
		{
			//remove edge associated with this node
			EdgeID e=getInEdges(nodeID);
			/*removing vertex cause the rearrange node index
			 * so make sure do it after get edge descriptor
			 */
			boost::remove_edge(e,tree);
		}

		boost::remove_vertex(nodeID,tree);
	}
	/**
	 * recursive version
	 * @param node
	 */
	void removeNode(string node)
	{
		auto nodeID = getNodeID(node);

		//recursive to its children first

		//get children name vec first before id is expired
		auto childrenIDs = getChildren(nodeID);
		auto n = childrenIDs.size();
		if(n>0)
		{
			vector<string> children(n);
			for(unsigned i = 0 ; i < n; i++)
			{
				children[i] = getNodePath(childrenIDs[i]);
			}
			for(auto child : children)
			{
				removeNode(child);
			}
		}

		//actually delete the node
		removeNode(nodeID);
	}

	/**
	 *
	 * It moves one node to the target parent.
	 *
	 * @param parent the target parent id
	 * @param child node id to be moved
	 */
	void moveNode(string node, string parent);

	/*
	 * Getter function for compensation member
	 * @return
	 */
	compensation get_compensation(){
		return comp;
	}
	void set_compensation(const compensation & _comp, bool is_update_prefix);
	void set_compensation(compensation && _comp, bool is_update_prefix);
	void printLocalTrans();
	/*
	 * the version without the need of loading data
	 * by supplying the extend_to value
	 */
	void extendGate(float extend_val, float extend_to);
	/*
	 * adjust gates by gains
	 */
	void adjustGate(map<string,float> &gains);

	/*
	 * transform gates
	 */
	void transform_gate();

	/*
	 * Apply post-transformation shifts to gates (like for magnetic gates)
	 */
	void shift_gate();

	void check_ungated_bool_node(VertexID u);
	/*
	 * traverse the tree to gate each pops
	 * assuming data have already been compensated and transformed
	 *
	 */
	void gating(MemCytoFrame & cytoframe, VertexID u,bool recompute=false, bool computeTerminalBool=true, bool skip_faulty_node = false);
	void gating(MemCytoFrame & cytoframe, VertexID u,bool recompute
			, bool computeTerminalBool, bool skip_faulty_node, INTINDICES &parentIndice);
	/*
	 * bool gating operates on the indices of reference nodes
	 * because they are global, thus needs to be combined with parent indices
	 * in cases of negated gate (i.e. !A & !B)
	 * @param u
	 * @return
	 */

	vector<bool> boolGating(MemCytoFrame & cytoframe, VertexID u, bool computeTerminalBool);
	/*
	 * external boolOpSpec can be provided .
	 * It is mainly used by openCyto rectRef gate
	 * (needs to be combined with parent indices)
	 *
	 * @param u
	 * @param boolOpSpec
	 * @return
	 */
	vector<bool> boolGating(MemCytoFrame & cytoframe, vector<BOOL_GATE_OP> boolOpSpec, bool computeTerminalBool);

	/*
	 * current output the graph in dot format
	 * and further covert it to gxl in order for Rgraphviz to read since it does not support dot directly
	 * right now the data exchange is through file system,it would be nice to do it in memory
	 */
	void drawGraph(string output);


	/**
	 * retrieve all the node IDs
	 *
	 * @param order accept 3 values: REGULAR(0) is the same original order by which nodes were added;
	 * 								 TSORT(1) topological order;
	 * 								 BFS(2) breadth first searching order
	 * @return a vector of node IDs
	 */


	VertexID_vec getVertices(unsigned short order = 0);
	/*
	 *  Unary predicate for checking whether a string is empty
	 * @param path
	 * @return
	 */
	static bool isEmpty(string path){
		return(path.empty());

	}
	/**
	 * retrieve the VertexID by the gating path
	 * @param gatePath single string containing full(or partial) gating path
	 *
	 * For example:
	 * \code
	 * gh->getNodeID("singlet");
	 * gh->getNodeID("CD3/CD4+");
	 * \endcode
	  */
	VertexID getNodeID(string gatePath);
	/*
	 * retrieve the VertexID by the gating path.
	 * this serves as a parser to convert generic gating path into internal node ID
	 * and it doesn't allow ambiguity (throw the exception when multiple nodes match)
	 * @param gatePath a string vector of full(or partial) gating path
	 * @return node id
	 */
	VertexID getNodeID(const deque<string> & gatePath);
	/*
	 *  find the most immediate common ancestor
	 *
	 * @param nodeIDs input node IDs
	 * @param nDepths the depths of ancestor node. It is used to as the measurement to determine which reference node to win when multiple matches occur.
	 * 											   The deeper it is, the nearer the reference is to the boolean node.
	 * @return ancestor ID
	 */
	VertexID getCommonAncestor(VertexID_vec nodeIDs, unsigned & nDepths);
	VertexID getRefNodeID(VertexID u, const deque<string> & refPath);

	VertexID_vec pathMatch(VertexID_vec leafIDs, const deque<string> & gatePath);
	/*
	 * retrieve the VertexIDs by the gating path.(bottom-up searching)
	 * This routine allows multiple matches
	 * @param ancestorID when gatePath is partial path, this node ID narrow the searching range.
	 * @param gatePath input
	 * @return node IDs that matches to the query path
	 */
	VertexID_vec queryByPath(VertexID ancestorID, const deque<string> & gatePath);

	/**
	 * check if v is the descendant of u
	 * @param u
	 * @param v
	 * @return
	 */
	bool isDescendant(VertexID u, VertexID v);

	/**
	 * search for all the nodes that matches the pop name given the ancestor node id
	 * @param u the ancestor node id to search from
	 * @param name the node name to search for
	 * @return the vector of node id that match
	 *
	 * For example:
	 * \code
	 * VertexID parentID = gh->getNodeID("CD3");
	 * //this may return two descendants: "CD3/CD4/CCR7+ 45RA+" and "CD3/CD8/CCR7+ 45RA+"
	 * gh->getDescendants(parentID, "CCR7+ 45RA+");
	 * \endcode

	 */
	VertexID_vec getDescendants(VertexID u,string name);



	/*
	 * retrieve VertexID that matches population name given an ancestor node
	 * It is used to search for the first node in the gating path (full or partial).
	 * This is different from getRefNodeID in the way that pop name must be uniquely identifiable in the tree.
	 * @param u the ancestor node id
	 * @param popName the population name to match
	 * @return node ID
	 */
	VertexID getDescendant(VertexID u,string popName);

	/**
	 * retrieve all the population paths
	 *
	 * The assumption is each node only has one parent.
	 *
	 * @param order passed to getVertices function
	 * @param fullPath flag indicates whether to return full path or partial path
	 * @param showHidden whether to include the hidden nodes
	 * @return
	 */
	vector<string> getNodePaths(unsigned short order,bool fullPath,bool showHidden);

	/**
	 * Compute the depth of the given node
	 *
	 * @param u node ID
	 */
	unsigned getNodeDepths(VertexID u);
	/**
	 * Convert node Id to abs path
	 * @param u
	 * @return
	 */
	string getNodePath(VertexID u,bool fullPath = true);
	/**
	 * Get ancestor node for the given node
	 *
	 * Assume getParent only returns one parent node per GatingHierarchy
	 *
	 * @param u the given node ID
	 * @param level specify the distance from the given node
	 */
	VertexID getAncestor(VertexID u,unsigned short level);
	/*
	 * using boost in_edges out_edges to retrieve adjacent vertices
	 * assuming only one parent for each node
	 */
	EdgeID getInEdges(VertexID target);

	/**
	 * Get parent node id for the given node
	 *
	 * @param target child ID
	 */
	VertexID getParent(VertexID target);
	/**
	 * retrieve all children nodes
	 *
	 * @param source parent node ID
	 */
	VertexID_vec getChildren(VertexID source);

	/*
	 * retrieve single child node by parent id and child name.
	 * @param source id of the source node
	 * @param childName the child node name
	 * @return the child node id if succeeds; otherwise return -1.
	 */
	int getChildren(VertexID source,string childName);

	/*
	 *
	 * make sure to use this API always since since it is safe way to access tree nodes due to the validity check
	 *
	 *since the vertex bundle should always exist as long as the  tree and node exist, thus it is safe
	 * to return the reference of it
	 */
	/**
	 * Retrieve the node properties
	 *
	 * It is the only way to access the gate, population indices and stats of the given node
	 * @param u node ID
	 * @return a reference to the nodeProperties object
	 */
	nodeProperties & getNodeProperty(VertexID u);
	populationTree & getTree(){return tree;};

	/*
	 * Retrieve a representation of the tree as a vector of directed edges represented as pairs of VertexIDs
	 * and a vector of VertexIDs of the leaf nodes (aka "tips").
	 *
	 * @param start Node ID of starting node to allow sub-graph extraction
	 */
	phylo getPhylo(VertexID start, bool fullPath = true);


	GatingHierarchyPtr  copy(bool is_copy_data, bool is_realize_data, const string & uri) const;
	/*
	 * It is mainly used by Rcpp API addTrans to propagate global trans map to each sample
	 * EDIT: But now also used by clone methods
	  * @param trans trans_map
	 */
	void addTransMap(trans_map tm){
		trans.setTransMap(tm);

	}
};
};

#endif /* GATINGHIERARCHY_HPP_ */
