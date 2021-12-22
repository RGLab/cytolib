/* Copyright 2019 Fred Hutchinson Cancer Research Center
 * See the included LICENSE file for details on the license that is granted to
 * the user of this software. tree.hpp
 *
 *  Created on: Mar 17, 2012
 *      Author: mike
 */

#ifndef TREE_HPP_
#define TREE_HPP_

#include "nodeProperties.hpp"
#include <boost/graph/adjacency_list.hpp>

namespace cytolib {
#define ROOTNODE 0

struct Edge {
  // nothing, probably. Or a weight, a distance, a direction, ...
};
typedef boost::adjacency_list< // adjacency_list is a template depending on :
    boost::vecS, //  The container used for egdes : here, std::list.
    boost::vecS, //  The container used for vertices: here, std::vector.
    boost::bidirectionalS, //  directed or undirected edges ?.
    nodeProperties *, Edge>
    populationTreeOld;

/*since we don't use pointer here
 * and has customized copy and assignment constructor defined for nodeProperties
 * class the entire graph adjacency_list is copiable now thus eliminate the need
 * for customized clone member functions
 */
/** populationTree represented as a boost graph (adjacency list) */
typedef boost::adjacency_list< // adjacency_list is a template depending on :
    boost::vecS, //  The container used for egdes : here, std::list.
    boost::vecS, //  The container used for vertices: here, std::vector.
    boost::bidirectionalS, //  directed or undirected edges ?.
    nodeProperties, boost::no_property>
    populationTree;
typedef populationTree::vertex_descriptor
    VertexID; //! ID for the node in the gating tree
typedef populationTree::vertex_iterator VertexIt;
typedef populationTree::edge_descriptor EdgeID;
typedef populationTree::edge_iterator EdgeIt;
}; // namespace cytolib

#endif /* TREE_HPP_ */
