/*
 * BFS.hpp
 *
 *  Created on: Feb 6, 2012
 *      Author: yasir
 */


#ifndef BFS_HPP_
#define BFS_HPP_

#include <boost/graph/adjacency_list.hpp>
//#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/lookup_edge.hpp>

#include <Eigen/Core>

namespace boost{
	enum vertex_transformation_t { vertex_transformation = 42 };
	enum edge_transformation_t { edge_transformation = 43 };
	enum edge_information_t { edge_information = 44 };
	enum vertex_cost_t { vertex_cost = 45 };

	BOOST_INSTALL_PROPERTY(vertex, transformation);
	BOOST_INSTALL_PROPERTY(edge, transformation);
	BOOST_INSTALL_PROPERTY(edge, information);
	BOOST_INSTALL_PROPERTY(vertex, cost);

};

using namespace boost;

typedef adjacency_list<
		vecS,
		vecS,
		bidirectionalS,
		property<vertex_cost_t, double>,
		property<edge_weight_t, double>,
		no_property,
		vecS
		>
		Graph;



typedef graph_traits<Graph>::vertex_iterator 		VertexIterator;
typedef graph_traits<Graph>::edge_iterator 			EdgeIterator;
typedef graph_traits<Graph>::vertex_descriptor		Vertex;
typedef graph_traits<Graph>::edge_descriptor		Edge;

typedef property_map<Graph, vertex_name_t>::type 			VertexNameType;
typedef property_map<Graph, vertex_index_t>::type 			VertexIndexType;
typedef property_map<Graph, vertex_transformation_t>::type 	VertexTransformationType;
typedef property_map<Graph, vertex_cost_t>::type 			VertexCostType;

typedef property_map<Graph, edge_name_t>::type   			EdgeNameType;
typedef property_map<Graph, edge_transformation_t>::type 	EdgeTransformationType;
typedef property_map<Graph, edge_information_t>::type 		EdgeInformationType;
typedef property_map<Graph, edge_weight_t>::type 			WeightMap;


#endif /* BFS_HPP_ */
