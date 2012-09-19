/*
 * Famus.hpp
 *
 *  Created on: Sep 19, 2012
 *      Author: yasir
 */

#ifndef FAMUS_HPP_
#define FAMUS_HPP_

#include <g2o/core/graph_optimizer_sparse.h>
#include <g2o/solvers/csparse/linear_solver_csparse.h>
#include <g2o/solvers/cholmod/linear_solver_cholmod.h>
#include <g2o/core/block_solver.h>
#include <g2o/types/slam2d/vertex_se2.h>
#include <g2o/types/slam2d/edge_se2.h>

#include "BFS.hpp"
#include "NeighbourSearch.hpp"


int expr = 0;
std::ofstream out("result.m");

//typedef std::pair<g2o::VertexSE2*, double> VertexCostPair;

using namespace g2o;

#define X_THRESH 1.50
#define Y_THRESH 1.50
#define THETA_THRESH 0.26

#define WRITE_MATLAB 0

typedef BlockSolver< BlockSolverTraits<-1, -1> >  SlamBlockSolver;
typedef LinearSolverCSparse<SlamBlockSolver::PoseMatrixType> SlamLinearSolver;

std::vector<int> leastUncertanityPath
(
		const std::string& filename,
		unsigned int start,
		unsigned int end
)
{
	SparseOptimizer optimizer;
	SlamLinearSolver * linearSolver = new SlamLinearSolver();
	linearSolver->setBlockOrdering(false);
	SlamBlockSolver * solver = new SlamBlockSolver(&optimizer,linearSolver);
	optimizer.setSolver(solver);

	if(!optimizer.load(filename.c_str()))
	{
		std::cerr<<" Could not load "<<filename<<std::endl;
		return std::vector<int>();
	}

	Graph g(optimizer.vertices().size());

	start = start%optimizer.vertices().size();
	end = end%optimizer.vertices().size();

	VertexIndexType indexMap;	indexMap = get(vertex_index, g);
	VertexCostType costMap;		costMap  = get(vertex_cost, g);


	optimizer.setVerbose(false);
	optimizer.findGauge()->setFixed(true);
	optimizer.initializeOptimization();
	optimizer.optimize(10);
	optimizer.solver()->computeMarginals();

	OptimizableGraph::VertexIDMap::const_iterator vIt = optimizer.vertices().begin(), vEnd = optimizer.vertices().end();
	OptimizableGraph::VertexIDMap::const_iterator vIt2 = optimizer.vertices().begin(); vIt2++;

	std::vector< double > vertexCostMap(optimizer.vertices().size(),0.);

	double logSum = 0 ; //, cost = 0;

	//std::cout<<"Calculating Cost .. "<<std::endl;

	clock_t now = clock();

	//std::ofstream out("result.m");

	if(WRITE_MATLAB and expr == 0 )
	{
		out<<" E_all {"<<expr+1<<"} = "<<optimizer.edges().size()<<";"<<std::endl;
		out<<" V_all {"<<expr+1<<"}= "<<optimizer.vertices().size()<<";"<<std::endl;

		out<<"cost {"<<expr+1<<"} = ["<<std::endl;
	}

	for (vIt = optimizer.vertices().begin() ; vIt != (vEnd) ; vIt++)
	{
		VertexSE2 * Vi = dynamic_cast<VertexSE2*>((vIt)->second);
		Eigen::MatrixXcd values = Vi->uncertainty().eigenvalues();
		logSum =  (log(values.real()(0,0))+log(values.real()(1,0))+log(values.real()(2,0)))/3.;
		vertexCostMap[Vi->id()]=exp(logSum);;
		if(WRITE_MATLAB and expr == 0)
		{
			out<<Vi->id()<<" "
					<<Vi->estimate().translation()[0]<<" "
					<<Vi->estimate().translation()[1]<<" "
					<<exp(logSum)<<" "<< std::setprecision(9)<<Vi->uncertainty().trace()<<std::endl;
		}
	}
	if(WRITE_MATLAB and expr == 0) out<<"];"<<std::endl;


	std::cerr<<(clock()-now)/1e3<<" ms"<<std::endl;
	//std::cout<<"Finding connected Nodes ... "<<std::endl;

	Points2D<double> pointCloud;

	pointCloud.pts.resize(optimizer.vertices().size());

	for(vIt = optimizer.vertices().begin(); vIt!=vEnd ; vIt++)
	{
		VertexSE2 * Vi = dynamic_cast<VertexSE2*>((vIt)->second);
		pointCloud.pts[Vi->id()].x = Vi->estimate().toVector()(0);
		pointCloud.pts[Vi->id()].y = Vi->estimate().toVector()(1);
	}

	my_kd_tree_t index(2, pointCloud, KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
	index.buildIndex();

	OptimizableGraph::EdgeSet::iterator eIt, eEnd;

	double query[2];

	for (vIt=optimizer.vertices().begin() ; vIt != (vEnd) ; vIt++)
	{
		VertexSE2 * Vi = dynamic_cast<VertexSE2*>((vIt)->second);
		query[0] = (double)Vi->estimate().toVector()(0);
		query[1] = (double)Vi->estimate().toVector()(1);

		//std::cout<<"Query : "<<Vi->id()<<std::endl;

		std::vector<std::pair<size_t,num_t> > ret_matches;

		nanoflann::SearchParams params;
		//params.sorted = false;

		double search_radius = 0.01*0.01;

		const size_t nMatches = index.radiusSearch(&query[0],search_radius, ret_matches, params);

		for (size_t i=1;i<nMatches;i++)
		{
			if((int)ret_matches[i].first< (int)Vi->id() and
					abs(static_cast<g2o::VertexSE2*>(optimizer.vertex(ret_matches[i].first))->estimate().toVector()(2) - Vi->estimate()[2])< THETA_THRESH and
					abs(ret_matches[i].first - Vi->id())>1
			)
			{
				//cout << ret_matches[i].first <<" ";//<< " dist["<< i << "]=" << ret_matches[i].second << endl;
				g2o::EdgeSE2* newEdge = new g2o::EdgeSE2;
				newEdge->vertices()[0] = Vi;
				newEdge->vertices()[1] = static_cast<g2o::VertexSE2*>(optimizer.vertex(ret_matches[i].first));

				newEdge->measurement().fromVector(Eigen::Vector3d(0,0,0));
				newEdge->setInformation(Eigen::Matrix3d::Identity());

				optimizer.addEdge(newEdge);
			}

		}
	}

	//std::cerr<<" Creating Boost Graph .. "<<std::endl;

	WeightMap weightMap; weightMap = get(edge_weight,g);

	for (eIt = optimizer.edges().begin() ; eIt!= optimizer.edges().end();
			eIt++)
	{
		g2o::EdgeSE2* _e = dynamic_cast<g2o::EdgeSE2*>(*eIt);

		Edge e ; bool inserted;
		boost::tie(e,inserted)= add_edge(_e->vertices()[0]->id(),_e->vertices()[1]->id(),g);
		weightMap[e] = vertexCostMap[_e->vertices()[1]->id()];

		boost::tie(e,inserted)= add_edge(_e->vertices()[1]->id(),_e->vertices()[0]->id(),g);
		weightMap[e] = vertexCostMap[_e->vertices()[0]->id()];

	}

	//std::cout<<"Reducing graph ..."<<std::endl;

	std::vector< std::pair<int,int> > points;
	std::vector<double> costs;

	for(vIt=optimizer.vertices().begin(); vIt!=vEnd; )
	{
		double cost=0; int members = 0;
		unsigned int gStart = vIt->second->id();

		double thisCost = 0;

		while(vIt->second->edges().size()==2 and !(vIt->second->id()==(int)start or vIt->second->id()==(int)end) and vIt!=vEnd)
		{

			//std::cout<<vIt->second->id()<<" ";
			cost+=vertexCostMap[vIt->second->id()];
			thisCost = vertexCostMap[vIt->second->id()];
			members++;
			vIt++;
		}

		if(members>0)
		{
			gStart = (gStart>0)? gStart-1:0;
			//std::cout<<" "<<gStart<<" -> "<<vIt->second->id()<<" : "<<cost<<std::endl;
			points.push_back(std::pair<int,int>(gStart,vIt->second->id()));
			costs.push_back(cost);
		}
		else
		{
			vIt++;
		}
	}

	//std::cout<<" Done with reducing!"<<std::endl;


	for(size_t i=0 ; i< points.size();i++)
	{
		eIt = optimizer.vertex(points[i].first)->edges().begin(); eEnd = optimizer.vertex(points[i].first)->edges().end();

		for( ; eIt!=eEnd ; eIt++)
		{
			if(points[i].second > points[i].first)
			{
				if((*eIt)->vertices()[0]->id() == points[i].first and (*eIt)->vertices()[1]->id() == points[i].first+1)
				{

					Edge e; bool exists;
					boost::tie(e,exists) = lookup_edge((*eIt)->vertices()[0]->id(),(*eIt)->vertices()[1]->id(),g);
					//std::cerr<<(*eIt)->vertices()[0]->id()<<" "<<(*eIt)->vertices()[1]->id()<<" "<<exists<<std::endl;

					optimizer.removeEdge(*eIt);
					if(exists) weightMap[e] = std::numeric_limits<double>::max();

				}

				if((*eIt)->vertices()[1]->id() == points[i].first and (*eIt)->vertices()[0]->id() ==  points[i].first+1)
				{
					Edge e; bool exists;
					tie(e,exists) = edge((*eIt)->vertices()[1]->id(),(*eIt)->vertices()[0]->id(),g);
					if(exists) weightMap[e] = std::numeric_limits<double>::max();
					//if(exists) remove_edge(e,g);
					optimizer.removeEdge(*eIt);
				}
			}
		}

		eIt = optimizer.vertex(points[i].second)->edges().begin(); eEnd = optimizer.vertex(points[i].second)->edges().end();

		for( ; eIt!=eEnd ; eIt++)
		{
			if(points[i].second > points[i].first) // this SHOULD always be true!
			{
				if((*eIt)->vertices()[0]->id() == points[i].second and (*eIt)->vertices()[1]->id() == points[i].second-1)
				{

					Edge e; bool exists;
					tie(e,exists) = edge((*eIt)->vertices()[0]->id(),(*eIt)->vertices()[1]->id(),g);
					if(exists) weightMap[e] = std::numeric_limits<double>::max();
					optimizer.removeEdge(*eIt);
				}
				if((*eIt)->vertices()[1]->id() == points[i].second and (*eIt)->vertices()[0]->id() == points[i].second-1)
				{
					Edge e; bool exists;
					tie(e,exists) = edge((*eIt)->vertices()[1]->id(),(*eIt)->vertices()[0]->id(),g);
					if(exists) weightMap[e] = std::numeric_limits<double>::max();
					optimizer.removeEdge(*eIt);
				}
			}
		}

		Edge e; bool added;

		tie(e,added) = add_edge(points[i].first,points[i].second,g);
		if(added) weightMap[e] = costs[i]+ vertexCostMap[points[i].second];

		tie(e,added) = add_edge(points[i].second,points[i].first,g);
		if(added) weightMap[e] = costs[i]+ vertexCostMap[points[i].first];

		g2o::EdgeSE2* newEdge = new EdgeSE2;

		newEdge->vertices()[0] = (VertexSE2*)optimizer.vertex(points[i].first);
		newEdge->vertices()[1] = (VertexSE2*)optimizer.vertex(points[i].second);

		optimizer.addEdge(newEdge);

	}


	for(size_t i=0 ; i< points.size();i++)
	{
		for(int j=points[i].first+1; j< points[i].second; j++)
		{
			optimizer.removeVertex(optimizer.vertex(j));
		}
	}


	Vertex s = vertex(start, g);

	std::vector<Vertex> p(num_vertices(g));
	std::vector<double> d(num_vertices(g),0.);

	d[start] = 0;

	dijkstra_shortest_paths(g, s, predecessor_map(&p[0]).distance_map(&d[0]));

	std::cout<<"Path Cost Minimum Uncertainty :"<< d[end]<<std::endl;

	if(WRITE_MATLAB) { out<<"PathCost{"<<expr+1<<"} = "<<d[end]<<";"<<std::endl; }

	std::set<std::pair<int,int> > pathSegments;

	pathSegments.insert(points.begin(), points.end());

	std::vector< int > Path;

	unsigned int currentnode = end;

	std::cout<<currentnode<<" ";

	if(WRITE_MATLAB)
	{
		out<<"minCostPath{"<<expr+1<<"} = [ "<<std::endl;
		out<<currentnode<<" ";
	}

	Path.push_back(currentnode);
	while(currentnode != start and p[currentnode]!=currentnode)
	{

		if(pathSegments.find(std::pair<int,int>(p[currentnode],currentnode))!=pathSegments.end() or
				pathSegments.find(std::pair<int,int>(currentnode,p[currentnode]))!=pathSegments.end())
		{
			if(currentnode > p[currentnode])
			{
				std::cout<<"[ ";
				for(size_t i=currentnode-1; i>p[currentnode];i--)
				{
					std::cout<<i<<" ";
					if(WRITE_MATLAB) out<<i<<" ";
					Path.push_back(i);
				}
				std::cout<<"] ";
			}
			else
			{
				std::cout<<"[ ";
				for(size_t i=currentnode+1; i<p[currentnode];i++)
				{
					std::cout<<i<<" ";
					if(WRITE_MATLAB) out<<i<<" ";
					Path.push_back(i);
				}
				std::cout<<"] ";
			}
		}
		std::cout<<p[currentnode]<<" ";
		out<<p[currentnode]<<" ";
		Path.push_back(p[currentnode]);
		currentnode = p[currentnode];

	}
	std::cout<<std::endl;

	if(WRITE_MATLAB )
	{
		out<<std::endl;
		out<<"];"<<std::endl;

		out<<" V_reduced  {"<<expr+1<<"}= "<< optimizer.vertices().size()<<";"<<std::endl;
		out<<" E_reduced  {"<<expr+1<<"}= "<<optimizer.edges().size()<<";"<<std::endl;
	}
	std::cout<<"# Decision Points "<<optimizer.vertices().size()<<std::endl;

	if(WRITE_MATLAB){ out<<"DecisionPoints {"<<expr+1<<"} = "<<optimizer.vertices().size()<<" ; "<<std::endl; }


	//optimizer.save("second.g2o");
	return Path;
}

std::vector<int> shortestPath
(
		const std::string& filename,
		unsigned int start,
		unsigned int end
)
{
	SparseOptimizer optimizer;
	SlamLinearSolver * linearSolver = new SlamLinearSolver();
	linearSolver->setBlockOrdering(false);
	SlamBlockSolver * solver = new SlamBlockSolver(&optimizer,linearSolver);
	optimizer.setSolver(solver);

	if(!optimizer.load(filename.c_str()))
	{
		std::cerr<<" Could not load "<<filename<<std::endl;
		return std::vector<int>();
	}

	Graph g(optimizer.vertices().size());

	start = start%optimizer.vertices().size();
	end = end%optimizer.vertices().size();

	VertexIndexType indexMap;	indexMap = get(vertex_index, g);
	VertexCostType costMap;		costMap = get(vertex_cost, g);


	optimizer.setVerbose(false);
	optimizer.vertex(0)->setFixed(true);
	optimizer.initializeOptimization();
	optimizer.optimize(10);

	OptimizableGraph::VertexIDMap::const_iterator vIt = optimizer.vertices().begin(), vEnd = optimizer.vertices().end();
	OptimizableGraph::VertexIDMap::const_iterator vIt2 = optimizer.vertices().begin(); vIt2++;

	optimizer.initializeOptimization();
	optimizer.optimize(1);

	optimizer.solver()->computeMarginals();

	std::vector< double > vertexCostMap(optimizer.vertices().size(),0.);

	//double logSum = 0 , cost = 0;

	//std::cout<<"Calculating Cost .. "<<std::endl;

	clock_t now = clock();

	//std::ofstream out("result.m");
	for (vIt = optimizer.vertices().begin() ; vIt != (vEnd) ; vIt++)
	{
		VertexSE2 * Vi = dynamic_cast<VertexSE2*>((vIt)->second);
		vertexCostMap[Vi->id()]=1.0;;
	}

	std::cerr<<(clock()-now)/1e3<<" ms"<<std::endl;
	//std::cout<<"Finding connected Nodes ... "<<std::endl;

	Points2D<double> pointCloud;

	pointCloud.pts.resize(optimizer.vertices().size());

	for(vIt = optimizer.vertices().begin(); vIt!=vEnd ; vIt++)
	{
		VertexSE2 * Vi = dynamic_cast<VertexSE2*>((vIt)->second);
		pointCloud.pts[Vi->id()].x = Vi->estimate().toVector()(0);
		pointCloud.pts[Vi->id()].y = Vi->estimate().toVector()(1);
	}

	my_kd_tree_t index(2, pointCloud, KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
	index.buildIndex();

	OptimizableGraph::EdgeSet::iterator eIt, eEnd;

	double query[2];

	for (vIt=optimizer.vertices().begin() ; vIt != (vEnd) ; vIt++)
	{
		VertexSE2 * Vi = dynamic_cast<VertexSE2*>((vIt)->second);
		query[0] = (double)Vi->estimate().toVector()(0);
		query[1] = (double)Vi->estimate().toVector()(1);


		std::vector<std::pair<size_t,num_t> > ret_matches;

		nanoflann::SearchParams params;

		double search_radius = 0.01*0.01;

		const size_t nMatches = index.radiusSearch(&query[0],search_radius, ret_matches, params);

		for (size_t i=1;i<nMatches;i++)
		{
			if((int)ret_matches[i].first< (int)Vi->id() and
					abs(static_cast<g2o::VertexSE2*>(optimizer.vertex(ret_matches[i].first))->estimate().toVector()(2) - Vi->estimate()[2])< THETA_THRESH and
					abs(ret_matches[i].first - Vi->id())>1
			)
			{
				//cout << ret_matches[i].first <<" ";//<< " dist["<< i << "]=" << ret_matches[i].second << endl;
				g2o::EdgeSE2* newEdge = new g2o::EdgeSE2;
				newEdge->vertices()[0] = Vi;
				newEdge->vertices()[1] = static_cast<g2o::VertexSE2*>(optimizer.vertex(ret_matches[i].first));

				newEdge->measurement().fromVector(Eigen::Vector3d(0,0,0));
				newEdge->setInformation(Eigen::Matrix3d::Identity());

				optimizer.addEdge(newEdge);
			}

		}
	}

	//std::cerr<<" Creating Boost Graph .. "<<std::endl;

	WeightMap weightMap; weightMap = get(edge_weight,g);

	for (eIt = optimizer.edges().begin() ; eIt!= optimizer.edges().end();
			eIt++)
	{
		g2o::EdgeSE2* _e = dynamic_cast<g2o::EdgeSE2*>(*eIt);

		//std::cerr<<_e->vertices()[0]->id()<<" "<<_e->vertices()[1]->id()<<std::endl;

		Edge e ; bool inserted;
		boost::tie(e,inserted)= add_edge(_e->vertices()[0]->id(),_e->vertices()[1]->id(),g);
		weightMap[e] = 1.0;

		boost::tie(e,inserted)= add_edge(_e->vertices()[1]->id(),_e->vertices()[0]->id(),g);
		weightMap[e] = 1.0;

	}


	Vertex s = vertex(start, g);

	std::vector<Vertex> p(num_vertices(g));
	std::vector<double> d(num_vertices(g),0.);

	d[start] = 0;

	dijkstra_shortest_paths(g, s, predecessor_map(&p[0]).distance_map(&d[0]));


	std::set<std::pair<int,int> > pathSegments;

	std::vector< int > path;

	unsigned int currentnode = end;
	std::cout<<currentnode<<" ";

	if(WRITE_MATLAB)
	{
		out<<"shortestPath{"<<expr+1<<"} = [ "<<std::endl;
		out<<currentnode<<" ";
	}

	path.push_back(currentnode);
	while(currentnode!= start and p[currentnode]!= currentnode)
	{
		std::cout<<p[currentnode]<<" ";
		if(WRITE_MATLAB) out<<p[currentnode]<<" ";
		path.push_back(p[currentnode]);
		currentnode = p[currentnode];

	}
	std::cout<<std::endl;

	if(WRITE_MATLAB)
	{
		out<<std::endl;
		out<<"];"<<std::endl;
	}
	return path;
}

#endif /* FAMUS_HPP_ */
