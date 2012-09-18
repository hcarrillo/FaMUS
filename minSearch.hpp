/*
 * minSearch.hpp
 *
 *  Created on: Feb 13, 2012
 *      Author: yasir
 */

#ifndef MINSEARCH_HPP_
#define MINSEARCH_HPP_

#include <Eigen/Core>

#define EXPLORE

bool minNode(std::set<int>& Q, std::vector<double>& d, int& minIndex, double& minCost)
{
	minCost = std::numeric_limits<double>::max();
	minIndex = -1;
	std::set<int>::iterator it = Q.begin(), end = Q.end();
	for( ; it!=end ; it++)
	{
		if(d[*it]< minCost)
		{
			minCost = d[*it];
			minIndex = *it;
		}
	}
	if(minIndex >= 0) return true;
	return false;
}

std::vector<int> explore(const g2o::SparseOptimizer& optimizer,
		g2o::VertexSE2* start, int direction, int EndID)
{
	std::vector<int> exploredNodes;

	g2o::VertexSE2* nextHop = start;
	g2o::VertexSE2* thisVertex = start;
	while(nextHop->edges().size()==2)
	{
		g2o::OptimizableGraph::EdgeSet::iterator eIt = nextHop->edges().begin(), eEnd = nextHop->edges().end();
		for ( ; eIt!=eEnd ; eIt++)
		{
			if(direction==1 and (*eIt)->vertices()[direction]->id() > thisVertex->id())
			{
				//std::cout<<nextHop->id()<<" ";
				exploredNodes.push_back(nextHop->id());
				if( nextHop->id()==EndID)
				{
					std::cout<<" Here "<<std::endl;
					return exploredNodes;
				}
				nextHop = (g2o::VertexSE2*)(*eIt)->vertices()[direction];
				thisVertex = nextHop;
				//std::cin.get();
				//explore(optimizer,(g2o::VertexSE2*)(*eIt)->vertices()[direction],direction);
			}
			if(direction == 0 and (*eIt)->vertices()[direction]->id() < thisVertex->id())
			{
				//std::cout<<nextHop->id()<<" ";
				exploredNodes.push_back(nextHop->id());
				if( nextHop->id()==EndID)
				{
					return exploredNodes;
				}
				nextHop = (g2o::VertexSE2*)(*eIt)->vertices()[direction];
				thisVertex = nextHop;
				//std::cin.get();
			}
		}
	}
	exploredNodes.push_back(nextHop->id()); //Include the last node as well.
	return exploredNodes;
}

std::vector<int>
minPath(
		const g2o::SparseOptimizer& optimizer,
		std::map<g2o::VertexSE2*, double>& costMap,
		int start,
		int end
		)
{
	int m = optimizer.vertices().size();
	std::vector<double> d(m,std::numeric_limits<double>::max());
	std::set<int> Q;

	std::vector<int> v(m,-1);

	for(int i=0 ; i< m ; i++) Q.insert(i);

	d[start] = 0;
	v[start] = start;

	int minIndex; double minCost;

	while(Q.find(end)!=Q.end())
	{
		if(!minNode(Q,d,minIndex,minCost)) std::cerr<<"No min Node"<<std::endl;
		//else std::cerr<<" Min node @ "<<minIndex<<" with val : "<< minCost << std::endl;

		if(minIndex == end)
		{
			std::cout<<"Goal Reached"<<std::endl;
			break;
		}

		Q.erase(minIndex);


		g2o::VertexSE2* Vi = static_cast<g2o::VertexSE2*>(optimizer.vertices().find(minIndex)->second);
		g2o::OptimizableGraph::EdgeSet::iterator eIt = Vi->edges().begin(), eEnd = Vi->edges().end();

		std::cout<<"Min Node : "<<Vi->id()<<" ";

		//TODO keep on looking till we find a decision node!

		for ( ; eIt!=eEnd ; eIt++)
		{
			g2o::EdgeSE2* thisEdge = dynamic_cast<g2o::EdgeSE2*>(*eIt);
			g2o::VertexSE2* source = dynamic_cast<g2o::VertexSE2*>(thisEdge->vertices()[0]);
			g2o::VertexSE2* target = dynamic_cast<g2o::VertexSE2*>(thisEdge->vertices()[1]);

			//std::cout<<"[ "<<target->edges().size();
			//std::cout<<" "<<source->edges().size()<<" ]" ;

			if(source->id() == Vi->id()) // This Vertex is the source
			{

				int targetID = target->id();
				double myDistance = minCost + costMap[target];
				//if(myDistance < d[targetID])
				//{
				//	d[targetID] = myDistance;
				//	v[targetID] = Vi->id();
				//}
#ifdef EXPLORE
					std::vector<int> reachable = explore(optimizer,source,1,end); // reachable[0] is the current Node
					double minDist = minCost ;
					std::cout<<" Looking Ahead !"<<std::endl;
					for(size_t i=0 ; reachable.size() and i< reachable.size()-1 ; i++)
					{
						std::cout<<reachable[i]<<" ";
						minDist += costMap[(g2o::VertexSE2*)optimizer.vertex(reachable[i])];
						Q.erase(reachable[i]);
					}
					std::cout<<std::endl;
					if(reachable.size()>0)
					{
						if(d[reachable[reachable.size()-1]] > minDist){
							d[reachable[reachable.size()-1]] = minDist;
							v[reachable[reachable.size()-1]] = Vi->id();
						}
					}
#endif

			}

			if(target->id() == Vi->id()) // This Vertex is the Target, lookup the source
			{
				int targetID = source->id();
				//double myDistance = minCost + costMap[source];
				//if(myDistance < d[targetID])
				//{
				//	d[targetID] = myDistance;
				//	v[targetID] = Vi->id();
				//}
#ifdef EXPLORE
				std::vector<int> reachable = explore(optimizer,source,0,end); // reachable[0] is the current Node
				double minDist = minCost ;
				std::cout<<"LOok back ";
				for(size_t i=0 ; reachable.size() and i< reachable.size()-1; i++)
				{
					std::cout<<reachable[i]<<" ";
					minDist += costMap[(g2o::VertexSE2*)optimizer.vertex(reachable[i])];
					Q.erase(reachable[i]);
				}
				std::cout<<std::endl;

				if(reachable.size() > 0)
				{
					if(d[reachable[reachable.size()-1]] > minDist)
					{
						d[reachable[reachable.size()-1]] = minDist;
						v[reachable[reachable.size()-1]] = Vi->id();
					}
				}
#endif
			}


		}
	}

	for(size_t i = 0 ; i<d.size() ; i++)
	{
		if(d[i] == std::numeric_limits<double>::max()) std::cout <<"X ";
		else std::cout<<d[i]<<" ";
	}
	std::cout<<std::endl;


	for(size_t i = 0 ; i<d.size() ; i++)
	{
		if(v[i] == -1) std::cout <<"X ";
		else std::cout<<v[i]<<" ";
	}
	std::cout<<std::endl;

	std::cout<<"Total Cost "<< d[end]<<std::endl;

	std::vector<int> path;

	int currentNode = end;

	while(currentNode!=start)
	{
		path.push_back(currentNode);
		//std::cout<<currentNode<< " "<< v[currentNode]<< std::endl;
		currentNode = v[currentNode];
	}
	path.push_back(start);

	return path;



}

#endif /* MINSEARCH_HPP_ */
