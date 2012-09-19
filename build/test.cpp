/*
 * test.cpp
 *
 *  Created on: Feb 6, 2012
 *      Author: yasir
 */
#include <iostream>
#include "Famus.hpp"


int main(int argc, char** argv)
{

	if(argc != 4)
	{
		std::cout<<"Usage "<<argv[0]<<" g2o_filename "<<" startVertexID endVertexID "<<std::endl;
		return 1;
	}

	int start = atoi(argv[2]), end = atoi(argv[3]);

	std::vector<int> path = leastUncertanityPath(argv[1],start,end);
	std::vector<int> path2 = shortestPath(argv[1],start,end);


	// Percentage Overlap between the two //

	std::set< int > minUncertanityPath(path.begin(), path.end());
	double matches = 0;
	for(std::vector<int>::iterator it = path2.begin(), end = path2.end(); it!=end ; it++)
	{
		if(minUncertanityPath.find(*it)!=minUncertanityPath.end())
		{
			matches = matches + 1;
		}
	}

	std::cout<<" % overlap: "<<(matches-2.0)/(double)(path2.size()-2)<<std::endl;

	return 0;

}
