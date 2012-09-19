/*
 * NeighbourSearch.hpp
 *
 *  Created on: Feb 20, 2012
 *      Author: yasir
 */

#ifndef NEIGHBOURSEARCH_HPP_
#define NEIGHBOURSEARCH_HPP_

#include <iostream>
#include <nanoflann.hpp>

using namespace nanoflann;
using std::cout;
using std::endl;

template<typename T>
struct Points2D
{
	struct Point
	{
		T  x,y;
	};

	std::vector<Point> pts;

	// Must return the number of data points
	inline size_t kdtree_get_point_count() const { return pts.size(); }

	// Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
	inline T kdtree_distance(const T *p1, const size_t idx_p2,size_t size) const
	{
		const T d0=p1[0]-pts[idx_p2].x;
		const T d1=p1[1]-pts[idx_p2].y;
		return d0*d0+d1*d1;
	}

	// Returns the dim'th component of the idx'th point in the class:
	// Since this is inlined and the "dim" argument is typically an immediate value, the
	//  "if/else's" are actually solved at compile time.
	inline T kdtree_get_pt(const size_t idx, int dim) const
	{
		if (dim==0) return pts[idx].x;
		else return pts[idx].y;
	}

	// Optional bounding-box computation: return false to default to a standard bbox computation loop.
	//   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
	//   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
	template <class BBOX>
	bool kdtree_get_bbox(BBOX &bb) const { return false; }

};

typedef double num_t;
typedef KDTreeSingleIndexAdaptor<
			L2_Simple_Adaptor<num_t, Points2D<num_t> > ,
			Points2D<num_t>,
			2 /* dim */
			> my_kd_tree_t;

#if 0
template<typename T, typename estimateType> // T : container Type , estimateType : type of estimate to cast to.
class NearestNeighbour2D
{
	typedef double num_t ;

	typedef KDTreeSingleIndexAdaptor<
			L2_Simple_Adaptor<num_t, Points2D<num_t> > ,
			Points2D<num_t>,
			2 /* dim */
			> my_kd_tree_t;

	my_kd_tree_t* index;

	typedef typename T::const_iterator Titerator;

public:

	bool init(T& points)
	{
		Points2D<num_t> pointCloud;
		pointCloud.pts.resize(points.size());

		std::cout<<pointCloud.pts.size()<<std::endl;

		int i =0 ;
		for(Titerator it = points.begin(); it!=points.end(); it++ , i++)
		{
			pointCloud.pts[it->first].x = (double)static_cast<estimateType*>(it->second)->estimate().toVector()(0);
			pointCloud.pts[it->first].y = (double)static_cast<estimateType*>(it->second)->estimate().toVector()(1);

			std::cout<<i<<" : "<<pointCloud.pts[i].x<<" , "<<pointCloud.pts[i].y<<std::endl;
		}

		index = new my_kd_tree_t(2 /*dim*/, pointCloud, KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
		index->buildIndex();
		return true;

	}

	bool query(const num_t* queryPoint, num_t radius)
	{
		const num_t search_radius = static_cast<num_t>(radius);
		std::vector<std::pair<size_t,num_t> > ret_matches;

		nanoflann::SearchParams params;
		//params.sorted = false;

		const size_t nMatches = index->radiusSearch(&queryPoint[0],search_radius, ret_matches, params);

		cout << "radiusSearch(): radius=" << search_radius << " -> " << nMatches << " matches\n";
		for (size_t i=0;i<nMatches;i++)
			cout << "idx["<< i << "]=" << ret_matches[i].first << " dist["<< i << "]=" << ret_matches[i].second << endl;
		cout << "\n";

		return true;
	}

};

#endif
#endif /* NEIGHBOURSEARCH_HPP_ */
