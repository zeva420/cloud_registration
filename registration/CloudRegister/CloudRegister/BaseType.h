#pragma once

#include <iostream>
#include <fstream>

#include <utility>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <list>
#include <queue>
#include <unordered_map>
#include <unordered_set>

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/StdVector>

//pcl
#include <pcl/common/common.h>
#include<pcl/point_types.h>
#include<pcl/io/pcd_io.h>



#ifndef EPS_FLOAT_DOUBLE
#define EPS_FLOAT_DOUBLE 0.00001
#endif



namespace Eigen
{
	template<class T, class Alloc = Eigen::aligned_allocator<T>>
	using vector = std::vector<T, Alloc>;

	template < class T, class Alloc = Eigen::aligned_allocator<T> >
	using list = std::list<T, Alloc>;

	template < class T, class Alloc = Eigen::aligned_allocator<T> >
	using deque = std::deque<T, Alloc>;

	template < class T,
		class Compare = std::less<T>,
		class Alloc = Eigen::aligned_allocator<T>
	>
		using set = std::set<T, Compare, Alloc>;

	template < class T,
		class Compare = std::less<T>,
		class Alloc = Eigen::aligned_allocator<T>
	>
		using multiset = std::multiset<T, Compare, Alloc>;

	template < class Key,
		class Hash = std::hash<Key>,
		class Pred = std::equal_to<Key>,
		class Alloc = Eigen::aligned_allocator<Key>
	>
		using unordered_set = std::unordered_set<Key, Hash, Pred, Alloc>;

	template < class Key,
		class Hash = std::hash<Key>,
		class Pred = std::equal_to<Key>,
		class Alloc = Eigen::aligned_allocator<Key>
	>
		using unordered_multiset = std::unordered_multiset<Key, Hash, Pred, Alloc>;

	template < class Key,
		class T,
		class Compare = std::less<Key>,
		class Alloc = Eigen::aligned_allocator<std::pair<const Key, T> >
	>
		using map = std::map<Key, T, Compare, Alloc>;

	template < class Key,
		class T,
		class Compare = std::less<Key>,
		class Alloc = Eigen::aligned_allocator<std::pair<const Key, T> >
	>
		using multimap = std::multimap<Key, T, Compare, Alloc>;

	template < class Key,
		class T,
		class Hash = std::hash<Key>,
		class Pred = std::equal_to<Key>,
		class Alloc = Eigen::aligned_allocator< std::pair<const Key, T> >
	>
		using unordered_map = std::unordered_map<Key, T, Hash, Pred, Alloc>;

	template < class Key,
		class T,
		class Hash = std::hash<Key>,
		class Pred = std::equal_to<Key>,
		class Alloc = Eigen::aligned_allocator< std::pair<const Key, T> >
	>
		using unordered_multimap = std::unordered_multimap<Key, T, Hash, Pred, Alloc>;

	template <class T, class Container = Eigen::deque<T> >
	using queue = std::queue<T, Container>;

	template <class T,
		class Container = Eigen::vector<T>,
		class Compare = std::less<typename Container::value_type> >
		using priority_queue = std::priority_queue<T, Container, Compare>;

	template < class T, size_t N >
	using array = std::array<T, N>;

	using MatrixXb_t = Matrix< bool, Eigen::Dynamic, Eigen::Dynamic >;
}
using Point = pcl::PointXYZ;
using PointCloud = pcl::PointCloud<Point>;

namespace CloudReg
{
	using seg_pair_t = std::pair<Eigen::Vector3d, Eigen::Vector3d>;
	using vec_seg_pair_t = std::vector<seg_pair_t>;

	struct calcMeassurment_t
	{
		double value = -1;
		std::vector<seg_pair_t> rangeSeg;
	};
}