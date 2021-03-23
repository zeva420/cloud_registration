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

	struct calcIdx2Meassurment_t
	{
		std::pair<std::size_t, std::size_t> idx;
		std::vector<calcMeassurment_t> vecCalcRet;
	};

	enum CloudItemType {
		CLOUD_BEAM_E,
		CLOUD_BOTTOM_E,
		CLOUD_WALL_E,
		CLOUD_TOP_E,
		CLOUD_MAX_E
	};

	struct CloudItem {
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
			CloudItem(pcl::PointCloud<pcl::PointXYZ>::Ptr pData)
			:pCloud_(pData)
		{
		}

		//origin data
		const pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud_;
		pcl::PointCloud<pcl::PointXYZ>::Ptr pCADCloud_;
		CloudItemType type_ = CLOUD_MAX_E;
		std::size_t parentIndex_ = 9999;
		Eigen::Vector4d cloudPlane_;
		Eigen::Vector4d cadPlane_;
		// if type is WALL include windows and door,the border.front belong to wall 
		std::vector<std::vector<seg_pair_t>> cloudBorder_;
		std::vector<std::vector<seg_pair_t>> cadBorder_;
	};

	using vecItems_t = Eigen::vector<CloudItem>;

	struct pairHash
	{
		size_t operator()(const std::pair<std::size_t, std::size_t> &pair_) const
		{
			return std::hash<std::size_t>()(pair_.first) ^ std::hash<std::size_t>()(pair_.second);
		}
	};

	// whitewash wall.
	struct Salient {
		float height_;
		float area_;
		Eigen::Vector3f boundingBoxMin_, boundingBoxMax_; // bounding box
	};

	struct Wall {
		// in & threshold
		double length_;
		double pos_;

		// out & temp
		double maxSalientHeight_{ 0. }; // the max height considered, not the actually height, 
		std::vector<Salient> salients_;

		double paintThickness_;
		double salientChipping_{ 0. }; // refer to maxSalientHeight.
		double wallChipping_{ 0. };
		std::vector<double> saliensChippingHeight_; // for each salients, 0 if no need to chip
	};

	enum TargetItemType {
		LEFT_TOP_E,
		LEFT_BOTTON_E,
		RIGHT_TOP_E,
		RIGHT_BOTTON_E
	};
}