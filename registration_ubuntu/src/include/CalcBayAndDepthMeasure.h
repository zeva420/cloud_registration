#pragma once


#include "CalcMeasureHelper.h"
namespace CloudReg 
{
	using vec_seg_pair_t = std::vector<seg_pair_t>;

	void calcDepthorBay(const std::vector<seg_pair_t>& rootBorder,
			const std::vector<vec_seg_pair_t>& allWallBorder,
			const std::map<std::size_t, std::vector<vec_seg_pair_t>>& holeBorder,
			const std::vector<PointCloud::Ptr>& vecCloud,
			const Eigen::vector<Eigen::Vector4d>& vecPlane,
			const int optType,
			const double calcLengthTh = 0.8);
}


