#pragma once

#include "CalcMeasureHelper.h"

namespace CloudReg 
{
	using vec_seg_pair_t = std::vector<seg_pair_t>;

	void calcWallFlatness(const vec_seg_pair_t& wallBorder, const std::vector<vec_seg_pair_t>& holeBorders,
			pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud, Eigen::Vector4d cadPlane, int index);
}