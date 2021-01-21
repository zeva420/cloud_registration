#pragma once

#include "CalcMeasureHelper.h"

namespace CloudReg 
{

	std::tuple<std::vector<calcMeassurment_t>, std::vector<seg_pair_t>>
	calWallFlatness(const vec_seg_pair_t& wallBorder, const std::vector<vec_seg_pair_t>& holeBorders,
			pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud, Eigen::Vector4d cadPlane, int index);
}