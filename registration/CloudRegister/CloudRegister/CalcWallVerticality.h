#pragma once

#include "CalcMeasureHelper.h"

namespace CloudReg 
{
	std::tuple<std::vector<calcMeassurment_t>, std::vector<seg_pair_t>> calcVerticality(
			const std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>& wallBorder, 
            const std::vector<std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>>& holeBorders,
			pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud, Eigen::Vector4d plane, int wallIndex);
}