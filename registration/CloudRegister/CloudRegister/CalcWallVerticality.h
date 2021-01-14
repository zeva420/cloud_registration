#pragma once

#include "CalcMeasureHelper.h"

namespace CloudReg 
{
	void calcVerticality(const std::pair<Eigen::Vector3d, Eigen::Vector3d>& horizen,
			const std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>& wallBorder, 
            const std::vector<std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>>& holeBorders,
			pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud, int index);
}