#pragma once


#include "CalcMeasureHelper.h"
namespace CloudReg 
{
	void calcHole(const std::pair<Eigen::Vector3d, Eigen::Vector3d>& horizen,
			std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>& holeBorder, 
			pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud);
}


