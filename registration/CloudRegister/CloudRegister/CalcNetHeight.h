#pragma once


#include "CalcMeasureHelper.h"

namespace CloudReg 
{

	std::tuple<std::vector<calcMeassurment_t>,std::vector<seg_pair_t>>
		CalcNetHeight(const std::vector<seg_pair_t>& roofBorder,
			const PointCloud::Ptr pCloud,
			const Eigen::Vector4d& plane,
			const std::string& name,
			const double calcLengthTh = 1.5);
	
	std::tuple<std::vector<calcMeassurment_t>, std::vector<seg_pair_t>>
		CalcPlaneRange(const std::vector<seg_pair_t>& roofBorder,
			const std::vector<seg_pair_t>& rootBorder,
			const std::vector<std::vector<seg_pair_t>>& allWallBorder,
			const PointCloud::Ptr pRoof,
			const PointCloud::Ptr pRoot,
			const double calcHeight = 1.0,
			const double calcLengthTh = 1.5);
}


