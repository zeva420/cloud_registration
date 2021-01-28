#pragma once


#include "CalcMeasureHelper.h"

namespace CloudReg 
{

	std::tuple<std::vector<calcIdx2Meassurment_t>,std::vector<seg_pair_t>>
		CalcNetHeight(const std::vector<seg_pair_t>& roofBorder,
			const PointCloud::Ptr pCloud,
			const Eigen::Vector4d& plane,
			const Eigen::Vector3d& center,
			const std::string& name,
			const double calcLengthTh = 1.5,
			const double moveRangeTh = 1.2);
	
	std::tuple<std::vector<calcIdx2Meassurment_t>, std::vector<calcIdx2Meassurment_t>,
		std::vector<seg_pair_t>, std::vector<seg_pair_t>>
		CalcHeightRange(const std::vector<seg_pair_t>& roofBorder,
			const std::vector<seg_pair_t>& rootBorder,
			const std::vector<std::vector<seg_pair_t>>& allWallBorder,
			const PointCloud::Ptr pRoof,
			const PointCloud::Ptr pRoot,
			const Eigen::Vector3d& center,
			const double calcHeight = 1.0,
			const double calcLengthTh = 1.5,
			const double moveRangeTh = 1.2);
}


