#pragma once


#include "CalcMeasureHelper.h"
namespace CloudReg 
{

	void CalcNetHeight(const std::vector<seg_pair_t>& roofBorder,
			const PointCloud::Ptr pCloud,
			const Eigen::Vector4d& plane,
			const double calcLengthTh = 1.5);
}


