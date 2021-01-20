#pragma once


#include "CalcMeasureHelper.h"
namespace CloudReg 
{
	std::map<std::pair<std::size_t, std::size_t>, std::vector<calcMeassurment_t>>
	CalcCorner(const std::vector<std::vector<seg_pair_t>>& allWallBorder,
			const std::map<std::size_t, std::vector<std::vector<seg_pair_t>>>& holeBorder,
			const std::vector<PointCloud::Ptr>& vecCloud,
			const double calcLengthTh = 0.13);
}


