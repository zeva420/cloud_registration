#pragma once

#include "CalcMeasureHelper.h"

namespace CloudReg 
{
	void calcSquareness(const std::vector<vec_seg_pair_t>& vecWall, const std::vector<vec_seg_pair_t>& holeBorders,
			std::vector<PointCloud::Ptr> pClouds, std::map<std::size_t, std::vector<std::size_t>> holeMap, const double calcLengthTh = 1.);
}