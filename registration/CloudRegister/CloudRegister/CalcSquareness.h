#pragma once

#include "CalcMeasureHelper.h"

namespace CloudReg 
{
	std::map<std::pair<int, int>,std::tuple<std::vector<calcMeassurment_t>, std::vector<seg_pair_t>>>
    calcSquareness(const std::vector<vec_seg_pair_t>& vecWall,std::vector<PointCloud::Ptr> pClouds, 
                        std::map<std::size_t, std::vector<vec_seg_pair_t>> holeMap, const double calcLengthTh);
}