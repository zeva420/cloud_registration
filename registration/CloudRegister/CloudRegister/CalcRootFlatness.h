#pragma once

#include "CalcMeasureHelper.h"

namespace CloudReg 
{

    std::vector<std::tuple<std::vector<calcMeassurment_t>, std::vector<seg_pair_t>>>
    calRootFlatness(const std::vector<seg_pair_t>& rootBorder,
                          Eigen::Vector4d plane, pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud,
                          const double calcLengthTh = 1.5);
}