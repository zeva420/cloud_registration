#pragma once

#include "CalcMeasureHelper.h"

namespace CloudReg 
{

    void calcRootFlatness(const std::vector<seg_pair_t>& rootBorder,
                          Eigen::Vector4d plane, pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud,
                          const double calcLengthTh = 1.5);
}