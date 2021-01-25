#pragma once


#include "CalcMeasureHelper.h"
namespace CloudReg 
{
	
	std::vector<std::vector<calcMeassurment_t>> calcHole(const seg_pair_t& horizen,
			const std::vector<seg_pair_t>& holeBorder,
			pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud,const std::string& name);
}


