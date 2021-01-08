#pragma once


#include "BaseType.h"
namespace CloudReg 
{
	using seg_pair_t = std::pair<Eigen::Vector3d, Eigen::Vector3d>;
	struct calcMeassurment_t 
	{
		double value;
		std::vector<seg_pair_t> rangeSeg;
	};
	
	void writePCDFile(const std::string& name, const std::vector<seg_pair_t>& segA, const std::vector<seg_pair_t>& segB);
	void writePCDFile(const std::string& name, const PointCloud::Ptr pCloud, std::vector<seg_pair_t>& border);
	
	Eigen::Vector3d calcPerpendicular(const Eigen::Vector3d &pt,
			const Eigen::Vector3d &start, const Eigen::Vector3d &end);
	
	void groupDirection(const Eigen::Vector3d& horizenSeg, const std::vector<seg_pair_t>& border, 
			std::vector<seg_pair_t>& vecVertical, std::vector<seg_pair_t>& vecHorizen);

	PointCloud::Ptr filerCloudByRange(pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud,
			const pcl::PointXYZ& min, const pcl::PointXYZ& max);
	
	Eigen::vector<Eigen::Vector3d> getNearestPt(const Eigen::vector<Eigen::Vector3d>& vecCalcPt, 
			const PointCloud::Ptr pCloud, const double maxDist);
	
	bool isRootInSeg(const seg_pair_t& seg, const Eigen::Vector3d& p);
	
	void groupDirectionIndex(const Eigen::Vector3d& horizenSeg, const std::vector<seg_pair_t>& border, 
			std::vector<std::size_t>& vecVertical, std::vector<std::size_t>& vecHorizen);
}


