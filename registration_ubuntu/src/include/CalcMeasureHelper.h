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
	
	std::tuple<std::size_t, std::size_t, int> getGrowAxisAndDir(const Eigen::Vector3d& sPt, const Eigen::Vector3d& ePt);

	std::tuple<bool, Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d>
		calcOverlap(seg_pair_t& toSeg, seg_pair_t& calcSeg);

	std::vector<seg_pair_t> calcBoxSegPair(std::vector<Eigen::Vector3d>& vecPt);

	//calculate ruler
	std::vector<Eigen::Vector3d> createRulerBox(seg_pair_t ruler, int thicknessDir, 
						double thickness, double width);
	bool calIntersection(seg_pair_t line1, seg_pair_t line2, Eigen::Vector3d& intersec);
	bool calRuler3d(const std::vector<seg_pair_t>& wallBorder, const std::vector<seg_pair_t>& holeBorder, 
					 const seg_pair_t& rotateLine, const Eigen::Vector3d& P0,
					 const float& theta, seg_pair_t& ruler);
}


