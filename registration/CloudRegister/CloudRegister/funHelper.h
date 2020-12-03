#pragma once

#include "BaseType.h"

//pcl
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>

#include "GeometryUtils.h"

namespace CloudReg 
{
	Eigen::vector<Eigen::Vector3d> ininterpolateSeg(const Eigen::Vector3d& sPoint, 
		const Eigen::Vector3d& ePoint, const double step);
	std::vector<std::string> splitByCharacter(const std::string& strtem, const char a);
	bool writePCDFile(const std::string& name, Eigen::vector<Eigen::Vector3d>& vecCloud);

	double calcArea(const Eigen::vector<Eigen::Vector2d>& vecPts);

	void uniformSampling(double radius,
							pcl::PointCloud<pcl::PointXYZ>::Ptr cloud,
							pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_filtered);

	double pointToPLaneDist(const Eigen::Vector4d &plane, const pcl::PointXYZ &p);

	pcl::PointXYZRGB getColorPtByDist(pcl::PointXYZ &p, double dist);
	
	void projectionToPlane(Eigen::Vector4d &plane,
							pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, 
							pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_projected);

	void searchBoundaries(pcl::PointCloud<pcl::PointXYZ>::Ptr input_cloud,
							std::vector<int> &boundIndices);

	pcl::PointCloud<pcl::PointXYZ>::Ptr calcCloudBorder(
			const std::string &name,
			pcl::PointCloud<pcl::PointXYZ>::Ptr cloud,
			Eigen::Vector4d &cloudPlane,
			std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> &cadBorder,
			std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> &cloudBorder);	

	double distToLine(const Eigen::Vector3d& p, 
						const Eigen::Vector3d& s, const Eigen::Vector3d& e);

	std::vector<int> clusterMainStructure(PointCloud::Ptr cloud, float distance) ;

	bool detectLineEndPoints(PointCloud::Ptr inliers, 
						Eigen::VectorXf &params,
						double radius, 
						std::pair<Eigen::Vector3d, Eigen::Vector3d> &segment);				
}


