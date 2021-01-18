#pragma once

#include "BaseType.h"

//pcl
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/common/common.h>

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
	
	void getWallColor(float dis, unsigned int & r, unsigned int & g, unsigned int & b);
	
	void projectionToPlane(Eigen::Vector4d &plane,
							pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, 
							pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_projected);

	void searchBoundaries(pcl::PointCloud<pcl::PointXYZ>::Ptr input_cloud,
							pcl::PointCloud<pcl::Normal>::Ptr normals,
							std::vector<int> &boundIndices);

	void searchBoundaries(pcl::PointCloud<pcl::PointXYZ>::Ptr input_cloud,
							std::vector<int> &boundIndices);

	std::vector<Eigen::Vector3d> convertCloudToEigenVec(
									const pcl::PointCloud<pcl::PointXYZ>::Ptr cloud);

	std::pair<double, std::pair<Eigen::Vector3d, Eigen::Vector3d>> findNearestSeg(
										const std::vector<Eigen::Vector3d> &vecPts, 
										const std::pair<Eigen::Vector3d, Eigen::Vector3d> &seg);

	std::pair<double, Eigen::Vector3d> findNearestPt(
			const std::vector<Eigen::Vector3d> &vecPts, const Eigen::Vector3d &point);

	bool groupPlanesBySamePt(const std::vector<std::vector<Eigen::Vector3d>> &segPtsOfPlanes,
						std::set<std::set<int>> &planeIdxGroup);

	bool interSectionOf3Planes(const std::vector<Eigen::Vector4d> &cloudPlaneVec,
							const std::set<std::set<int>> &idGroups, 
							std::vector<Eigen::Vector3d> &focalPointVec);

	std::vector<Eigen::Vector3d> calcWallNodes(const std::string &name, 
			pcl::PointCloud<pcl::PointXYZ>::Ptr cloud,
			Eigen::Vector4d &cloudPlane,
			const std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> &outerSegs);

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

	void planeFitting(float distTh, 
					pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, 
					Eigen::VectorXf &coeff, std::vector<int> &inlierIdxs);
					
	Eigen::Vector4d calcPlaneParam(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud);

	bool isPointInPolygon2D(const Eigen::Vector2d &point, 
			const Eigen::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>> &segments);

	bool isPointInPolygon3D(const Eigen::Vector3d &point, 
			const Eigen::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> &segments);
			
	bool interSectionOfLineToLine(const Eigen::VectorXd &line1, 
					const Eigen::VectorXd &line2, Eigen::Vector3d &interSectionPt);

	bool interSectionOfLineToPlane(const Eigen::VectorXd &line,
		const Eigen::Vector4d &plane, Eigen::Vector3d &interSectionPt);

	bool interSectionOfPlaneToPlane(const Eigen::Vector4d &plane1,
		const Eigen::Vector4d &plane2, Eigen::VectorXd &interSectionLine);

	double calcCorner(const Eigen::Vector3d &n1, const Eigen::Vector3d &n2);

	double calcCloudPairCorner(const std::string &name,
						const pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud1,
						const pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud2,
						const Eigen::Vector3d &floorPt, double height,
						Eigen::Vector3d center, Eigen::Vector4d bottomPlane);

	// cornerPoint(2, 0) == z
	//double calcCorner_beta(PointCloud::Ptr cloud1, PointCloud::Ptr cloud2, const Eigen::Vector3f& cornerPoint, float z);

}


