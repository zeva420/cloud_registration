#pragma once

#include "BaseType.h"



#include "GeometryUtils.h"

namespace CloudReg 
{
	Eigen::vector<Eigen::Vector3d> ininterpolateSeg(const Eigen::Vector3d& sPoint, 
		const Eigen::Vector3d& ePoint, const double step);
	std::vector<std::string> splitByCharacter(const std::string& strtem, const char a);

	bool writePCDFile(const std::string& name, Eigen::vector<Eigen::Vector3d>& vecCloud);

	double calcArea(const Eigen::vector<Eigen::Vector2d>& vecPts);

	void uniformSampling(const double radius,
							const pcl::PointCloud<pcl::PointXYZ>::Ptr cloud,
							pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_filtered);

	double pointToPLaneDist(const Eigen::Vector4d &plane, const pcl::PointXYZ &p);
	double pointToPLaneDist(const Eigen::VectorXf &plane, const pcl::PointXYZ &p);

	Eigen::Vector3d pointToPlaneRoot(const Eigen::Vector4d &plane, const Eigen::Vector3d &point);

	pcl::PointXYZRGB getColorPtByDist(pcl::PointXYZ &p, const double dist);
	
	void getWallColor(const float dis, unsigned int & r, unsigned int & g, unsigned int & b);
	
	void projectionToPlane(const Eigen::Vector4d &plane,
							const pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, 
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
										const std::pair<Eigen::Vector3d, Eigen::Vector3d> &seg,
										const double searchDist);

	std::pair<double, Eigen::Vector3d> findNearestPt(
			const std::vector<Eigen::Vector3d> &vecPts, const Eigen::Vector3d &point);

	bool groupPlanesBySamePt(const std::vector<std::vector<Eigen::Vector3d>> &segPtsOfPlanes,
						std::set<std::set<int>> &planeIdxGroup);

	bool interSectionOf3Planes(const std::vector<Eigen::Vector4d> &cloudPlaneVec,
							const std::set<std::set<int>> &idGroups, 
							std::vector<Eigen::Vector3d> &focalPointVec);

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


}


