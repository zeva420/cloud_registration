#pragma once

#include "BaseType.h"
#include "glog/logging.h"
#include "CADModel.h"

//pcl
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>


namespace CloudReg
{

class Segmentation
{
public:
    Segmentation()
	{
	}

    ~Segmentation()
    {

    }

    bool run(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, CADModel &cadModel);

private:

	bool downSampling(double radius, pcl::PointCloud<pcl::PointXYZ>::Ptr cloud,
                        pcl::PointCloud<pcl::PointXYZ>::Ptr samplingCloud);

	bool statisticsForPointZ(float binSizeTh,
                pcl::PointCloud<pcl::PointXYZ>::Ptr cloud,
                std::vector<std::pair<int, std::vector<int>>> &zToNumVec);

	bool calibrateDirectionToAxisZ(
            			pcl::PointCloud<pcl::PointXYZ>::Ptr inputCloud, 
						CADModel &cad, Eigen::Matrix4d &T);

	bool calibrateCloudToCadPose(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, 
						CADModel &cad, Eigen::Matrix4d &bestT);

	bool calibrateCloudByUpOrDown(CADModel &cad, pcl::PointCloud<pcl::PointXYZ>::Ptr ptSubSet1, 
    							pcl::PointCloud<pcl::Normal>::Ptr normalSubSet1, 
								double &moveZ);

	bool calibrateCloudByWall(CADModel &cad, pcl::PointCloud<pcl::PointXYZ>::Ptr ptSubSet2, 
    							pcl::PointCloud<pcl::Normal>::Ptr normalSubSet2, 
								Eigen::Matrix<double, 2, 3> &T);

	bool doSegmentation(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, 
						CADModel &cad);

	pcl::PointCloud<pcl::PointXYZ>::Ptr segmentateSingleModelItem(
						pcl::PointCloud<pcl::PointXYZ>::Ptr inputCloud, 
						const ModelItem &item, const std::string &name);

	void estimateNormals(float radius, 
						pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, 
						pcl::PointCloud<pcl::Normal>::Ptr normals);

	void regionGrowingSegmentation(
						pcl::PointCloud<pcl::PointXYZ>::Ptr cloud,
						pcl::PointCloud <pcl::Normal>::Ptr normals,
						std::vector <pcl::PointIndices> &clusters);

	pcl::PointCloud<pcl::Normal>::Ptr getSubSet(pcl::PointCloud<pcl::Normal>::Ptr cloud, 
                const std::vector<int>& indices, bool negative);

	Eigen::Matrix<double, 2, 3> estimateTransform(
                const Eigen::Vector2d& s1, const Eigen::Vector2d& e1,
	            const Eigen::Vector2d& s2, const Eigen::Vector2d& e2) ;

	Eigen::Vector2d to_point_xy(const Eigen::Vector3d &point);

	bool planeFitAndCheck(pcl::PointCloud<pcl::PointXYZ>::Ptr nearPoints, 
                        Eigen::Vector4d &fitPlane, std::vector<int> &inlierIndices);

	pcl::PointCloud<pcl::PointXYZRGB>::Ptr to_rgp_cloud(
                            const pcl::PointCloud<pcl::PointXYZ>::Ptr cloud,
                            int r, int g, int b, double sampleRadius = 0.0);

};

} //namespace CloudReg

