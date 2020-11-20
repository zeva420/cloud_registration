#pragma once

#include <string>
#include <vector>

#include<pcl/point_types.h>
#include<pcl/io/pcd_io.h>

namespace CloudReg
{
	struct CloudItem {
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
			CloudItem(pcl::PointCloud<pcl::PointXYZ>::Ptr pData)
			:pCloud_(pData)
		{
		}

		const pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud_;
		Eigen::Vector3d cloudPlane_;
		Eigen::Vector3d cadPlane_;
		std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> cloudBorder_;
		std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> cadBorder_;
	};

	enum CloudItemType {
		CLOUD_BEAM_E,
		CLOUD_BOTTOM_E,
		CLOUD_WALL_E,
		CLOUD_TOP_E,
		CLOUD_MAX_E
	};
	using vecItems_t = std::vector<CloudItem>;
	using pairCloud_t = std::pair<CloudItem*, CloudItem*>;

	class  __declspec(dllexport) CloudRegister
	{
	public:
		CloudRegister();
		~CloudRegister();

		
		bool run(std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr>& vecCloudPtr,
			const std::string& CAD_File);

		const std::map<CloudItem, vecItems_t>& 
			getAllCloudPlane() const;

		const std::map<pairCloud_t, std::pair<double, double>>&
			getAllCorner() const;

		pcl::PointCloud<pcl::PointXYZRGB>::Ptr 
			calcDistError(const pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud_,
			const Eigen::Vector3d& plane, const double downRatio= 0.1) const;

		pcl::PointCloud<pcl::PointXYZRGB>::Ptr
			genCloudByModel(const Eigen::Vector3d& planePara,
				const std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>& border)const;

	private:
		
		std::map<CloudItem, vecItems_t> mapCloudItem_;
		std::map<pairCloud_t, std::pair<double, double>> mapCorner_;

	};
}

