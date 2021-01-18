#pragma once

#include <string>
#include <vector>

#include<pcl/point_types.h>
#include<pcl/io/pcd_io.h>

namespace CloudReg
{
	enum CloudItemType {
		CLOUD_BEAM_E,
		CLOUD_BOTTOM_E,
		CLOUD_WALL_E,
		CLOUD_TOP_E,
		CLOUD_MAX_E
	};

	struct CloudItem {
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
			CloudItem(pcl::PointCloud<pcl::PointXYZ>::Ptr pData)
			:pCloud_(pData)
		{
		}

		//origin data
		const pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud_;
		pcl::PointCloud<pcl::PointXYZ>::Ptr pCADCloud_;
		CloudItemType type_ = CLOUD_MAX_E;
		std::size_t parentIndex_ = 9999;
		Eigen::Vector4d cloudPlane_;
		Eigen::Vector4d cadPlane_;
		// if type is WALL include windows and door
		std::vector<std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>> cloudBorder_; 
		std::vector<std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>> cadBorder_;
	};

	
	using vecItems_t = std::vector<CloudItem, Eigen::aligned_allocator<CloudItem>>;
	using pairCloud_t = std::pair<CloudItem*, CloudItem*>;

	class CADModel;
	class TransformOptimize;
	class  __declspec(dllexport) CloudRegister
	{
	public:
		CloudRegister();
		~CloudRegister();

		
		bool run(std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr>& vecCloudPtr,
			const std::string& CAD_File);

		const std::map<CloudItemType, vecItems_t>&
			getAllCloudPlane() const;

		//value first is 0.3 second 1.5
		const std::map<pairCloud_t, std::pair<double, double>>&
			getAllCorner() const;

		//p_rgb.r = distError; radius = 0.05m
		pcl::PointCloud<pcl::PointXYZI>::Ptr 
			calcDistError(const pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud_,
			const Eigen::Vector4d& plane, const double radius = 0.05) const;


	private:
		
		void fillRet(CADModel& cad, TransformOptimize& optimitor);
		void calcAllCloudBorder(CADModel& cad);
		int findMatchCloud(const Eigen::Vector4d &plane,
			std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &vecOrigCloud);

		std::map<CloudItemType, vecItems_t> mapCloudItem_;
		std::map<pairCloud_t, std::pair<double, double>> mapCorner_;

	};
}

