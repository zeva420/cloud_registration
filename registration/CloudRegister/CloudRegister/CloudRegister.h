#pragma once

#include <string>
#include <vector>

#include "TransformOptimize.h"

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
		std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> cloudBorder_; 
		std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> cadBorder_;
	};

	
	using vecItems_t = std::vector<CloudItem>;
	using pairCloud_t = std::pair<CloudItem*, CloudItem*>;

	class CADModel;
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

		//p_rgb.r = distError;
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr 
			calcDistError(const pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud_,
			const Eigen::Vector4d& plane, const double downRatio= 0.1) const;


	private:
		
		void fillRet(CADModel& cad, TransformOptimize::optCloudRets &optRets);
		std::map<CloudItemType, vecItems_t> mapCloudItem_;
		std::map<pairCloud_t, std::pair<double, double>> mapCorner_;

	};
}

