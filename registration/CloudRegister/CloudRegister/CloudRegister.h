#pragma once

#include "BaseType.h"

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
		// if type is WALL include windows and door,the border.front belong to wall 
		std::vector<std::vector<seg_pair_t>> cloudBorder_;
		std::vector<std::vector<seg_pair_t>> cadBorder_;
	};

	

	using vecItems_t = Eigen::vector<CloudItem>;
	
	class CADModel;
	class TransformOptimize;
	#ifdef UBUNTU_SWITCH
	class CloudRegister
	#else
	class  __declspec(dllexport) CloudRegister
	#endif
	{
	public:
		CloudRegister();
		~CloudRegister();

		//bNeedOptimize = false not run optimize
		bool run(std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr>& vecCloudPtr,
			const std::string& CAD_File, const bool bNeedOptimize = true);

		const std::map<CloudItemType, vecItems_t>&
			getAllCloudPlane() const;

		//p_rgb.r = distError; radius = 0.05m
		pcl::PointCloud<pcl::PointXYZI>::Ptr 
			calcDistError(const pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud_,
			const Eigen::Vector4d& plane, const double radius = 0.05) const;
		
		//calcAllCorner: corner of two neighbouring walls (in height 0.3m and 1.5m)
		std::map<std::pair<std::size_t, std::size_t>, std::vector<calcMeassurment_t>>
			calcAllCorner();

		//calcLengthTh: the shorest wall length, 
		std::tuple<std::vector<std::vector<calcMeassurment_t>>, std::vector<seg_pair_t>> 
			calcRoofNetHeight(const double calcLengthTh = 1.5);

		//first roof second root
		//calcHeight: the height from bottom 
		//calcLengthTh: the shorest wall length
		std::tuple<std::vector<std::vector<calcMeassurment_t>>, std::vector<std::vector<calcMeassurment_t>>,std::vector<seg_pair_t>>
			calcPlaneRange(const double calcHeight = 1.0,const double calcLengthTh = 1.5);

		//calcLengthTh: the shorest wall length
		std::tuple<std::map<std::pair<std::size_t, std::size_t>,
			std::vector<calcMeassurment_t>>, std::vector<seg_pair_t>>
			calcDepth(const double calcLengthTh = 0.8);

		//calcLengthTh: the shorest wall length
		std::tuple<std::map<std::pair<std::size_t, std::size_t>,
			std::vector<calcMeassurment_t>>, std::vector<seg_pair_t>>
			calcBay(const double calcLengthTh = 0.8);

		//key.first = wall index  key.second = hole border index in cloudBorder_
		//calcMeassurment_t order height width cross root
		std::map<std::pair<std::size_t, std::size_t>, std::vector<std::vector<calcMeassurment_t>>>
			calcAllHole();

		//planeType = cloud, use cloud plane, else use cad plane
		std::map<int, std::tuple<std::vector<calcMeassurment_t>, std::vector<seg_pair_t>>>
			calcWallVerticality(const std::string& planeType = "cad");

		std::map<int, std::tuple<std::vector<calcMeassurment_t>, std::vector<seg_pair_t>>>
			calcWallFlatness(const std::string& planeType = "cad");

		std::map<std::pair<int, int>,std::tuple<std::vector<calcMeassurment_t>, std::vector<seg_pair_t>>>
    		calcAllSquareness(const double calcLengthTh = 1.);

		std::vector<std::tuple<std::vector<calcMeassurment_t>, std::vector<seg_pair_t>>> 
    		calcRootFlatness(const std::string& planeType = "cad", const double calcLengthTh = 1.5);

	private:
		
		void fillRet(CADModel& cad, TransformOptimize& optimitor);
		void calcAllCloudBorder(CADModel& cad);

		std::map<CloudItemType, vecItems_t> mapCloudItem_;

	};
}

