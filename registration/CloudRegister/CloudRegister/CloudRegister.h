#pragma once

#include "BaseType.h"

namespace CloudReg
{
	
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
			const std::string& CAD_File, const bool bNeedOptimize = true, const bool bOriginCloud = false, const bool changeCADOrder = false);

		const std::map<CloudItemType, vecItems_t>&
			getAllCloudPlane() const;

		//p_rgb.r = distError; radius = 0.05m
		pcl::PointCloud<pcl::PointXYZI>::Ptr 
			calcDistError(const pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud_,
			const Eigen::Vector4d& plane, const double radius = 0.05) const;
		
		//calcAllCorner: corner of two neighbouring walls (in height 0.3m and 1.5m)
		std::map<std::pair<std::size_t, std::size_t>, std::vector<calcMeassurment_t>>
			calcAllCorner(const double calcLengthTh = 0.09);

		//calcLengthTh: the shorest wall length, 
		std::tuple<std::vector<calcIdx2Meassurment_t>,std::vector<seg_pair_t>>
			calcRoofNetHeight(const double calcLengthTh = 1.5);

		//first roof second root
		//calcHeight: the height from bottom 
		//calcLengthTh: the shorest wall length
		std::tuple<std::vector<calcIdx2Meassurment_t>,std::vector<calcIdx2Meassurment_t>,
			std::vector<seg_pair_t>, std::vector<seg_pair_t>>
			calcPlaneRange(const double calcHeight = 1.0,const double calcLengthTh = 1.5, const double moveRangeTh = 1.0);


		//calcLengthTh: the shorest wall length
		std::tuple<std::vector<calcIdx2Meassurment_t>, std::vector<seg_pair_t>>
			calcDepth(const double calcLengthTh = 0.8);

		//calcLengthTh: the shorest wall length
		std::tuple<std::vector<calcIdx2Meassurment_t>, std::vector<seg_pair_t>>
			calcBay(const double calcLengthTh = 0.8);

		//key.first = wall index  key.second = hole border index in cloudBorder_
		//calcMeassurment_t order height width cross root
		std::map<std::pair<std::size_t, std::size_t>, std::vector<std::vector<calcMeassurment_t>>>
			calcAllHole();

		//planeType = cloud, use cloud plane, else use cad plane
		std::map<int, std::tuple<std::vector<calcMeassurment_t>, std::vector<seg_pair_t>>>
			calcWallVerticality(const std::string& planeType = "local"); // "cloud" "cad" "local"

		std::map<int, std::tuple<std::vector<calcMeassurment_t>, std::vector<seg_pair_t>>>
			calcWallFlatness(const std::string& planeType = "local"); // "cloud" "cad" "local"

		std::map<std::pair<int, int>,std::tuple<std::vector<calcMeassurment_t>, std::vector<seg_pair_t>>>
    		calcAllSquareness(const double calcLengthTh = 1.);

		std::vector<std::tuple<std::vector<calcMeassurment_t>, std::vector<seg_pair_t>>> 
    		calcRootFlatness(const std::string& planeType = "local", const double calcLengthTh = 1.5); // "cloud" "cad" "local"

	
		std::vector<Wall> whitewashPaint(double minSalientArea = 0., double maxSalientHeight = 1.,
			double minWallPaintThickness = 0.005, double minSalientPaintThickness = 0.002, 
			double designedPaintThickness = 0.01, double lowDeviation = -0.01, double highDeviation = 0.01, double deviationCompensation = 0.002);
		
		calcMeassurment_t getTargetPoint(const TargetItemType ptType, const std::size_t wallIndex,
			double hDis, double vDis, double radius);
	private:
		
		void fillRet(CADModel& cad, TransformOptimize& optimitor);
		void calcAllCloudBorder(CADModel& cad);

		std::map<CloudItemType, vecItems_t> mapCloudItem_;
		std::vector<Wall> walls_;
		Eigen::Vector3d centerPt_;

	};
}

