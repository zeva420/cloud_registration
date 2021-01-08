#include "CloudRegister.h"


#include "BaseType.h"
#include "CADModel.h"
#include "CoarseMatching.h"
#include "TransformOptimize.h"
#include "funHelper.h"
#include "Segmentation.h"
#include "CalcHoleMeasure.h"
#include "CalcBayAndDepthMeasure.h"
#include "CalcWallVerticality.h"

#include <pcl/common/transforms.h>

namespace CloudReg {
CloudRegister::CloudRegister() {
	google::InitGoogleLogging("Cloud");
	FLAGS_log_dir = "./";

#define VISUALIZATION_ENABLED
#ifdef VISUALIZATION_ENABLED
	google::LogToStderr();
#endif
}

CloudRegister::~CloudRegister() {
	google::ShutdownGoogleLogging();
}

bool CloudRegister::run(std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr>& vecCloudPtr,
	const std::string& CAD_File) {

	if (vecCloudPtr.empty())
	{
		LOG(ERROR) << "empty cloud input";
		//return false;
	}
#if 0
	//for calc corner
	std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> vecOrigCloud;
	for (auto cloud1 : vecCloudPtr)
	{
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud2(new pcl::PointCloud<pcl::PointXYZ>());
		pcl::copyPointCloud(*cloud1, *cloud2);
		vecOrigCloud.push_back(cloud2);
	}
#endif
	CADModel model;
	model.initCAD(CAD_File);

	auto cadCloud = model.genTestFragCloud();
	auto& wall  = cadCloud[ITEM_WALL_E];
	auto vecHole = model.getTypedModelItems(ITEM_HOLE_E);
	auto vecWall = model.getTypedModelItems(ITEM_WALL_E);
	auto vecRoot = model.getTypedModelItems(ITEM_BOTTOM_E);

	std::vector<vec_seg_pair_t> allWallBorder;
	for(auto& wall : vecWall)
		allWallBorder.push_back(wall.segments_);
	std::map<std::size_t, std::vector<vec_seg_pair_t>> holeBorder;

	calcDepth(vecRoot.front().segments_, allWallBorder, holeBorder ,wall,0);
	calcDepth(vecRoot.front().segments_, allWallBorder, holeBorder ,wall,1);
	/*for(std::size_t i  = 0 ; i< vecHole.size(); i++)

	// //
	// for (std::size_t i = 0; i < vecWall.size(); i++)
	// {
	// 	auto &wal = vecWall[i];
	// 	std::vector<std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>> holeBorders;
	// 	for (std::size_t j = 0; j < vecHole.size(); j++)
	// 	{
	// 		auto& hole = vecHole[j];
	// 		if (hole.parentIndex_ == i)
	// 			holeBorders.emplace_back(hole.segments_);
	// 	}
	// 	calcVerticality(wal.segments_.back(), wal.segments_, holeBorders, wall[i], i);
	// 	// break;
	// }
	// return true;

	
	for (std::size_t i = 0; i < vecWall.size(); i++)
	{
		auto wal = vecWall[i];
		std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> holeBorders;
		std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> wallBorders;
		wallBorders = wal.segments_;
		for (std::size_t j = 0; j < vecHole.size(); j++)
		{
			auto& hole = vecHole[j];
			if (hole.parentIndex_ == i)
				holeBorders.insert(holeBorders.end(), hole.segments_.begin(), hole.segments_.end());
		}
		testRuler(wallBorders, holeBorders, wall[i], i);
		
	}
	
	
	return true;

	for(std::size_t i  = 0 ; i< vecHole.size(); i++)
	{
		//if ( i < 2) continue;
		auto& hole = vecHole[i];
		if (hole.parentIndex_ < wall.size()-1)
		{
			std::cout << "index:" << i << std::endl;
			calcHole(vecWall[hole.parentIndex_].segments_.back(), hole.segments_, wall[hole.parentIndex_]);
			break;
		}
	}*/
	return true;
	LOG(INFO) << "cad model loaded: " << model.toString() << ". from: " << CAD_File;
	
	// wall segmentation: (PointCloud, CADModel)-> [PointCloud]
	Segmentation SegmentationObj;
	if(!SegmentationObj.run(vecCloudPtr.front(), model))
	{
		LOG(INFO) << "Segmentation failed.";
		return false;
	}
	return true;

	// coarse match: ([PointCloud], CADModel)-> ([transformed & filtered PointCloud])
	CoarseMatching cm;
	auto re = cm.run(vecCloudPtr, model);
	if (!re.isValid()) {
		LOG(INFO) << "coarse matching failed.";
		return false;
	}

	// registration
	std::string logStr = "";
	TransformOptimize obj("refined Transform Opt", logStr);
	auto cloud = re.getAllPieces();
    Eigen::Vector3f tmpPt = re.T_.block<3,3>(0,0) * Eigen::Vector3f(0,0,0) + re.T_.block<3,1>(0,3);
    Eigen::Vector3d center(tmpPt(0), tmpPt(1), tmpPt(2));
	if(!obj.run(cloud, model, center))
	{
		LOG(INFO) << "transform opt failed.";
		return false;
	}

	//fill return value
	fillRet(model, obj);

	//calc corner
#if 0
	auto optRets = obj.getRet();
	if (optRets.count(TransformOptimize::CloudType::BOTTOM_E))
	{
		auto &ret = optRets[TransformOptimize::BOTTOM_E];
		for (auto cloud : vecOrigCloud)
		{
			pcl::PointCloud<pcl::PointXYZ>::Ptr transformed_cloud(new pcl::PointCloud<pcl::PointXYZ>());
			pcl::transformPointCloud(*cloud, *transformed_cloud, re.T_);
			cloud->swap(*transformed_cloud);
			transformed_cloud->clear();
			pcl::transformPointCloud(*cloud, *transformed_cloud, ret.T_);
			cloud->swap(*transformed_cloud);
		}

		center = ret.T_.block<3, 3>(0, 0) * center + ret.T_.block<3, 1>(0, 3);
		calcAllCorner(model, center, vecOrigCloud);
	}
#endif
	return true;
}

const std::map<CloudItemType, vecItems_t>& CloudRegister::getAllCloudPlane() const
{
	return mapCloudItem_;
}

const std::map<pairCloud_t, std::pair<double, double>>& CloudRegister::getAllCorner() const
{
	return mapCorner_;
}

bool CloudRegister::findSameSegment(
				const std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> &segments1,
				const std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> &segments2,
				std::pair<int, int> &idxPair)
{
	for (std::size_t i = 0; i < segments1.size(); i++)
	{
		const Eigen::Vector3d &s1 = segments1[i].first;
		const Eigen::Vector3d &e1 = segments1[i].second;
		for (std::size_t j = 0; j < segments2.size(); j++)
		{
			const Eigen::Vector3d &s2 = segments2[j].first;
			const Eigen::Vector3d &e2 = segments2[j].second;
			if ((s1 - s2).norm() < 0.0001 && (e1 - e2).norm() < 0.0001)
			{
				idxPair = std::make_pair(i, j);
				return true;
			}
			if ((s1 - e2).norm() < 0.0001 && (e1 - s2).norm() < 0.0001)
			{
				idxPair = std::make_pair(i, j);
				return true;
			}
		}
	}

	return false;
}

int CloudRegister::findMatchCloud(const Eigen::Vector4d &plane,
		std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &vecOrigCloud)
{
	double min = 20000;
	int bestIdx = -1;
	for (int i = 0; i < vecOrigCloud.size(); i++)
	{
		auto origCloud = vecOrigCloud[i];
		double aveDist = 0.0;
		for (auto &p : origCloud->points)
		{
			double dist = pointToPLaneDist(plane, p);
			aveDist += std::fabs(dist);
		}
		if (origCloud->size() > 0)
		{
			aveDist /= origCloud->size();
		}
		if (aveDist < min)
		{
			min = aveDist;
			bestIdx = i;
		}
	}
	return bestIdx;
}

void CloudRegister::calcAllCorner(CADModel& cad, Eigen::Vector3d center,
	std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &vecOrigCloud)
{
	const auto &mapCloud = getAllCloudPlane();
	const auto &it1 = mapCloud.find(CLOUD_BOTTOM_E);
	const vecItems_t &cloudbottoms = it1->second;
	double floorZ = cloudbottoms.front().pCloud_->front().z;
	Eigen::Vector4d bottomPlane = cloudbottoms.front().cloudPlane_;

	auto cadWalls = cad.getTypedModelItems(ITEM_WALL_E);
	const auto &it2 = mapCloud.find(CLOUD_WALL_E);
	const vecItems_t &cloudWalls = it2->second;
	std::vector<std::pair<int, int>> wallIdxPairs;
	for (std::size_t i = 0; i < cloudWalls.size()-1; i++)
	{
		wallIdxPairs.push_back(std::make_pair(i, i+1));
	}
	wallIdxPairs.push_back(std::make_pair(cloudWalls.size()-1, 0));

	std::vector<std::pair<double, double>> cornerPairs;
	std::vector<std::pair<double, double>> cornerPairs2;
	for (auto &it : wallIdxPairs)
	{
		int i = it.first;
		int j = it.second;
		LOG(INFO) << "---------wall pair:" << i << "-" << j;
		Eigen::VectorXd interSectionLine(6);
		if (false == interSectionOfPlaneToPlane(cloudWalls[i].cloudPlane_,
			cloudWalls[j].cloudPlane_, interSectionLine)) continue;

		Eigen::Vector3d floorPt;
		if (false == interSectionOfLineToPlane(interSectionLine,
			bottomPlane, floorPt)) continue;

		/*double v1 = calcCloudPairCorner(std::to_string(i) + "-" + std::to_string(j) + "-0.3m",
			cloudWalls[i].pCloud_, cloudWalls[j].pCloud_, floorPt, 0.3, center, bottomPlane);
		double v2 = calcCloudPairCorner(std::to_string(i) + "-" + std::to_string(j) + "-1.5m",
			cloudWalls[i].pCloud_, cloudWalls[j].pCloud_, floorPt, 1.5, center, bottomPlane);
		cornerPairs.push_back(std::make_pair(v1, v2));
		*/
		int matchIdx1 = findMatchCloud(cloudWalls[i].cadPlane_, vecOrigCloud);
		int matchIdx2 = findMatchCloud(cloudWalls[j].cadPlane_, vecOrigCloud);
		double v3 = calcCloudPairCorner(std::to_string(i) + "-" + std::to_string(j) + "-0.3m-orig",
			vecOrigCloud[matchIdx1], vecOrigCloud[matchIdx2], floorPt, 0.3, center, bottomPlane);
		double v4 = calcCloudPairCorner(std::to_string(i) + "-" + std::to_string(j) + "-1.5m-orig",
			vecOrigCloud[matchIdx1], vecOrigCloud[matchIdx2], floorPt, 1.5, center, bottomPlane);
		cornerPairs2.push_back(std::make_pair(v3, v4));
	}

	/*for (std::size_t i = 0; i < cornerPairs.size(); i++)
	{
		int idx1 = wallIdxPairs[i].first;
		int idx2 = wallIdxPairs[i].second;
		LOG(INFO) << "*****wall pair:" << idx1 << "-" << idx2
			<< ", 0.3m-1.5m: " << cornerPairs[i].first << " " << cornerPairs[i].second;
	}*/
	for (std::size_t i = 0; i < cornerPairs2.size(); i++)
	{
		int idx1 = wallIdxPairs[i].first;
		int idx2 = wallIdxPairs[i].second;
		LOG(INFO) << "*****wall pair:" << idx1 << "-" << idx2
			<< ", 0.3m-1.5m: " << cornerPairs2[i].first << " " << cornerPairs2[i].second;
	}
}

pcl::PointCloud<pcl::PointXYZI>::Ptr
CloudRegister::calcDistError(const pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud_,
	const Eigen::Vector4d& plane, const double radius) const
{
	pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud_filtered(new pcl::PointCloud<pcl::PointXYZ>());
	uniformSampling(radius, pCloud_, pCloud_filtered);

	pcl::PointCloud<pcl::PointXYZI>::Ptr pCloud_rgb(new pcl::PointCloud<pcl::PointXYZI>());
	for (auto &p : pCloud_filtered->points)
	{
		double dist = pointToPLaneDist(plane, p);

		pcl::PointXYZI p_rgb;
		p_rgb.x = p.x;
		p_rgb.y = p.y;
		p_rgb.z = p.z;
		p_rgb.intensity = dist;

		pCloud_rgb->push_back(p_rgb);
	}

	return pCloud_rgb;
}


void CloudRegister::fillRet(CADModel& cad, TransformOptimize& optimitor)
{
	mapCloudItem_.clear();
	auto optRets = optimitor.getRet();
	auto cadCloud = cad.genTestFragCloud();

	if (optRets.count(TransformOptimize::CloudType::BOTTOM_E))
	{
		auto &ret = optRets[TransformOptimize::BOTTOM_E];
		auto Botton = cad.getTypedModelItems(ITEM_BOTTOM_E).front();
		pcl::PointCloud<pcl::PointXYZ>::Ptr pData = ret.vecCloud_.front();
		CloudItem item(pData);
		item.pCADCloud_ = cadCloud[ITEM_BOTTOM_E].front();
		item.type_ = CLOUD_BOTTOM_E;
		item.cloudPlane_ = ret.vecCloudPlane_.front();
		item.cadPlane_ = ret.vecCadPlane_.front();
		item.cadBorder_.insert(item.cadBorder_.end(), Botton.segments_.begin(), Botton.segments_.end());
		LOG(INFO) << "*********************bottom************************";
		auto boundPoints = calcCloudBorder("bottom",
				pData, item.cloudPlane_, item.cadBorder_, item.cloudBorder_);
		mapCloudItem_[CLOUD_BOTTOM_E].emplace_back(item);
#ifdef VISUALIZATION_ENABLED
		pcl::io::savePCDFile("boundPoints-bottom.pcd", *boundPoints);
#endif
	}

	if (optRets.count(TransformOptimize::CloudType::TOP_E))
	{
		auto &ret = optRets[TransformOptimize::CloudType::TOP_E];
		auto Top = cad.getTypedModelItems(ITEM_TOP_E).front();
		pcl::PointCloud<pcl::PointXYZ>::Ptr pData = ret.vecCloud_.front();
		CloudItem item(pData);
		item.type_ = CLOUD_TOP_E;
		item.pCADCloud_ = cadCloud[ITEM_TOP_E].front();
		item.cloudPlane_ = ret.vecCloudPlane_.front();
		item.cadPlane_ = ret.vecCadPlane_.front();
		item.cadBorder_.insert(item.cadBorder_.end(), Top.segments_.begin(), Top.segments_.end());
		LOG(INFO) << "*********************top************************";
		auto boundPoints = calcCloudBorder("top",
				pData, item.cloudPlane_, item.cadBorder_, item.cloudBorder_);
		mapCloudItem_[CLOUD_TOP_E].emplace_back(item);
#ifdef VISUALIZATION_ENABLED
		pcl::io::savePCDFile("boundPoints-top.pcd", *boundPoints);
#endif
	}

	auto vecWall = cad.getTypedModelItems(ITEM_WALL_E);
	auto vecHole = cad.getTypedModelItems(ITEM_HOLE_E);
	if (optRets.count(TransformOptimize::CloudType::WALL_E) 
			&& vecWall.size() == optRets[TransformOptimize::CloudType::WALL_E].vecCloud_.size())
	{
		auto &ret = optRets[TransformOptimize::CloudType::WALL_E];
		for (int i = 0; i < vecWall.size(); i++)
		{
			auto& wall = vecWall[i];
			pcl::PointCloud<pcl::PointXYZ>::Ptr pData = ret.vecCloud_[i];;
			CloudItem item(pData);
			item.type_ = CLOUD_WALL_E;
			item.pCADCloud_ = cadCloud[ITEM_WALL_E][i];
			item.cloudPlane_ = ret.vecCloudPlane_[i];
			item.cadPlane_ = ret.vecCadPlane_[i];
			item.cadBorder_.insert(item.cadBorder_.end(), wall.segments_.begin(), wall.segments_.end());
			for (auto& hole : vecHole)
			{
				if (i != hole.parentIndex_) continue;
				item.cadBorder_.insert(item.cadBorder_.end(), hole.segments_.begin(), hole.segments_.end());
			}	
			LOG(INFO) << "*********************wall:" << i << "************************";
			auto boundPoints = calcCloudBorder("wall-" + std::to_string(i),
					pData, item.cloudPlane_, item.cadBorder_, item.cloudBorder_);
			mapCloudItem_[CLOUD_WALL_E].emplace_back(item);		
#ifdef VISUALIZATION_ENABLED
			pcl::io::savePCDFile("boundPoints-wall-" + std::to_string(i) + ".pcd", *boundPoints);
#endif
		}
	}
	
	auto vecBeam = cad.getTypedModelItems(ITEM_BEAM_E);
	if (optRets.count(TransformOptimize::CloudType::BEAM_E) 
			&& vecWall.size() == optRets[TransformOptimize::CloudType::BEAM_E].vecCloud_.size())
	{
		auto &ret = optRets[TransformOptimize::CloudType::BEAM_E];
		for (int i = 0; i < vecBeam.size(); i++)
		{
			auto& beam = vecBeam[i];
			pcl::PointCloud<pcl::PointXYZ>::Ptr pData = ret.vecCloud_[i];
			CloudItem item(pData);
			item.type_ = CLOUD_BEAM_E;
			item.pCADCloud_ = cadCloud[ITEM_BEAM_E][i];
			item.parentIndex_ = item.parentIndex_;
			item.cloudPlane_ = ret.vecCloudPlane_[i];
			item.cadPlane_ = ret.vecCadPlane_[i];
			item.cadBorder_.insert(item.cadBorder_.end(), beam.segments_.begin(), beam.segments_.end());
			LOG(INFO) << "*********************beam:" << i << "************************";
			auto boundPoints = calcCloudBorder("beam-" + std::to_string(i),
					pData, item.cloudPlane_, item.cadBorder_, item.cloudBorder_);
			mapCloudItem_[CLOUD_BEAM_E].emplace_back(item);
#ifdef VISUALIZATION_ENABLED
			pcl::io::savePCDFile("boundPoints-beam-" + std::to_string(i) + ".pcd", *boundPoints);
#endif
		}
	}
}

}
