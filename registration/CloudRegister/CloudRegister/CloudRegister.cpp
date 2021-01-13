#include "CloudRegister.h"


#include "BaseType.h"
#include "CADModel.h"
#include "CoarseMatching.h"
#include "TransformOptimize.h"
#include "funHelper.h"
#include "CloudSegment.h"

#include <pcl/common/transforms.h>

namespace CloudReg {
CloudRegister::CloudRegister() {
	google::InitGoogleLogging("Cloud");
	FLAGS_log_dir = "./";
	
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
		return false;
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

	LOG(INFO) << "cad model loaded: " << model.toString() << ". from: " << CAD_File;
	
	// wall segmentation: (PointCloud, CADModel)-> [PointCloud]

	CloudSegment cs(vecCloudPtr.front(), model);
	auto sr = cs.run();
	if(!sr.valid()){
		LOG(WARNING)<< "failed to segment cloud.";
		return false;
	}
#if 0
	return true;

	// coarse match: ([PointCloud], CADModel)-> ([transformed & filtered PointCloud])
	CoarseMatching cm;

	cm.segment(vecCloudPtr.front(), model);
	return true;

	auto re = cm.run(vecCloudPtr, model);
	if (!re.isValid()) {
		LOG(INFO) << "coarse matching failed.";
		return false;
	}
#endif
	// registration
	std::string logStr = "";
	TransformOptimize obj("refined Transform Opt", logStr);
	std::map<ModelItemType, std::vector<PointCloud::Ptr> > clouds;
	for (auto &it : sr.clouds_)
	{
		clouds[it.first] = std::vector<PointCloud::Ptr>();
		for (auto &piece : it.second)
		{
			clouds[it.first].push_back(piece.cloud_);
		}
	}
    Eigen::Vector3f tmpPt = sr.T_.block<3,3>(0,0) * Eigen::Vector3f(0,0,0) + sr.T_.block<3,1>(0,3);
    Eigen::Vector3d center(tmpPt(0), tmpPt(1), tmpPt(2));
	if(!obj.run(clouds, model, center))
	{
		LOG(INFO) << "transform opt failed.";
		return false;
	}

	//fill return value
	fillRet(model, obj);

	//calc corner
#if 0
	auto optRets = obj.getRet();
	if (optRets.valid())
	{
		for (auto cloud : vecOrigCloud)
		{
			pcl::PointCloud<pcl::PointXYZ>::Ptr transforCloud(new pcl::PointCloud<pcl::PointXYZ>());
			pcl::transformPointCloud(*cloud, *transforCloud, sr.T_);
			cloud->swap(*transforCloud);
			transforCloud->clear();
			pcl::transformPointCloud(*cloud, *transforCloud, optRets.T_);
			cloud->swap(*transforCloud);
		}

		center = optRets.T_.block<3, 3>(0, 0) * center + optRets.T_.block<3, 1>(0, 3);
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

void CloudRegister::calcAllCloudBorder(CADModel& cad)
{
	LOG(INFO) << "*******************calcAllCloudBorder*********************";
	std::vector<Eigen::Vector4d> cloudPlaneVec;
	std::vector<std::vector<Eigen::Vector3d>> cadSegPtsOfPlanes;
	auto pushNeedData = [&](ModelItemType cadType, CloudItemType cloudType)->void {
		if (mapCloudItem_.count(cloudType))
		{
			const vecItems_t &vecCloudItems = mapCloudItem_[cloudType];
			auto vecCadItems = cad.getTypedModelItems(cadType);
			if (vecCadItems.size() !=  vecCloudItems.size()) return;

			for (int i = 0; i < vecCloudItems.size(); i++)
			{
				std::vector<Eigen::Vector3d> segPoints;
				for (auto &seg : vecCadItems[i].segments_)
				{
					segPoints.push_back(seg.first);
				}				
				cloudPlaneVec.push_back(vecCloudItems[i].cloudPlane_);
				cadSegPtsOfPlanes.push_back(segPoints);
			}
		}
	};
	pushNeedData(ITEM_WALL_E, CLOUD_WALL_E);
	pushNeedData(ITEM_BOTTOM_E, CLOUD_BOTTOM_E);
	pushNeedData(ITEM_TOP_E, CLOUD_TOP_E);
	pushNeedData(ITEM_BEAM_E, CLOUD_BEAM_E);
	LOG(INFO) << "cloudPlaneVec:" << cloudPlaneVec.size() << " cadSegPtsOfPlanes:" << cadSegPtsOfPlanes.size();

	//get 3 planes' idx which share same cad seg point
	std::set<std::set<int>> planeIdxGroup;
	if (false == groupPlanesBySamePt(cadSegPtsOfPlanes, planeIdxGroup))
	{
		LOG(ERROR) << "groupPlanesBySamePt Failed.";
		return;
	}
	LOG(INFO) << "planeIdxGroup:" << planeIdxGroup.size();

	//calc the interSection point of each 3 planes
	std::vector<Eigen::Vector3d> focalPointVec;
	if (false == interSectionOf3Planes(cloudPlaneVec, planeIdxGroup, focalPointVec))
	{
		LOG(ERROR) << "interSectionOf3Planes Failed.";
		return;
	}


	auto findNearestSegOfNodeOrCloudPt = [&](const std::vector<Eigen::Vector3d> &pCloud, 
			const std::vector<Eigen::Vector3d> &vecNodes, 
			const std::pair<Eigen::Vector3d, Eigen::Vector3d> &seg)
			->std::pair<Eigen::Vector3d, Eigen::Vector3d> {
		auto ret1 = findNearestSeg(vecNodes, seg);
		if (ret1.first < 0.3)
		{
			LOG(INFO) << "find node seg, dist1:" <<ret1.first;
			return ret1.second;
		}

		auto ret2 = findNearestSeg(pCloud, seg);
		LOG(INFO) << "empty, to find cloud seg, dist1:" <<ret2.first;
		return ret2.second;
	};

	//find nearest_pt (in intersection pts, or cloud pts) for each cadBorder
	std::vector<std::string> typeNames{"Beam","Bottom","Wall","Top","Unknow"};
	for (auto &it : mapCloudItem_)
	{
		std::string name = typeNames[it.first];
		auto &vecItems = it.second;
		for (int i = 0; i < vecItems.size(); i++)
		{
			LOG(INFO) << "------get cloud border for " << name << "-" << i << "------";
			auto &item = vecItems[i];
			std::vector<Eigen::Vector3d> vecNodes;
			if (CLOUD_WALL_E == it.first && item.cadBorder_.size() > 1)
			{
				LOG(INFO) << "wall-" << i << " has holes";
				vecNodes = calcWallNodes(name + std::to_string(i), item.pCloud_, item.cloudPlane_);
			}

			auto cloudPts = convertCloudToEigenVec(item.pCloud_);
			for (int k = 0; k < item.cadBorder_.size(); k++)
			{
				const auto &vecSegs = item.cadBorder_[k];
				std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> cloudSegs;
				if (0 == k)
				{
					for (const auto &seg : vecSegs)
					{
						std::pair<Eigen::Vector3d, Eigen::Vector3d> bestSeg;
						LOG(INFO) << "+++is outer contour seg";
						bestSeg = findNearestSegOfNodeOrCloudPt(cloudPts, focalPointVec, seg);
						cloudSegs.push_back(bestSeg);
					}
				}
				else
				{
					for (const auto &seg : vecSegs)
					{
						std::pair<Eigen::Vector3d, Eigen::Vector3d> bestSeg;
						LOG(INFO) << "---is hole seg";
						bestSeg = findNearestSegOfNodeOrCloudPt(cloudPts, vecNodes, seg);
						cloudSegs.push_back(bestSeg);
					}
				}
				item.cloudBorder_.push_back(cloudSegs);
			}
			if (item.cadBorder_.size() != item.cloudBorder_.size())
			{
				LOG(ERROR) << "current item cadBorder_.size != cloudBorder_.size";
			}
#ifdef VISUALIZATION_ENABLED
			if (!vecNodes.empty())
			{
				pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud(new pcl::PointCloud<pcl::PointXYZ>());
				for (size_t i = 0; i < vecNodes.size(); ++i)
				{
					pcl::PointXYZ p(vecNodes[i](0), vecNodes[i](1), vecNodes[i](2));
					pCloud->push_back(p);
				}			
				std::string file_name = "vecNodes-" + name + "-" + std::to_string(i)  + ".pcd";
				pcl::io::savePCDFile(file_name, *pCloud);	
			}

			{
				std::vector<Eigen::Vector3d> vecPoints;
				for (auto& segVec : item.cadBorder_) {
					for (auto& pt_pair : segVec) {
						auto vec_tmp = ininterpolateSeg(pt_pair.first, pt_pair.second, 0.05);
						vecPoints.insert(vecPoints.end(), vec_tmp.begin(), vec_tmp.end());
					}
				}
				pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud(new pcl::PointCloud<pcl::PointXYZ>());
				for (size_t i = 0; i < vecPoints.size(); ++i)
				{
					pcl::PointXYZ p(vecPoints[i](0), vecPoints[i](1), vecPoints[i](2));
					pCloud->push_back(p);
				}			
				std::string file_name = "cadSegs-plane-" + name + "-" + std::to_string(i) + ".pcd";
				pcl::io::savePCDFile(file_name, *pCloud);	
			}
			{
				std::vector<Eigen::Vector3d> vecPoints;
				for (auto& segVec : item.cloudBorder_) {
					for (auto& pt_pair : segVec) {
						auto vec_tmp = ininterpolateSeg(pt_pair.first, pt_pair.second, 0.05);
						vecPoints.insert(vecPoints.end(), vec_tmp.begin(), vec_tmp.end());
					}
				}
				pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud(new pcl::PointCloud<pcl::PointXYZ>());
				for (size_t i = 0; i < vecPoints.size(); ++i)
				{
					pcl::PointXYZ p(vecPoints[i](0), vecPoints[i](1), vecPoints[i](2));
					pCloud->push_back(p);
				}			
				std::string file_name = "cloudSegs-plane-" + name + "-" + std::to_string(i) + ".pcd";
				pcl::io::savePCDFile(file_name, *pCloud);	
			}
#endif			
		}		
	}

	LOG(INFO) << "***************************************";
	return;
}

void CloudRegister::fillRet(CADModel& cad, TransformOptimize& optimitor)
{
	mapCloudItem_.clear();
	auto optRets = optimitor.getRet();
	auto cadCloud = cad.genTestFragCloud();

	if (optRets.mapClouds_.count(ITEM_BOTTOM_E))
	{
		auto &ret = optRets.mapClouds_[ITEM_BOTTOM_E];
		auto Botton = cad.getTypedModelItems(ITEM_BOTTOM_E).front();
		pcl::PointCloud<pcl::PointXYZ>::Ptr pData = ret.front().cloud_;
		CloudItem item(pData);
		item.pCADCloud_ = cadCloud[ITEM_BOTTOM_E].front();
		item.type_ = CLOUD_BOTTOM_E;
		item.cloudPlane_ = ret.front().cloudPlane_;
		item.cadPlane_ = ret.front().cadPlane_;
		item.cadBorder_.push_back(Botton.segments_);
		mapCloudItem_[CLOUD_BOTTOM_E].emplace_back(item);
	}

	if (optRets.mapClouds_.count(ITEM_TOP_E))
	{
		auto &ret = optRets.mapClouds_[ITEM_TOP_E];
		auto Top = cad.getTypedModelItems(ITEM_TOP_E).front();
		pcl::PointCloud<pcl::PointXYZ>::Ptr pData = ret.front().cloud_;
		CloudItem item(pData);
		item.type_ = CLOUD_TOP_E;
		item.pCADCloud_ = cadCloud[ITEM_TOP_E].front();
		item.cloudPlane_ = ret.front().cloudPlane_;
		item.cadPlane_ = ret.front().cadPlane_;
		item.cadBorder_.push_back(Top.segments_);
		mapCloudItem_[CLOUD_TOP_E].emplace_back(item);
	}

	auto vecWall = cad.getTypedModelItems(ITEM_WALL_E);
	auto vecHole = cad.getTypedModelItems(ITEM_HOLE_E);
	if (optRets.mapClouds_.count(ITEM_WALL_E) 
			&& vecWall.size() == optRets.mapClouds_[ITEM_WALL_E].size())
	{
		auto &ret = optRets.mapClouds_[ITEM_WALL_E];
		for (int i = 0; i < vecWall.size(); i++)
		{
			auto& wall = vecWall[i];
			pcl::PointCloud<pcl::PointXYZ>::Ptr pData = ret[i].cloud_;
			CloudItem item(pData);
			item.type_ = CLOUD_WALL_E;
			item.pCADCloud_ = cadCloud[ITEM_WALL_E][i];
			item.cloudPlane_ = ret[i].cloudPlane_;
			item.cadPlane_ = ret[i].cadPlane_;
			item.cadBorder_.push_back(wall.segments_);
			for (auto& hole : vecHole)
			{
				if (i != hole.parentIndex_) continue;
				item.cadBorder_.push_back(hole.segments_);
			}	
			mapCloudItem_[CLOUD_WALL_E].emplace_back(item);		
		}
	}
	
	auto vecBeam = cad.getTypedModelItems(ITEM_BEAM_E);
	if (optRets.mapClouds_.count(ITEM_BEAM_E) 
			&& vecWall.size() == optRets.mapClouds_[ITEM_BEAM_E].size())
	{
		auto &ret = optRets.mapClouds_[ITEM_BEAM_E];
		for (int i = 0; i < vecBeam.size(); i++)
		{
			auto& beam = vecBeam[i];
			pcl::PointCloud<pcl::PointXYZ>::Ptr pData = ret[i].cloud_;
			CloudItem item(pData);
			item.type_ = CLOUD_BEAM_E;
			item.pCADCloud_ = cadCloud[ITEM_BEAM_E][i];
			item.parentIndex_ = item.parentIndex_;
			item.cloudPlane_ = ret[i].cloudPlane_;
			item.cadPlane_ = ret[i].cadPlane_;
			item.cadBorder_.push_back(beam.segments_);
			mapCloudItem_[CLOUD_BEAM_E].emplace_back(item);
		}
	}

	calcAllCloudBorder(cad);
}

}