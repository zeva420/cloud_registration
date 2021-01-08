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

	auto isHoleSeg = [](const std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> &holeSegs, 
				const std::pair<Eigen::Vector3d, Eigen::Vector3d> &seg)->bool {
		for (int i = 0; i < holeSegs.size(); i++)
		{
			auto &holeSeg = holeSegs[i];
			bool flag1 = ((seg.first - holeSeg.first).norm() < 1e-6) 
							&& ((seg.second - holeSeg.second).norm() < 1e-6);
			bool flag2 = ((seg.first - holeSeg.second).norm() < 1e-6) 
							&& ((seg.second - holeSeg.first).norm() < 1e-6);
			if (flag1 || flag2)
			{
				return true;
			}
		}
		return false;
	};
	std::set<int> wallIdxWithHole;
	std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> holeSegs;
	auto vecHole = cad.getTypedModelItems(ITEM_HOLE_E);
	for (auto& hole : vecHole)
	{
		if (!hole.segments_.empty())
		{
			wallIdxWithHole.insert(hole.parentIndex_);
			holeSegs.insert(holeSegs.end(), hole.segments_.begin(), hole.segments_.end());
		}
	}

	auto findNearestNodeOrCloudPt = [&](const std::vector<Eigen::Vector3d> &pCloud, 
			const std::vector<Eigen::Vector3d> &vecNodes, Eigen::Vector3d &pt)->Eigen::Vector3d {
		LOG(INFO) << "+++is hole seg";
		auto ret1 = findNearestPt(vecNodes, pt);
		if (ret1.first < 0.3)
		{
			LOG(INFO) << "find node, dist1:" <<ret1.first;
			return ret1.second;
		}
		auto ret2 = findNearestPt(pCloud, pt);
		LOG(INFO) << "empty, to find nearest pt, dist1:" <<ret2.first;
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
			if (CLOUD_WALL_E == it.first && wallIdxWithHole.count(i))
			{
				LOG(INFO) << "wall-" << i << " has holes";
				vecNodes = calcWallNodes(name + std::to_string(i), item.pCloud_, item.cloudPlane_);
			}

			auto cloudPts = convertCloudToEigenVec(item.pCloud_);
			for (auto seg : item.cadBorder_)
			{
				Eigen::Vector3d p1, p2;
				if (isHoleSeg(holeSegs, seg) && !vecNodes.empty()) 
				{
					LOG(INFO) << "+++is hole seg";
					p1 = findNearestNodeOrCloudPt(cloudPts, vecNodes, seg.first);
					p2 = findNearestNodeOrCloudPt(cloudPts, vecNodes, seg.second);
				}
				else
				{
					LOG(INFO) << "+++is outer contour seg";
					p1 = findNearestNodeOrCloudPt(cloudPts, focalPointVec, seg.first);
					p2 = findNearestNodeOrCloudPt(cloudPts, focalPointVec, seg.second);
				}
				item.cloudBorder_.push_back(std::make_pair(p1, p2));
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
				for (auto& pt_pair : item.cadBorder_) {
					auto vec_tmp = ininterpolateSeg(pt_pair.first, pt_pair.second, 0.05);
					vecPoints.insert(vecPoints.end(), vec_tmp.begin(), vec_tmp.end());
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
				for (auto& pt_pair : item.cloudBorder_) {
					auto vec_tmp = ininterpolateSeg(pt_pair.first, pt_pair.second, 0.05);
					vecPoints.insert(vecPoints.end(), vec_tmp.begin(), vec_tmp.end());
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
		// LOG(INFO) << "*********************bottom************************";
		// auto boundPoints = calcCloudBorder("bottom",
		// 		pData, item.cloudPlane_, item.cadBorder_, item.cloudBorder_);
		mapCloudItem_[CLOUD_BOTTOM_E].emplace_back(item);
// #ifdef VISUALIZATION_ENABLED
// 		pcl::io::savePCDFile("boundPoints-bottom.pcd", *boundPoints);
// #endif
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
		// LOG(INFO) << "*********************top************************";
		// auto boundPoints = calcCloudBorder("top",
		// 		pData, item.cloudPlane_, item.cadBorder_, item.cloudBorder_);
		mapCloudItem_[CLOUD_TOP_E].emplace_back(item);
// #ifdef VISUALIZATION_ENABLED
// 		pcl::io::savePCDFile("boundPoints-top.pcd", *boundPoints);
// #endif
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
			// LOG(INFO) << "*********************wall:" << i << "************************";
			// auto boundPoints = calcCloudBorder("wall-" + std::to_string(i),
			// 		pData, item.cloudPlane_, item.cadBorder_, item.cloudBorder_);
			mapCloudItem_[CLOUD_WALL_E].emplace_back(item);		
// #ifdef VISUALIZATION_ENABLED
// 			pcl::io::savePCDFile("boundPoints-wall-" + std::to_string(i) + ".pcd", *boundPoints);
// #endif
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
// 			LOG(INFO) << "*********************beam:" << i << "************************";
// 			auto boundPoints = calcCloudBorder("beam-" + std::to_string(i),
// 					pData, item.cloudPlane_, item.cadBorder_, item.cloudBorder_);
			mapCloudItem_[CLOUD_BEAM_E].emplace_back(item);
// #ifdef VISUALIZATION_ENABLED
// 			pcl::io::savePCDFile("boundPoints-beam-" + std::to_string(i) + ".pcd", *boundPoints);
// #endif
		}
	}

	calcAllCloudBorder(cad);
}

}