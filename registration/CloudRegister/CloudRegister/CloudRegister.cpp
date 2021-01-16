#include "CloudRegister.h"


#include "BaseType.h"
#include "funHelper.h"
#include "CalcMeasureHelper.h"
#include "CADModel.h"
#include "TransformOptimize.h"
#include "CloudSegment.h"
#include "CalcNetHeight.h"
#include "CalcBayAndDepthMeasure.h"
#include "CalcHoleMeasure.h"

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
		return false;
	}

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

	return true;
}

const std::map<CloudItemType, vecItems_t>& CloudRegister::getAllCloudPlane() const
{
	return mapCloudItem_;
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

	const double distTh = 0.3;
	auto findNearestSegOfNodeOrCloudPt = [&](const std::vector<Eigen::Vector3d> &pCloud, 
			const std::vector<Eigen::Vector3d> &vecNodes, 
			const seg_pair_t &seg)
			->std::pair<double, seg_pair_t> {
		auto ret1 = findNearestSeg(vecNodes, seg);
		if (ret1.first < distTh)
		{
			LOG(INFO) << "find node seg, dist1:" <<ret1.first;
			return ret1;
		}

		auto ret2 = findNearestSeg(pCloud, seg);
		LOG(INFO) << "empty, to find cloud seg, dist1:" <<ret2.first;
		return ret2;
	};

	//find nearest_pt (in intersection pts, or cloud pts) for each outer cadBorder
	std::vector<std::string> typeNames{"Beam","Bottom","Wall","Top","Unknow"};
	for (auto &it : mapCloudItem_)
	{
		std::string name = typeNames[it.first];
		auto &vecItems = it.second;
		for (int i = 0; i < vecItems.size(); i++)
		{
			LOG(INFO) << "------get outer border for " << name << "-" << i << "------";
			auto &item = vecItems[i];
			auto cloudPts = convertCloudToEigenVec(item.pCloud_);

			//front is the wall outer border
			const auto &vecSegs = item.cadBorder_.front();
			std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> cloudSegs;
			for (const auto &seg : vecSegs)
			{
				std::pair<double, std::pair<Eigen::Vector3d, Eigen::Vector3d>> bestSeg;
				bestSeg = findNearestSegOfNodeOrCloudPt(cloudPts, focalPointVec, seg);
				if (bestSeg.first < distTh) cloudSegs.push_back(bestSeg.second);
			}
			if (cloudSegs.size() != vecSegs.size())
				LOG(WARNING) << "the border size miss match: " << name;


			item.cloudBorder_.push_back(cloudSegs);

			//refine cut
			{	
				//if (it.first == CLOUD_WALL_E)
				{
					std::vector<Eigen::Vector3d> vecPts;
					for (auto& seg : item.cloudBorder_.front())
					{
						vecPts.emplace_back(seg.first);												
					}
					
					auto pNew = filerCloudByConvexHull(item.pCloud_, vecPts);
					item.pCloud_->swap(*pNew);
				}
			}
		}		
	}

	//find nearest_pt (in vecNodes pts, or cloud pts) for each hole cadBorder
	for (auto &it : mapCloudItem_)
	{
		std::string name = typeNames[it.first];
		auto &vecItems = it.second;
		for (int i = 0; i < vecItems.size(); i++)
		{
			auto &item = vecItems[i];
			//if no holes, skip
			if (CLOUD_WALL_E != it.first || item.cadBorder_.size() <= 1) continue;

			LOG(INFO) << "------get hole border for " << name << "-" << i << "------";
			const auto &outerSegs = item.cloudBorder_.front();

			std::vector<Eigen::Vector3d> vecNodes;
			vecNodes = calcWallNodes(name + std::to_string(i), item.pCloud_, item.cloudPlane_, outerSegs);
		
			//find hole segs
			auto cloudPts = convertCloudToEigenVec(item.pCloud_);
			for (int k = 0; k < item.cadBorder_.size(); k++)
			{
				if (0 == k) continue;
				const auto &vecSegs = item.cadBorder_[k];
				std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> cloudSegs;
				for (const auto &seg : vecSegs)
				{
					std::pair<double, std::pair<Eigen::Vector3d, Eigen::Vector3d>> bestSeg;
					bestSeg = findNearestSegOfNodeOrCloudPt(cloudPts, vecNodes, seg);
					if (bestSeg.first < distTh) cloudSegs.push_back(bestSeg.second);
				}
				
				if (cloudSegs.size() != vecSegs.size())
					LOG(WARNING) << "the border size miss match: " << name;

				item.cloudBorder_.push_back(cloudSegs);
			}
			if (item.cadBorder_.size() != item.cloudBorder_.size())
			{
				LOG(ERROR) << "current item cadBorder_.size:" << item.cadBorder_.size()
					 << " != cloudBorder_.size:" << item.cloudBorder_.size();
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
#endif			
		}		
	}

#ifdef VISUALIZATION_ENABLED
	for (const auto &it : mapCloudItem_)
	{
		std::string name = typeNames[it.first];
		const auto &vecItems = it.second;
		for (int i = 0; i < vecItems.size(); i++)
		{
			const auto &item = vecItems[i];
			//cad border
			{
				std::vector<Eigen::Vector3d> vecPoints;
				for (auto& segVec : item.cadBorder_) {
					for (auto& pt_pair : segVec) {
						auto vec_tmp = ininterpolateSeg(pt_pair.first, pt_pair.second, 0.05);
						vecPoints.insert(vecPoints.end(), vec_tmp.begin(), vec_tmp.end());
					}
				}
				pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud(new pcl::PointCloud<pcl::PointXYZ>());
				for (size_t j = 0; j < vecPoints.size(); ++j)
				{
					pcl::PointXYZ p(vecPoints[j](0), vecPoints[j](1), vecPoints[j](2));
					pCloud->push_back(p);
				}			
				std::string file_name = "cadSegs-plane-" + name + "-" + std::to_string(i) + ".pcd";
				pcl::io::savePCDFile(file_name, *pCloud);	
			}
			//cloud outer border
			{
				std::vector<Eigen::Vector3d> vecPoints;
				for (auto& segVec : item.cloudBorder_) {
					for (auto& pt_pair : segVec) {
						auto vec_tmp = ininterpolateSeg(pt_pair.first, pt_pair.second, 0.05);
						vecPoints.insert(vecPoints.end(), vec_tmp.begin(), vec_tmp.end());
					}
				}
				pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud(new pcl::PointCloud<pcl::PointXYZ>());
				for (size_t j = 0; j < vecPoints.size(); ++j)
				{
					pcl::PointXYZ p(vecPoints[j](0), vecPoints[j](1), vecPoints[j](2));
					pCloud->push_back(p);
				}			
				std::string file_name = "cloudSegs-plane-" + name + "-" + std::to_string(i) + ".pcd";
				pcl::io::savePCDFile(file_name, *pCloud);	
			}
			//hole border
			if (CLOUD_WALL_E == it.first && item.cadBorder_.size() > 1)
			{
				std::vector<Eigen::Vector3d> vecPoints;
				for (int k = 0; k < item.cloudBorder_.size(); k++) {
					if (0 == k) continue;
					const auto& segVec = item.cloudBorder_[k]; 
					for (auto& pt_pair : segVec) {
						auto vec_tmp = ininterpolateSeg(pt_pair.first, pt_pair.second, 0.05);
						vecPoints.insert(vecPoints.end(), vec_tmp.begin(), vec_tmp.end());
					}
				}
				pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud(new pcl::PointCloud<pcl::PointXYZ>());
				for (size_t j = 0; j < vecPoints.size(); ++j)
				{
					pcl::PointXYZ p(vecPoints[j](0), vecPoints[j](1), vecPoints[j](2));
					pCloud->push_back(p);
				}			
				std::string file_name = "holeSegs-plane-" + name + "-" + std::to_string(i) + ".pcd";
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
	auto cadCloud = cad.genFragCloud();

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
	
	if (optRets.mapClouds_.count(ITEM_BEAM_E))
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

 			//LOG(INFO) << "*********************beam:" << i << "************************";
// 			auto boundPoints = calcCloudBorder("beam-" + std::to_string(i),
// 					pData, item.cloudPlane_, item.cadBorder_, item.cloudBorder_);
			mapCloudItem_[CLOUD_BEAM_E].emplace_back(item);
		}
	}

	calcAllCloudBorder(cad);
}

std::tuple<std::vector<calcMeassurment_t>, std::vector<seg_pair_t>>
CloudRegister::calcRoofNetHeight(const double calcLengthTh)
{
	const auto& itemRoof = mapCloudItem_[CLOUD_TOP_E].front();
	const auto& itemRoot = mapCloudItem_[CLOUD_BOTTOM_E].front();

	auto vecRet = CalcNetHeight(itemRoof.cloudBorder_.front(),itemRoof.pCloud_,
		itemRoot.cloudPlane_,"roof_net_height.pcd", calcLengthTh);

	return vecRet;
}

//first roof second root
std::tuple<std::vector<calcMeassurment_t>, std::vector<calcMeassurment_t>, std::vector<seg_pair_t>>
CloudRegister::calcPlaneRange(const double calcHeight, const double calcLengthTh)
{
	const auto& itemRoof = mapCloudItem_[CLOUD_TOP_E].front();
	const auto& itemRoot = mapCloudItem_[CLOUD_BOTTOM_E].front();
	const auto& itemWall = mapCloudItem_[CLOUD_WALL_E];
	std::vector<std::vector<seg_pair_t>> allWallBorder;
	for (auto& item : itemWall)
	{
		allWallBorder.emplace_back(item.cloudBorder_.front());
	}

	auto vecRet = CalcHeightRange(itemRoof.cloudBorder_.front(), itemRoot.cloudBorder_.front(),
		allWallBorder,itemRoof.pCloud_, itemRoot.pCloud_, calcHeight, calcLengthTh);

	return vecRet;
}

std::tuple<std::map<std::pair<std::size_t, std::size_t>,
	std::vector<calcMeassurment_t>>, std::vector<seg_pair_t>>
CloudRegister::calcDepth(const double calcLengthTh)
{

	const auto& itemRoot = mapCloudItem_[CLOUD_BOTTOM_E].front();
	const auto& itemWall = mapCloudItem_[CLOUD_WALL_E];
	std::vector<std::vector<seg_pair_t>> allWallBorder;
	std::map<std::size_t, std::vector<vec_seg_pair_t>> holeBorder;
	std::vector<PointCloud::Ptr> vecCloud;
	Eigen::vector<Eigen::Vector4d> vecPlane;
	for(std::size_t i = 0; i < itemWall.size(); i++)
	{
		const auto& item = itemWall[i];
		allWallBorder.emplace_back(item.cloudBorder_.front());
		vecCloud.emplace_back(item.pCloud_);
		vecPlane.emplace_back(item.cloudPlane_);

		for (std::size_t j = 1; j < item.cloudBorder_.size(); j++)
		{
			if (item.cloudBorder_[j].size() != item.cadBorder_[j].size())
				LOG(ERROR) << "cloudBorder error, need check";

			holeBorder[i].emplace_back(item.cloudBorder_[j]);
		}
	}
	
	auto ret = calcDepthorBay(itemRoot.cloudBorder_.front(), allWallBorder, holeBorder,vecCloud, vecPlane,0, calcLengthTh);

	return ret;

}

std::tuple<std::map<std::pair<std::size_t, std::size_t>,
	std::vector<calcMeassurment_t>>, std::vector<seg_pair_t>>
	CloudRegister::calcBay(const double calcLengthTh)
{
	const auto& itemRoot = mapCloudItem_[CLOUD_BOTTOM_E].front();
	const auto& itemWall = mapCloudItem_[CLOUD_WALL_E];

	std::vector<std::vector<seg_pair_t>> allWallBorder;
	std::map<std::size_t, std::vector<vec_seg_pair_t>> holeBorder;
	std::vector<PointCloud::Ptr> vecCloud;
	Eigen::vector<Eigen::Vector4d> vecPlane;
	for (std::size_t i = 0; i < itemWall.size(); i++)
	{
		const auto& item = itemWall[i];
		allWallBorder.emplace_back(item.cloudBorder_.front());
		vecCloud.emplace_back(item.pCloud_);
		vecPlane.emplace_back(item.cloudPlane_);

		for (std::size_t j = 1; j < item.cloudBorder_.size(); j++)
		{
			if (item.cloudBorder_[j].size() != item.cadBorder_[j].size())
				LOG(ERROR) << "cloudBorder error, need check";

			holeBorder[i].emplace_back(item.cloudBorder_[j]);


		}
	}

	auto ret = calcDepthorBay(itemRoot.cloudBorder_.front(), allWallBorder, holeBorder, vecCloud, vecPlane, 1, calcLengthTh);

	return ret;
}

std::map<std::pair<std::size_t, std::size_t>,std::vector<calcMeassurment_t>> 
CloudRegister::calcAllHole()
{
	std::map<std::pair<std::size_t, std::size_t>, std::vector<calcMeassurment_t>> mapRet;
	const auto& itemWall = mapCloudItem_[CLOUD_WALL_E];
	for (std::size_t i = 0; i < itemWall.size(); i++)
	{
		//if (i != 1) continue;

		const auto& item = itemWall[i];
		for (std::size_t j = 1; j < item.cloudBorder_.size(); j++)
		{
			//if (j != 3) continue;

			if (item.cloudBorder_[j].size() != item.cadBorder_[j].size())
				LOG(ERROR) << "cloudBorder error, need check";

			LOG(INFO) << "calc hole:" << j << " in wall " << i;
			std::string name = "hole_" + std::to_string(i) + "_" + std::to_string(j) + ".pcd";
			auto ret = calcHole(item.cloudBorder_.front().back(),
				item.cloudBorder_[j], item.pCloud_,name);

			if (!ret.empty())
			{
				mapRet[std::make_pair(i, j)] = ret;
			}
		}

	}
	return mapRet;
}
}