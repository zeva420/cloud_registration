#include "CloudRegister.h"


#include "BaseType.h"
#include "CADModel.h"
#include "CoarseMatching.h"

#include "funHelper.h"

#include <pcl/visualization/pcl_visualizer.h>

namespace CloudReg {
CloudRegister::CloudRegister() {
	google::InitGoogleLogging("Cloud");
	FLAGS_log_dir = "./log";

	google::LogToStderr();
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
	// scale model to meters
	model.scaleModel(0.001);

	// wall segmentation: (PointCloud, CADModel)-> [PointCloud]

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
	auto optRets = obj.run(cloud, model);
	if (optRets.empty()) {
		LOG(INFO) << "transform opt failed.";
		return false;
	}

	//rescale model
	//model.scaleModel(1000.0);

	//fill return value
	fillRet(model, optRets);

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

pcl::PointCloud<pcl::PointXYZRGB>::Ptr
CloudRegister::calcDistError(const pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud_,
	const Eigen::Vector4d& plane, const double downRatio) const
{
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr pCloud_rgb(new pcl::PointCloud<pcl::PointXYZRGB>());
	for (auto &p : pCloud_->points)
	{
		double dist = pointToPLaneDist(plane, p);

		pcl::PointXYZRGB p_rgb;
		p_rgb.x = p.x;
		p_rgb.y = p.y;
		p_rgb.z = p.z;
		p_rgb.r = dist;
		p_rgb.g = 0;
		p_rgb.b = 0;

		pCloud_rgb->push_back(p_rgb);
	}

	return pCloud_rgb;
}


void CloudRegister::fillRet(CADModel& cad, TransformOptimize::optCloudRets &optRets)
{
	pcl::visualization::PCLVisualizer viewer("dist");
	mapCloudItem_.clear();
	

	if (optRets.count(TransformOptimize::CloudType::BOTTOM_E))
	{
		LOG(INFO) << "************BOTTOM_E*****";
		auto &ret = optRets[TransformOptimize::BOTTOM_E];
		auto Botton = cad.getTypedModelItems(ITEM_BOTTOM_E).front();
		pcl::PointCloud<pcl::PointXYZ>::Ptr pData = ret.vecCloud_.front();
		CloudItem item(pData);
		item.type_ = CLOUD_BOTTOM_E;
		item.cloudPlane_ = ret.vecCloudPlane_.front();
		item.cadPlane_ = ret.vecCadPlane_.front();
		item.cadBorder_.insert(item.cadBorder_.begin(), Botton.segments_.begin(), Botton.segments_.end());
		auto boundPoints = calcCloudBorder(pData, item.cloudPlane_, item.cadBorder_, item.cloudBorder_);
		mapCloudItem_[CLOUD_BOTTOM_E].emplace_back(item);

		auto cloud_rgb = calcDistError(pData, item.cloudPlane_, 1.0);
		viewer.addPointCloud(boundPoints, "bottom-boundPoints");
		viewer.addPointCloud(cloud_rgb, "bottom-dist-cloud");
		for(int idx = 0; idx < item.cloudBorder_.size(); idx++)
		{
			std::pair<Eigen::Vector3d, Eigen::Vector3d> &it = item.cloudBorder_[idx];
			pcl::PointXYZ s(it.first(0), it.first(1), it.first(2));
			pcl::PointXYZ e(it.second(0), it.second(1), it.second(2));
			LOG(INFO) << "-----s:" << s.x << "," << s.y << "," << s.z;
			LOG(INFO) << "-----e:" << e.x << "," << e.y << "," << e.z;
			viewer.addSphere(s, 0.2f, "bpttom-" + std::to_string(idx));
		}
	}

	if (optRets.count(TransformOptimize::CloudType::TOP_E))
	{
		LOG(INFO) << "************TOP_E*****";
		auto &ret = optRets[TransformOptimize::CloudType::TOP_E];
		auto Top = cad.getTypedModelItems(ITEM_TOP_E).front();
		pcl::PointCloud<pcl::PointXYZ>::Ptr pData = ret.vecCloud_.front();
		CloudItem item(pData);
		item.type_ = CLOUD_TOP_E;
		item.cloudPlane_ = ret.vecCloudPlane_.front();
		item.cadPlane_ = ret.vecCadPlane_.front();
		item.cadBorder_.insert(item.cadBorder_.begin(), Top.segments_.begin(), Top.segments_.end());
		auto boundPoints = calcCloudBorder(pData, item.cloudPlane_, item.cadBorder_, item.cloudBorder_);
		mapCloudItem_[CLOUD_TOP_E].emplace_back(item);

		auto cloud_rgb = calcDistError(pData, item.cloudPlane_, 1.0);
		viewer.addPointCloud(boundPoints, "top-boundPoints");
		viewer.addPointCloud(cloud_rgb, "top-dist-cloud");
		for(int idx = 0; idx < item.cloudBorder_.size(); idx++)
		{
			std::pair<Eigen::Vector3d, Eigen::Vector3d> &it = item.cloudBorder_[idx];
			pcl::PointXYZ s(it.first(0), it.first(1), it.first(2));
			pcl::PointXYZ e(it.second(0), it.second(1), it.second(2));
			LOG(INFO) << "-----s:" << s.x << "," << s.y << "," << s.z;
			LOG(INFO) << "-----e:" << e.x << "," << e.y << "," << e.z;
			viewer.addSphere(s, 0.2f, "top-" + std::to_string(idx));
		}
	}

	
	auto vecWall = cad.getTypedModelItems(ITEM_WALL_E);
	if (optRets.count(TransformOptimize::CloudType::WALL_E) 
			&& vecWall.size() == optRets[TransformOptimize::CloudType::WALL_E].vecCloud_.size())
	{
		auto &ret = optRets[TransformOptimize::CloudType::WALL_E];
		for (int i = 0; i < vecWall.size(); i++)
		{
			LOG(INFO) << "************WALL_E: " << i << "**************";
			auto& wall = vecWall[i];
			pcl::PointCloud<pcl::PointXYZ>::Ptr pData = ret.vecCloud_[i];;
			CloudItem item(pData);
			item.type_ = CLOUD_WALL_E;
			item.cloudPlane_ = ret.vecCloudPlane_[i];
			item.cadPlane_ = ret.vecCadPlane_[i];
			item.cadBorder_.insert(item.cadBorder_.begin(), wall.segments_.begin(), wall.segments_.end());
			auto boundPoints = calcCloudBorder(pData, item.cloudPlane_, item.cadBorder_, item.cloudBorder_);
			mapCloudItem_[CLOUD_WALL_E].emplace_back(item);

			auto cloud_rgb = calcDistError(pData, item.cloudPlane_, 1.0);
			viewer.addPointCloud(boundPoints, "wall-" + std::to_string(i) + "-boundPoints");
			viewer.addPointCloud(cloud_rgb, "wall-" + std::to_string(i) + "dist-cloud");
			for(int idx = 0; idx < item.cloudBorder_.size(); idx++)
			{
				std::pair<Eigen::Vector3d, Eigen::Vector3d> &it = item.cloudBorder_[idx];
				pcl::PointXYZ s(it.first(0), it.first(1), it.first(2));
				pcl::PointXYZ e(it.second(0), it.second(1), it.second(2));
				LOG(INFO) << "-----s:" << s.x << "," << s.y << "," << s.z;
				LOG(INFO) << "-----e:" << e.x << "," << e.y << "," << e.z;
				viewer.addSphere(s, 0.2f, "wall-" + std::to_string(i) + std::to_string(idx));
			}
		}
	}

	auto vecHole = cad.getTypedModelItems(ITEM_HOLE_E);
	for (auto& hole : vecHole)
	{
		auto& item = mapCloudItem_[CLOUD_WALL_E][hole.parentIndex_];
		item.cadBorder_.insert(item.cadBorder_.begin(), hole.segments_.begin(), hole.segments_.end());
	}
	
	auto vecBeam = cad.getTypedModelItems(ITEM_BEAM_E);
	if (optRets.count(TransformOptimize::CloudType::BEAM_E) 
			&& vecWall.size() == optRets[TransformOptimize::CloudType::BEAM_E].vecCloud_.size())
	{
		auto &ret = optRets[TransformOptimize::CloudType::BEAM_E];
		for (int i = 0; i < vecBeam.size(); i++)
		{
			LOG(INFO) << "************BEAM_E: " << i << "**************";
			auto& beam = vecBeam[i];
			pcl::PointCloud<pcl::PointXYZ>::Ptr pData = ret.vecCloud_[i];
			CloudItem item(pData);
			item.type_ = CLOUD_BEAM_E;
			item.parentIndex_ = item.parentIndex_;
			item.cloudPlane_ = ret.vecCloudPlane_[i];
			item.cadPlane_ = ret.vecCadPlane_[i];
			item.cadBorder_.insert(item.cadBorder_.begin(), beam.segments_.begin(), beam.segments_.end());
			auto boundPoints = calcCloudBorder(pData, item.cloudPlane_, item.cadBorder_, item.cloudBorder_);
			mapCloudItem_[CLOUD_BEAM_E].emplace_back(item);

			auto cloud_rgb = calcDistError(pData, item.cloudPlane_, 1.0);
			viewer.addPointCloud(boundPoints, "beam-" + std::to_string(i) + "-boundPoints");
			viewer.addPointCloud(cloud_rgb, "beam-" + std::to_string(i) + "dist-cloud");
			for(int idx = 0; idx < item.cloudBorder_.size(); idx++)
			{
				std::pair<Eigen::Vector3d, Eigen::Vector3d> &it = item.cloudBorder_[idx];
				pcl::PointXYZ s(it.first(0), it.first(1), it.first(2));
				pcl::PointXYZ e(it.second(0), it.second(1), it.second(2));
				LOG(INFO) << "-----s:" << s.x << "," << s.y << "," << s.z;
				LOG(INFO) << "-----e:" << e.x << "," << e.y << "," << e.z;
				viewer.addSphere(s, 0.2f, "beam-" + std::to_string(i) + std::to_string(idx));
			}
		}
	}	

	while (!viewer.wasStopped())
	{
		viewer.spinOnce();
	}
}

}