#include "CloudRegister.h"


#include "BaseType.h"
#include "CADModel.h"
#include "CoarseMatching.h"
#include "TransformOptimize.h"
#include "funHelper.h"


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

	CADModel model;
	model.initCAD(CAD_File);

	LOG(INFO) << "cad model loaded: " << model.toString() << ". from: " << CAD_File;
	

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
	if(!obj.run(cloud, model))
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

const std::map<pairCloud_t, std::pair<double, double>>& CloudRegister::getAllCorner() const
{
	return mapCorner_;
}

pcl::PointCloud<pcl::PointXYZRGB>::Ptr
CloudRegister::calcDistError(const pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud_,
	const Eigen::Vector4d& plane, const double radius) const
{
	pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud_filtered(new pcl::PointCloud<pcl::PointXYZ>());
	uniformSampling(radius, pCloud_, pCloud_filtered);

	pcl::PointCloud<pcl::PointXYZRGB>::Ptr pCloud_rgb(new pcl::PointCloud<pcl::PointXYZRGB>());
	for (auto &p : pCloud_filtered->points)
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