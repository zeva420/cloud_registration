#include "CloudRegister.h"


#include "BaseType.h"
#include "CADModel.h"
#include "CoarseMatching.h"
#include "TransformOptimize.h"

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
	auto optRet = obj.run(cloud, model);
	if (!re.isValid()) {
		LOG(INFO) << "transform opt failed.";
		return false;
	}

	//rescale model
	model.scaleModel(1000.0);

	//fill return value
	fillRet(model);

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
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr pCloud = nullptr;
	return pCloud;
}


void CloudRegister::fillRet(CADModel& cad)
{
	mapCloudItem_.clear();
	

	{
		auto Botton = cad.getTypedModelItems(ITEM_BOTTOM_E).front();
		pcl::PointCloud<pcl::PointXYZ>::Ptr pData = nullptr;
		CloudItem item(pData);
		item.type_ = CLOUD_BOTTOM_E;
		item.cadBorder_.insert(item.cadBorder_.begin(), Botton.segments_.begin(), Botton.segments_.end());
		mapCloudItem_[CLOUD_BOTTOM_E].emplace_back(item);
	}

	{
		auto Top = cad.getTypedModelItems(ITEM_TOP_E).front();
		pcl::PointCloud<pcl::PointXYZ>::Ptr pData = nullptr;
		CloudItem item(pData);
		item.type_ = CLOUD_TOP_E;
		item.cadBorder_.insert(item.cadBorder_.begin(), Top.segments_.begin(), Top.segments_.end());
		mapCloudItem_[CLOUD_TOP_E].emplace_back(item);
	}

	
	auto vecWall = cad.getTypedModelItems(ITEM_WALL_E);
	for (auto& wall : vecWall)
	{
		pcl::PointCloud<pcl::PointXYZ>::Ptr pData = nullptr;
		CloudItem item(pData);
		item.type_ = CLOUD_WALL_E;
		item.cadBorder_.insert(item.cadBorder_.begin(), wall.segments_.begin(), wall.segments_.end());
		mapCloudItem_[CLOUD_WALL_E].emplace_back(item);
	}

	auto vecHole = cad.getTypedModelItems(ITEM_HOLE_E);
	for (auto& hole : vecHole)
	{
		auto& item = mapCloudItem_[CLOUD_WALL_E][hole.parentIndex_];
		item.cadBorder_.insert(item.cadBorder_.begin(), hole.segments_.begin(), hole.segments_.end());
	}
	
	auto vecBeam = cad.getTypedModelItems(ITEM_BEAM_E);
	for (auto& beam : vecWall)
	{
		pcl::PointCloud<pcl::PointXYZ>::Ptr pData = nullptr;
		CloudItem item(pData);
		item.type_ = CLOUD_BEAM_E;
		item.parentIndex_ = item.parentIndex_;
		item.cadBorder_.insert(item.cadBorder_.begin(), beam.segments_.begin(), beam.segments_.end());
		mapCloudItem_[CLOUD_BEAM_E].emplace_back(item);
	}


	
}

}