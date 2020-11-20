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
	obj.run(cloud, model);

	return true;
}

const std::map<CloudItem, vecItems_t>& CloudRegister::getAllCloudPlane() const
{
	return mapCloudItem_;
}

const std::map<pairCloud_t, std::pair<double, double>>& CloudRegister::getAllCorner() const
{
	return mapCorner_;
}

pcl::PointCloud<pcl::PointXYZRGB>::Ptr
CloudRegister::calcDistError(const pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud_,
	const Eigen::Vector3d& plane, const double downRatio) const
{
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr pCloud = nullptr;
	return pCloud;
}

pcl::PointCloud<pcl::PointXYZRGB>::Ptr
CloudRegister::genCloudByModel(const Eigen::Vector3d& planePara,
	const std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>& border) const
{
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr pCloud = nullptr;
	return pCloud;
}


}