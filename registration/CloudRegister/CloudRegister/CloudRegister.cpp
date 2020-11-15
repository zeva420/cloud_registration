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
	if (!cm.run(vecCloudPtr, model)) {
		LOG(INFO) << "coarse matching failed.";
		return false;
	}

	// registration
    std::string logStr = "";
    TransformOptimize obj("refined Transform Opt", logStr);
    obj.run(vecCloudPtr, model);

	return true;
}
}