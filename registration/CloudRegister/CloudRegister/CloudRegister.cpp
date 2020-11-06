#include "CloudRegister.h"


#include "BaseType.h"
#include "CADModel.h"

namespace CloudReg
{
	CloudRegister::CloudRegister()
	{
		google::InitGoogleLogging("Cloud");
		FLAGS_log_dir = "./log";
	}

	CloudRegister::~CloudRegister()
	{
		google::ShutdownGoogleLogging();
	}

	bool CloudRegister::run(std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr>& vecCloudPtr,
		const std::string& CAD_File)
	{

		LOG(INFO) << "file: " << CAD_File;

		CADModel model;
		model.initCAD(CAD_File);
		return true;
	}
}