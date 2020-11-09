#include "CloudRegister.h"


#include "BaseType.h"
#include "CADModel.h"
#include "CoarseMatching.h"

namespace CloudReg
{
	CloudRegister::CloudRegister()
	{
		google::InitGoogleLogging("Cloud");
		FLAGS_log_dir = "./log";

		google::LogToStderr();
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

		// wall segmentation: (PointCloud, CADModel)-> [PointCloud]

		// coarse match: ([PointCloud], CADModel)-> ([filtered PointCloud], Mat4d)
		CoarseMatching cm;
		if(!cm.run(vecCloudPtr, model)){
			LOG(INFO)<< "coarse matching failed.";
			return false;
		}

		// registration

		return true;
	}
}