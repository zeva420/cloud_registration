#pragma once

#include <string>
#include <vector>

#include<pcl/point_types.h>
#include<pcl/io/pcd_io.h>

class  __declspec(dllexport) CloudRegister
{
public:
	CloudRegister();
	~CloudRegister();

	bool run(std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr>& vecCloudPtr,
		const std::string& CAD_File);
private:

};

