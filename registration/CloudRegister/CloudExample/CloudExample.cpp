// CloudExample.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>

#include "FileHelper.h"
#include "CloudRegister.h"



int main()
{
	if (__argc < 3)
	{
		std::cout << "CloudExample PCD_Dir CAD_File"<< std::endl;
		return -1;

	}

	std::string pcd_dir = __argv[1];
	std::string cad_file = __argv[2];
	std::vector<std::string> pcd_list;
	FileHelper::getFilenamesFromdir(pcd_dir,pcd_list);

	std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> vecCloudPtr;
	for (auto& file : pcd_list)
	{
		pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud(new pcl::PointCloud<pcl::PointXYZ>);
		if (pcl::io::loadPCDFile<pcl::PointXYZ>(file, *pCloud) == -1) {
			std::cout << "Couldn't read file " << file << std::endl;
			return -1;
		}
		vecCloudPtr.emplace_back(pCloud);
	}

	CloudReg::CloudRegister obj;
	obj.run(vecCloudPtr, cad_file);
   
	return 0;
}

