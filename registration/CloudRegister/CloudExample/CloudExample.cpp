// CloudExample.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include "FileHelper.h"
#include "CloudRegister.h"

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/StdVector>

std::vector<Eigen::Vector3d> ininterpolateSeg(const Eigen::Vector3d& sPoint, const Eigen::Vector3d& ePoint, const double step)
{
	std::vector<Eigen::Vector3d> value;


	const std::size_t number = (sPoint - ePoint).norm() / step;

	const double r_N = 1.0 / number;
	const Eigen::Vector3d step_AB = (ePoint - sPoint) * r_N;

	value.emplace_back(sPoint);
	for (std::size_t i = 1; i < number; i++)
	{
		value.emplace_back(Eigen::Vector3d(step_AB * i + sPoint));
	}
	value.emplace_back(ePoint);
	return value;
}

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
   
	auto mapCloud = obj.getAllCloudPlane();
	std::vector<std::string> itemName{"Beam","Bottom","Wall","Top","Unknow"};
	for (auto& value : mapCloud)
	{
		const std::string name = itemName[value.first];

		for(std::size_t index = 0; index < value.second.size(); index++)
		{
			auto& item = value.second[index];
			std::vector<Eigen::Vector3d> cadPoints;
			for (auto& vecSegs : item.cadBorder_) {
				for (auto& pt_pair : vecSegs) {
					auto vec_tmp = ininterpolateSeg(pt_pair.first, pt_pair.second, .1);
					cadPoints.insert(cadPoints.end(), vec_tmp.begin(), vec_tmp.end());
				}
			}

			std::string file_name = name + "_" + std::to_string(index) + ".pcd";
			pcl::PointCloud<pcl::PointXYZ> cloud;
			cloud.width = cadPoints.size() + item.pCloud_->size();
			cloud.height = 1;
			cloud.is_dense = false;
			cloud.points.resize(cloud.width * cloud.height);

			for (size_t i = 0; i < item.pCloud_->size(); ++i)
				cloud.points[i] = item.pCloud_->points[i];
			
			for (size_t i = 0; i < cadPoints.size(); ++i)
			{
				cloud.points[item.pCloud_->size()+i].x = cadPoints[i][0];
				cloud.points[item.pCloud_->size()+i].y = cadPoints[i][1];
				cloud.points[item.pCloud_->size()+i].z = cadPoints[i][2];
			}

			pcl::io::savePCDFile(file_name, cloud);
			
			file_name = "cad_"+ file_name;
			pcl::io::savePCDFile(file_name, *item.pCADCloud_);

			std::vector<Eigen::Vector3d> cloudBorder;
			for (auto& vecSegs : item.cloudBorder_) {
				for (auto& pt_pair : vecSegs) {
					auto vec_tmp = ininterpolateSeg(pt_pair.first, pt_pair.second, 0.05);
					cloudBorder.insert(cloudBorder.end(), vec_tmp.begin(), vec_tmp.end());
				}
			}
			pcl::PointCloud<pcl::PointXYZ>::Ptr pCloudBorder(new pcl::PointCloud<pcl::PointXYZ>());
			for (size_t i = 0; i < cloudBorder.size(); ++i)
			{
				pcl::PointXYZ p(cloudBorder[i](0), cloudBorder[i](1), cloudBorder[i](2));
				pCloudBorder->push_back(p);
			}			
			file_name = "cloudBorder-" + name + "_" + std::to_string(index) + ".pcd";
			pcl::io::savePCDFile(file_name, *pCloudBorder);
		}
	}
	return 0;
}

