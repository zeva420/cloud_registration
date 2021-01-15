// CloudExample.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include "FileHelper.h"
#include "CloudRegister.h"

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/StdVector>
#ifdef UBUNTU_SWITCH
#include <pcl/keypoints/impl/uniform_sampling.hpp>
#else
#include <pcl/filters/uniform_sampling.h>
#endif

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

void writePCDFile(const std::string& name, const pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud,
	std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>& border)
{
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr pCloudRGB(new pcl::PointCloud<pcl::PointXYZRGB>());

	for (auto& pt : pCloud->points)
	{
		pcl::PointXYZRGB p2;
		p2.x = pt.x;
		p2.y = pt.y;
		p2.z = pt.z;
		p2.r = 0;
		p2.g = 255;
		p2.b = 0;

		pCloudRGB->push_back(p2);
	}

	for (auto& seg : border)
	{
		auto vecPts = ininterpolateSeg(seg.first, seg.second, 0.01);
		for (auto& pt : vecPts)
		{
			pcl::PointXYZRGB p2;
			p2.x = pt[0];
			p2.y = pt[1];
			p2.z = pt[2];
			p2.r = 255;
			p2.g = 0;
			p2.b = 0;

			pCloudRGB->push_back(p2);
		}
	}
	pcl::io::savePCDFile(name, *pCloudRGB);

}

void uniformSampling(double radius,
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud,
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_filtered)
{
	pcl::UniformSampling<pcl::PointXYZ> filter;
	filter.setInputCloud(cloud);
	filter.setRadiusSearch(radius);

	#ifdef UBUNTU_SWITCH
	pcl::PointCloud<int> keypointIndices;
	filter.compute(keypointIndices);
	pcl::copyPointCloud(*cloud, keypointIndices.points, *cloud_filtered);
	#else
	filter.filter(*cloud_filtered);
	#endif
}

#ifdef UBUNTU_SWITCH
int main(int argc, char** argv)
#else
int main()
#endif
{
	#ifdef UBUNTU_SWITCH
	if (argc < 3)
	#else
	if (__argc < 3)
	#endif
	{
		std::cout << "CloudExample PCD_Dir CAD_File"<< std::endl;
		return -1;

	}

	#ifdef UBUNTU_SWITCH
	std::string pcd_dir = argv[1];
	std::string cad_file = argv[2];
	#else
	std::string pcd_dir = __argv[1];
	std::string cad_file = __argv[2];
	#endif

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
#if 1
	for (auto& value : mapCloud)
	{
		const std::string name = itemName[value.first];

		for(std::size_t index = 0; index < value.second.size(); index++)
		{
			auto& item = value.second[index];
			std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> vecBorder;
			
			
			for (auto& vecSegs : item.cadBorder_) {
				vecBorder.insert(vecBorder.end(), vecSegs.begin(), vecSegs.end());
			}
			pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud_filtered(new pcl::PointCloud<pcl::PointXYZ>());
			uniformSampling(0.01, item.pCloud_, pCloud_filtered);

			std::string file_name = name + "_" + std::to_string(index) + ".pcd";
			writePCDFile(file_name, pCloud_filtered, vecBorder);

			
			file_name = "cad_cloud"+ file_name;
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
#endif
	obj.calcRoofNetHeight();
	return 0;
}

