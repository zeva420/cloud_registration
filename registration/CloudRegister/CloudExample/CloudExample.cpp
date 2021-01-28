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

double pointToPLaneDist(const Eigen::Vector4d &plane, const pcl::PointXYZ &p)
{
	Eigen::Vector3d n = plane.block<3,1>(0,0);
	double d  = plane(3);
	Eigen::Vector3d point(p.x, p.y, p.z);
	double dist = (n.dot(point) + d) / n.norm();    
	return dist;
}

void getWallColor(float dis, unsigned int & r, unsigned int & g, unsigned int & b)
{
	float rf = 0.0f;
	float gf = 0.0f;
	float bf = 0.0f;
	float factor = 0.004f;//颜色刻度因子

	#pragma region MyRegion 新色谱方案(以蓝色为主色调)
	//色系组成为红色（255,0,0）—凹、紫色（255,0,255）—凹、蓝色（0,0,255）—平、绿色（0,255,0）—凸、黄色（255,255,0）—凸
	if (dis < 0)
	{
		if (dis > -factor) //介于紫色和蓝色之间
		{
			rf = 1.0f - (dis + factor) / factor;
			gf = 0.0f;
			bf = 1.0f;
		}
		else if (dis > (-factor * 1.5))
		{
			rf = 1.0f;
			gf = 0.0f;
			bf = 1.0f - (dis + factor) / dis;
		}
		else //介于红色和紫色之间
		{
			rf = 1.0f;
			gf = 0.0f;
			bf = 1.0f - (dis + factor) / dis - 0.15;  //减去0.15作为颜色补偿
		}
	}
	else if (dis == 0) //蓝色
	{
		rf = 0.0f;
		gf = 0.0f;
		bf = 1.0f;
	}
	else //介于蓝色和紫色之间
	{
		if (dis <= factor) //介于蓝色和绿色之间
		{
			rf = 0.0f;
			gf = dis / factor;
			bf = (factor - dis) / factor;
		}
		else //介于绿色和黄色之间
		{
			rf = (dis - factor) / dis + 0.4;  //增加0.4作为颜色补偿
			gf = 1.0f;
			bf = 0.0f;
		}
	}

	#pragma endregion

	rf = rf > 1.0 ? 1.0 : rf;
	gf = gf > 1.0 ? 1.0 : gf;
	bf = bf > 1.0 ? 1.0 : bf;

	rf = rf < 0.0 ? 0.0 : rf;
	gf = gf < 0.0 ? 0.0 : gf;
	bf = bf < 0.0 ? 0.0 : bf;

	r = (unsigned int)(rf*255.0);
	g = (unsigned int)(gf*255.0);
	b = (unsigned int)(bf*255.0);
}

pcl::PointXYZRGB getColorPtByDist(pcl::PointXYZ &p, double dist)
{
	unsigned int r = 255;
	unsigned int g = 255; 
	unsigned int b = 255;
	getWallColor(dist, r, g, b);

	pcl::PointXYZRGB p_rgb;
	p_rgb.x = p.x;
	p_rgb.y = p.y;
	p_rgb.z = p.z;
	p_rgb.r = r;
	p_rgb.g = g;
	p_rgb.b = b;
	return p_rgb;
}

void writePCDFile(const std::string& name, const pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud,
	std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>& border)
{
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr pCloudRGB(new pcl::PointCloud<pcl::PointXYZRGB>());

	if (pCloud != nullptr)
	{
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

			//dist rgb
            {
                pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_rgb(new pcl::PointCloud<pcl::PointXYZRGB>);
                for (auto &p : pCloud_filtered->points)
                {
                    double dist = pointToPLaneDist(item.cloudPlane_, p);
                    pcl::PointXYZRGB p_rgb = getColorPtByDist(p, dist);
                    cloud_rgb->push_back(p_rgb);
                }

				std::string fileName = "dist-to-cloudPlane-" + name + "_" + std::to_string(index) + ".pcd";
				pcl::io::savePCDFile(fileName, *cloud_rgb);
            }			
		}
	}
#endif

#if 0
	using namespace CloudReg;
	{
		std::vector<seg_pair_t> vecSeg;
		std::vector<calcIdx2Meassurment_t> vecRet;
		std::tie(vecRet, vecSeg) = obj.calcRoofNetHeight();

		for (auto& oneVec : vecRet)
		{
			for (auto& one : oneVec.vecCalcRet)
			{
				std::cout << "calcRoofNetHeight: " << std::to_string(oneVec.idx.first)
					<< " " << std::to_string(oneVec.idx.second)
					<< " " << one.value << std::endl;
				vecSeg.insert(vecSeg.end(), one.rangeSeg.begin(), one.rangeSeg.end());
			}
		}
		writePCDFile("roof_net_height.pcd", nullptr,vecSeg);

	}
	
	{
		std::vector<seg_pair_t> vecSeg1;
		std::vector<seg_pair_t> vecSeg2;
		std::vector<calcIdx2Meassurment_t> vecRoof;
		std::vector<calcIdx2Meassurment_t> vecRoot;
		std::tie(vecRoof, vecRoot, vecSeg1, vecSeg2) = obj.calcPlaneRange();
		std::vector<seg_pair_t> vecSegTmp = vecSeg1;
		for (auto& oneVec : vecRoof)
		{
			for (auto& one : oneVec.vecCalcRet)
			{
				std::cout << "calcPlaneRange roof: " << std::to_string(oneVec.idx.first)
					<< " " << std::to_string(oneVec.idx.second)
					<< " " << one.value << std::endl;
				vecSegTmp.insert(vecSegTmp.end(), one.rangeSeg.begin(), one.rangeSeg.end());
			}
		}

		writePCDFile("roof_range.pcd", nullptr, vecSegTmp);

		vecSegTmp.clear();
		vecSegTmp = vecSeg2;
		for (auto& oneVec : vecRoot)
		{
			for (auto& one : oneVec.vecCalcRet)
			{
				std::cout << "calcPlaneRange root: " << std::to_string(oneVec.idx.first)
					<< " " << std::to_string(oneVec.idx.second)
					<< " " << one.value << std::endl;
				vecSegTmp.insert(vecSegTmp.end(), one.rangeSeg.begin(), one.rangeSeg.end());
			}
		}
		writePCDFile("root_range.pcd", nullptr, vecSegTmp);
	}
	
	
	{
		std::vector<seg_pair_t> vecSeg;
		std::vector<calcIdx2Meassurment_t> vecRet;
		std::tie(vecRet, vecSeg) = obj.calcDepth();
		for (auto& item : vecRet)
		{
			for (auto& value : item.vecCalcRet)
			{
				std::cout << "calcDepth: " << item.idx.first << "- " << item.idx.second 
					<< " value:" << value.value<< std::endl;
				vecSeg.insert(vecSeg.end(), value.rangeSeg.begin(), value.rangeSeg.end());
			}
		}
		writePCDFile("room_depth.pcd", nullptr, vecSeg);
	}

	{
		std::vector<seg_pair_t> vecSeg;
		std::vector<calcIdx2Meassurment_t> vecRet;
		std::tie(vecRet, vecSeg) = obj.calcBay();
		for (auto& item : vecRet)
		{
			for (auto& value : item.vecCalcRet)
			{
				std::cout << "calcBay: " << item.idx.first << "- " << item.idx.second
					<< " value:" << value.value << std::endl;
				vecSeg.insert(vecSeg.end(), value.rangeSeg.begin(), value.rangeSeg.end());
			}
		}
		writePCDFile("room_bay.pcd", nullptr, vecSeg);
	}
	
	{
		std::vector<seg_pair_t> vecSeg;
		std::map<std::pair<std::size_t, std::size_t>, std::vector<std::vector<calcMeassurment_t>>> vecRet;
		vecRet= obj.calcAllHole();
		for (auto& item : vecRet)
		{
			for (auto& vectype : item.second)
			{
				for (auto& value : vectype) {
					std::cout << "calcAllHole: " << item.first.first << " - " << item.first.second
						<< " value:" << value.value << std::endl;
					vecSeg.insert(vecSeg.end(), value.rangeSeg.begin(), value.rangeSeg.end());
				}
			}
		}
		//writePCDFile("room_hole.pcd", nullptr, vecSeg);
	}
#endif
	 //obj.calcWallVerticality();
	//obj.calcWallFlatness();
	 //obj.calcAllSquareness();
	//obj.calcRootFlatness();
	//obj.calcAllCorner();
	return 0;
}

