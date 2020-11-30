#include "funHelper.h"

#include <pcl/features/normal_3d.h>
#include <pcl/features/boundary.h>
#include <pcl/filters/project_inliers.h>

namespace CloudReg
{
	Eigen::vector<Eigen::Vector3d> ininterpolateSeg(const Eigen::Vector3d& sPoint, const Eigen::Vector3d& ePoint, const double step)
	{
		Eigen::vector<Eigen::Vector3d> value;


		const std::size_t number = (sPoint - ePoint).norm() / step;

		const double r_N = 1.0 / number;
		const Eigen::Vector3d step_AB = (ePoint - sPoint) * r_N;

		value.emplace_back(sPoint);
		for (std::size_t i = 1; i < number; i++)
		{
			value.emplace_back(Eigen::Vector3d(step_AB * i + sPoint));
		}
		value.emplace_back(ePoint);
		// LOG(INFO) << "ininterpolateSeg number:" << value.size();
		return value;
	}

	std::vector<std::string> splitByCharacter(const std::string& strtem, const char a)
	{
		std::vector<std::string> strvec;

		std::string::size_type pos1, pos2;
		pos2 = strtem.find(a);
		pos1 = 0;
		while (std::string::npos != pos2)
		{
			strvec.push_back(strtem.substr(pos1, pos2 - pos1));

			pos1 = pos2 + 1;
			pos2 = strtem.find(a, pos1);
		}
		strvec.push_back(strtem.substr(pos1));
		return strvec;
	}

	bool writePCDFile(const std::string& name, Eigen::vector<Eigen::Vector3d>& vecCloud)
	{
		if (vecCloud.empty()) return false;

		pcl::PointCloud<pcl::PointXYZ> cloud;
		cloud.width = vecCloud.size();
		cloud.height = 1;
		cloud.is_dense = false;
		cloud.points.resize(cloud.width * cloud.height);

		for (size_t i = 0; i < cloud.points.size(); ++i)
		{
			cloud.points[i].x = vecCloud[i][0];
			cloud.points[i].y = vecCloud[i][1];
			cloud.points[i].z = vecCloud[i][2];
		}

		pcl::io::savePCDFile(name, cloud);

		return true;
	}

	double calcArea(const Eigen::vector<Eigen::Vector2d>& vecPts)
	{
		if (vecPts.size() < 3) return 0.0;

		double area = 0.0;
		for (std::size_t i = 2; i < vecPts.size();i++)
		{
			Eigen::Vector2d OB = vecPts[0] - vecPts[i - 1];
			Eigen::Vector2d OC = vecPts[0] - vecPts[i];
			area += OB[0] * OC[1] - OC[0] * OB[1];
			
		}
		return fabs(0.5* area);
	}

	double pointToPLaneDist(const Eigen::Vector4d &plane, const pcl::PointXYZ &p)
	{
		Eigen::Vector3d n = plane.block<3,1>(0,0);
		double d  = plane(3);
		Eigen::Vector3d point(p.x, p.y, p.z);
		double dist = std::abs(n.dot(point) + d) / n.norm();    
		return dist;
	}

	pcl::PointXYZRGB getColorPtByDist(pcl::PointXYZ &p, double dist)
	{
		std::vector<std::tuple<int, int, int>> color;
		color.push_back(std::make_tuple(0,0,255)); //rgb_blue, [-inf, -0.01]
		color.push_back(std::make_tuple(0,255,255)); //rgb_cyan, [-0.01, 0]
		color.push_back(std::make_tuple(0,255,0)); //rgb_green, [0, 0.01]
		color.push_back(std::make_tuple(255,255,0)); //rgb_yellow, [0.01, 0.02]
		color.push_back(std::make_tuple(238,134,149)); //rgb_pink, [0.02, 0.03]
		color.push_back(std::make_tuple(255,0,0)); //rgb_red, [0.03, +inf]

		int size = color.size() - 1;
		double min = -0.01;
		double max = 0.04;
		dist = std::max(min, dist);
		dist = std::min(max, dist);

		int index = std::ceil(size * (dist - min) / (max - min));
		std::tuple<int, int, int> rgb = color[index];

		pcl::PointXYZRGB p_rgb;
		p_rgb.x = p.x;
		p_rgb.y = p.y;
		p_rgb.z = p.z;
		p_rgb.r = std::get<0>(rgb);
		p_rgb.g = std::get<1>(rgb);
		p_rgb.b = std::get<2>(rgb);
		return p_rgb;
	}

	void projectionToPlane(Eigen::Vector4d &plane,
							pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, 
							pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_projected)
	{
		pcl::ModelCoefficients::Ptr coefficients(new pcl::ModelCoefficients());
		coefficients->values.resize (4);
		coefficients->values[0] = plane(0);
		coefficients->values[1] = plane(1);
		coefficients->values[2] = plane(2);
		coefficients->values[3] = plane(3);

		// Create the filtering object
		pcl::ProjectInliers<pcl::PointXYZ> proj;
		proj.setModelType (pcl::SACMODEL_PLANE);
		proj.setInputCloud (cloud);
		proj.setModelCoefficients (coefficients);
		proj.filter (*cloud_projected);
	}

	void searchBoundaries(pcl::PointCloud<pcl::PointXYZ>::Ptr input_cloud,
							std::vector<int> &boundIndices)
	{
		int arg_kNearest = 100;
		int arg_Angle = 144;

		pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);
		pcl::PointCloud<pcl::Boundary> boundaries;
		pcl::BoundaryEstimation<pcl::PointXYZ, pcl::Normal, pcl::Boundary> est; 
		pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>());

		pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> normEst;
		normEst.setInputCloud(input_cloud);
		normEst.setSearchMethod(tree);
		// normEst.setRadiusSearch(2); 
		normEst.setKSearch(50); 
		normEst.compute(*normals);

		est.setInputCloud(input_cloud);
		est.setInputNormals(normals);
		est.setAngleThreshold(M_PI*(arg_Angle % 360) / 180);
		est.setSearchMethod(tree);
		est.setKSearch(arg_kNearest);  
		est.compute(boundaries);

		for (int i = 0; i<input_cloud->size(); i++)
		{
			uint8_t x = (boundaries.points[i].boundary_point);
			int a = static_cast<int>(x);
			if (a == 1)
			{
				boundIndices.push_back(i);
			}
		}
	}

	pcl::PointCloud<pcl::PointXYZ>::Ptr calcCloudBorder(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud,
			Eigen::Vector4d &cloudPlane,
			std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> &cadBorder,
			std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> &cloudBorder)
	{
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_filter(new pcl::PointCloud<pcl::PointXYZ>);
        projectionToPlane(cloudPlane, cloud, cloud_filter);

		std::vector<int> boundIndices;
		searchBoundaries(cloud, boundIndices);
		auto boundPoints = geo::getSubSet(cloud, boundIndices, false);

		for (auto &seg1 : cadBorder)
		{
			Eigen::Vector3d nearestS;
			Eigen::Vector3d nearestE;
			double minDistS = 2000.0;
			double minDistE = 2000.0;
			for (auto &v : boundPoints->points)
			{	
				Eigen::Vector3d p(v.x, v.y, v.z);
				double distS = (p-seg1.first).norm();
				double distE = (p-seg1.second).norm();
				if (distS < minDistS)
				{
					minDistS = distS;
					nearestS = p;
				}
				if (distE < minDistE)
				{
					minDistE = distE;
					nearestE = p;
				}
			}
			cloudBorder.push_back(std::make_pair(nearestS, nearestE));
		}

		return boundPoints;
	}
	

}