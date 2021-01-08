#include "funHelper.h"

#include <pcl/features/normal_3d.h>
#include <pcl/features/boundary.h>
#include <pcl/filters/project_inliers.h>
#include <pcl/filters/uniform_sampling.h>
#include <pcl/filters/radius_outlier_removal.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/segmentation/extract_clusters.h>

#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/sample_consensus/ransac.h>
#include <pcl/sample_consensus/sac_model_plane.h>

#include <pcl/filters/passthrough.h>

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

	void uniformSampling(double radius,
						pcl::PointCloud<pcl::PointXYZ>::Ptr cloud,
						pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_filtered)
	{
		pcl::UniformSampling<pcl::PointXYZ> filter;
		filter.setInputCloud(cloud);
		filter.setRadiusSearch(radius);
		filter.filter(*cloud_filtered);
	}

	double pointToPLaneDist(const Eigen::Vector4d &plane, const pcl::PointXYZ &p)
	{
		Eigen::Vector3d n = plane.block<3,1>(0,0);
		double d  = plane(3);
		Eigen::Vector3d point(p.x, p.y, p.z);
		double dist = (n.dot(point) + d) / n.norm();    
		return dist;
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
	
	void searchBoundaries(pcl::PointCloud<pcl::PointXYZ>::Ptr input_cloud,
							pcl::PointCloud<pcl::Normal>::Ptr normals,
							std::vector<int> &boundIndices)
	{
		int arg_kNearest = 100;
		int arg_Angle = 144;

		pcl::PointCloud<pcl::Boundary> boundaries;
		pcl::BoundaryEstimation<pcl::PointXYZ, pcl::Normal, pcl::Boundary> est; 
		pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>());

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

	std::vector<Eigen::Vector3d> convertCloudToEigenVec(
			const pcl::PointCloud<pcl::PointXYZ>::Ptr cloud)
	{
		std::vector<Eigen::Vector3d> vecPts;	
		for (const auto &p : cloud->points)
		{
			Eigen::Vector3d pt(p.x, p.y, p.z);
			vecPts.push_back(pt);
		}	
		return vecPts;
	}

	std::pair<double, Eigen::Vector3d> findNearestPt(
			const std::vector<Eigen::Vector3d> &vecPts, const Eigen::Vector3d &point)
	{
		double minDist = double(RAND_MAX);
		int idx = 0;
		for (int i = 0; i < vecPts.size(); i++)
		{
			auto &pt = vecPts[i];
			double dist = (point - pt).norm();
			if (dist < minDist)
			{	
				minDist = dist;
				idx = i;
			}
		}
		return std::make_pair(minDist, vecPts[idx]);
	}

	bool groupPlanesBySamePt(const std::vector<std::vector<Eigen::Vector3d>> &segPtsOfPlanes,
						std::set<std::set<int>> &planeIdxGroup)
	{
		for (int i = 0; i < segPtsOfPlanes.size(); i++)
		{
			for (auto &pt1 : segPtsOfPlanes[i])
			{
				std::set<int> idxSet;
				idxSet.insert(i);
				for (int j = 0; j < segPtsOfPlanes.size(); j++)
				{
					if (i == j) continue;
					bool findSame = false;
					for (auto &pt2 : segPtsOfPlanes[j])
					{
						if ((pt1 - pt2).norm() < 1e-6)
						{
							findSame = true;
							break;
						}
					}
					if (findSame)
					{
						idxSet.insert(j);
					}
				}
				planeIdxGroup.insert(idxSet);
			}
		}

		LOG(INFO) << "planeIdxGroup:" << planeIdxGroup.size();
		if (planeIdxGroup.empty())
		{
			return false;
		}
		return true;
	}

	bool interSectionOf3Planes(const std::vector<Eigen::Vector4d> &cloudPlaneVec,
							const std::set<std::set<int>> &idGroups, 
							std::vector<Eigen::Vector3d> &focalPointVec)
	{
		for (const auto &currSet : idGroups)
		{
			if (currSet.size() != 3)
			{
				LOG(ERROR) << "neighbour currSet size:" << currSet.size() << " is not 3";
				return false;
			}
			std::vector<int> vecIds(currSet.begin(), currSet.end());
			auto idx1 = vecIds[0];
			auto idx2 = vecIds[1];
			auto idx3 = vecIds[2];
			LOG(INFO) << "idxSet:" << idx1 << "," << idx2 << "," << idx3;
			const Eigen::Vector4d &plane1 = cloudPlaneVec[idx1];
			const Eigen::Vector4d &plane2 = cloudPlaneVec[idx2];
			const Eigen::Vector4d &plane3 = cloudPlaneVec[idx3];
			
			Eigen::VectorXd interSectionLine(6);
			if (false == interSectionOfPlaneToPlane(plane1, plane2, interSectionLine)) continue;
			LOG(INFO) << "interSectionLine:\n" << interSectionLine;

			Eigen::Vector3d interSectionPt;
			if (false == interSectionOfLineToPlane(interSectionLine, plane3, interSectionPt)) continue;
			LOG(INFO) << "interSectionPt:\n" << interSectionPt;

			focalPointVec.push_back(interSectionPt);
		}
		LOG(INFO) << "focalPointVec:" << focalPointVec.size();
		if (focalPointVec.empty())
		{
			return false;
		}
		return true;
	}

	std::vector<Eigen::Vector3d> calcWallNodes(const std::string &name, 
			pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, Eigen::Vector4d &cloudPlane)
	{
		LOG(INFO) << "*******************calcWallNodes*********************";
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_filter(new pcl::PointCloud<pcl::PointXYZ>());
        projectionToPlane(cloudPlane, cloud, cloud_filter);

		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_sampling(new pcl::PointCloud<pcl::PointXYZ>);
        uniformSampling(0.01, cloud_filter, cloud_sampling);

		pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);
		pcl::Normal ptNormal(cloudPlane(0), cloudPlane(1), cloudPlane(2));
		for (int i = 0; i < cloud_sampling->size(); i++)
		{
			normals->push_back(ptNormal);
		}
		std::vector<int> boundIndices;
		searchBoundaries(cloud_sampling, normals, boundIndices);
		auto boundPoints = geo::getSubSet(cloud_sampling, boundIndices, false);

		pcl::PointCloud<pcl::PointXYZ>::Ptr inputPoints(new pcl::PointCloud<pcl::PointXYZ>());
		pcl::copyPointCloud(*boundPoints, *inputPoints);

		float disthresh = 0.01;
		std::vector<Eigen::VectorXf> lineCoeffs;
		std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> linePoints;
		while (inputPoints->size() > 10)
		{
			std::vector<int> indices;
			Eigen::VectorXf params;
			std::tie(indices, params) = geo::detectOneLineRansac(inputPoints, disthresh);
			auto inliers = geo::getSubSet(inputPoints, indices, false);
			lineCoeffs.push_back(params);
			linePoints.push_back(inliers);       

			auto tmpLeft = geo::getSubSet(inputPoints, indices, true);
			inputPoints->swap(*tmpLeft);
		}
		LOG(INFO) << "lineCoeffs:" << lineCoeffs.size() << ", linePoints:" << linePoints.size();
 
		//calc intersection node between two lines
		std::vector<Eigen::Vector3d> vecNodes;
		for (int i = 0; i < lineCoeffs.size(); i++)
		{
			auto tmp1 = lineCoeffs[i];
			Eigen::VectorXd line1(6);
			line1 << tmp1(0), tmp1(1), tmp1(2), tmp1(3), tmp1(4), tmp1(5);
			for (int j = i+1; j < lineCoeffs.size(); j++)
			{
				auto tmp2 = lineCoeffs[j];
				Eigen::VectorXd line2(6);
				line2 << tmp2(0), tmp2(1), tmp2(2), tmp2(3), tmp2(4), tmp2(5);
				Eigen::Vector3d interSectionPt;
				if (false == interSectionOfLineToLine(line1, line2, interSectionPt)) continue;
				auto vecPts1 = convertCloudToEigenVec(linePoints[i]);
				auto vecPts2 = convertCloudToEigenVec(linePoints[j]);

				auto ret1 = findNearestPt(vecPts1, interSectionPt);
				auto ret2 = findNearestPt(vecPts2, interSectionPt);
				if (ret1.first > 1.5 || ret2.first > 1.5) continue;

				vecNodes.push_back(interSectionPt);
			}
		}
		LOG(INFO) << "vecNodes:" << vecNodes.size();

//debug files		
#ifdef VISUALIZATION_ENABLED
		{	
			std::default_random_engine e;
    		std::uniform_real_distribution<double> random(0,1);
			pcl::PointCloud<pcl::PointXYZRGB>::Ptr pCloud(new pcl::PointCloud<pcl::PointXYZRGB>());
			for (size_t i = 0; i < linePoints.size(); ++i)
			{
				auto cloud = linePoints[i];
				int r = int(random(e)*255);
				int g = int(random(e)*255);
				int b = int(random(e)*255);
				for (auto &p1 : cloud->points)
				{
					pcl::PointXYZRGB p2;
					p2.x = p1.x;
					p2.y = p1.y;
					p2.z = p1.z;
					p2.r = r;
					p2.g = g;
					p2.b = b;
					pCloud->push_back(p2);
				}	
			}
			std::string file_name = "findLine-" + name + ".pcd";
			pcl::io::savePCDFile(file_name, *pCloud);				
		}
#endif

		return vecNodes;
	}
	
	pcl::PointCloud<pcl::PointXYZ>::Ptr calcCloudBorder(const std::string &name, 
			pcl::PointCloud<pcl::PointXYZ>::Ptr cloud,
			Eigen::Vector4d &cloudPlane,
			std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> &cadBorder,
			std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> &cloudBorder)
	{
		LOG(INFO) << "*******************calcCloudBorder*********************";
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_filter(new pcl::PointCloud<pcl::PointXYZ>());
        projectionToPlane(cloudPlane, cloud, cloud_filter);

		std::vector<int> boundIndices;
		searchBoundaries(cloud_filter, boundIndices);
		auto boundPoints = geo::getSubSet(cloud_filter, boundIndices, false);
		pcl::PointCloud<pcl::PointXYZ>::Ptr inputPoints(new pcl::PointCloud<pcl::PointXYZ>());
		pcl::copyPointCloud(*boundPoints, *inputPoints);

		float disthresh = 0.03;
		float connect_thresh = 0.1;
		float radius = 0.1;
		std::vector<Eigen::VectorXf> lineCoeffs;
		std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> linePoints;
		std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> lineSegs;
		while (inputPoints->size() > 10)
		{
			std::vector<int> indices;
			Eigen::VectorXf params;
			std::tie(indices, params) = geo::detectOneLineRansac(inputPoints, disthresh);

			auto inliers = geo::getSubSet(inputPoints, indices, false);
			pcl::PointCloud<pcl::PointXYZ>::Ptr inlierPoints(new pcl::PointCloud<pcl::PointXYZ>());
			pcl::copyPointCloud(*inliers, *inlierPoints);

			while (inlierPoints->size() > 10)
			{
				std::vector<int> segIndex = clusterMainStructure(inlierPoints, connect_thresh);
				auto segInliers = geo::getSubSet(inlierPoints, segIndex, false);

				std::pair<Eigen::Vector3d, Eigen::Vector3d> segment;
				if (true == detectLineEndPoints(segInliers, params, radius, segment))
				{
					lineSegs.push_back(segment);
					lineCoeffs.push_back(params);
					linePoints.push_back(inliers);
				}

				auto tmpLeft = geo::getSubSet(inlierPoints, segIndex, true);
				inlierPoints->swap(*tmpLeft);				
			}
			
			auto tmpLeft = geo::getSubSet(inputPoints, indices, true);
			inputPoints->swap(*tmpLeft);
		}
		LOG(INFO) << "cadBorder:" << cadBorder.size() << ", lineSegs:" << lineSegs.size();

		for (auto &seg1 : cadBorder)
		{
			if ((seg1.second - seg1.first).norm() < 0.000001) continue;
			LOG(INFO) << "*****************seg1 len:" 
					<< (seg1.second - seg1.first).norm() 
					<<"***********************";
			std::map<double, int> mapDist2Idx;
			for (int i = 0; i < lineSegs.size(); i++)
			{
				auto &seg2 = lineSegs[i];
				double dotProduct = (seg1.second - seg1.first).normalized().dot((seg2.second - seg2.first).normalized());
				double lenDiff = std::fabs((seg1.second - seg1.first).norm() - (seg2.second - seg2.first).norm());
				if (std::fabs(dotProduct) < 0.98) continue;

				double dist1 = distToLine(seg2.first, seg1.first, seg1.second);
				double dist2 = distToLine(seg2.second, seg1.first, seg1.second);
				double dist3 = distToLine(seg1.first, seg2.first, seg2.second);
				double dist4 = distToLine(seg1.second, seg2.first, seg2.second);
				double aveDist = (dist1 + dist2 + dist3 + dist4) / 4.0;
				LOG(INFO) << "idx:" << i << " seg2:" << (seg2.second - seg2.first).norm()
					/*<< " dot:" << dotProduct*/ << " lenDiff:" << lenDiff
					<< " dist1:" << dist1 << " dist2:" << dist2 
					<< " dist3:" << dist3 << " dist4:" << dist4
					<< " aveDist:" << aveDist;
				if (aveDist > 0.3) continue;

				mapDist2Idx[aveDist] = i;		
			}

			LOG(INFO) << "mapDist2Idx:" << mapDist2Idx.size();
			if (!mapDist2Idx.empty())
			{
				LOG(INFO) << "+++first one, aveDist:" << mapDist2Idx.begin()->first 
					<< ", idx:" << mapDist2Idx.begin()->second;
				auto &seg2 = lineSegs[mapDist2Idx.begin()->second];
				double dotProduct = (seg1.second - seg1.first).normalized().dot((seg2.second - seg2.first).normalized());
				LOG(INFO) << "dotProduct:" << dotProduct;
				if (dotProduct > 0)
				{
					cloudBorder.push_back(std::make_pair(seg2.first, seg2.second));
				}
				else
				{
					cloudBorder.push_back(std::make_pair(seg2.second, seg2.first));
				}
				
				continue;
			}
			else
			{
				LOG(INFO) << "+++empty, to find nearest pt"; 
			}

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

//debug files		
#ifdef VISUALIZATION_ENABLED
		{	
			std::default_random_engine e;
    		std::uniform_real_distribution<double> random(0,1);
			pcl::PointCloud<pcl::PointXYZRGB>::Ptr pCloudBorder(new pcl::PointCloud<pcl::PointXYZRGB>());
			for (size_t i = 0; i < linePoints.size(); ++i)
			{
				auto cloud = linePoints[i];
				int r = int(random(e)*255);
				int g = int(random(e)*255);
				int b = int(random(e)*255);
				for (auto &p1 : cloud->points)
				{
					pcl::PointXYZRGB p2;
					p2.x = p1.x;
					p2.y = p1.y;
					p2.z = p1.z;
					p2.r = r;
					p2.g = g;
					p2.b = b;
					pCloudBorder->push_back(p2);
				}	
			}
			std::string file_name = "findLine-" + name + ".pcd";
			pcl::io::savePCDFile(file_name, *pCloudBorder);				
		}
		{
			std::default_random_engine e;
    		std::uniform_real_distribution<double> random(0,1);

			pcl::PointCloud<pcl::PointXYZRGB>::Ptr pCloudBorder(new pcl::PointCloud<pcl::PointXYZRGB>());
			for (auto &seg : lineSegs)
			{
				int r = int(random(e)*255);
				int g = int(random(e)*255);
				int b = int(random(e)*255);
				auto vec_tmp = ininterpolateSeg(seg.first, seg.second, 0.05);
				for (auto &p1 : vec_tmp)
				{
					pcl::PointXYZRGB p2;
					p2.x = p1(0);
					p2.y = p1(1);
					p2.z = p1(2);
					p2.r = r;
					p2.g = g;
					p2.b = b;
					pCloudBorder->push_back(p2);
				}
			}
		
			std::string file_name = "lineSegs-" + name + ".pcd";
			pcl::io::savePCDFile(file_name, *pCloudBorder);	
		}
		{
			std::vector<Eigen::Vector3d> cloudBorder;
			for (auto& pt_pair : cadBorder) {
				auto vec_tmp = ininterpolateSeg(pt_pair.first, pt_pair.second, 0.05);
				cloudBorder.insert(cloudBorder.end(), vec_tmp.begin(), vec_tmp.end());
			}
			pcl::PointCloud<pcl::PointXYZ>::Ptr pCloudBorder(new pcl::PointCloud<pcl::PointXYZ>());
			for (size_t i = 0; i < cloudBorder.size(); ++i)
			{
				pcl::PointXYZ p(cloudBorder[i](0), cloudBorder[i](1), cloudBorder[i](2));
				pCloudBorder->push_back(p);
			}			
			std::string file_name = "cadSegs-" + name + ".pcd";
			pcl::io::savePCDFile(file_name, *pCloudBorder);	
		}
#endif

		return boundPoints;
	}

	double distToLine(const Eigen::Vector3d& p, const Eigen::Vector3d& s, const Eigen::Vector3d& e) 
	{
		Eigen::Vector3d sp = p - s;
		Eigen::Vector3d se = e - s;
		double se2 = se.squaredNorm();
		if (se2 < 1e-8) return sp.norm();

		double t = sp.dot(se) / se2;
		if (t < 0.f || t > 1.f) return std::min(sp.norm(), (p - e).norm());
		// if (t > 1.f) return (p - e).norm();

		return std::fabs(sp.cross(se).norm() / std::sqrt(se2));
	}		

	std::vector<int> clusterMainStructure(PointCloud::Ptr cloud, float distance) 
	{
		pcl::search::KdTree<Point>::Ptr tree(new pcl::search::KdTree<Point>);
		tree->setInputCloud(cloud);

		std::vector<pcl::PointIndices> clusters;
		pcl::EuclideanClusterExtraction<Point> ece;
		ece.setClusterTolerance(distance);
		ece.setSearchMethod(tree);
		ece.setInputCloud(cloud);
		ece.extract(clusters);

		pcl::PointIndices mainIdx = clusters.front();
		std::vector<int> index;
		for (auto &idx : mainIdx.indices)
		{
			index.push_back(idx);
		}
		return index;
	}

	bool detectLineEndPoints(PointCloud::Ptr inliers, Eigen::VectorXf &params,
						double radius, 
						std::pair<Eigen::Vector3d, Eigen::Vector3d> &segment) 
	{
		PointCloud::Ptr body(new PointCloud());
		{
			pcl::RadiusOutlierRemoval<Point> ror;
			ror.setInputCloud(inliers);
			ror.setRadiusSearch(radius);
			ror.setMinNeighborsInRadius(10);
			ror.filter(*body);
		}
		//LOG(INFO) << "get segment cloud: " << inliers->size() << " -> " << body->size();
		if (body->size() < 10)
		{
			return false;
		}

		// we now get ends
		Point p = geo::P_(params.block<3, 1>(0, 0));
		Point n = geo::P_(params.block<3, 1>(3, 0));

		auto projlens = ll::mapf([p, n](const Point& q) {return geo::dot(q - p, n); }, body->points);
		auto pr = std::minmax_element(projlens.begin(), projlens.end());
		std::size_t i = std::distance(projlens.begin(), pr.first);
		std::size_t j = std::distance(projlens.begin(), pr.second);

		Eigen::Vector3d s(body->points[i].x, body->points[i].y, body->points[i].z);
		Eigen::Vector3d e(body->points[j].x, body->points[j].y, body->points[j].z);
		segment = std::make_pair(s, e);
		return true;
	}

	double calcCorner(const Eigen::Vector3d &n1, const Eigen::Vector3d &n2)
	{
		const double weightTh = 130.0;//mm
		double angle = std::acos(n1.dot(n2));
		double diffAngle = std::fabs(3.1415/2.0 - angle);
		double v = std::fabs(n1.dot(n2)) * weightTh;

		//LOG(INFO) << "angle:" << angle << "->" << (angle * 180.0 / 3.14);
		LOG(INFO) << "diffAngle:" << diffAngle << "->" << (diffAngle * 180.0 / 3.14)
			<< " sin:" << std::sin(diffAngle) <<" v:" << v;

		return v;
	}

	void planeFitting(float distTh, 
					pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, 
					Eigen::VectorXf &coeff, std::vector<int> &inlierIdxs)
	{
		pcl::SampleConsensusModelPlane<pcl::PointXYZ>::Ptr model(
							new pcl::SampleConsensusModelPlane<pcl::PointXYZ>(cloud));
		pcl::RandomSampleConsensus<pcl::PointXYZ> ransac(model);
		ransac.setDistanceThreshold(distTh);
		ransac.computeModel();
		ransac.getInliers(inlierIdxs);
		
		//ax+by+cz+d=0，coeff: a,b,c,d
		ransac.getModelCoefficients(coeff);
	}

	Eigen::Vector4d calcPlaneParam(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud)
	{
		Eigen::MatrixXd A(cloud->size(), 4);
		for (std::size_t r = 0; r < cloud->size(); r++) {
			const auto& p = cloud->points[r];
			A(r, 0) = p.x;
			A(r, 1) = p.y;
			A(r, 2) = p.z;
			A(r, 3) = 1.;
		}

		// then solve Ax = 0
		Eigen::MatrixXd ATA = A.transpose()* A;
		//todo: assert that ATA must be full rank
		Eigen::Vector4d x = ATA.jacobiSvd(Eigen::ComputeFullV).matrixV().col(ATA.rows()-1);
		Eigen::Vector3d n = x.block<3, 1>(0, 0);
		double len = n.norm();
		x *= 1./len;

		return x;
	}

	bool isPointInPolygon2D(const Eigen::Vector2d &point, 
			const Eigen::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>> &segments)
	{
		const double PI = 3.1415;
		if (segments.size() < 2)
		{
			return false;
		}

		auto dist_to_seg = [](const Eigen::Vector2d &point, 
				const std::pair<Eigen::Vector2d, Eigen::Vector2d> &seg)->double {
			Eigen::Vector2d s = seg.first;
			Eigen::Vector2d e = seg.second;
			Eigen::Vector2d n = e - s;
			Eigen::Vector2d sp = point - s;
			if (n.norm() < 1e-7 || sp.norm() < 1e-7) return double(RAND_MAX);

			double dot = sp.dot(n) / n.squaredNorm();
			if (dot < 0.0 || dot > 1.0) return double(RAND_MAX);

			double dist = std::fabs(sp(0)*n(1) - sp(1)*n(0)) / n.norm();	
			return dist;	
		};

		for (auto &seg : segments)
		{
			if (dist_to_seg(point, seg) < 1e-7)
			{
				return true;
			}
		}

		double angleSum = 0;
		for (auto &seg : segments)
		{
			Eigen::Vector2d a = seg.first - point;
			Eigen::Vector2d b = seg.second - point;
			double angle = std::atan2(a(1), a(0)) - atan2(b(1), b(0));
        
			if (angle >= PI)
				angle = angle - PI * 2;
			else if (angle <= -PI)
				angle = angle + PI * 2;
			angleSum += angle;
		}
		angleSum = std::fabs(angleSum);
		if (std::fabs(angleSum - PI * 2) < 1e-7)
		{
			return true;
		}

		return false;
	}

	bool isPointInPolygon3D(const Eigen::Vector3d &point, 
			const Eigen::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> &segments)
	{

		auto to_pcl_point = [](const Eigen::Vector3d &p)->pcl::PointXYZ {
			return pcl::PointXYZ(p(0), p(1), p(2));
		};
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>());
		for (auto &seg : segments)
		{
			cloud->push_back(to_pcl_point(seg.first));
			cloud->push_back(to_pcl_point(seg.second));
		}
        Eigen::VectorXf coeff;
        std::vector<int> inlierIdxs;
        planeFitting(0.003, cloud, coeff, inlierIdxs);
		Eigen::Vector4d plane(coeff(0), coeff(1), coeff(2), coeff(3));	
		Eigen::Vector3d n = plane.block<3,1>(0,0);
		if (std::fabs(n.dot(point) + plane(3)) > 1e-5) return false;

		auto findMaxAxis = [](const Eigen::Vector3d &n)->int
		{
			int maxIdx = -1;
			if (n.norm() > 1e-5)
			{
				std::vector<double> valueVec;
				valueVec.push_back(std::fabs(n(0)));
				valueVec.push_back(std::fabs(n(1)));
				valueVec.push_back(std::fabs(n(2)));
				auto maxItor = std::max_element(valueVec.begin(), valueVec.end());
				maxIdx = std::distance(valueVec.begin(), maxItor);
			}
			return maxIdx;
		};
		int maxAxis = findMaxAxis(n);
		if (-1 == maxAxis) return false;
	
		auto to_2D_point = [](int deletCoordinate, const Eigen::Vector3d &p)->Eigen::Vector2d {
			Eigen::Vector2d pt2D;
			int j = 0;
			for (int i = 0; i < p.size(); i++)
			{
				if (i == deletCoordinate) continue;
				pt2D(j) = p(i);
				j++;
			}
			return pt2D;
		};
		Eigen::Vector2d point2D = to_2D_point(maxAxis, point);
		Eigen::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>> segments2D;
		for (auto &seg : segments)
		{
			auto pair2D = std::make_pair(to_2D_point(maxAxis, seg.first), to_2D_point(maxAxis, seg.second));
			segments2D.push_back(pair2D);
		}
		
		if (false == isPointInPolygon2D(point2D, segments2D))
		{
			return false;
		}
		return true;
	}

	bool interSectionOfLineToLine(const Eigen::VectorXd &line1, 
					const Eigen::VectorXd &line2, Eigen::Vector3d &interSectionPt)
	{
		Eigen::Vector3d p1 = line1.block<3, 1>(0, 0);
		Eigen::Vector3d n1 = line1.block<3, 1>(3, 0);	
		Eigen::Vector3d p2 = line2.block<3, 1>(0, 0);
		Eigen::Vector3d n2 = line2.block<3, 1>(3, 0);
		if (std::fabs(n1.dot(n2)) > 0.9997)// parallel
		{
			return false;
		}

		auto dist_to_line = [](const Eigen::Vector3d &point, const Eigen::VectorXd &line)->double {
			Eigen::Vector3d p = line.block<3, 1>(0, 0);
			Eigen::Vector3d n = line.block<3, 1>(3, 0);	
			double dist = (point - p).cross(n).norm();	
			return dist;	
		};

		if (dist_to_line(p1, line2) < 1e-7)
		{
			interSectionPt = p1;
			return true;
		}
		Eigen::Vector3d A = p1;
		Eigen::Vector3d B = A + 2*n1;
		if (dist_to_line(B, line2) < 1e-7)
		{
			interSectionPt = B;
			return true;
		}
		
		auto verticalPt_to_line = [](const Eigen::Vector3d &point, 
					const Eigen::VectorXd &line)->Eigen::Vector3d {
			Eigen::Vector3d p = line.block<3, 1>(0, 0);
			Eigen::Vector3d n = line.block<3, 1>(3, 0);	
			Eigen::Vector3d verticalPt = (point - p).dot(n) * n + p;
			return verticalPt;	
		};
		
		// AO / BO = AC / BD = s, O = (s*B - A) / (s-1)
		Eigen::Vector3d C = verticalPt_to_line(A, line2);
		Eigen::Vector3d D = verticalPt_to_line(B, line2);
		Eigen::Vector3d AC = C - A;
		Eigen::Vector3d BD = D - B;
		double flag = (AC.dot(BD) > 0) ? 1.0 : -1.0;
		double s = flag * AC.norm() / BD.norm();
		interSectionPt = (s*B - A) / (s-1);
		double dist1 = dist_to_line(interSectionPt, line1);
		double dist2 = dist_to_line(interSectionPt, line2);
		//LOG(INFO) << "interSectionPt dist to line1:" << dist1 << ", dist to line2:" << dist2;
		return true;
	}

	bool interSectionOfLineToPlane(const Eigen::VectorXd &line, 
					const Eigen::Vector4d &plane, Eigen::Vector3d &interSectionPt)
	{
		Eigen::Vector3d p1 = line.block<3, 1>(0, 0);
		Eigen::Vector3d n1 = line.block<3, 1>(3, 0);
		Eigen::Vector3d n2 = plane.block<3, 1>(0, 0);
		double d2 = plane(3);
		if (n1.norm() < 1e-6 || n2.norm() < 1e-6)
		{
			LOG(WARNING) << "n1-norm:" << n1.norm() << ", n2-norm:" << n2.norm() << " are too small";
			return false;
		}
		
		double dot = n1.dot(n2);
		if (std::fabs(dot) < 1e-6)
		{
			LOG(WARNING) << "std::fabs(dot):" << std::fabs(dot) << " is too small";
			return false;
		}
		double dist = n2.dot(p1) + d2;
		double t = (-dist / dot);
		
		interSectionPt = p1 + t * n1;
		double distToPlane = n2.dot(interSectionPt) + d2;
		//LOG(INFO) << "interSectionPt:" << interSectionPt;
		//LOG(INFO) << "p1 to plane dist:" << dist << ", interSectionPt To plane dist:" << distToPlane;
		if (std::fabs(distToPlane) > 1e-6)
		{
			LOG(WARNING) << "std::fabs(distToPlane):" << std::fabs(distToPlane) << " is not zero";
			return false;
		}
		return true;
	}

	bool interSectionOfPlaneToPlane(const Eigen::Vector4d &plane1, 
					const Eigen::Vector4d &plane2, Eigen::VectorXd &interSectionLine)
	{
		Eigen::Vector3d n1 = plane1.block<3, 1>(0, 0);
		Eigen::Vector3d n2 = plane2.block<3, 1>(0, 0);
		double d1 = plane1(3);
		Eigen::Vector3d n = n1.cross(n2);
		if (n1.norm() < 1e-6 || n2.norm() <1e-6 || n.norm() < 1e-6)
		{
			LOG(WARNING) << "n1-norm:" << n1.norm() << ", n2-norm:" << n2.norm() 
				<< "n-norm:" << n.norm() << " are too small";
			return false;
		}
		
		std::vector<double> valueVec;
		valueVec.push_back(std::fabs(n1(0)));
		valueVec.push_back(std::fabs(n1(1)));
		valueVec.push_back(std::fabs(n1(2)));
		auto maxItor = std::max_element(valueVec.begin(), valueVec.end());
		auto maxIdx = std::distance(valueVec.begin(), maxItor);

		Eigen::Vector3d pt_random; // take a random point in plane1
		double tmp = 0;
		for (int i = 0; i < pt_random.size(); i++)
		{
			if (i == maxIdx) continue;
			pt_random(i) = 1.0;
			tmp += n1(i) * pt_random(i);
		}
		pt_random(maxIdx) = (-d1 - tmp) / n1(maxIdx);

		double dist = n1.dot(pt_random) + d1;
		//LOG(INFO) << "pt_random:" << pt_random;
		//LOG(INFO) << "pt_random to plane1 dist:" << dist;
		if (std::fabs(dist) > 1e-6)
		{
			LOG(WARNING) << "std::fabs(dist):" << std::fabs(dist) << " is not zero";
			return false;
		}

		//lineA is in plane1, go through pt_random,and is perpendicular to interSectionLine
		Eigen::VectorXd lineA(6);
		Eigen::Vector3d dir_lineA = n.cross(n1);
		lineA << pt_random(0), pt_random(1), pt_random(2), dir_lineA(0), dir_lineA(1), dir_lineA(2);

		Eigen::Vector3d pt;
		if (false == interSectionOfLineToPlane(lineA, plane2, pt))
		{
			LOG(WARNING) << "interSectionOfLineToPlane Failed";
			return false;
		}
		interSectionLine << pt(0), pt(1), pt(2), n(0), n(1), n(2);
		return true;
	}

	//not finish
	double calcCloudPairCorner(const std::string &name,
						const pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud1,
						const pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud2,
						const Eigen::Vector3d &floorPt, double height,
						Eigen::Vector3d center, Eigen::Vector4d bottomPlane)
	{
		LOG(INFO) << "===height:" << height << " floorPt.z:" << floorPt(2);
		const double distTh = 0.2;//m
		const double halfLen = 0.025;
		Eigen::Vector3d centerPt = floorPt;
		centerPt(2) += height;

		pcl::PointCloud<pcl::PointXYZ>::Ptr subSet1(new pcl::PointCloud<pcl::PointXYZ>());
		for (auto &p : pCloud1->points)
		{
			Eigen::Vector3d v(p.x, p.y, p.z);
			if ((v - centerPt).norm() < distTh)
			{
				subSet1->push_back(p);
			}
		}
		pcl::PointCloud<pcl::PointXYZ>::Ptr filteredCloud(new pcl::PointCloud<pcl::PointXYZ>());
		pcl::PassThrough<pcl::PointXYZ> pass;
		pass.setInputCloud (subSet1);
		pass.setFilterFieldName ("z");
		pass.setFilterLimits (centerPt(2)- halfLen, centerPt(2)+ halfLen);
		pass.filter (*filteredCloud);
		subSet1->swap(*filteredCloud);

		pcl::PointCloud<pcl::PointXYZ>::Ptr subSet2(new pcl::PointCloud<pcl::PointXYZ>());
		for (auto &p : pCloud2->points)
		{
			Eigen::Vector3d v(p.x, p.y, p.z);
			if ((v - centerPt).norm() < distTh)
			{
				subSet2->push_back(p);
			}
		}
		filteredCloud->clear();
		pass.setInputCloud (subSet2);
		pass.setFilterFieldName ("z");
		pass.setFilterLimits (centerPt(2)- halfLen, centerPt(2)+ halfLen);
		pass.filter (*filteredCloud);
		subSet2->swap(*filteredCloud);

		LOG(INFO) << "subSet1:" << subSet1->size() << " subSet2:" << subSet2->size();
		if (subSet1->size() < 50 || subSet2->size() < 50)
		{
			return 0.0;
		}

		//重新分割两个平面
		pcl::PointCloud<pcl::PointXYZ>::Ptr subSetMerged(new pcl::PointCloud<pcl::PointXYZ>());
		subSetMerged->insert(subSetMerged->end(), subSet1->begin(), subSet1->end());
		subSetMerged->insert(subSetMerged->end(), subSet2->begin(), subSet2->end());
		
		Eigen::VectorXf coeff1;
		std::vector<int> inlierIdxs1;
		planeFitting(0.002, subSetMerged, coeff1, inlierIdxs1);
		auto left1 = geo::getSubSet(subSetMerged, inlierIdxs1, true);

        Eigen::VectorXf coeff2;
        std::vector<int> inlierIdxs2;
        planeFitting(0.002, left1, coeff2, inlierIdxs2);
		auto left2 = geo::getSubSet(left1, inlierIdxs2, true);

		auto orig_inliers1 = geo::getSubSet(subSetMerged, inlierIdxs1, false);
		auto orig_inliers2 = geo::getSubSet(left1, inlierIdxs2, false);

		auto dist_to_plane = [](Eigen::Vector3f &p, Eigen::VectorXf &plane)->float {
			Eigen::Vector3f n = plane.block<3,1>(0,0);
			float d  = plane(3);
			float dist = std::abs(n.dot(p) + d) / n.norm(); 
			return dist;
		};
		for (auto &it : left2->points)
		{
			Eigen::Vector3f p(it.x, it.y, it.z);
			float dist1 = dist_to_plane(p, coeff1);
			float dist2 = dist_to_plane(p, coeff2);
			auto tmpPtr = (dist1 < dist2) ? orig_inliers1 : orig_inliers2;
			tmpPtr->push_back(it);
		}

		pcl::PointCloud<pcl::PointXYZ>::Ptr inliers1(new pcl::PointCloud<pcl::PointXYZ>());
		pcl::PointCloud<pcl::PointXYZ>::Ptr inliers2(new pcl::PointCloud<pcl::PointXYZ>());
#if 1
		//平面法向求夹角
		inliers1 = orig_inliers1;
		inliers2 = orig_inliers2;
		Eigen::Vector4d plane1 = calcPlaneParam(inliers1);
		Eigen::Vector4d plane2 = calcPlaneParam(inliers2);
		Eigen::Vector3d n1 = plane1.block<3,1>(0,0);
		Eigen::Vector3d n2 = plane2.block<3,1>(0,0);
		Eigen::Vector3d p1(inliers1->front().x, inliers1->front().y, inliers1->front().z);
		Eigen::Vector3d p2(inliers2->front().x, inliers2->front().y, inliers2->front().z);
		p1 = p1 - center;
		p2 = p2 - center;
		n1 = (n1.dot(p1) > 0) ? n1 : (-1.0 * n1);
		n2 = (n2.dot(p2) > 0) ? n2 : (-1.0 * n2);
		//LOG(INFO) << "n1:" << n1(0) << "," << n1(1) << "," << n1(2);
		//LOG(INFO) << "n2:" << n2(0) << "," << n2(1) << "," << n2(2);
		double v = calcCorner(n1, n2);
#else
		//投影到ＸＯＹ平面，通过检测线来获取方向，线方向求夹角
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_filter(new pcl::PointCloud<pcl::PointXYZ>());
		projectionToPlane(bottomPlane, orig_inliers1, cloud_filter);
		orig_inliers1->swap(*cloud_filter);

		cloud_filter->clear();
		projectionToPlane(bottomPlane, orig_inliers2, cloud_filter);
		orig_inliers2->swap(*cloud_filter);

		Eigen::Vector3d n1;
		{
			std::vector<int> indices;
			Eigen::VectorXf params;
			std::tie(indices, params) = geo::detectOneLineRansac(orig_inliers1, 0.002);
			Eigen::Vector3f tmp = params.block<3, 1>(3, 0);
			n1 = Eigen::Vector3d(tmp(0), tmp(1), tmp(2));
			inliers1 = geo::getSubSet(orig_inliers1, indices, false);
		}
		Eigen::Vector3d n2;
		{
			std::vector<int> indices;
			Eigen::VectorXf params;
			std::tie(indices, params) = geo::detectOneLineRansac(orig_inliers2, 0.002);
			Eigen::Vector3f tmp = params.block<3, 1>(3, 0);
			n2 = Eigen::Vector3d(tmp(0), tmp(1), tmp(2));
			inliers2 = geo::getSubSet(orig_inliers2, indices, false);
		}
		LOG(INFO) << "line1 inlierRate:" << 100.0 * (float)inliers1->size() / (float)orig_inliers2->size()
			<< ", line2 inlierRate:" << 100.0 * (float)inliers2->size() / (float)orig_inliers2->size();

		Eigen::Vector3d p1(inliers1->front().x, inliers1->front().y, inliers1->front().z);
		Eigen::Vector3d p2(inliers2->front().x, inliers2->front().y, inliers2->front().z);
		p1 = p1 - Eigen::Vector3d(center(0), center(1), floorPt(2));
		p2 = p2 - Eigen::Vector3d(center(0), center(1), floorPt(2));
		n1 = (n1.dot(p1) > 0) ? n1 : (-1.0 * n1);
		n2 = (n2.dot(p2) > 0) ? n2 : (-1.0 * n2);
		//LOG(INFO) << "n1:" << n1(0) << "," << n1(1) << "," << n1(2);
		//LOG(INFO) << "n2:" << n2(0) << "," << n2(1) << "," << n2(2);
		double v = calcCorner(n1, n2);
#endif

#ifdef VISUALIZATION_ENABLED
		{
			pcl::PointCloud<pcl::PointXYZRGB>::Ptr pCloudBorder(new pcl::PointCloud<pcl::PointXYZRGB>());
			for (auto &p1 : subSet1->points)
			{
				pcl::PointXYZRGB p2;
				p2.x = p1.x;
				p2.y = p1.y;
				p2.z = p1.z;
				p2.r = 255;
				p2.g = 0;
				p2.b = 0;
				pCloudBorder->push_back(p2);
			}

			for (auto &p1 : subSet2->points)
			{
				pcl::PointXYZRGB p2;
				p2.x = p1.x;
				p2.y = p1.y;
				p2.z = p1.z;
				p2.r = 0;
				p2.g = 255;
				p2.b = 0;
				pCloudBorder->push_back(p2);
			}
			std::string file_name = "corner-" + name + ".pcd";
			pcl::io::savePCDFile(file_name, *pCloudBorder);	
		}

		{
			pcl::PointCloud<pcl::PointXYZRGB>::Ptr pCloudBorder(new pcl::PointCloud<pcl::PointXYZRGB>());
			for (auto &p1 : inliers1->points)
			{
				pcl::PointXYZRGB p2;
				p2.x = p1.x;
				p2.y = p1.y;
				p2.z = p1.z;
				p2.r = 255;
				p2.g = 0;
				p2.b = 0;
				pCloudBorder->push_back(p2);
			}

			for (auto &p1 : inliers2->points)
			{
				pcl::PointXYZRGB p2;
				p2.x = p1.x;
				p2.y = p1.y;
				p2.z = p1.z;
				p2.r = 0;
				p2.g = 255;
				p2.b = 0;
				pCloudBorder->push_back(p2);
			}
			std::string file_name = "inliers-" + name + ".pcd";
			pcl::io::savePCDFile(file_name, *pCloudBorder);	
		}
#endif
		return v;
	}

}

// #include "g2o/core/base_edge.h"
#include "g2o/core/base_vertex.h"
#include "g2o/core/base_unary_edge.h"
#include "g2o/core/block_solver.h"
#include "g2o/solvers/eigen/linear_solver_eigen.h"
#include "g2o/core/optimization_algorithm_levenberg.h"
#include "g2o/core/robust_kernel_impl.h"
#include "g2o/types/sba/types_sba.h"

using Vec4 = Eigen::Matrix<double, 4, 1>;

// [x, y, theta, phi]
namespace g2o {
class VAngleMeasurer : public BaseVertex<4, Vec4> {
public:
	static constexpr double LENGTH = 0.2;

	Eigen::Vector2d xy() const { return _estimate.block<2, 1>(0, 0); }
	double theta() const { return _estimate(2, 0); }
	double phi() const { return _estimate(3, 0); }

	void setToOriginImpl() override { _estimate = Vec4::Zero(); }
	void oplusImpl(const double* v) override {
		Eigen::Map<const Vec4> dx(v);
		_estimate += dx;
	}

	bool read(std::istream& is) override { return true; }
	bool write(std::ostream& os) const override { return true; }

	std::string to_string() const {
		return ll::unsafe_format("AngleMeasure: (%.3f, %.3f), %.3f, %.3f",
			_estimate(0), _estimate(1), 180. / 3.1415926* _estimate(2), 180. / 3.1415926 * _estimate(3));
	}
};

class EDistanceToMeasurer2D : public BaseUnaryEdge<1, double, VAngleMeasurer> {
public:
	EDistanceToMeasurer2D(const Point& p) :
		g2o::BaseUnaryEdge<1, double, VAngleMeasurer>(), point_(p.x, p.y) {}

	void computeError() override {
		// distance to each ray.
		auto vam = static_cast<VAngleMeasurer*>(_vertices[0]);
		Eigen::Vector2d op = point_ - vam->xy();
		Eigen::Vector2d dir1 = dir(vam->theta());
		Eigen::Vector2d dir2 = dir(vam->phi());

		double t1 = op.dot(dir1);
		double t2 = op.dot(dir2);

		_error(0) = cross(op, std::abs(t1) > std::abs(t2) ? dir1 : dir2);
	}

	bool read(std::istream& is) override { return true; }
	bool write(std::ostream& os) const override { return true; }

	using g2o::BaseEdge<1, double>::information;

private:
	Eigen::Vector2d point_;

	Eigen::Vector2d dir(double angle) const { return Eigen::Vector2d(std::cos(angle), std::sin(angle)); }
	double cross(const Eigen::Vector2d& v1, const Eigen::Vector2d& v2) { return v1(0) * v2(1) - v1(1) * v2(0); }
};
}


namespace CloudReg {
double calcCorner_beta(PointCloud::Ptr cloud1, PointCloud::Ptr cloud2, const Eigen::Vector3f& cornerPoint, float z) {
	constexpr float SLICE_THICKNESS = 0.05f;
	constexpr float PIECE_MAX_LENGTH = 0.3f;
	constexpr float PIECE_MAX_LENGTH_SQUARED = PIECE_MAX_LENGTH * PIECE_MAX_LENGTH;

	const float zLow = z - SLICE_THICKNESS;
	const float zHigh = z + SLICE_THICKNESS;

	auto gain_sub_cloud = [&](PointCloud::Ptr cloud) {
		auto sub = geo::passThrough(cloud, "z", zLow, zHigh);
		sub = geo::filterPoints(sub, [&cornerPoint, PIECE_MAX_LENGTH_SQUARED](const Point& p) {
			return (geo::V_(p) - cornerPoint).squaredNorm() < PIECE_MAX_LENGTH_SQUARED;
		});

		if (sub->size() < 50) sub = cloud;

		return sub;
	};

	auto sub1 = gain_sub_cloud(cloud1);
	auto sub2 = gain_sub_cloud(cloud2);

	// now optimize
	g2o::SparseOptimizer optimier; // actually no need to be sparse.
	{
		auto ls = std::make_unique<g2o::LinearSolverEigen<g2o::BlockSolverX::PoseMatrixType>>();
		auto bs = std::make_unique<g2o::BlockSolverX>(std::move(ls));
		auto algo = new g2o::OptimizationAlgorithmLevenberg(std::move(bs));
		optimier.setAlgorithm(algo);
	}

	int id=0;

	auto v = new g2o::VAngleMeasurer();
	{
		// simply use middle point to initialize
		auto dir_angle = [](PointCloud::Ptr cloud, const Eigen::Vector3f& o)-> double{
			double cenx{ 0. }, ceny{ 0. }, dn { static_cast<double>(cloud->size()) };
			for(const auto& p: cloud->points){
				cenx += p.x; ceny += p.y;
			}

			return std::atan2(ceny/ dn- o(1), cenx/ dn- o(0));
		};

		Vec4 estimate;
		estimate(0, 0) = cornerPoint(0);
		estimate(1, 0) = cornerPoint(1);
		estimate(2, 0) = dir_angle(sub1, cornerPoint);
		estimate(3, 0) = dir_angle(sub2, cornerPoint);
		v->setEstimate(estimate);
	}
	v->setId(id++);
	v->setFixed(false);
	optimier.addVertex(v);

	Eigen::Matrix<double, 1, 1> info = Eigen::Matrix<double, 1, 1>::Identity();
	for(const auto& p: sub1->points){
		auto e = new g2o::EDistanceToMeasurer2D(p);
		e->setVertex(0, v);
		e->setMeasurement(0.);
		e->information() = info;
		optimier.addEdge(e);
	}
	for (const auto& p : sub2->points) {
		auto e = new g2o::EDistanceToMeasurer2D(p);
		e->setVertex(0, v);
		e->setMeasurement(0.);
		e->information() = info;
		optimier.addEdge(e);
	}

	LOG(INFO) << "before: "<< v->to_string();

	optimier.setVerbose(true);
	optimier.initializeOptimization();
	optimier.optimize(10);

	LOG(INFO) << "after: " << v->to_string();

	return std::cos(v->theta()- v->phi())* 130.;
}

}