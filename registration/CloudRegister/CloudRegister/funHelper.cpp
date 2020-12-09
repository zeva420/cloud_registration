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

		std::size_t size = color.size() - 1;
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
		searchBoundaries(cloud, boundIndices);
		auto boundPoints = geo::getSubSet(cloud, boundIndices, false);
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
#if 1
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
		LOG(INFO) << "get segment cloud: " << inliers->size() << " -> " << body->size();
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

	// double calcCorner(const pairCloud_t &cloudPair)
	double calcCorner(const Eigen::Vector3d &n1, const Eigen::Vector3d &n2)
	{
		const double weightTh = 130;//mm
		double angle = std::acos(n1.dot(n2));
		double v1 = std::sin(angle) * weightTh;
		
		double diffAngle = std::fabs(3.1415/2.0 - angle);
		double v2 = std::sin(diffAngle) * weightTh;

		LOG(INFO) << "angle:" << angle << "->" << (angle * 180.0 / 3.14) 
			<< " sin:" << std::sin(angle) <<" v1:" << v1;
		LOG(INFO) << "diffAngle:" << diffAngle << "->" << (diffAngle * 180.0 / 3.14)
			<< " sin:" << std::sin(diffAngle) <<" v2:" << v2;
		LOG(INFO) << "sin90:"<< std::sin(3.1415/2.0);

		return v2;
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
		
		//ax+by_cz_d=0，coeff: a,b,c,d
		ransac.getModelCoefficients(coeff);
		LOG(INFO) << "plane coeff " << coeff[0] << " " <<coeff[1] 
			<< " " << coeff[2] << " " << coeff[3] 
			<< ", inlierRate:" << 100.0 * double(inlierIdxs.size()) / double(cloud->size()) << "%";
	}

	double calcCloudPairCorner(const std::string &name,
						const pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud1,
						const pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud2,
						const Eigen::Vector3d &floorPt, double floorZ, double height)
	{
		LOG(INFO) << "===height:" << height << " floorZ:" << floorZ << " floorPt.z:" << floorPt(2);
		const double distTh = 0.2;//m
		Eigen::Vector3d centerPt = floorPt;
		centerPt(2) += height;

		pcl::PointCloud<pcl::PointXYZ>::Ptr subSet1(new pcl::PointCloud<pcl::PointXYZ>());
		// for (auto &p : pCloud1->points)
		// {
		// 	Eigen::Vector3d v(p.x, p.y, p.z);
		// 	if ((v - centerPt).norm() < distTh)
		// 	{
		// 		subSet1->push_back(p);
		// 	}
		// }
		pcl::PassThrough<pcl::PointXYZ> pass;
		pass.setInputCloud (pCloud1);
		pass.setFilterFieldName ("z");
		pass.setFilterLimits (floorZ+height-distTh, floorZ+height+distTh);
		pass.filter (*subSet1);

		pcl::PointCloud<pcl::PointXYZ>::Ptr subSet2(new pcl::PointCloud<pcl::PointXYZ>());
		// for (auto &p : pCloud2->points)
		// {
		// 	Eigen::Vector3d v(p.x, p.y, p.z);
		// 	if ((v - centerPt).norm() < distTh)
		// 	{
		// 		subSet2->push_back(p);
		// 	}
		// }
		pass.setInputCloud (pCloud2);
		pass.setFilterFieldName ("z");
		pass.setFilterLimits (floorZ+height-distTh, floorZ+height+distTh);
		pass.filter (*subSet2);

		LOG(INFO) << "subSet1:" << subSet1->size() << " subSet2:" << subSet2->size();
		if (subSet1->size() < 50)
		{
			subSet1->clear();
			subSet1 = pCloud1;
		}
		if (subSet2->size() < 50)
		{
			subSet2->clear();
			subSet2 = pCloud2;
		}

		pcl::PointCloud<pcl::PointXYZ>::Ptr filterdCLoud(new pcl::PointCloud<pcl::PointXYZ>());
		uniformSampling(0.002, subSet1, filterdCLoud);
		LOG(INFO) << "subSet1:" << subSet1->size() << " after uniform:" << filterdCLoud->size();
		subSet1->swap(*filterdCLoud);
		filterdCLoud->clear();

		uniformSampling(0.002, subSet2, filterdCLoud);
		LOG(INFO) << "subSet2:" << subSet2->size() << " after uniform:" << filterdCLoud->size();
		subSet2->swap(*filterdCLoud);

        Eigen::VectorXf coeff1;
        std::vector<int> inlierIdxs1;
        planeFitting(0.02, subSet1, coeff1, inlierIdxs1);

        Eigen::VectorXf coeff2;
        std::vector<int> inlierIdxs2;
        planeFitting(0.02, subSet2, coeff2, inlierIdxs2);

		Eigen::Vector4d plane1(coeff1(0), coeff1(1), coeff1(2), coeff1(3));
		Eigen::Vector4d plane2(coeff2(0), coeff2(1), coeff2(2), coeff2(3));
		Eigen::Vector3d n1 = plane1.block<3,1>(0,0);
		Eigen::Vector3d n2 = plane2.block<3,1>(0,0);
		Eigen::Vector3d p1(subSet1->front().x, subSet1->front().y, subSet1->front().z);
		Eigen::Vector3d p2(subSet2->front().x, subSet2->front().y, subSet2->front().z);
		n1 = (n1.dot(p1) > 0) ? n1 : (-1.0 * n1);
		n2 = (n2.dot(p2) > 0) ? n2 : (-1.0 * n2);
		LOG(INFO) << "n1:" << n1(0) << "," << n1(1) << "," << n1(2);
		LOG(INFO) << "n2:" << n2(0) << "," << n2(1) << "," << n2(2);
		double v = calcCorner(n1, n2);

//debug file
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

		return v;
	}

}