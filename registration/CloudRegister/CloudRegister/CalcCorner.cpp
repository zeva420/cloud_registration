#include "CalcCorner.h"

#include "funHelper.h"
#include <pcl/common/common.h>


namespace CloudReg
{
	std::tuple<PointCloud::Ptr, std::vector<seg_pair_t>> 
		calcCornerArea(const seg_pair_t& seg, const PointCloud::Ptr pWall, double height, bool bLeft)
	{
		pcl::PointXYZ min;
		pcl::PointXYZ max;		
		pcl::getMinMax3D(*pWall, min, max);

		const double calcWidth = 0.005;
		const double calcLength = 0.13; 
		std::size_t optIndex,indexOther; 
		int dir;
		std::tie(optIndex,indexOther,dir) = getWallGrowAxisAndDir(seg.first, seg.second);

		Eigen::Vector3d pt1,pt2;
		if(bLeft)
		{
			pt2 = pt1 = seg.first;
			pt1[optIndex] += dir * calcLength;

			pt1[2] = height;
			pt2[2] = height;

		}else{
			pt2 = pt1 = seg.second;
			pt2[optIndex] -= dir * calcLength;
		
			pt1[2] = height;
			pt2[2] = height;
		}

		double thickness = indexOther == 0 ? (max.x - min.x) : (max.y - min.y);	
		auto vecPt = createRulerBox(std::make_pair(pt1,pt2),indexOther,thickness,calcWidth*2);
		auto rangeSeg = calcBoxSegPair(vecPt);
		auto filerPt = getRulerCorners(vecPt); 
		auto pCloud = filerCloudByConvexHull(pWall,filerPt); 
		//writePCDFile("test.pcd", pCloud, rangeSeg);
	
		return std::make_tuple(pCloud,rangeSeg);
	}

	std::map<std::pair<std::size_t, std::size_t>, std::vector<calcMeassurment_t>>
	CalcCorner(const std::vector<std::vector<seg_pair_t>>& allWallBorder,
			const std::map<std::size_t, std::vector<std::vector<seg_pair_t>>>& holeBorder,
			const std::vector<PointCloud::Ptr>& vecCloud,
			const double calcLengthTh)
	{
		std::map<std::pair<std::size_t, std::size_t>, std::vector<calcMeassurment_t>> result;
		if (allWallBorder.size() != vecCloud.size())
		{
			LOG(WARNING) << "size not match";
			return result;
		}

		std::vector<std::pair<std::size_t, std::size_t>> calcIdx;
		for(std::size_t i = 1; i < vecCloud.size(); i++)
		{
			calcIdx.emplace_back(std::make_pair(i-1,i));
			if (i == vecCloud.size() -1)
				calcIdx.emplace_back(std::make_pair(i,0));
		}

		auto calcLength = [](seg_pair_t& value)->double { return (value.first - value.second).norm();};
		for(auto& idx : calcIdx)
		{
			auto leftWall = allWallBorder[idx.first];
			auto rightWall = allWallBorder[idx.second];

			if (calcLength(leftWall.back()) < calcLengthTh ||
					calcLength(rightWall.back()) < calcLengthTh)
				continue;

			//check hole
			{
				auto iter = holeBorder.find(idx.first);
				if(iter != holeBorder.end())
				{
					auto& vecHole = iter->second;
					auto ptA = leftWall.back().second;
				
					std::vector<Eigen::Vector3d> vecPt;
					for(auto& hole : vecHole)
					vecPt.emplace_back(hole.back().first);

					std::size_t optIndex,indexOther; 
					int dir;
					std::tie(optIndex,indexOther,dir) = getWallGrowAxisAndDir(leftWall.back().first, leftWall.back().second);
			
					if(dir > 0)
					{
						std::sort(vecPt.begin(), vecPt.end(),[&](const Eigen::Vector3d& left, const Eigen::Vector3d& right){
							 return left[optIndex] < right[optIndex];});
					}else{
						std::sort(vecPt.begin(), vecPt.end(),[&](const Eigen::Vector3d& left, const Eigen::Vector3d& right){
							 return left[optIndex] > right[optIndex];});
					}
					auto ptB = vecPt.front();
					if ((ptA-ptB).norm() < calcLengthTh)
						continue;
				}
			}
			
			{
				
				auto iter = holeBorder.find(idx.second);
				if(iter != holeBorder.end())
				{
					auto& vecHole = iter->second;
					auto ptA = leftWall.back().first;
				
					std::vector<Eigen::Vector3d> vecPt;
					for(auto& hole : vecHole)
						vecPt.emplace_back(hole.back().second);

					std::size_t optIndex,indexOther; 
					int dir;
					std::tie(optIndex,indexOther,dir) = getWallGrowAxisAndDir(rightWall.back().first, rightWall.back().second);
				
					if(dir > 0)
					{
						std::sort(vecPt.begin(), vecPt.end(),[&](const Eigen::Vector3d& left, const Eigen::Vector3d& right){
							 return left[optIndex] < right[optIndex];});
					}else{
						std::sort(vecPt.begin(), vecPt.end(),[&](const Eigen::Vector3d& left, const Eigen::Vector3d& right){
							 return left[optIndex] > right[optIndex];});
					}
					auto ptB = vecPt.back();
					if ((ptA-ptB).norm() < calcLengthTh)
						continue;
				}
			}

			auto getCornerByPlaneNorm = [](PointCloud::Ptr left, PointCloud::Ptr right,
				const Eigen::Vector3d &center)->double {
				LOG(INFO) << "left size:" << left->size() << ", right size:" << right->size();		
				if (left->size() < 50 || right->size() < 50) return -1;

				Eigen::Vector4d plane1 = calcPlaneParam(left);
				Eigen::Vector4d plane2 = calcPlaneParam(right);
				Eigen::Vector3d n1 = plane1.block<3,1>(0,0);
				Eigen::Vector3d n2 = plane2.block<3,1>(0,0);
				Eigen::Vector3d p1(left->front().x, left->front().y, left->front().z);
				Eigen::Vector3d p2(right->front().x, right->front().y, right->front().z);
				p1 = p1 - center;
				p2 = p2 - center;
				n1 = (n1.dot(p1) > 0) ? n1 : (-1.0 * n1);
				n2 = (n2.dot(p2) > 0) ? n2 : (-1.0 * n2);
				// LOG(INFO) << "n1:" << n1(0) << "," << n1(1) << "," << n1(2);
				// LOG(INFO) << "n2:" << n2(0) << "," << n2(1) << "," << n2(2);
				double v = calcCorner(n1, n2);	
				return v;				
			};

			// getCornerByPlaneNorm
			LOG(INFO) << "calc between: " << std::to_string(idx.first) << " - " << std::to_string(idx.second);
			Eigen::Vector3d center = ((leftWall.back().first - leftWall.back().second) 
				+ (rightWall.back().second - rightWall.back().first)) / 2.0;
			// LOG(INFO) << "center:" << center;
			//0.3
			{
				auto left = calcCornerArea(leftWall.back(), vecCloud[idx.first],0.3,true);
				auto right = calcCornerArea(rightWall.back(), vecCloud[idx.second],0.3,false);
				double v = getCornerByPlaneNorm(std::get<0>(left), std::get<0>(right), center);	
				calcMeassurment_t meassurment;
				meassurment.value = v;
				meassurment.rangeSeg = std::get<1>(left);
				meassurment.rangeSeg.insert(meassurment.rangeSeg.end(), std::get<1>(right).begin(), std::get<1>(right).end());
				result[std::make_pair(idx.first, idx.second)].push_back(meassurment);

#ifdef VISUALIZATION_ENABLED
				{
					pcl::PointCloud<pcl::PointXYZRGB>::Ptr pCloud(new pcl::PointCloud<pcl::PointXYZRGB>());
					for (auto &p1 : std::get<0>(left)->points)
					{
						pcl::PointXYZRGB p2;
						p2.x = p1.x;
						p2.y = p1.y;
						p2.z = p1.z;
						p2.r = 255;
						p2.g = 0;
						p2.b = 0;
						pCloud->push_back(p2);
					}

					for (auto &p1 : std::get<0>(right)->points)
					{
						pcl::PointXYZRGB p2;
						p2.x = p1.x;
						p2.y = p1.y;
						p2.z = p1.z;
						p2.r = 0;
						p2.g = 255;
						p2.b = 0;
						pCloud->push_back(p2);
					}
					std::string file_name = "corner-" + std::to_string(idx.first) 
						+ "-" + std::to_string(idx.second) + "-0.3m" + ".pcd";
					if (!pCloud->empty()) pcl::io::savePCDFile(file_name, *pCloud);	
				}				
#endif
			}
			//1.5
			{
				auto left = calcCornerArea(leftWall.back(), vecCloud[idx.first],1.5,true);
				auto right = calcCornerArea(rightWall.back(), vecCloud[idx.second],1.5,false);
				double v = getCornerByPlaneNorm(std::get<0>(left), std::get<0>(right), center);	
				calcMeassurment_t meassurment;
				meassurment.value = v;		
				meassurment.rangeSeg = std::get<1>(left);
				meassurment.rangeSeg.insert(meassurment.rangeSeg.end(), std::get<1>(right).begin(), std::get<1>(right).end());
				result[std::make_pair(idx.first, idx.second)].push_back(meassurment);	
#ifdef VISUALIZATION_ENABLED
				{
					pcl::PointCloud<pcl::PointXYZRGB>::Ptr pCloud(new pcl::PointCloud<pcl::PointXYZRGB>());
					for (auto &p1 : std::get<0>(left)->points)
					{
						pcl::PointXYZRGB p2;
						p2.x = p1.x;
						p2.y = p1.y;
						p2.z = p1.z;
						p2.r = 255;
						p2.g = 0;
						p2.b = 0;
						pCloud->push_back(p2);
					}

					for (auto &p1 : std::get<0>(right)->points)
					{
						pcl::PointXYZRGB p2;
						p2.x = p1.x;
						p2.y = p1.y;
						p2.z = p1.z;
						p2.r = 0;
						p2.g = 255;
						p2.b = 0;
						pCloud->push_back(p2);
					}
					std::string file_name = "corner-" + std::to_string(idx.first) 
						+ "-" + std::to_string(idx.second) + "-1.5m" + ".pcd";
					if (!pCloud->empty()) pcl::io::savePCDFile(file_name, *pCloud);	
				}				
#endif			
			}

		}

		return result;
	}

}//namespace
