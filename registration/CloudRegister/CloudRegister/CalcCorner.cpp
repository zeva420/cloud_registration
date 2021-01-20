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

	void saveTwoPieceCloud(const std::string fileName, 
			const PointCloud::Ptr left, const PointCloud::Ptr right)
	{
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr pCloud(new pcl::PointCloud<pcl::PointXYZRGB>());
		for (auto &p1 : left->points)
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

		for (auto &p1 : right->points)
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

		if (!pCloud->empty()) pcl::io::savePCDFile(fileName, *pCloud);		
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
					const Eigen::Vector3d &center)
					->std::tuple<double, PointCloud::Ptr, PointCloud::Ptr> {
				LOG(INFO) << "left size:" << left->size() << ", right size:" << right->size();		
				if (left->size() < 18 || right->size() < 18) 
					return std::make_tuple(-1, left, right);

				Eigen::VectorXf coeff1;
				std::vector<int> inlierIdxs1;
				planeFitting(0.005, left, coeff1, inlierIdxs1);
				auto inliers1 = geo::getSubSet(left, inlierIdxs1, false);

				Eigen::VectorXf coeff2;
				std::vector<int> inlierIdxs2;
				planeFitting(0.005, right, coeff2, inlierIdxs2);
				auto inliers2 = geo::getSubSet(right, inlierIdxs2, false);
				if (inliers1->size() < 15 || inliers2->size() < 15)
					return std::make_tuple(-1, left, right);

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
				// LOG(INFO) << "n1:" << n1(0) << "," << n1(1) << "," << n1(2);
				// LOG(INFO) << "n2:" << n2(0) << "," << n2(1) << "," << n2(2);
				double v = calcCorner(n1, n2);

				return std::make_tuple(v, inliers1, inliers2);				
			};

			// getCornerByPlaneNorm
			LOG(INFO) << "calc between: " << std::to_string(idx.first) << " - " << std::to_string(idx.second);
			Eigen::Vector3d center = ((leftWall.back().second - leftWall.back().first) 
				+ (rightWall.back().first - rightWall.back().second)) / 2.0;
			// LOG(INFO) << "center:" << center;
			//0.3
			{
				auto left = calcCornerArea(leftWall.back(), vecCloud[idx.first],0.3,true);
				auto right = calcCornerArea(rightWall.back(), vecCloud[idx.second],0.3,false);
				auto ret = getCornerByPlaneNorm(std::get<0>(left), std::get<0>(right), center);	
				calcMeassurment_t meassurment;
				meassurment.value = std::get<0>(ret);
				meassurment.rangeSeg = std::get<1>(left);
				meassurment.rangeSeg.insert(meassurment.rangeSeg.end(), std::get<1>(right).begin(), std::get<1>(right).end());
				result[std::make_pair(idx.first, idx.second)].push_back(meassurment);

#ifdef VISUALIZATION_ENABLED
				{
					std::string file_name = "corner-" + std::to_string(idx.first) 
							+ "-" + std::to_string(idx.second) + "-0.3m" + ".pcd";
					saveTwoPieceCloud(file_name, std::get<0>(left), std::get<0>(right));
				}
				{
					std::string file_name = "inliers-" + std::to_string(idx.first) 
							+ "-" + std::to_string(idx.second) + "-0.3m" + ".pcd";
					saveTwoPieceCloud(file_name, std::get<1>(ret), std::get<2>(ret));
				}			
#endif
			}
			//1.5
			{
				auto left = calcCornerArea(leftWall.back(), vecCloud[idx.first],1.5,true);
				auto right = calcCornerArea(rightWall.back(), vecCloud[idx.second],1.5,false);
				auto ret = getCornerByPlaneNorm(std::get<0>(left), std::get<0>(right), center);	
				calcMeassurment_t meassurment;
				meassurment.value = std::get<0>(ret);		
				meassurment.rangeSeg = std::get<1>(left);
				meassurment.rangeSeg.insert(meassurment.rangeSeg.end(), std::get<1>(right).begin(), std::get<1>(right).end());
				result[std::make_pair(idx.first, idx.second)].push_back(meassurment);	
#ifdef VISUALIZATION_ENABLED
				{
					std::string file_name = "corner-" + std::to_string(idx.first) 
							+ "-" + std::to_string(idx.second) + "-1.5m" + ".pcd";
					saveTwoPieceCloud(file_name, std::get<0>(left), std::get<0>(right));
				}
				{
					std::string file_name = "inliers-" + std::to_string(idx.first) 
							+ "-" + std::to_string(idx.second) + "-1.5m" + ".pcd";
					saveTwoPieceCloud(file_name, std::get<1>(ret), std::get<2>(ret));
				}		
#endif			
			}

		}

		return result;
	}

}//namespace
