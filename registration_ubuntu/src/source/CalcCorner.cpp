#include "CalcCorner.h"

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
		std::tie(optIndex,indexOther,dir) = getGrowAxisAndDir(seg.first, seg.second);

		Eigen::Vector3d pt1,pt2;
		if(bLeft)
		{
			pt2 = pt1 = seg.second;
			pt1[optIndex] -= dir * calcLength;

			pt1[2] = height;
			pt2[2] = height;

		}else{
			pt2 = pt1 = seg.first;
			pt2[optIndex] += dir * calcLength;
		
			pt1[2] = height;
			pt2[2] = height;
		}

		double thickness = indexOther == 0 ? (max.x - min.x) : (max.y - min.y);	
		auto vecPt = createRulerBox(std::make_pair(pt1,pt2),indexOther,thickness,calcWidth*2);
		auto rangeSeg = calcBoxSegPair(vecPt);

		
		auto pCloud = filerCloudByConvexHull(pWall, vecPt);
		//writePCDFile("test.pcd", pCloud, rangeSeg);
	
		return std::make_tuple(pCloud,rangeSeg);
	}

	void CalcCorner(const std::vector<std::vector<seg_pair_t>>& allWallBorder,
			const std::map<std::size_t, std::vector<std::vector<seg_pair_t>>>& holeBorder,
			const std::vector<PointCloud::Ptr>& vecCloud,
			const double calcLengthTh)
	{
		if (allWallBorder.size() != vecCloud.size())
		{
			LOG(WARNING) << "size not match";
			return;
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
					std::tie(optIndex,indexOther,dir) = getGrowAxisAndDir(leftWall.back().first, leftWall.back().second);
			
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
					std::tie(optIndex,indexOther,dir) = getGrowAxisAndDir(rightWall.back().first, rightWall.back().second);
				
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

			LOG(INFO) << "calc between: " << std::to_string(idx.first) << " - " << std::to_string(idx.second);
			//0.3
			{
				auto left = calcCornerArea(leftWall.back(), vecCloud[idx.first],0.3,true);
				auto right = calcCornerArea(rightWall.back(), vecCloud[idx.second],0.3,false);
			}
			
			//1.5
			{
				auto left = calcCornerArea(leftWall.back(), vecCloud[idx.first],1.5,true);
				auto right = calcCornerArea(rightWall.back(), vecCloud[idx.second],1.5,false);
			}

		}
	}

}//namespace
