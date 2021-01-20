#include "CalcCorner.h"

#include "funHelper.h"
#include <pcl/common/common.h>
#define VISUALIZATION_ENABLED

namespace CloudReg
{
	std::tuple<PointCloud::Ptr, std::vector<seg_pair_t>> 
		calcCornerArea(const seg_pair_t& seg, const PointCloud::Ptr pWall, double height, bool bLeft,
			double calcWidth, double calcLength)
	{
		pcl::PointXYZ min;
		pcl::PointXYZ max;		
		pcl::getMinMax3D(*pWall, min, max);

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

	calcMeassurment_t meassurmentProcess(
				const std::pair<std::size_t, std::size_t> &idxPair,
				const std::pair<seg_pair_t, seg_pair_t>& segPair, 
				const std::pair<const PointCloud::Ptr, const PointCloud::Ptr> &cloudPair, 
				double height)
	{
		const auto &leftSeg = segPair.first;
		const auto &rightSeg = segPair.second;
		const auto pLeftCloud = cloudPair.first;
		const auto pRightCloud = cloudPair.second;

		Eigen::Vector3d center = ((leftSeg.second - leftSeg.first) 
									+ (rightSeg.first - rightSeg.second)) / 2.0;
		// LOG(INFO) << "center:" << center;
		const double calcWidth_first = 0.02;
		const double calcLength_first = 0.25; 
		const double calcWidth_second = 0.005;
		const double calcLength_second = 0.13; 

		const double cloudSize_first = 20;
		const double cloudSize_second = 10;
		const double planeFitDistTh = 0.003;
		
		//get large-scale area
		auto roughLeft = calcCornerArea(leftSeg, pLeftCloud, height, true, 
											calcWidth_first, calcLength_first);
		auto roughRight = calcCornerArea(rightSeg, pRightCloud, height, false, 
											calcWidth_first, calcLength_first);
		if (std::get<0>(roughLeft)->size() < cloudSize_first 
					|| std::get<0>(roughRight)->size() < cloudSize_first)
		{
			LOG(WARNING) << "roughLeft size:" << std::get<0>(roughLeft)->size()
				<< "or roughRight:" << std::get<0>(roughRight)->size() << " < " << cloudSize_first;
			calcMeassurment_t meassurment;
			meassurment.value = -1;
			meassurment.rangeSeg = std::get<1>(roughLeft);
			meassurment.rangeSeg.insert(meassurment.rangeSeg.end(), 
						std::get<1>(roughRight).begin(), std::get<1>(roughRight).end());
			return meassurment;
		}

		//plane fitting
		Eigen::VectorXf coeff1;
		std::vector<int> inlierIdxs1;
		planeFitting(planeFitDistTh, std::get<0>(roughLeft), coeff1, inlierIdxs1);
		auto inliers1 = geo::getSubSet(std::get<0>(roughLeft), inlierIdxs1, false);

		Eigen::VectorXf coeff2;
		std::vector<int> inlierIdxs2;
		planeFitting(planeFitDistTh, std::get<0>(roughRight), coeff2, inlierIdxs2);
		auto inliers2 = geo::getSubSet(std::get<0>(roughRight), inlierIdxs2, false);

		//get small-scale area
		auto left = calcCornerArea(leftSeg, inliers1, height, true, calcWidth_second, calcLength_second);
		auto right = calcCornerArea(rightSeg, inliers2, height, false, calcWidth_second, calcLength_second);
		if (std::get<0>(left)->size() < cloudSize_second 
					|| std::get<0>(right)->size() < cloudSize_second)
		{
			LOG(WARNING) << "left size:" << std::get<0>(left)->size()
				<< "or right:" << std::get<0>(right)->size() << " < " << cloudSize_second;
			calcMeassurment_t meassurment;
			meassurment.value = -1;
			meassurment.rangeSeg = std::get<1>(left);
			meassurment.rangeSeg.insert(meassurment.rangeSeg.end(), 
						std::get<1>(right).begin(), std::get<1>(right).end());
			return meassurment;
		}

		//get corner value
		auto getCornerByPlaneNorm = [](PointCloud::Ptr left, PointCloud::Ptr right,
				const Eigen::Vector3d &center)
				->std::tuple<double, PointCloud::Ptr, PointCloud::Ptr> {

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

			return std::make_tuple(v, left, right);				
		};

		auto ret = getCornerByPlaneNorm(std::get<0>(left), std::get<0>(right), center);	
		calcMeassurment_t meassurment;
		meassurment.value = std::get<0>(ret);
		meassurment.rangeSeg = std::get<1>(left);
		meassurment.rangeSeg.insert(meassurment.rangeSeg.end(), 
					std::get<1>(right).begin(), std::get<1>(right).end());
		
#ifdef VISUALIZATION_ENABLED
		{
			std::string file_name = "corner-" + std::to_string(idxPair.first) 
					+ "-" + std::to_string(idxPair.second) + "-0.3m" + ".pcd";
			saveTwoPieceCloud(file_name, std::get<0>(roughLeft), std::get<0>(roughRight));
		}
		{
			std::string file_name = "inliers-" + std::to_string(idxPair.first) 
					+ "-" + std::to_string(idxPair.second) + "-0.3m" + ".pcd";
			saveTwoPieceCloud(file_name, std::get<1>(ret), std::get<2>(ret));
		}			
#endif

		return meassurment;
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

			// getCornerByPlaneNorm
			LOG(INFO) << "calc between: " << std::to_string(idx.first) << " - " << std::to_string(idx.second);
			//0.3
			{
				calcMeassurment_t meassurment 
					= meassurmentProcess(idx, std::make_pair(leftWall.back(), rightWall.back()), 
								std::make_pair(vecCloud[idx.first], vecCloud[idx.second]), 0.3);
				result[std::make_pair(idx.first, idx.second)].push_back(meassurment);
			}
			//1.5
			{
				calcMeassurment_t meassurment 
					= meassurmentProcess(idx, std::make_pair(leftWall.back(), rightWall.back()), 
								std::make_pair(vecCloud[idx.first], vecCloud[idx.second]), 1.5);
				result[std::make_pair(idx.first, idx.second)].push_back(meassurment);		
			}
		}

		return result;
	}

}//namespace
