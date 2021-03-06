#include "CalcCorner.h"

#include "funHelper.h"
#include "Threshold.h"

//#define VISUALIZATION_ENABLED

namespace CloudReg
{
	std::tuple<PointCloud::Ptr, std::vector<seg_pair_t>> 
		calcCornerArea(const seg_pair_t& seg, const seg_pair_t& highSeg, const PointCloud::Ptr pWall, double height, bool bLeft,
			double calcWidth, double calcLength, double dropLength)
	{
		std::size_t optIndex,indexOther; 
		int dir;
		std::tie(optIndex,indexOther,dir) = getWallGrowAxisAndDir(seg.first, seg.second);
		auto allPt = ininterpolateSeg(seg.first, seg.second, 0.001);
		auto allHighPt = ininterpolateSeg(highSeg.first, highSeg.second, 0.001);
		std::size_t moveStep = (calcLength * 1000);
		std::size_t moveHighStep = (height * 1000);
		std::size_t dropStep = (dropLength * 1000);

		Eigen::Vector3d pt1,pt2;
		if(bLeft)
		{
			//pt2 = pt1 = seg.first;
			pt2 = allPt[dropStep];
			pt1 = allPt[moveStep];
			pt1[2] = allHighPt[moveHighStep][2];
			pt2[2] = allHighPt[moveHighStep][2];


		}else{
			//pt2 = pt1 = seg.second;
			pt1 = allPt[allPt.size() - dropStep -1];
			pt2 = allPt[allPt.size() - moveStep - 1];
			pt1[2] = allHighPt[moveHighStep][2];
			pt2[2] = allHighPt[moveHighStep][2];
		}


		auto vecPt = createRulerBox(std::make_pair(pt1,pt2),indexOther,0.02,calcWidth*2);
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
				double height, const seg_pair_t& highSeg)
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
		const double calcWidth_second = TheThreshold::instance()->get_corner_rulerWidth();
		const double calcLength_second = TheThreshold::instance()->get_corner_rulerLength();

		const double cloudSize_first = 20;
		const double cloudSize_second = 10;
		const double planeFitDistTh = 0.003;
		const double dropLength = TheThreshold::instance()->get_corner_dropLength();

		//get large-scale area
		auto roughLeft = calcCornerArea(leftSeg, highSeg,pLeftCloud, height, true,
											calcWidth_first, calcLength_first,0);
		auto roughRight = calcCornerArea(rightSeg, highSeg,pRightCloud, height, false,
											calcWidth_first, calcLength_first,0);
		if (std::get<0>(roughLeft)->size() < cloudSize_first 
					|| std::get<0>(roughRight)->size() < cloudSize_first)
		{
			LOG(WARNING) << "roughLeft size:" << std::get<0>(roughLeft)->size()
				<< " or roughRight:" << std::get<0>(roughRight)->size() << " < " << cloudSize_first
				<< " height:" << height;
			
			auto left = calcCornerArea(leftSeg, highSeg, pLeftCloud, height, true, calcWidth_second, calcLength_second, dropLength);
			auto right = calcCornerArea(rightSeg, highSeg, pRightCloud, height, false, calcWidth_second, calcLength_second, dropLength);
			calcMeassurment_t meassurment;
			meassurment.value = -1;
			meassurment.rangeSeg = std::get<1>(left);
			meassurment.rangeSeg.insert(meassurment.rangeSeg.end(),
				std::get<1>(right).begin(), std::get<1>(right).end());
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
		

		PointCloud::Ptr inliers1_new(new PointCloud());
		PointCloud::Ptr inliers2_new(new PointCloud());
		auto checkPt = [&](const PointCloud::Ptr pCheck) {
			for (auto& pt : pCheck->points)
			{
				double dist1 = std::fabs(pointToPLaneDist(coeff1,pt));
				double dist2 = std::fabs(pointToPLaneDist(coeff2,pt));
				if (dist1 < dist2) inliers1_new->push_back(pt);
				else	inliers2_new->push_back(pt);			
			}			
		};

			
		checkPt(inliers1);
		checkPt(inliers2);
		//std::cout << "before:"<< inliers1->points.size() << " -- " << inliers2->points.size() << std::endl;
		//std::cout << "after:" <<inliers1_new->points.size() << " -- " << inliers2_new->points.size() << std::endl;
		
		inliers1->swap(*inliers1_new);
		inliers2->swap(*inliers2_new);

		//get small-scale area
		auto left = calcCornerArea(leftSeg, highSeg,inliers1, height, true, calcWidth_second, calcLength_second, dropLength);
		auto right = calcCornerArea(rightSeg, highSeg, inliers2, height, false, calcWidth_second, calcLength_second, dropLength);
		if (std::get<0>(left)->size() < cloudSize_second 
					|| std::get<0>(right)->size() < cloudSize_second)
		{
			LOG(WARNING) << "left size:" << std::get<0>(left)->size()
				<< "or right:" << std::get<0>(right)->size() << " < " << cloudSize_second
				<< " height:" << height;
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
					+ "-" + std::to_string(idxPair.second) + "-" + std::to_string(height) + ".pcd";
			saveTwoPieceCloud(file_name, std::get<0>(roughLeft), std::get<0>(roughRight));
		}
		{
			std::string file_name = "inliers-" + std::to_string(idxPair.first) 
					+ "-" + std::to_string(idxPair.second) + "-" + std::to_string(height) + ".pcd";
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
			{
				LOG(INFO) << "CalcCorner: between wall: " << std::to_string(idx.first) << " - " 
					<< std::to_string(idx.second) << " too short:";
				continue;
			}
			
			//check hole
			{
				//the dir is first = right
				auto iter = holeBorder.find(idx.first);
				if(iter != holeBorder.end())
				{
					auto& vecHole = iter->second;
					auto ptA = leftWall.back().first;

					std::vector<Eigen::Vector3d> vecPt;
					for (auto& hole : vecHole)
					{
						vecPt.emplace_back(hole.back().first);
						vecPt.emplace_back(hole.back().second);
					}
					

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
					
					if (ptB[2] < 0.3 && (ptA - ptB).norm() < calcLengthTh)
					{						
						LOG(INFO) << "CalcCorner: check hole in wall " << std::to_string(idx.first) 
							<< " too short: " << (ptA - ptB).norm();
						continue;
					}
				}
			}
			
			{
				
				auto iter = holeBorder.find(idx.second);
				if(iter != holeBorder.end())
				{
					auto& vecHole = iter->second;
					auto ptA = rightWall.back().second;
					
				
					std::vector<Eigen::Vector3d> vecPt;
					for (auto& hole : vecHole) {
						vecPt.emplace_back(hole.back().first);
						vecPt.emplace_back(hole.back().second);
					}

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
					if (ptB[2] < 0.3 && (ptA - ptB).norm() < calcLengthTh)
					{
						LOG(INFO) << "CalcCorner: check hole in wall " << std::to_string(idx.second) 
							<< " too short: " << (ptA - ptB).norm();
						continue;
					}
				}
			}

			auto highSeg = leftWall[leftWall.size()-2];
			std::swap(highSeg.first,highSeg.second);
			// getCornerByPlaneNorm
			LOG(INFO) << "CalcCorner:calc between: " << std::to_string(idx.first) << " - " << std::to_string(idx.second);
			//0.3
			{
				calcMeassurment_t meassurment 
					= meassurmentProcess(idx, std::make_pair(leftWall.back(), rightWall.back()), 
								std::make_pair(vecCloud[idx.first], vecCloud[idx.second]), 0.3, highSeg);
				result[std::make_pair(idx.first, idx.second)].push_back(meassurment);
			}
			//1.5
			{
				calcMeassurment_t meassurment 
					= meassurmentProcess(idx, std::make_pair(leftWall.back(), rightWall.back()), 
								std::make_pair(vecCloud[idx.first], vecCloud[idx.second]), 1.5, highSeg);
				result[std::make_pair(idx.first, idx.second)].push_back(meassurment);		
			}
		}

		return result;
	}

}//namespace
