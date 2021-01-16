#include "CalcBayAndDepthMeasure.h"

#include "funHelper.h"
#include <pcl/common/common.h>

namespace CloudReg
{


	calcMeassurment_t calcArea(PointCloud::Ptr pWall, const Eigen::Vector4d& plane, 
			const Eigen::Vector3d& sPt, const Eigen::Vector3d& ePt, bool bStart)
	{
		std::size_t optIndex,indexOther; 
		int dir;

		const double calcPara = 0.2;
		const double calcParaZ = 0.5;
		const double calcHalfPara = 0.05;
		std::tie(optIndex,indexOther,dir) = getWallGrowAxisAndDir(sPt,ePt);
		Eigen::Vector3d calcPt;
		if (bStart)
		{
			calcPt = sPt;
			calcPt[optIndex] += dir * calcPara ;
		}else{
			calcPt = ePt;
			calcPt[optIndex] -= dir * calcPara;
		}

		pcl::PointXYZ min;
		pcl::PointXYZ max;		
		pcl::getMinMax3D(*pWall, min, max);
		min.z = calcParaZ - calcHalfPara;
		max.z = calcParaZ + calcHalfPara;
		
		seg_pair_t seg(calcPt,calcPt);
		double thickness = 0.0;	
		seg.first[2] = min.z;
		seg.second[2] = max.z;

		if (optIndex == 0)
		{
			min.x = calcPt[optIndex] - dir * calcHalfPara;
			max.x = calcPt[optIndex] + dir * calcHalfPara;
			thickness = std::fabs(max.y - min.y);
		}else{
			min.y = calcPt[optIndex] - dir * calcHalfPara;
			max.y = calcPt[optIndex] + dir * calcHalfPara;
			thickness = std::fabs(max.x - min.x);
		}

		calcMeassurment_t item;
		auto vecPt = createRulerBox(seg,indexOther,thickness,calcHalfPara*2);
		item.rangeSeg = calcBoxSegPair(vecPt);
		
		auto filerPt = getRulerCorners(vecPt);	
		auto pCloud = filerCloudByConvexHull(pWall, filerPt);
		if (pCloud->points.empty())
		{
			LOG(WARNING) << "filerCloudByConvexHull failed";
			return item;
		}
		item.value = 0.0;
		for(auto& pt : pCloud->points)
		{
			item.value += fabs(pointToPLaneDist(plane,pt));
		}
		item.value = item.value/pCloud->points.size();
		

		LOG(INFO) << "calcDepthorBay avgDist:" << item.value;
		//writePCDFile("test.pcd", pCloud, item.rangeSeg);

		return item;
	}

	std::vector<calcMeassurment_t> calcOverlapWithHole(const Eigen::Vector3d& sPt, const Eigen::Vector3d ePt,
			const vec_seg_pair_t& vecHole, 
			const double calcLengthTh, PointCloud::Ptr pWall, const Eigen::Vector4d& plane)
	{
		std::vector<calcMeassurment_t> vecRet;
		std::vector<seg_pair_t> vecHorizen;	
		std::vector<seg_pair_t> vecVertical;
		const auto horizenSeg = sPt - ePt;
		groupDirection(horizenSeg, vecHole, vecVertical, vecHorizen);
		
		if (vecHorizen.size() < 2 || vecVertical.size() < 2)
		{
			LOG(ERROR) << "groupDirection failed: " << vecHorizen.size() << " -- " << vecVertical.size();
			return vecRet;
		}

		seg_pair_t toSeg = std::make_pair(sPt,ePt);
		seg_pair_t calcSeg = std::make_pair(Eigen::Vector3d(0,0,0), Eigen::Vector3d(0,0,0));
		for(auto& value : vecHorizen)
		{
			if((value.first - value.second).squaredNorm() 
					> (calcSeg.first- calcSeg.second).squaredNorm())
				calcSeg = value;
		}

		bool hasOverlap;
		Eigen::Vector3d s1Pt, e1Pt, s2Pt,e2Pt;
		std::tie(hasOverlap, s1Pt, e1Pt, s2Pt,e2Pt) = calcOverlap(toSeg,calcSeg);

		if(!hasOverlap ||  (s1Pt-e1Pt).norm() < EPS_FLOAT_DOUBLE || 
				(s2Pt-e2Pt).norm() < EPS_FLOAT_DOUBLE)
			return vecRet;

		bool hasLeft = true;
		bool hasRight = true;
		Eigen::Vector3d left_s = toSeg.first;
		Eigen::Vector3d left_e = s1Pt;
		const double left_length = (toSeg.first - s1Pt).norm();
		const double right_length = (toSeg.second - e1Pt).norm();

		if (left_length < calcLengthTh) {

			if (right_length < calcLengthTh)
			{
				hasLeft = false;
			}
			else 
			{
				LOG(INFO) << "recalc left side pt because of hole";
				left_s = e1Pt;
				left_e = toSeg.second;
			}
			
		}
		
		Eigen::Vector3d right_s = e1Pt;
		Eigen::Vector3d right_e = toSeg.second;
		if (right_length < calcLengthTh)
		{
			if (left_length < calcLengthTh)
			{
				hasRight = false;
			}
			else
			{
				LOG(INFO) << "re calc right side pt because of hole";
				right_s = toSeg.first;
				right_e = s1Pt;
			}
			
		}

		LOG(INFO) << "has enough length: "<< calcLengthTh << " leftFlag :" << hasLeft << " rightFlag:" << hasRight;

		if(hasLeft)
		{
			auto value = calcArea(pWall, plane, left_s, left_e, true);
			if(!value.rangeSeg.empty())
				vecRet.emplace_back(value);
		}
		if(hasRight)
		{
			auto value = calcArea(pWall, plane, right_s, right_e, false);
			if (!value.rangeSeg.empty())
				vecRet.emplace_back(value);
		}

		return vecRet;
	}

	std::tuple<std::map<std::pair<std::size_t, std::size_t>, 
		std::vector<calcMeassurment_t>>,std::vector<seg_pair_t>>
	calcDepthorBay(const std::vector<seg_pair_t>& rootBorder,
			const std::vector<vec_seg_pair_t>& allWallBorder,
			const std::map<std::size_t, std::vector<vec_seg_pair_t>>& holeBorder,
			const std::vector<PointCloud::Ptr>& vecCloud,
			const Eigen::vector<Eigen::Vector4d>& vecPlane,
			const int optType,
			const double calcLengthTh)
	{
		const std::string optName = optType == 0 ? "depth" : "bay";
		std::vector<seg_pair_t> vecCutSeg;
		std::map<std::pair<std::size_t, std::size_t>, std::vector<calcMeassurment_t>> mapCalcRet;

		if (rootBorder.size() != allWallBorder.size())
		{
			LOG(ERROR) << "the root and wall size not match: " << rootBorder.size() << "-" << allWallBorder.size();
			return std::make_tuple(mapCalcRet, vecCutSeg);
		}
	
		const Eigen::Vector3d& horizenSeg = rootBorder.front().first - rootBorder.front().second;
		std::vector<std::size_t> vecVerticalIndex;
		std::vector<std::size_t> vecHorizenIndex;
		groupDirectionIndex(horizenSeg, rootBorder, vecVerticalIndex, vecHorizenIndex);

		
		const auto& calcIndex = optType == 0 ? vecHorizenIndex : vecVerticalIndex;

		
		for(std::size_t i = 0; i< calcIndex.size(); i++)
		{
			seg_pair_t toSeg = rootBorder[calcIndex[i]];
			if ((toSeg.first - toSeg.second).norm() < calcLengthTh)
				continue;
			
			//if(calcIndex[i] != 3) continue;

			for(std::size_t j = i+1 ; j < calcIndex.size(); j++)
			{
				//if(calcIndex[j] != 5) continue;
				seg_pair_t calcSeg = rootBorder[calcIndex[j]];
				if ((calcSeg.first - calcSeg.second).norm() < calcLengthTh)
					continue;

				//LOG(INFO)<< "calcSeg: " << vecHorizenIndex[i] <<" " << vecHorizenIndex[j];

				bool hasOverlap;
				Eigen::Vector3d s1Pt, e1Pt, s2Pt,e2Pt;
				std::tie(hasOverlap, s1Pt, e1Pt, s2Pt,e2Pt) = calcOverlap(toSeg,calcSeg);
				
				if (!hasOverlap) continue;

				if((s1Pt-e1Pt).norm() < calcLengthTh || (s2Pt-e2Pt).norm() < calcLengthTh)
					continue;
		
				vecCutSeg.emplace_back(std::make_pair(s1Pt,s2Pt));
				vecCutSeg.emplace_back(std::make_pair(e1Pt,e2Pt));
				LOG(INFO)<< "type: "<< optName  << " findSeg:" << calcIndex[i] <<" " << calcIndex[j];

				auto save_key = std::make_pair(calcIndex[j], calcIndex[i]);				
				PointCloud::Ptr pWall = vecCloud[calcIndex[j]];
				auto& plane = vecPlane[calcIndex[i]];
				
				auto iterHole = holeBorder.find(calcIndex[j]);
				if (iterHole != holeBorder.end())
				{
					auto curHole = iterHole->second;
					for(auto& hole : curHole)
					{
						LOG(INFO) << "calcOverlapWithHole:" << calcIndex[j];
						auto vecRet = calcOverlapWithHole(s2Pt,e2Pt,hole,0.2, pWall,plane);
						if (!vecRet.empty())
						{
							mapCalcRet[save_key].insert(mapCalcRet[save_key].end(), vecRet.begin(), vecRet.end());
						}
					}
				}
				else
				{
					auto left = calcArea(pWall, plane, s2Pt, e2Pt, true);
					auto right = calcArea(pWall, plane, s2Pt, e2Pt, false);
					mapCalcRet[save_key].emplace_back(left);
					mapCalcRet[save_key].emplace_back(right);
				}
			}
		}

		for(auto& value : mapCalcRet)
		{
			
			std::vector<seg_pair_t> vecRange;
			for(auto& item : value.second)
			{
				vecRange.insert(vecRange.end(), item.rangeSeg.begin(), item.rangeSeg.end());
				LOG(INFO) << value.first.first << " - " << value.first.second << " :avgDist:" << item.value;
			}
#ifdef VISUALIZATION_ENABLED
			std::string name = optName + "_wall_" +std::to_string(value.first.first) + "_" + std::to_string(value.first.second) + ".pcd";
			pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud_filtered(new pcl::PointCloud<pcl::PointXYZ>());
			uniformSampling(0.01, vecCloud[value.first.first], pCloud_filtered);
			writePCDFile(name, pCloud_filtered, vecRange);
#endif // DEBUG
		}
		//std::string name = "root_" + optName + ".pcd";
		//writePCDFile(name,rootBorder, vecCutSeg);
		return std::make_tuple(mapCalcRet, vecCutSeg);
	}

	

}//namespace
