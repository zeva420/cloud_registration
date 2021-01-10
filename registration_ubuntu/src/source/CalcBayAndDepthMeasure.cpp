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
		std::tie(optIndex,indexOther,dir) = getGrowAxisAndDir(sPt,ePt);
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

		auto pCloud = filerCloudByRange(pWall,min,max);
		calcMeassurment_t item;
		item.value = 0.0;
		for(auto& pt : pCloud->points)
		{
			item.value += fabs(pointToPLaneDist(plane,pt));
		}
		item.value = item.value/pCloud->points.size();
		
		auto vecPt = createRulerBox(seg,indexOther,thickness,calcHalfPara*2);
		item.rangeSeg = calcBoxSegPair(vecPt);

		//LOG(INFO) << "avgDist:" << item.value;
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
		if((toSeg.first - s1Pt).norm() < calcLengthTh)
			hasLeft = false;
			
		if((toSeg.second - e1Pt).norm() < calcLengthTh)
			hasRight = false;

		LOG(INFO) << "has enough length: "<< calcLengthTh << " leftFlag :" << hasLeft << " rightFlag:" << hasRight;

		if(hasLeft)
		{
			auto value = calcArea(pWall, plane, toSeg.first, s1Pt, true);
			vecRet.emplace_back(value);
		}
		if(hasRight)
		{
			auto value = calcArea(pWall, plane, e1Pt, toSeg.second, false);
			vecRet.emplace_back(value);
		}

		return vecRet;
	}

	void calcDepthorBay(const std::vector<seg_pair_t>& rootBorder,
			const std::vector<vec_seg_pair_t>& allWallBorder,
			const std::map<std::size_t, std::vector<vec_seg_pair_t>>& holeBorder,
			const std::vector<PointCloud::Ptr>& vecCloud,
			const Eigen::vector<Eigen::Vector4d>& vecPlane,
			const int optType,
			const double calcLengthTh)
	{
		const std::string optName = optType == 0 ? "depth" : "bay";

		if (rootBorder.size() != allWallBorder.size())
		{
			LOG(ERROR) << "the root and wall size not match: " << rootBorder.size() << "-" << allWallBorder.size();
			return;
		}
	
		const Eigen::Vector3d& horizenSeg = rootBorder.front().first - rootBorder.front().second;
		std::vector<std::size_t> vecVerticalIndex;
		std::vector<std::size_t> vecHorizenIndex;
		groupDirectionIndex(horizenSeg, rootBorder, vecVerticalIndex, vecHorizenIndex);

		std::vector<seg_pair_t> vecCutSeg;
		const auto& calcIndex = optType == 0 ? vecHorizenIndex : vecVerticalIndex;

		std::map<std::size_t, std::vector<calcMeassurment_t>> mapCalcRet;
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

				auto& calcvecRet = mapCalcRet[calcIndex[j]];
				
				PointCloud::Ptr pWall = vecCloud[calcIndex[j]];
				auto& plane = vecPlane[calcIndex[i]];
				
				auto iterHole = holeBorder.find(calcIndex[j]);
				if (iterHole != holeBorder.end())
				{
					auto curHole = iterHole->second;
					for(auto& hole : curHole)
					{
						LOG(INFO) << "calcOverlapWithHole";
						auto vecRet = calcOverlapWithHole(s2Pt,e2Pt,hole,0.2, pWall,plane);
						if (!vecRet.empty())
						{
							calcvecRet.insert(calcvecRet.end(), vecRet.begin(), vecRet.end());
						}
					}
				}
				else
				{
					auto left = calcArea(pWall, plane, s2Pt, e2Pt, true);
					auto right = calcArea(pWall, plane, s2Pt, e2Pt, false);
					calcvecRet.emplace_back(left);
					calcvecRet.emplace_back(right);
				}
			}
		}

		for(auto& value : mapCalcRet)
		{
			std::string name = "wall_"+ std::to_string(value.first) +"_" + optName + ".pcd";
			std::vector<seg_pair_t> vecRange;
			for(auto& item : value.second)
			{
				vecRange.insert(vecRange.end(), item.rangeSeg.begin(), item.rangeSeg.end());
				LOG(INFO) << value.first << ":avgDist:" << item.value;
			}
			writePCDFile(name,vecCloud[value.first], vecRange);
		}
		std::string name = "root_" + optName + ".pcd";
		writePCDFile(name,rootBorder, vecCutSeg);
	}

}//namespace
