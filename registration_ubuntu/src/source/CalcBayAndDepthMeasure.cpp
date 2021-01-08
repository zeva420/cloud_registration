#include "CalcBayAndDepthMeasure.h"

#include "funHelper.h"

namespace CloudReg
{



	std::tuple<int,Eigen::Vector3d, Eigen::Vector3d, std::vector<seg_pair_t>> calcOverlap(seg_pair_t& toSeg, 
			seg_pair_t& calcSeg)
	{
		if((toSeg.first - toSeg.second).dot(calcSeg.first - calcSeg.second) < EPS_FLOAT_DOUBLE)
		{
			auto tmp = calcSeg;
			calcSeg.first = tmp.second;
			calcSeg.second = tmp.first;
		}

		bool findS1 = isRootInSeg(toSeg, calcSeg.first);
		bool findE1 = isRootInSeg(toSeg, calcSeg.second);
		bool findS2 = isRootInSeg(calcSeg, toSeg.first);
		bool findE2 = isRootInSeg(calcSeg, toSeg.second);
				
		LOG(INFO)<< findS1 << " " << findE1 << " " << findS2 << " " << findE2;
		
		Eigen::Vector3d sPt, ePt;
		int ret = -1;
		std::vector<seg_pair_t> vecSeg;
		if (findE1 && findS1)
		{
			sPt = calcPerpendicular(calcSeg.first, toSeg.first,toSeg.second);
			ePt = calcPerpendicular(calcSeg.second, toSeg.first,toSeg.second);
			ret = 0;
			
			vecSeg.emplace_back(std::make_pair(calcSeg.first, sPt));
			vecSeg.emplace_back(std::make_pair(calcSeg.second, ePt));

		}else if (findE2 && findS2)
		{
			sPt = calcPerpendicular(toSeg.first, calcSeg.first,calcSeg.second);
			ePt = calcPerpendicular(toSeg.second, calcSeg.first,calcSeg.second);
			ret = 1;
			
			vecSeg.emplace_back(std::make_pair(toSeg.first, sPt));
			vecSeg.emplace_back(std::make_pair(toSeg.second, ePt));
		}else if(findS1 && findE2)
		{
			sPt = calcPerpendicular(calcSeg.first ,toSeg.first,toSeg.second);
			ePt = toSeg.second;
			ret = 0;

			auto pt = calcPerpendicular(toSeg.second,calcSeg.first,calcSeg.second);
			vecSeg.emplace_back(std::make_pair(calcSeg.first, sPt));
			vecSeg.emplace_back(std::make_pair(toSeg.second, pt));
		}
		else if (findS2 && findE1)
		{
			sPt = toSeg.first;
			ePt = calcPerpendicular(calcSeg.second, toSeg.first,toSeg.second);
			ret = 0;
			
			auto pt = calcPerpendicular(toSeg.first,calcSeg.first,calcSeg.second);
			vecSeg.emplace_back(std::make_pair(calcSeg.second, ePt));
			vecSeg.emplace_back(std::make_pair(toSeg.first, pt));
		}

		return std::tuple<int, Eigen::Vector3d, Eigen::Vector3d, std::vector<seg_pair_t>>(ret, sPt, ePt,vecSeg);
	}
	
	void calcDepth(const std::vector<seg_pair_t>& rootBorder,
			const std::vector<vec_seg_pair_t>& allWallBorder,
			const std::map<std::size_t, std::vector<vec_seg_pair_t>>& holeBorder,
			const std::vector<PointCloud::Ptr>& vecCloud,
			const int optType,
			const double calcLengthTh)
	{
		if (rootBorder.size() != allWallBorder.size())
		{
			LOG(ERROR) << "the root and wall size not match: " << rootBorder.size() << "-" << allWallBorder.size();
		}
	
		const Eigen::Vector3d& horizenSeg = rootBorder.front().first - rootBorder.front().second;
		std::vector<std::size_t> vecVerticalIndex;
		std::vector<std::size_t> vecHorizenIndex;
		groupDirectionIndex(horizenSeg, rootBorder, vecVerticalIndex, vecHorizenIndex);

		std::vector<seg_pair_t> vecCutSeg;
		const auto& calcIndex = optType == 0 ? vecHorizenIndex : vecVerticalIndex;

		for(std::size_t i = 0; i< calcIndex.size(); i++)
		{
			seg_pair_t toSeg = rootBorder[calcIndex[i]];
			if ((toSeg.first - toSeg.second).norm() < calcLengthTh)
				continue;

			for(std::size_t j = i+1 ; j < calcIndex.size(); j++)
			{
				seg_pair_t calcSeg = rootBorder[calcIndex[j]];
				if ((calcSeg.first - calcSeg.second).norm() < calcLengthTh)
					continue;

				//LOG(INFO)<< "calcSeg: " << vecHorizenIndex[i] <<" " << vecHorizenIndex[j];

				Eigen::Vector3d sPt, ePt;
				int hasOverlap;
				std::vector<seg_pair_t> vecSeg;
				std::tie(hasOverlap, sPt, ePt, vecSeg) = calcOverlap(toSeg,calcSeg);
				
				if (hasOverlap < 0) continue;

				if((sPt-ePt).norm() < calcLengthTh)
					continue;
			
				vecCutSeg.insert(vecCutSeg.end(), vecSeg.begin(), vecSeg.end());
				LOG(INFO)<< "findSeg:" << vecHorizenIndex[i] <<" " << vecHorizenIndex[j];
			}
		}

		if (optType == 0) writePCDFile("root_depth.pcd",rootBorder, vecCutSeg);
		else writePCDFile("root_bay.pcd",rootBorder, vecCutSeg);
	}

}//namespace
