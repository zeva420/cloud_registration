#include "CalcHoleMeasure.h"

#include "funHelper.h"
#include <pcl/common/common.h>
#include <pcl/filters/crop_box.h>
#include <pcl/filters/passthrough.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>

namespace CloudReg
{


	void calcHoleMaxAndMin(std::vector<seg_pair_t> vecHorizen, std::vector<seg_pair_t> vecVertical,
			const double extendRange, std::vector<Eigen::Vector3d>& vecFilerPt)
	{

		std::sort(vecHorizen.begin(), vecHorizen.end(), [&](const seg_pair_t& left, const seg_pair_t& right){
				double distLeft = (left.first - left.second).squaredNorm();	
				double distRight = (right.first - right.second).squaredNorm();	
				return distLeft < distRight;});
		
		
		std::sort(vecVertical.begin(), vecVertical.end(), [&](const seg_pair_t& left, const seg_pair_t& right){
				double distLeft = (left.first - left.second).squaredNorm();	
				double distRight = (right.first - right.second).squaredNorm();	
				return distLeft < distRight;});

		seg_pair_t segA, segB;
		segA = vecHorizen.back();
		segB = vecVertical.back();
	
		segA.first[2] = segB.first[2] > segB.second[2] ? segB.first[2]/2 : segB.second[2]/2;
		segA.second[2] = segA.first[2];
		std::size_t optIndex,indexOther; 
		int dir;
		std::tie(optIndex,indexOther,dir) = getWallGrowAxisAndDir(segA.first,segA.second);
		segA.first[optIndex] -= extendRange;
		segA.second[optIndex] += extendRange;
		
		double width = (segB.first - segB.second).norm() + extendRange * 2; 
		auto vecPt = createRulerBox(segA,indexOther,0.1,width);
		vecFilerPt = getRulerCorners(vecPt);	

	}



	std::vector<calcMeassurment_t> calcHoleCross(const std::vector<seg_pair_t>& holeBorder, const PointCloud::Ptr pCloud)
	{
		std::vector<calcMeassurment_t> vecRet;
		Eigen::vector<Eigen::Vector3d> vecCalcPt;
		vecCalcPt.emplace_back(holeBorder[0].second);
		vecCalcPt.emplace_back(holeBorder[2].second);
		vecCalcPt.emplace_back(holeBorder[1].second);
		vecCalcPt.emplace_back(holeBorder[3].second);
		
		auto vecRawPt = getNearestPt(vecCalcPt, pCloud,0.1*0.1);
		for(std::size_t i = 1; i < vecRawPt.size(); i+=2)
		{
			calcMeassurment_t item;
			item.value = (vecRawPt[i-1] - vecRawPt[i]).norm();
			item.rangeSeg.emplace_back(std::make_pair(vecRawPt[i-1],vecRawPt[i]));
			vecRet.emplace_back(item);
			LOG(INFO) << "cross dist :" << item.value;
		}

		return vecRet;
	}

	std::vector<calcMeassurment_t> calcHoleDist(const seg_pair_t& baseSeg, const Eigen::vector<Eigen::Vector3d>& vecPts, 
			const PointCloud::Ptr pCloud, const std::size_t distTh, bool bHeight)
	{
	 	
		std::vector<calcMeassurment_t> vecRet;
		double length = (vecPts.back() - vecPts.front()).norm();
		if (length < 0.5f)
		{
			LOG(ERROR) << "too short length: " << length;
			return vecRet;
		}

		Eigen::vector<Eigen::Vector3d> vecCalcPt;
		const std::size_t begin = distTh;//150mm
		const std::size_t end = vecPts.size()- begin -1;//150mm
		vecCalcPt.emplace_back(vecPts[begin]);
		vecCalcPt.emplace_back(vecPts[end]);
		if (length > 1.5 && bHeight)
		{
			std::size_t mid = (vecPts.size()-1)/2;
			vecCalcPt.emplace_back(vecPts[mid]);
			
		}

		auto vecRawPt = getNearestPt(vecCalcPt, pCloud, 0.1*0.1);
		for(auto& pt : vecRawPt)
		{
			auto root = calcPerpendicular(pt, baseSeg.first, baseSeg.second);
			
			calcMeassurment_t item;
			item.value = (pt -root).norm();
			item.rangeSeg.emplace_back(std::make_pair(pt,root));
			vecRet.emplace_back(item);
			LOG(INFO) << bHeight << " dist :" << item.value;
		}

		return vecRet;
	}

	void calcHole(const std::pair<Eigen::Vector3d, Eigen::Vector3d>& horizen,
			std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>& holeBorder, 
			pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud)
	{

		std::vector<seg_pair_t> vecHorizen;	
		std::vector<seg_pair_t> vecVertical;
		const auto horizenSeg = horizen.first - horizen.second;
		groupDirection(horizenSeg,holeBorder, vecVertical, vecHorizen);

		if (vecHorizen.size() < 2 || vecVertical.size() < 2)
		{
			LOG(ERROR) << "groupDirection failed: " << vecHorizen.size() << " -- " << vecVertical.size();
			return;
		}

		std::vector<Eigen::Vector3d> vecFilerPt;
		calcHoleMaxAndMin(vecHorizen,vecVertical, 0.05, vecFilerPt);
		auto rangeCloud = filerCloudByConvexHull(pCloud,vecFilerPt);
		if (rangeCloud->points.empty()) 
		{
			LOG(ERROR) << "filerCloudByRange failed";
			return;
		}
		//pcl::io::savePCDFile("cloud.pcd", *pCloud);
		//pcl::io::savePCDFile("rangeCloud.pcd", *rangeCloud);

		std::vector<seg_pair_t> vecRange;

		//height
		{
			Eigen::Vector3d A1 = vecHorizen.front().first - vecHorizen.front().second;
			Eigen::Vector3d A2 = vecHorizen.back().first - vecHorizen.back().second;

			std::size_t baseIndex =  A1.norm() >= A2.norm() ? 0 : vecHorizen.size()-1;
			const seg_pair_t& baseSeg = vecHorizen[baseIndex];
		
			for(std::size_t i = 0; i < vecHorizen.size(); i++)
			{
				if (i == baseIndex) continue;

				const auto& seg = vecHorizen[i];
				auto vecPts = ininterpolateSeg(seg.first, seg.second, 0.01f);
				auto vecHeight = calcHoleDist(baseSeg, vecPts, rangeCloud, 15, true);

				for(auto& item : vecHeight)
					vecRange.insert(vecRange.end(), item.rangeSeg.begin(), item.rangeSeg.end());
			}
		}

		//width
		{
			Eigen::Vector3d A1 = vecVertical.front().first - vecVertical.front().second;
			Eigen::Vector3d A2 = vecVertical.back().first - vecVertical.back().second;

			std::size_t baseIndex =  A1.norm() >= A2.norm() ? 0 : vecVertical.size()-1;
			
			const seg_pair_t& baseSeg = vecVertical[baseIndex];

			for(std::size_t i = 0; i < vecVertical.size(); i++)
			{
				if (i == baseIndex) continue;

				const auto& seg = vecVertical[i];
				auto vecPts = ininterpolateSeg(seg.first, seg.second, 0.01f);
				auto vecWidth = calcHoleDist(baseSeg, vecPts, rangeCloud, 15, false);
				
				for(auto& item : vecWidth)
					vecRange.insert(vecRange.end(), item.rangeSeg.begin(), item.rangeSeg.end());
			}
		}

		//cross
		if (holeBorder.size() == 4)
		{
			auto vecCross = calcHoleCross(holeBorder,rangeCloud);
			
			for(auto& item : vecCross)
				vecRange.insert(vecRange.end(), item.rangeSeg.begin(), item.rangeSeg.end());
		}

		//to root
		if (vecHorizen.front().first[2] > 0.1)
		{
			const auto& seg = vecHorizen.front();
			auto vecPts = ininterpolateSeg(seg.first, seg.second, 0.01f);
			auto vecRoot = calcHoleDist(horizen, vecPts, rangeCloud, 15, true);
			
			for(auto& item : vecRoot)
				vecRange.insert(vecRange.end(), item.rangeSeg.begin(), item.rangeSeg.end());
		}
		
		writePCDFile("test.pcd", pCloud, vecRange);
	}
}
