#include "CalcHoleMeasure.h"

#include "funHelper.h"

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
	
		const double halfHeight = fabs(segB.first[2] - segB.second[2])/2;
		segA.first[2] = segB.first[2] < segB.second[2] ? segB.first[2] + halfHeight : segB.second[2] + halfHeight;
		segA.second[2] = segA.first[2];
		std::size_t optIndex,indexOther; 
		int dir;
		std::tie(optIndex,indexOther,dir) = getWallGrowAxisAndDir(segA.first,segA.second);
		segA.first[optIndex] -= extendRange* dir;
		segA.second[optIndex] += extendRange* dir;
		
		double width = (segB.first - segB.second).norm() + extendRange * 2; 
		auto vecPt = createRulerBox(segA,indexOther,0.1,width);
		vecFilerPt = getRulerCorners(vecPt);	

		//writePCDFile("test.pcd",std::vector<seg_pair_t>(), calcBoxSegPair(vecPt));
	}



	std::vector<calcMeassurment_t> calcHoleCross(const std::vector<seg_pair_t>& holeBorder, const PointCloud::Ptr pCloud)
	{
		std::vector<calcMeassurment_t> vecRet;
				
		for (std::size_t i = 0; i < holeBorder.size(); i++)
		{
			seg_pair_t toSeg = holeBorder[i];

			for (std::size_t j = i + 1; j < holeBorder.size(); j++)
			{
				seg_pair_t calcSeg = holeBorder[j];
				
				bool hasOverlap;
				Eigen::Vector3d s1Pt, e1Pt, s2Pt, e2Pt;
				std::tie(hasOverlap, s1Pt, e1Pt, s2Pt, e2Pt) = calcOverlap(toSeg, calcSeg);

				if (!hasOverlap) continue;
				

				if ((s1Pt - e1Pt).norm() < 0.01 || (s2Pt - e2Pt).norm() < 0.01) continue;
				LOG(INFO) << "calc betweenn: " << i << " " << j;

				Eigen::vector<Eigen::Vector3d> vecCalcPt;
				vecCalcPt.emplace_back(s1Pt);
				vecCalcPt.emplace_back(e1Pt);
				vecCalcPt.emplace_back(s2Pt);
				vecCalcPt.emplace_back(e2Pt);
				
				auto vecRawPt = getNearestPt(vecCalcPt, pCloud, 0.1*0.1);

				{
					calcMeassurment_t item;
					item.value = (vecRawPt[0] - vecRawPt[3]).norm();
					item.rangeSeg.emplace_back(std::make_pair(vecRawPt[0], vecRawPt[3]));
					vecRet.emplace_back(item);
					LOG(INFO) << "cross dist :" << item.value;
				}

				{
					calcMeassurment_t item;
					item.value = (vecRawPt[1] - vecRawPt[2]).norm();
					item.rangeSeg.emplace_back(std::make_pair(vecRawPt[1], vecRawPt[2]));
					vecRet.emplace_back(item);
					LOG(INFO) << "cross dist :" << item.value;
				}

			}
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

	std::vector<std::vector<calcMeassurment_t>> calcHole(const seg_pair_t& horizen,
			const std::vector<seg_pair_t>& holeBorder,
			pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud,
			const std::string& name)
	{

		std::vector<std::vector<calcMeassurment_t>> vecRet;
		std::vector<seg_pair_t> vecHorizen;	
		std::vector<seg_pair_t> vecVertical;
		const auto horizenSeg = horizen.first - horizen.second;
		groupDirection(horizenSeg,holeBorder, vecVertical, vecHorizen);

		if (vecHorizen.size() < 2 || vecVertical.size() < 2)
		{
			LOG(ERROR) << "groupDirection failed: " << vecHorizen.size() << " -- " << vecVertical.size();
			return vecRet;
		}

		std::vector<Eigen::Vector3d> vecFilerPt;
		calcHoleMaxAndMin(vecHorizen,vecVertical, 0.05, vecFilerPt);
		auto rangeCloud = filerCloudByConvexHull(pCloud,vecFilerPt);
		if (rangeCloud->points.empty()) 
		{
			LOG(ERROR) << "filerCloudByRange failed";
			return vecRet;
		}

		//pcl::io::savePCDFile("cloud.pcd", *pCloud);
		//pcl::io::savePCDFile("rangeCloud.pcd", *rangeCloud);

		std::vector<seg_pair_t> vecRange;
		
		//height
		{
			LOG(INFO) << "calc hole height";
			Eigen::Vector3d A1 = vecHorizen.front().first - vecHorizen.front().second;
			Eigen::Vector3d A2 = vecHorizen.back().first - vecHorizen.back().second;

			std::size_t baseIndex =  (A2.norm() - A1.norm()) > 0.1 ?  vecHorizen.size()-1 : 0;
			const seg_pair_t& baseSeg = vecHorizen[baseIndex];
		
			for(std::size_t i = 0; i < vecHorizen.size(); i++)
			{
				if (i == baseIndex) continue;

				const auto& seg = vecHorizen[i];
				auto vecPts = ininterpolateSeg(seg.first, seg.second, 0.01f);
				auto vecHeight = calcHoleDist(baseSeg, vecPts, rangeCloud, 15, true);
				vecRet.emplace_back(vecHeight);
				for(auto& item : vecHeight)
					vecRange.insert(vecRange.end(), item.rangeSeg.begin(), item.rangeSeg.end());
			}
		}

		//width
		{
			LOG(INFO) << "calc hole width";
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
				vecRet.emplace_back(vecWidth);
				
				for(auto& item : vecWidth)
					vecRange.insert(vecRange.end(), item.rangeSeg.begin(), item.rangeSeg.end());
			}
		}

		//cross
		{
			LOG(INFO) << "calc hole cross";
			auto vecCross = calcHoleCross(vecVertical,rangeCloud);
			vecRet.emplace_back(vecCross);
			for(auto& item : vecCross)
				vecRange.insert(vecRange.end(), item.rangeSeg.begin(), item.rangeSeg.end());
		}

		//to root
		if (vecHorizen.front().first[2] > 0.1)
		{
			LOG(INFO) << "calc hole root";
			const auto& seg = vecHorizen.front();
			auto vecPts = ininterpolateSeg(seg.first, seg.second, 0.01f);
			auto vecRoot = calcHoleDist(horizen, vecPts, rangeCloud, 15, true);
			vecRet.emplace_back(vecRoot);
		
			for(auto& item : vecRoot)
				vecRange.insert(vecRange.end(), item.rangeSeg.begin(), item.rangeSeg.end());
		}
#ifdef VISUALIZATION_ENABLED		
		pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud_filtered(new pcl::PointCloud<pcl::PointXYZ>());
		uniformSampling(0.01, pCloud, pCloud_filtered);
		writePCDFile(name, pCloud_filtered, vecRange);
#endif

		return vecRet;
	}
}
