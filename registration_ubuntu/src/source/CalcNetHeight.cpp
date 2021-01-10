#include "CalcBayAndDepthMeasure.h"

#include "funHelper.h"

namespace CloudReg
{

	calcMeassurment_t calcArea(PointCloud::Ptr pRoof, const Eigen::Vector4d& plane, 
			const Eigen::Vector3d& sPt, const Eigen::Vector3d& ePt, const Eigen::Vector3d& pt)
	{
		const double calcHalfPara = 0.05;
		
		pcl::PointXYZ min;
		pcl::PointXYZ max;		
		pcl::getMinMax3D(*pRoof, min, max);

		min.x = pt[0] - calcHalfPara;
		min.y = pt[1] - calcHalfPara;
		max.x = pt[0] + calcHalfPara;
		max.y = pt[1] + calcHalfPara;
		
		auto pCloud = filerCloudByRange(pRoof,min,max);
		calcMeassurment_t item;
		item.value = 0.0;
		for(auto& pt : pCloud->points)
		{
			item.value += fabs(pointToPLaneDist(plane,pt));
		}
		item.value = item.value/pCloud->points.size();

		double thickness = std::fabs(max.z - min.z);	
		seg_pair_t seg;
		seg.first << min.x , min.y + (max.y - min.y)/2 ,min.z;
		seg.second << max.x, min.y + (max.y - min.y)/2 ,min.z;
		auto vecPt = createRulerBox(seg,2,thickness,calcHalfPara*2);
		item.rangeSeg = calcBoxSegPair(vecPt);
		
		//LOG(INFO) << "avgDist:" << item.value;
		//writePCDFile("test.pcd", pCloud, item.rangeSeg);
		return item;
	}
	Eigen::vector<Eigen::Vector3d> calcRangBox(const Eigen::Vector3d& s1Pt, const Eigen::Vector3d& e1Pt, 
			const Eigen::Vector3d& s2Pt, const Eigen::Vector3d& e2Pt, const double calcLength)
	{
		std::size_t g_optIndex,g_indexOther; 
		int g_dir;
		std::tie(g_optIndex,g_indexOther,g_dir) = getGrowAxisAndDir(s1Pt,s2Pt);
		
		Eigen::vector<Eigen::Vector3d> vecCalcPt;
		{
			std::size_t optIndex,indexOther; 
			int dir;
			std::tie(optIndex,indexOther,dir) = getGrowAxisAndDir(s1Pt,e1Pt);

			Eigen::Vector3d pt1(s1Pt),pt2(e1Pt);
			pt1[optIndex] += dir* calcLength;
			pt1[indexOther] += g_dir *calcLength;	
			
			pt2[optIndex] -= dir* calcLength;
			pt2[indexOther] += g_dir *calcLength;	

			vecCalcPt.emplace_back(pt1);
			vecCalcPt.emplace_back(pt2);
		}
		
		{
			std::size_t optIndex,indexOther; 
			int dir;
			std::tie(optIndex,indexOther,dir) = getGrowAxisAndDir(s2Pt,e2Pt);

			Eigen::Vector3d pt3(s2Pt),pt4(e2Pt);
			pt3[optIndex] += dir* calcLength;
			pt3[indexOther] -= g_dir *calcLength;	
			
			pt4[optIndex] -= dir* calcLength;
			pt4[indexOther] -= g_dir *calcLength;	

			vecCalcPt.emplace_back(pt3);
			vecCalcPt.emplace_back(pt4);
		}

		Eigen::Vector3d  segA = (vecCalcPt[0] - vecCalcPt[3]).normalized();
		Eigen::Vector3d  segB = (vecCalcPt[1] - vecCalcPt[2]).normalized();
		Eigen::VectorXd line1(6);
		line1 << vecCalcPt[0][0],vecCalcPt[0][1],vecCalcPt[0][2],segA[0],segA[1],segA[2]; 
		
		Eigen::VectorXd line2(6);
		line2 << vecCalcPt[1][0],vecCalcPt[1][1],vecCalcPt[1][2],segB[0],segB[1],segB[2];

		Eigen::Vector3d interSectionPt;
		if (!interSectionOfLineToLine(line1, line2,interSectionPt))
			LOG(ERROR) << "calc interSectionOfLineToLine failed";
		else	
			vecCalcPt.emplace_back(interSectionPt);

		return vecCalcPt;
	}

	void CalcNetHeight(const std::vector<seg_pair_t>& roofBorder,
			const PointCloud::Ptr pCloud,
			const Eigen::Vector4d& plane,
			const double calcLengthTh)
	{
		const Eigen::Vector3d& horizenSeg = roofBorder.front().first - roofBorder.front().second;
		std::vector<std::size_t> vecVerticalIndex;
		std::vector<std::size_t> vecHorizenIndex;
		groupDirectionIndex(horizenSeg, roofBorder, vecVerticalIndex, vecHorizenIndex);
		
		std::vector<seg_pair_t> vecCutSeg;
		std::vector<calcMeassurment_t> vecRet;
		for(std::size_t i = 0; i< vecVerticalIndex.size(); i++)
		{
			seg_pair_t toSeg = roofBorder[vecVerticalIndex[i]];
			if ((toSeg.first - toSeg.second).norm() < calcLengthTh)
				continue;

			for(std::size_t j = i+1; j< vecVerticalIndex.size(); j++)
			{
				seg_pair_t calcSeg = roofBorder[vecVerticalIndex[j]];
				if ((calcSeg.first - calcSeg.second).norm() < calcLengthTh)
					continue;
				
				bool hasOverlap;
				Eigen::Vector3d s1Pt, e1Pt, s2Pt,e2Pt;
				std::tie(hasOverlap, s1Pt, e1Pt, s2Pt,e2Pt) = calcOverlap(toSeg,calcSeg);
				
				if (!hasOverlap) continue;
				
			
				if((s1Pt-e1Pt).norm() < calcLengthTh 
					|| (s2Pt-e2Pt).norm() < calcLengthTh)
				{
					LOG(INFO)<< "too short seg:" << vecVerticalIndex[i] <<" " << vecVerticalIndex[j];
					continue;
				}
			
				vecCutSeg.emplace_back(std::make_pair(s1Pt,s2Pt));
				vecCutSeg.emplace_back(std::make_pair(e1Pt,e2Pt));

				LOG(INFO)<< "findSeg:" << vecVerticalIndex[i] <<" " << vecVerticalIndex[j];
			
				auto vecCalcPt = calcRangBox(s1Pt,e1Pt,s2Pt,e2Pt,0.3);
				for(auto& pt : vecCalcPt)
				{
					auto ret = calcArea(pCloud, plane, s1Pt, e1Pt, pt);
					vecRet.emplace_back(ret);
					vecCutSeg.insert(vecCutSeg.end(),ret.rangeSeg.begin(), ret.rangeSeg.end());
				}
			}

		}
		std::string name = "roof_net_height.pcd";
		writePCDFile(name,roofBorder, vecCutSeg);
	}
}//namespace
