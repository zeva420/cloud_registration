#include "CalcNetHeight.h"
#include "funHelper.h"

namespace CloudReg
{

	void calcPlane(PointCloud::Ptr pCloud, const Eigen::Vector3d& sPt, const Eigen::Vector3d& ePt, const Eigen::Vector3d& pt, Eigen::Vector4d& plane)
	{
		auto calcPtOther = pointToPlaneRoot(plane, pt);

		const double calcHalfPara = 0.05;
		const double planeFitDistTh = 0.003;
		pcl::PointXYZ min;
		pcl::PointXYZ max;
		pcl::getMinMax3D(*pCloud, min, max);

		min.x = calcPtOther[0] - calcHalfPara;
		min.y = calcPtOther[1] - calcHalfPara;
		max.x = calcPtOther[0] + calcHalfPara;
		max.y = calcPtOther[1] + calcHalfPara;

		seg_pair_t seg;
		seg.first << min.x, min.y + (max.y - min.y) / 2, min.z;
		seg.second << max.x, min.y + (max.y - min.y) / 2, min.z;
		auto vecPt = createRulerBox(seg, 2, 0.02, calcHalfPara * 2);


		auto filerPt = getRulerCorners(vecPt);
		auto pLeft = filerCloudByConvexHull(pCloud, filerPt);

		if (!pLeft->points.empty())
		{
			Eigen::VectorXf coeff;
			std::vector<int> inlierIdxs;

			planeFitting(planeFitDistTh, pLeft, coeff, inlierIdxs);
			if (inlierIdxs.size() > 50)
			{
				plane[0] = coeff[0];
				plane[1] = coeff[1];
				plane[2] = coeff[2];
				plane[3] = coeff[3];
				LOG(INFO) << "use cloud get new plane";
			}
		}
	}

	calcMeassurment_t calcArea(PointCloud::Ptr pRoof, const Eigen::Vector4d& plane, 
			const Eigen::Vector3d& sPt, const Eigen::Vector3d& ePt, const Eigen::Vector3d& pt,
			bool hasMoreLine)
	{
		const double calcHalfPara = 0.005;
		
		pcl::PointXYZ min;
		pcl::PointXYZ max;		
		pcl::getMinMax3D(*pRoof, min, max);

		min.x = pt[0] - calcHalfPara;
		min.y = pt[1] - calcHalfPara;
		max.x = pt[0] + calcHalfPara;
		max.y = pt[1] + calcHalfPara;
		
		calcMeassurment_t item;

		double thickness = std::fabs(max.z - min.z);	
		seg_pair_t seg;
		seg.first << min.x , min.y + (max.y - min.y)/2 ,min.z;
		seg.second << max.x, min.y + (max.y - min.y)/2 ,min.z;
		auto vecPt = createRulerBox(seg,2,thickness,calcHalfPara*2);
		

		auto filerPt = getRulerCorners(vecPt);	
		auto pCloud = filerCloudByConvexHull(pRoof, filerPt);
		item.rangeSeg = calcBoxSegPair(vecPt);

		if (hasMoreLine)
		{
			auto vecTmp = getRulerCorners(vecPt);
			for (auto& pt : vecTmp)
			{
				auto root = pointToPlaneRoot(plane, pt);
				item.rangeSeg.emplace_back(std::make_pair(pt, root));
			}
			LOG(INFO) << "add more line: " << vecTmp.size();
		}

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
		LOG(INFO) << "clacHeight avgDist:" << item.value;

		//writePCDFile("test.pcd", pCloud, item.rangeSeg);
		return item;
	}

	void moveCheck(Eigen::vector<Eigen::Vector3d>& vecCalcPt, const Eigen::Vector3d& s1Pt, const Eigen::Vector3d& e1Pt,
		const Eigen::Vector3d& s2Pt, const Eigen::Vector3d& e2Pt, const Eigen::Vector3d& centerPt, const double moveRangeTh)
	{
		const double xCenterMax = centerPt[0] + moveRangeTh;
		const double xCenterMin = centerPt[0] - moveRangeTh;
		const double yCenterMax = centerPt[1] + moveRangeTh;
		const double yCenterMin = centerPt[1] - moveRangeTh;
		const double moveExtend = 0.2;

		for (auto& pt : vecCalcPt)
		{
			if (pt[2] - centerPt[2] > 1.0) continue; //not check roof

			double dist = std::sqrtf(std::pow(pt[0] - centerPt[0], 2) + std::pow(pt[1] - centerPt[1], 2));
			if(dist <= moveRangeTh)
			{
				
				double xMax, yMax, xMin, yMin;

				Eigen::Vector3d segA = s1Pt - e1Pt;
				Eigen::Vector3d segB = s1Pt - s2Pt;


				if (std::fabs(segA[0]) > 0.1)
				{
					xMax = s1Pt[0] > e1Pt[0] ? s1Pt[0] : e1Pt[0];
					xMin = s1Pt[0] < e1Pt[0] ? s1Pt[0] : e1Pt[0];
					yMax = s1Pt[1] > s2Pt[1] ? s1Pt[1] : s2Pt[1];
					yMin = s1Pt[1] < s2Pt[1] ? s1Pt[1] : s2Pt[1];

				}
				else
				{
					yMax = s1Pt[1] > e1Pt[1] ? s1Pt[1] : e1Pt[1];
					yMin = s1Pt[1] < e1Pt[1] ? s1Pt[1] : e1Pt[1];
					xMax = s1Pt[0] > s2Pt[0] ? s1Pt[0] : s2Pt[0];
					xMin = s1Pt[0] < s2Pt[0] ? s1Pt[0] : s2Pt[0];
				}

				//std::cout << pt << std::endl;
				//std::cout << centerPt << std::endl;

				if ((xCenterMax + moveExtend) < xMax)
				{
					pt[0] = xCenterMax + moveExtend;
					LOG(INFO) << "moveCheck ok, get new pt";
				}
				else if ((xCenterMin - moveExtend) > xMin)
				{
					pt[0] = xCenterMin - moveExtend;
					LOG(INFO) << "moveCheck ok, get new pt";
				}
				else if ((yCenterMax + moveExtend) < yMax)
				{
					pt[1] = yCenterMax + moveExtend;
					LOG(INFO) << "moveCheck ok, get new pt";
				}
				else if ((yCenterMin - moveExtend) > yMin)
				{
					pt[1] = yCenterMin - moveExtend;
					LOG(INFO) << "moveCheck ok, get new pt";
				}
				else
				{
					LOG(INFO) << "moveCheck failed, use old pt";
				}

			}
		}
		

	}

	Eigen::vector<Eigen::Vector3d> calcRangBox(const Eigen::Vector3d& s1Pt, const Eigen::Vector3d& e1Pt, 
			const Eigen::Vector3d& s2Pt, const Eigen::Vector3d& e2Pt, const Eigen::Vector3d& centerPt,
			const double calcLength, const double moveRangeTh)
	{
		std::size_t g_optIndex,g_indexOther; 
		int g_dir;
		std::tie(g_optIndex,g_indexOther,g_dir) = getWallGrowAxisAndDir(s1Pt,s2Pt);
		
		Eigen::vector<Eigen::Vector3d> vecCalcPt;
		Eigen::Vector3d pt1(s1Pt), pt2(e1Pt);
		{
			std::size_t optIndex,indexOther; 
			int dir;
			std::tie(optIndex,indexOther,dir) = getWallGrowAxisAndDir(s1Pt,e1Pt);

			
			pt1[optIndex] += dir* calcLength;
			pt1[indexOther] += g_dir *calcLength;	
			
			pt2[optIndex] -= dir* calcLength;
			pt2[indexOther] += g_dir *calcLength;	

			
		}
		
		Eigen::Vector3d pt3(s2Pt), pt4(e2Pt);
		{
			std::size_t optIndex,indexOther; 
			int dir;
			std::tie(optIndex,indexOther,dir) = getWallGrowAxisAndDir(s2Pt,e2Pt);

			
			pt3[optIndex] += dir* calcLength;
			pt3[indexOther] -= g_dir *calcLength;	
			
			pt4[optIndex] -= dir* calcLength;
			pt4[indexOther] -= g_dir *calcLength;	

			
		}

		
		vecCalcPt.emplace_back(pt1);
		vecCalcPt.emplace_back(pt2);
		vecCalcPt.emplace_back(pt4);
		vecCalcPt.emplace_back(pt3);


		Eigen::Vector3d  segA = (pt1 - pt4).normalized();
		Eigen::Vector3d  segB = (pt2 - pt3).normalized();
		Eigen::VectorXd line1(6);
		line1 << pt1[0],pt1[1], pt1[2],segA[0],segA[1],segA[2];
		
		Eigen::VectorXd line2(6);
		line2 << pt2[0], pt2[1], pt2[2],segB[0],segB[1],segB[2];

		Eigen::Vector3d interSectionPt;
		if (!interSectionOfLineToLine(line1, line2, interSectionPt))
		{
			LOG(ERROR) << "calc interSectionOfLineToLine failed";
			return vecCalcPt;
		}
		vecCalcPt.emplace_back(interSectionPt);

		moveCheck(vecCalcPt, s1Pt, e1Pt, s2Pt, e2Pt, centerPt, moveRangeTh);
		

		return vecCalcPt;
	}

	void calcAvgDiff(std::vector<calcIdx2Meassurment_t>& vecRet)
	{
	
		for (auto& value : vecRet)
		{
			double minValue = 9999.0;
			for (std::size_t i = 0; i < value.vecCalcRet.size(); i++)
			{
				if (value.vecCalcRet[i].value > 0.0 && value.vecCalcRet[i].value < minValue)
					minValue = value.vecCalcRet[i].value;

				//LOG(INFO) << "orig value:" << vecRet.vecCalcRet[i].value;
			}

			if (minValue > 0.0 && minValue < 9999.0) {
				std::for_each(value.vecCalcRet.begin(), value.vecCalcRet.end(), [&](calcMeassurment_t& item) {
					if(item.value > 0) item.value -= minValue; });
			}

		}

	}

	std::tuple<std::vector<calcIdx2Meassurment_t>,std::vector<seg_pair_t >>
		CalcNetHeight(const std::vector<seg_pair_t>& roofBorder,
			const PointCloud::Ptr pCloud,
			const PointCloud::Ptr pOther,
			const Eigen::Vector4d& plane,
			const Eigen::Vector3d& center,
			const std::string& name,
			const double calcLengthTh,
			const double moveRangeTh,
			bool hasMoreLine,
			bool getNewPlane)
	{
		const Eigen::Vector3d& horizenSeg = roofBorder.front().first - roofBorder.front().second;
		std::vector<std::size_t> vecVerticalIndex;
		std::vector<std::size_t> vecHorizenIndex;
		groupDirectionIndex(horizenSeg, roofBorder, vecVerticalIndex, vecHorizenIndex);
		
		LOG(INFO) << "vecVerticalIndex: " << vecVerticalIndex.size() << " vecHorizenIndex:" << vecHorizenIndex.size()
			<< " " << hasMoreLine;

		std::vector<seg_pair_t> vecCutSeg;
		std::vector<calcIdx2Meassurment_t> vecRet;
		for(std::size_t i = 0; i< vecVerticalIndex.size(); i++)
		{
			seg_pair_t toSeg = roofBorder[vecVerticalIndex[i]];
			if ((toSeg.first - toSeg.second).norm() < calcLengthTh)
				continue;

			for (int j = vecVerticalIndex.size() - 1; j >= i + 1; j--)
			{
				seg_pair_t calcSeg = roofBorder[vecVerticalIndex[j]];
				if ((calcSeg.first - calcSeg.second).norm() < calcLengthTh)
					continue;

				//LOG(INFO) << "begin calc seg:" << vecVerticalIndex[i] << " "
				//	<< vecVerticalIndex[j];

				bool hasOverlap;
				Eigen::Vector3d s1Pt, e1Pt, s2Pt,e2Pt;
				std::tie(hasOverlap, s1Pt, e1Pt, s2Pt,e2Pt) = calcOverlap(toSeg,calcSeg);
				
				if (!hasOverlap) continue;
				
			
				if((s1Pt-e1Pt).norm() < calcLengthTh 
					|| (s2Pt-e2Pt).norm() < calcLengthTh)
				{
					LOG(INFO)<< "too short seg:" << vecVerticalIndex[i] <<" " 
						<< vecVerticalIndex[j] << " :" << (s1Pt -e1Pt).norm() 
						<<" " << (s2Pt-e2Pt).norm();
					continue;
				}
			
				vecCutSeg.emplace_back(std::make_pair(s1Pt,s2Pt));
				vecCutSeg.emplace_back(std::make_pair(e1Pt,e2Pt));

				LOG(INFO) << "findSeg:" << vecVerticalIndex[i] <<" " << vecVerticalIndex[j];
				auto vecCalcPt = calcRangBox(s1Pt,e1Pt,s2Pt,e2Pt, center,0.3, moveRangeTh);
				std::vector<calcMeassurment_t> tmp;
				for (auto& pt : vecCalcPt)
				{
					Eigen::Vector4d otherPlane = plane;
					if (getNewPlane)
						calcPlane(pOther, s1Pt, e1Pt, pt, otherPlane);
					
					auto ret = calcArea(pCloud, otherPlane, s1Pt, e1Pt, pt,hasMoreLine);
					if (!ret.rangeSeg.empty())
					{
						tmp.emplace_back(ret);
						//vecCutSeg.insert(vecCutSeg.end(), ret.rangeSeg.begin(), ret.rangeSeg.end());
					}
				}
				if (!tmp.empty())
				{
					calcIdx2Meassurment_t ret;
					ret.vecCalcRet.swap(tmp);
					ret.idx = std::make_pair(vecVerticalIndex[i], vecVerticalIndex[j]);
					vecRet.emplace_back(ret);
				}
					

			}

		}
		//writePCDFile(name,roofBorder, vecCutSeg);
		
		return std::make_tuple(vecRet, vecCutSeg);
	}
	
	std::tuple<std::vector<calcIdx2Meassurment_t>, std::vector<calcIdx2Meassurment_t>,
		std::vector<seg_pair_t>, std::vector<seg_pair_t>>
		CalcHeightRange(const std::vector<seg_pair_t>& roofBorder,
			const std::vector<seg_pair_t>& rootBorder,
			const std::vector<std::vector<seg_pair_t>>& allWallBorder,
			const PointCloud::Ptr pRoof,
			const PointCloud::Ptr pRoot,
			const Eigen::Vector3d& center,
			const double calcHeight,
			const double calcLengthTh,
			const double moveRangeTh)
	{
		
		std::vector<Eigen::Vector3d> vecPt;
		Eigen::Vector3d calcPt = allWallBorder.front().front().first;
		calcPt[2] += calcHeight;
		vecPt.emplace_back(calcPt);	
		for(std::size_t i = 1; i< allWallBorder.size(); i++)
		{
			auto& seg = allWallBorder[i].front();
			auto rootPt = calcPerpendicular(calcPt, seg.first, seg.second);
			vecPt.emplace_back(rootPt);
			calcPt = rootPt;
		}
		
		pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud(new pcl::PointCloud<pcl::PointXYZ>());
		for (size_t i = 0; i < vecPt.size(); ++i)
		{
			pcl::PointXYZ p2;
			p2.x = vecPt[i][0];
			p2.y = vecPt[i][1];
			p2.z = vecPt[i][2];
			pCloud->push_back(p2);
		}
		Eigen::Vector4d calcPlane = calcPlaneParam(pCloud);
	
		std::vector<seg_pair_t> vecSeg1, vecSeg2;
		for (std::size_t i = 1; i < vecPt.size(); i++) {
			vecSeg1.emplace_back(seg_pair_t(vecPt[i - 1], vecPt[i]));

		}
		vecSeg1.emplace_back(seg_pair_t(vecPt.back(), vecPt.front()));
		vecSeg2 = vecSeg1;
		//writePCDFile("rangePlane.pcd",rootBorder, vecSeg);
		
		
		std::vector<calcIdx2Meassurment_t> vecroofRet;
		{
			const std::string name = "roof_height_range.pcd";			
			std::vector<seg_pair_t> cutSeg;
			std::tie(vecroofRet, cutSeg) = CalcNetHeight(roofBorder,pRoof, pRoot, calcPlane, center, name,calcLengthTh, moveRangeTh);
			calcAvgDiff(vecroofRet);
			vecSeg1.insert(vecSeg1.end(),cutSeg.begin(), cutSeg.end());
		}
		
		std::vector<calcIdx2Meassurment_t> vecrootRet;
		{
			const std::string name = "root_height_range.pcd";
			std::vector<seg_pair_t> cutSeg;
			std::tie(vecrootRet, cutSeg) = CalcNetHeight(rootBorder,pRoot, pRoof, calcPlane, center, name, calcLengthTh, moveRangeTh);
			calcAvgDiff(vecrootRet);
			vecSeg2.insert(vecSeg2.end(), cutSeg.begin(), cutSeg.end());
		}
		return std::make_tuple(vecroofRet,vecrootRet, vecSeg1, vecSeg2);

	}
}//namespace
