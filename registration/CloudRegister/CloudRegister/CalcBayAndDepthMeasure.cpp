#include "CalcBayAndDepthMeasure.h"

#include "funHelper.h"
#include "Threshold.h"

namespace CloudReg
{
	calcMeassurment_t calcArea(PointCloud::Ptr pWall, const Eigen::Vector4d& plane,
		const Eigen::Vector3d& pt, const std::size_t index)
	{
		const double calcParaZ = TheThreshold::instance()->get_bay_depth_calcParaZ();
		const double calcHalfPara = TheThreshold::instance()->get_bay_depth_calcHalfPara();
		pcl::PointXYZ min;
		pcl::PointXYZ max;
		pcl::getMinMax3D(*pWall, min, max);
		min.z = calcParaZ - calcHalfPara;
		max.z = calcParaZ + calcHalfPara;

		seg_pair_t seg(pt, pt);
		double thickness = 0.05;
		seg.first[2] = min.z;
		seg.second[2] = max.z;

		calcMeassurment_t item;
		auto vecPt = createRulerBox(seg, index, thickness, calcHalfPara * 2);
		item.rangeSeg = calcBoxSegPair(vecPt);
		auto vecTmp = getRulerCorners(vecPt);
		for (auto& pt : vecTmp)
		{
			auto root = pointToPlaneRoot(plane, pt);
			item.rangeSeg.emplace_back(std::make_pair(pt, root));
		}


		auto filerPt = getRulerCorners(vecPt);
		auto pCloud = filerCloudByConvexHull(pWall, filerPt);
		if (pCloud->points.empty())
		{
			LOG(WARNING) << "filerCloudByConvexHull failed";
			return item;
		}
		item.value = 0.0;
		for (auto& pt : pCloud->points)
		{
			item.value += fabs(pointToPLaneDist(plane, pt));
		}
		item.value = item.value / pCloud->points.size();
			
		return item;
	}


	Eigen::Vector3d getCalcPt(const Eigen::Vector3d& sPt, const Eigen::Vector3d& ePt, bool bStart)
	{
		std::size_t optIndex, indexOther;
		int dir;
		
		std::tie(optIndex, indexOther, dir) = getWallGrowAxisAndDir(sPt, ePt);
		Eigen::Vector3d calcPt(0,0,0);
		const double calcHalfPara = 0.005;

		const double calcPara = 0.2;
		if (bStart)
		{
			calcPt = sPt;
			calcPt[optIndex] += dir * calcPara;
		}
		else {
			calcPt = ePt;
			calcPt[optIndex] -= dir * calcPara;
		}

		return calcPt;
	}

	std::vector<seg_pair_t> getOverlapWithHole(const Eigen::Vector3d& sPt, const Eigen::Vector3d ePt,
		const std::vector<vec_seg_pair_t>& vecHole)
	{
		std::vector<seg_pair_t> vecSeg;
		const auto horizenSeg = sPt - ePt;
		Eigen::Vector3d ptA = sPt;
		Eigen::Vector3d ptB = ePt;


		for(std::size_t i = 0; i< vecHole.size(); i++)
		{
			auto& hole = vecHole[i];
			std::vector<seg_pair_t> vecHorizen;
			std::vector<seg_pair_t> vecVertical;
			
			groupDirection(horizenSeg, hole, vecVertical, vecHorizen);

			if (vecHorizen.size() < 2 || vecVertical.size() < 2)
			{
				LOG(WARNING) << "groupDirection failed: " << vecHorizen.size() << " -- " << vecVertical.size();
				continue;
			}
			seg_pair_t toSeg = std::make_pair(ptA, ptB);
			double length1 = (vecHorizen.front().first - vecHorizen.front().second).norm();
			double length2 = (vecHorizen.back().first - vecHorizen.back().second).norm();
			seg_pair_t calcSeg = length1 >= length2 ? vecHorizen.front() : vecHorizen.back();

			if (vecHorizen.front().first[2] < 0.5)
			{

				bool hasOverlap;
				Eigen::Vector3d s1Pt, e1Pt, s2Pt, e2Pt;
				std::tie(hasOverlap, s1Pt, e1Pt, s2Pt, e2Pt) = calcOverlap(toSeg, calcSeg);

				if (hasOverlap)
				{
					if (ptA.norm() != s1Pt.norm())
					{
						LOG(INFO) << "getOverlapWithHole seg length: " << (ptA - s1Pt).norm();
						vecSeg.emplace_back(std::make_pair(ptA, s1Pt));
					}
					ptA = e1Pt;
				}
			}

			if (i == vecHole.size() - 1)
			{
				if(ptB.norm() != ptA.norm())
				{
					LOG(INFO) << "getOverlapWithHole seg length: " << (ptA - ptB).norm();
					vecSeg.emplace_back(std::make_pair(ptA, ptB));
					break;
				}
				
			}
		}

		return vecSeg;
	}

	std::vector<seg_pair_t> mergeOverlap(const std::vector<seg_pair_t>& segA, const std::vector<seg_pair_t>& segB)
	{
		std::vector<seg_pair_t> mergeSeg;
		for (std::size_t i = 0; i < segA.size(); i++)
		{
			seg_pair_t toSeg = segA[i];
		
			
			for (int j = segB.size()-1; j >=0; j--)
			{
				seg_pair_t calcSeg = segB[j];

				bool hasOverlap;
				Eigen::Vector3d s1Pt, e1Pt, s2Pt, e2Pt;
				std::tie(hasOverlap, s1Pt, e1Pt, s2Pt, e2Pt) = calcOverlap(toSeg, calcSeg);

				if (!hasOverlap) continue;
				
				if ((s1Pt - e1Pt).norm() > 0.02)
				{
					mergeSeg.emplace_back(std::make_pair(s1Pt, e1Pt));
					LOG(INFO) << "seg length: " << (s1Pt - e1Pt).norm();

				}
								
			}
		}

		return mergeSeg;
	}

	std::vector<Eigen::Vector3d> getMoveRange(const Eigen::Vector3d& sPt, const Eigen::Vector3d ePt, bool bLeft)
	{
		std::vector<Eigen::Vector3d> vecPts;

		const double MaxMoveTh = 0.2;
		Eigen::Vector3d pt(0, 0, 0);
		Eigen::Vector3d ptLeft(0, 0, 0); //move to edge
		Eigen::Vector3d ptRight(0, 0, 0);//move to center

		std::size_t optIndex, indexOther;
		int dir;
		std::tie(optIndex, indexOther, dir) = getWallGrowAxisAndDir(sPt, ePt);
		
		if (bLeft)
		{
			pt = getCalcPt(sPt, ePt, true);
			ptRight = ptLeft = pt;			
			ptLeft[optIndex] = pt[optIndex] - MaxMoveTh * dir;
			ptRight[optIndex] = pt[optIndex] + MaxMoveTh * dir;

			
		}
		else
		{
			pt = getCalcPt(sPt, ePt, false);
			ptRight = ptLeft = pt;
			ptLeft[optIndex] = pt[optIndex] + MaxMoveTh * dir;
			ptRight[optIndex] = pt[optIndex] - MaxMoveTh * dir;

			
			
		}

		vecPts.emplace_back(ptLeft);
		vecPts.emplace_back(ptRight);
		vecPts.emplace_back(pt);

		return vecPts;
	}

	std::tuple<std::vector<calcIdx2Meassurment_t>,std::vector<seg_pair_t>>
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
		std::vector<calcIdx2Meassurment_t> mapCalcRet;

		if (rootBorder.size() != allWallBorder.size())
		{
			LOG(WARNING) << "the root and wall size not match: " << rootBorder.size() << "-" << allWallBorder.size();
			return std::make_tuple(mapCalcRet, vecCutSeg);
		}
	
		const Eigen::Vector3d& horizenSeg = rootBorder.front().first - rootBorder.front().second;
		std::vector<std::size_t> vecVerticalIndex;
		std::vector<std::size_t> vecHorizenIndex;
		groupDirectionIndex(horizenSeg, rootBorder, vecVerticalIndex, vecHorizenIndex);

		auto rangeSegCalc = [](std::vector<seg_pair_t> OverlapSeg, seg_pair_t& rangeSeg, bool& moveOk) {
			for (auto& seg : OverlapSeg)
			{
				bool hasOverlap;
				Eigen::Vector3d s1Pt, e1Pt, s2Pt, e2Pt;
				std::tie(hasOverlap, s1Pt, e1Pt, s2Pt, e2Pt) = calcOverlap(rangeSeg, seg);

				if (hasOverlap)
				{
					double length = (s1Pt - e1Pt).norm();
					if (length < 0.2)
					{
						auto vecPt = ininterpolateSeg(s1Pt,e1Pt,0.001);
						std::size_t moveStep = (length*1000)/4;
						if (vecPt.size() > moveStep*2)
						{
							s1Pt = vecPt[moveStep];
							e1Pt = vecPt[vecPt.size() - moveStep - 1];
						}
						LOG(INFO) << "rangeSegCalc: add overlap:" << length;

					}
					rangeSeg.first = s1Pt;
					rangeSeg.second = e1Pt;
					moveOk = true;
					LOG(INFO) << "rangeSegCalc: " << (rangeSeg.first - rangeSeg.second).norm();
					break;
				}

			}
		};

		//parallel plane
		auto plane2planeDist = [](const Eigen::Vector4d& planeA, const Eigen::Vector4d& planeB)->double
		{
			auto n = planeA.block<3, 1>(0, 0);
			return fabs(planeA[3] - planeB[3]) / n.norm();
		};
		
		const auto& calcIndex = optType == 0 ? vecHorizenIndex : vecVerticalIndex;

		
		for(std::size_t i = 0; i< calcIndex.size(); i++)
		{
			seg_pair_t toSeg = rootBorder[calcIndex[i]];
			if ((toSeg.first - toSeg.second).norm() < calcLengthTh)
			{
				LOG(INFO) << "type: " << optName << " findSeg:" << calcIndex[i] 
					<< " too short:" << (toSeg.first - toSeg.second).norm();
				continue;
			}

			//if(calcIndex[i] != 0) continue;

			for (int j = calcIndex.size()-1; j >= i+1; j--)
			{
				//if(calcIndex[j] != 6) continue;
				
				seg_pair_t calcSeg = rootBorder[calcIndex[j]];
				if ((calcSeg.first - calcSeg.second).norm() < calcLengthTh)
				{
					LOG(INFO) << "type: " << optName << " findSeg:" << calcIndex[j]
						<< " too short:" << (calcSeg.first - calcSeg.second).norm();
					continue;
				}

				//LOG(INFO) << "type: " << optName << " calcSeg: " << vecHorizenIndex[i] <<" " << vecHorizenIndex[j];
				bool hasOverlap;
				Eigen::Vector3d s1Pt, e1Pt, s2Pt,e2Pt;
				//toSeg and calcSeg has same dir
				std::tie(hasOverlap, s1Pt, e1Pt, s2Pt,e2Pt) = calcOverlap(toSeg,calcSeg);
				
				if (!hasOverlap) continue;

				if ((s1Pt - e1Pt).norm() < calcLengthTh || (s2Pt - e2Pt).norm() < calcLengthTh)
				{
					LOG(INFO) << "type: " << optName << " findSeg:" << calcIndex[i] << " " << calcIndex[j]
						<< " too short: " << (s1Pt - e1Pt).norm() << " " << (s2Pt - e2Pt).norm();
					continue;
				}
				vecCutSeg.emplace_back(std::make_pair(s1Pt,s2Pt));
				vecCutSeg.emplace_back(std::make_pair(e1Pt,e2Pt));
				LOG(INFO)<< "type: "<< optName  << " findSeg:" << calcIndex[i] <<" " << calcIndex[j];

			
				std::vector<seg_pair_t> OverlapSegJ;
				{
					auto iterHole = holeBorder.find(calcIndex[j]);
					if (iterHole != holeBorder.end())
					{
						auto curAllHole = iterHole->second;
						OverlapSegJ = getOverlapWithHole(e2Pt, s2Pt, curAllHole);
						LOG(INFO) << "OverlapSegJ: " << OverlapSegJ.size();
					}
					else
					{
						OverlapSegJ.emplace_back(std::make_pair(e2Pt, s2Pt));
					}
				}

				std::vector<seg_pair_t> OverlapSegI;
				{
					auto iterHole = holeBorder.find(calcIndex[i]);
					if (iterHole != holeBorder.end())
					{
						auto curAllHole = iterHole->second;
						OverlapSegI = getOverlapWithHole(s1Pt, e1Pt, curAllHole);
						LOG(INFO) << "OverlapSegI: " << OverlapSegI.size();
					}
					else
					{
						OverlapSegI.emplace_back(std::make_pair(s1Pt, e1Pt));
					}
				}

				//update OverlapSegJ by OverlapSegI
				std::vector<seg_pair_t> OverlapSeg = mergeOverlap(OverlapSegJ, OverlapSegI);
				LOG(INFO) << "mergeOverlap: " << OverlapSeg.size();

				calcIdx2Meassurment_t save_value;
				save_value.idx = std::make_pair(calcIndex[i], calcIndex[j]);
				
				PointCloud::Ptr pWallI = vecCloud[calcIndex[i]];
				PointCloud::Ptr pWallJ = vecCloud[calcIndex[j]];
				auto& planeI = vecPlane[calcIndex[i]];
				auto& planeJ = vecPlane[calcIndex[j]];

								
				auto movePtandClac = [&](std::vector<seg_pair_t> OverlapSeg,seg_pair_t& rangeSegEdge, 
					seg_pair_t& rangeSegCenter,bool bLeft)->calcMeassurment_t {
					
					const std::string sideName = bLeft == true ? "left side" : "right side";
					std::size_t optIndex, indexOther;
					int dir;
					std::tie(optIndex, indexOther, dir) = getWallGrowAxisAndDir(e2Pt, s2Pt);

					bool moveOk = false;
					Eigen::Vector3d calcPt(0, 0, 0);
					if (!OverlapSeg.empty())
					{
						//try move to edge first
						
						rangeSegCalc(OverlapSeg, rangeSegEdge, moveOk);

						if (!moveOk)
						{	
							//try move to center
							rangeSegCalc(OverlapSeg, rangeSegCenter, moveOk);
							if (moveOk) {
								//center side has overlap 
								if (bLeft) calcPt = rangeSegCenter.first;
								else calcPt = rangeSegCenter.second;
								LOG(INFO) << sideName << " move to center ok";
							}

						}
						else
						{
							//edge side has overlap 
							if(bLeft) calcPt = rangeSegEdge.second;
							else calcPt = rangeSegEdge.first;
							LOG(INFO) << sideName << " use ori pos or move to edge ok";
						}

					}

					if (!moveOk) {
						//use ori pos 0.2
						if (bLeft) calcPt = rangeSegEdge.second;
						else calcPt = rangeSegEdge.first;
						LOG(INFO) << sideName <<" use ori position";
					}


					auto calcPtOther = pointToPlaneRoot(planeI, calcPt);
					auto calcRet = calcArea(pWallJ, planeI, calcPt, indexOther);
					if (calcRet.value < EPS_FLOAT_DOUBLE)
					{
						LOG(INFO) << sideName <<" try to use other";
						calcRet = calcArea(pWallI, planeJ, calcPtOther, indexOther);
						if (calcRet.value < EPS_FLOAT_DOUBLE)
						{
							LOG(INFO) << sideName <<" use plane to plane";
							calcRet.value = plane2planeDist(planeI, planeJ);
						}

					}

					return calcRet;
				};

				//left
				{
					auto moveRange = getMoveRange(e2Pt, s2Pt, true);
					seg_pair_t rangeSegEdge = std::make_pair(moveRange[0], moveRange[2]);
					seg_pair_t rangeSegCenter = std::make_pair(moveRange[2], moveRange[1]);
					auto value = movePtandClac(OverlapSeg, rangeSegEdge, rangeSegCenter,true);
					save_value.vecCalcRet.emplace_back(value);
				}

				//right
				{
					auto moveRange = getMoveRange(e2Pt, s2Pt, false);
					seg_pair_t rangeSegEdge = std::make_pair(moveRange[2], moveRange[0]);
					seg_pair_t rangeSegCenter = std::make_pair(moveRange[1], moveRange[2]);
					auto value = movePtandClac(OverlapSeg, rangeSegEdge, rangeSegCenter, false);
					save_value.vecCalcRet.emplace_back(value);
				}
			
				mapCalcRet.emplace_back(save_value);


			}
		}

		for(auto& value : mapCalcRet)
		{
			
			std::vector<seg_pair_t> vecRange;
			for(auto& item : value.vecCalcRet)
			{
				vecRange.insert(vecRange.end(), item.rangeSeg.begin(), item.rangeSeg.end());
				LOG(INFO) << value.idx.first << " - " << value.idx.second << " :avgDist:" << item.value;
			}

#ifdef VISUALIZATION_ENABLED
			//std::string name = optName + "_wall_" +std::to_string(value.first.first) + "_" + std::to_string(value.first.second) + ".pcd";
			//pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud_filtered(new pcl::PointCloud<pcl::PointXYZ>());
			//uniformSampling(0.01, vecCloud[value.first.first], pCloud_filtered);
			//writePCDFile(name, pCloud_filtered, vecRange);
#endif // DEBUG
		}
		//std::string name = "root_" + optName + ".pcd";
		//writePCDFile(name,rootBorder, vecCutSeg);
		return std::make_tuple(mapCalcRet, vecCutSeg);
	}

	

}//namespace
