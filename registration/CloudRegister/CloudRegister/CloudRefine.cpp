#include "CloudRefine.h"
#include "funHelper.h"
#include "CalcMeasureHelper.h"

namespace CloudReg
{
	void splitePt(PointCloud::Ptr pCloudA, const Eigen::Vector4d& planeA,
		PointCloud::Ptr pCloudB, const Eigen::Vector4d& planeB)
	{
		LOG(INFO) << "before:" << pCloudA->points.size() << " -- " << pCloudB->points.size();

		pcl::PointCloud<pcl::PointXYZ>::Ptr subSetMerged(new pcl::PointCloud<pcl::PointXYZ>());
		subSetMerged->insert(subSetMerged->end(), pCloudA->begin(), pCloudA->end());
		subSetMerged->insert(subSetMerged->end(), pCloudB->begin(), pCloudB->end());

		Eigen::VectorXf coeff1;
		std::vector<int> inlierIdxs1;
		planeFitting(0.002, subSetMerged, coeff1, inlierIdxs1);
		auto left1 = geo::getSubSet(subSetMerged, inlierIdxs1, true);

		Eigen::VectorXf coeff2;
		std::vector<int> inlierIdxs2;
		planeFitting(0.002, left1, coeff2, inlierIdxs2);
		auto left2 = geo::getSubSet(left1, inlierIdxs2, true);

		auto orig_inliers1 = geo::getSubSet(subSetMerged, inlierIdxs1, false);
		auto orig_inliers2 = geo::getSubSet(left1, inlierIdxs2, false);
		LOG(INFO) << "after plane " << orig_inliers1->points.size() << " -- "
			<< orig_inliers2->points.size();

		auto dist_to_plane = [](Eigen::Vector3f &p, Eigen::VectorXf &plane)->float {
			Eigen::Vector3f n = plane.block<3, 1>(0, 0);
			float d = plane(3);
			float dist = std::abs(n.dot(p) + d) / n.norm();
			return dist;
		};

		for (auto &it : left2->points)
		{
			Eigen::Vector3f p(it.x, it.y, it.z);
			float dist1 = dist_to_plane(p, coeff1);
			float dist2 = dist_to_plane(p, coeff2);
			auto tmpPtr = (dist1 < dist2) ? orig_inliers1 : orig_inliers2;
			tmpPtr->push_back(it);
		}

		{
			Eigen::Vector3d nA(std::fabs(planeA[0]), std::fabs(planeA[1]), std::fabs(planeA[2]));
			Eigen::Vector3d nANew(std::fabs(coeff1[0]), std::fabs(coeff1[1]), std::fabs(coeff1[2]));
			
			
			Eigen::Vector3d c12 = nA.cross(nANew);
			if (c12.norm() <= 0.1)
			{
				pCloudA->swap(*orig_inliers1);
				pCloudB->swap(*orig_inliers2);
			}
			else
			{
				pCloudA->swap(*orig_inliers2);
				pCloudB->swap(*orig_inliers1);
			}
		}
		LOG(INFO) << "after:" << pCloudA->points.size() << " -- " << pCloudB->points.size();
	}

	void refineTop(std::map<CloudItemType, vecItems_t>& ret)
	{
				
		//b is roof
		auto refineOperate = [](const std::vector<seg_pair_t>& calcBorderA,const std::vector<seg_pair_t>& calcBorderB, 
			PointCloud::Ptr pCloudA, PointCloud::Ptr pCloudB, const Eigen::Vector4d& planeA, const Eigen::Vector4d& planeB, 
			const seg_pair_t& calcSeg, bool bTop = true){

			const double calcLength = 0.1;

			PointCloud::Ptr pRefineA, pLeftA;
			{
				std::vector<seg_pair_t> vecVertical;
				std::vector<seg_pair_t> vecHorizen; 
				groupDirection((calcSeg.first - calcSeg.second), calcBorderA, vecVertical, vecHorizen);

				Eigen::Vector3d sPt, ePt;
				if (bTop)
				{
					sPt = vecHorizen.back().first;
					ePt = vecHorizen.back().second;
				}
				else
				{
					sPt = vecHorizen.front().first;
					ePt = vecHorizen.front().second;
				}
				
				std::size_t optIndex, indexOther;
				int dir;
				std::tie(optIndex, indexOther, dir) = getWallGrowAxisAndDir(sPt, ePt);

				auto vecPt = createRulerBox(vecHorizen.back(), indexOther, 0.05, calcLength * 2);
				auto filerPt = getRulerCorners(vecPt);
				pRefineA = filerCloudByConvexHull(pCloudA, filerPt);
				pLeftA = filerCloudByConvexHull(pCloudA, filerPt, false);
				
				
				/*auto rangeSeg = calcBoxSegPair(vecPt);
				pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud_filtered(new pcl::PointCloud<pcl::PointXYZ>());
				uniformSampling(0.01, pRefineA, pCloud_filtered);
				writePCDFile("testA1.pcd", pCloud_filtered, rangeSeg);
				uniformSampling(0.01, pLeftA, pCloud_filtered);
				writePCDFile("testA2.pcd", pCloud_filtered, rangeSeg);*/

			}
			
			PointCloud::Ptr pRefineB, pLeftB;
			{
				
				auto vecPt = createRulerBox(calcSeg, 2, 0.05, calcLength * 2);
				auto filerPt = getRulerCorners(vecPt);
				pRefineB = filerCloudByConvexHull(pCloudB, filerPt);
				pLeftB = filerCloudByConvexHull(pCloudB, filerPt, false);

				/*auto rangeSeg = calcBoxSegPair(vecPt);
				pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud_filtered(new pcl::PointCloud<pcl::PointXYZ>());
				uniformSampling(0.01, pRefineB, pCloud_filtered);
				writePCDFile("testB1.pcd", pCloud_filtered, rangeSeg);
				uniformSampling(0.01, pLeftB, pCloud_filtered);
				writePCDFile("testB2.pcd", pCloud_filtered, rangeSeg);*/

			}

			splitePt(pRefineA, planeA, pRefineB, planeB);
			pRefineA->insert(pRefineA->end(), pLeftA->begin(), pLeftA->end());
			pRefineB->insert(pRefineB->end(), pLeftB->begin(), pLeftB->end());

			pCloudA->swap(*pRefineA);
			pCloudB->swap(*pRefineB);

		};

		auto& roof = ret[CLOUD_TOP_E].front();
		auto& roofBorder = roof.cloudBorder_.front();
		auto& vecBeam = ret[CLOUD_BEAM_E];
		auto& vecWall = ret[CLOUD_WALL_E];

		std::set<std::size_t> checked;


		//beam refine
		for (std::size_t i = 0; i < vecBeam.size(); i+=2)
		{
			auto& beamA = vecBeam[i];
			auto& beamBorderA = beamA.cloudBorder_.front();
			

			//beam to top
			{
				const seg_pair_t& calcSeg = roofBorder[beamA.parentIndex_];
				refineOperate(beamBorderA, roofBorder, beamA.pCloud_, roof.pCloud_,
					beamA.cloudPlane_, roof.cloudPlane_, calcSeg);
				checked.insert(beamA.parentIndex_);
			}


			//beam to wall
			auto& wall = vecWall[beamA.parentIndex_];
			auto& wallBorder = wall.cloudBorder_.front();
			auto& beamB = vecBeam[i + 1];
			auto& beamBorderB = beamB.cloudBorder_.front();
			{
				const seg_pair_t& calcSeg = beamBorderB.back();				
				refineOperate(wallBorder, beamBorderB, wall.pCloud_, beamB.pCloud_,
					wall.cloudPlane_, beamB.cloudPlane_, calcSeg);
				
			}
						
		}

		//wall refine
		for (std::size_t i = 0; i < vecWall.size(); i ++)
		{
			if (checked.count(i)) continue;

			auto& wall = vecWall[i];
			auto& wallBorder = wall.cloudBorder_.front();

			//wall to top
			{
				const seg_pair_t& calcSeg = roofBorder[i];
				refineOperate(wallBorder, roofBorder, wall.pCloud_, roof.pCloud_,
					wall.cloudPlane_, roof.cloudPlane_, calcSeg);
			}
			
		}
		
		
		auto pNew = refineBySegment(roof.cloudBorder_.front(), roof.pCloud_);
		roof.pCloud_->swap(*pNew);
	}

	void refindWall(std::map<CloudItemType, vecItems_t>& ret)
	{
		auto& vecWall = ret[CLOUD_WALL_E];
		std::vector<std::pair<std::size_t, std::size_t>> calcIdx;
		for (std::size_t i = 1; i < vecWall.size(); i++)
		{
			calcIdx.emplace_back(std::make_pair(i - 1, i));
			if (i == vecWall.size() - 1)
				calcIdx.emplace_back(std::make_pair(i, 0));
		}

		auto getRangeCloud = [](const std::vector<seg_pair_t>& boder, PointCloud::Ptr pCloud, bool bLeft)->std::tuple<PointCloud::Ptr, PointCloud::Ptr> {
			const double calcLength = 0.1;
			PointCloud::Ptr pRefine, pLeft;

			seg_pair_t seg = boder.back();
			std::vector<seg_pair_t> vecVertical;
			std::vector<seg_pair_t> vecHorizen;
			groupDirection(seg.first - seg.second, boder, vecVertical, vecHorizen);

			seg_pair_t calcSeg;

			if (bLeft)
				calcSeg = vecVertical.back();
			else
				calcSeg = vecVertical.front();

			
			std::size_t optIndex, indexOther;
			int dir;
			std::tie(optIndex, indexOther, dir) = getWallGrowAxisAndDir(seg.first, seg.second);

			auto vecPt = createRulerBox(calcSeg, indexOther, 0.05, calcLength * 2);
			auto filerPt = getRulerCorners(vecPt);
			pRefine = filerCloudByConvexHull(pCloud, filerPt);
			pLeft = filerCloudByConvexHull(pCloud, filerPt, false);

			LOG(INFO)<<pCloud->points.size() << " -- " << pRefine->points.size() << " -- " << pLeft->points.size();
			return std::make_tuple(pRefine, pLeft);
		};
		
		for (auto& idx : calcIdx)
		{
			auto& leftWall = vecWall[idx.first];
			auto& rightWall = vecWall[idx.second];
		
			//left
			PointCloud::Ptr pRefineA, pLeftA;
			std::tie(pRefineA, pLeftA) = getRangeCloud(leftWall.cloudBorder_.front(), 
					leftWall.pCloud_,true);		
			

			//right
			PointCloud::Ptr pRefineB, pLeftB;
			std::tie(pRefineB, pLeftB) = getRangeCloud(rightWall.cloudBorder_.front(),
				rightWall.pCloud_, false);

			splitePt(pRefineA, leftWall.cloudPlane_, pRefineB, rightWall.cloudPlane_);
			pRefineA->insert(pRefineA->end(), pLeftA->begin(), pLeftA->end());
			pRefineB->insert(pRefineB->end(), pLeftB->begin(), pLeftB->end());

			leftWall.pCloud_->swap(*pRefineA);
			rightWall.pCloud_->swap(*pRefineB);
			

		}

		for (auto& wall : vecWall)
		{
			std::vector<Eigen::Vector3d> vecPts;
			for (auto& seg : wall.cloudBorder_.front())
			{
				vecPts.emplace_back(seg.first);
			}
			auto pNew = filerCloudByConvexHull(wall.pCloud_, vecPts);
			wall.pCloud_->swap(*pNew);
		}
		
	}
	bool CloudRefine::run(std::map<CloudItemType, vecItems_t>& ret)
	{
		LOG(INFO) << "begin refine cloud";

		//refine cut
		for(auto& value: ret)
		{
			auto vecItem = value.second;

			if (value.first == CLOUD_WALL_E)
			{
				for (auto& wall : vecItem)
				{
					std::vector<Eigen::Vector3d> vecPts;
					for (auto& seg : wall.cloudBorder_.front())
					{
						vecPts.emplace_back(seg.first);
					}
					auto pNew = filerCloudByConvexHull(wall.pCloud_, vecPts);
					wall.pCloud_->swap(*pNew);

					for (std::size_t i = 1; i < wall.cloudBorder_.size(); i++)
					{
						if (wall.cloudBorder_[i].size() == 4)
						{
							std::vector<Eigen::Vector3d> vecPts;
							for (auto& seg : wall.cloudBorder_[i])
							{
								vecPts.emplace_back(seg.first);
							}
							auto pNew = filerCloudByConvexHull(wall.pCloud_, vecPts, false);
							wall.pCloud_->swap(*pNew);
						}
						else {
							refineByHole(wall.cloudBorder_[i], wall.pCloud_);
						}
					}
				}

								
			}
			else if (value.first == CLOUD_BOTTOM_E || value.first == CLOUD_TOP_E)
			{
				auto pNew = refineBySegment(vecItem.front().cloudBorder_.front(), vecItem.front().pCloud_);
				vecItem.front().pCloud_->swap(*pNew);
			}
		}

		//refineTop(ret);
		//refindWall(ret);

		LOG(INFO) << "finished refine cloud";
		return true;
	}
}