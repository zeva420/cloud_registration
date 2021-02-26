#include "CloudBorder.h"
#include "funHelper.h"
#include "CalcMeasureHelper.h"
#include "glog/logging.h"

#define VISUALIZATION_ENABLED
namespace CloudReg
{

std::map<size_t, std::map<size_t, std::vector<Eigen::Vector3d>>> CloudBorder::calcWallNodes(
		const std::string &name, 
		const PointCloud::Ptr cloud, 
		const Eigen::Vector4d &cloudPlane,
		const std::vector<std::vector<seg_pair_t>> &cadBorder,
		const Eigen::Vector4d &cadPlane,
		const std::vector<seg_pair_t> &cloudOuterSegs)
{
	//be careful with cloudPlane and cadPlane, when to use the former, when to use the latter
	prepareData(name, cloud, cloudPlane, cadBorder, cadPlane, cloudOuterSegs);

	//find cloud outerSeg candidate which has overlap with holeSeg
	std::vector<seg_pair_t> candidateSegs;
	getOverlapOfOuterSegAndHoleSeg(candidateSegs);

	//detect lines in cloud wall
	std::vector<Eigen::VectorXf> detectLineCoeffs;
	std::vector<PointCloud::Ptr> detectLinePoints;
	detectLinesInCloudWall(detectLineCoeffs, detectLinePoints);

	//check inside line, delete outerside line
	std::vector<Eigen::VectorXf> leftLineCoeffs;
	std::vector<PointCloud::Ptr> leftLinePoints;
	removeOuterSideLines(detectLineCoeffs, detectLinePoints,
				leftLineCoeffs, leftLinePoints);

	std::vector<Eigen::VectorXf> lineCoeffs = leftLineCoeffs;
	std::vector<PointCloud::Ptr> linePoints = leftLinePoints;
	std::set<int> candidateSegIdxs;
	addCandidateSegs(candidateSegs, lineCoeffs, linePoints, candidateSegIdxs);
	LOG(INFO) << "detectLineCoeffs:" << detectLineCoeffs.size()
		<< " leftLineCoeffs:" << leftLineCoeffs.size()
		<< " candidateSegs:" << candidateSegs.size()
		<< " lineCoeffs:" << lineCoeffs.size()
		<< " cloudOuterSegs:" << cloudOuterSegs_.size();

	std::map<size_t, std::map<size_t, std::vector<size_t>>> mapHole2Lines;
	matchLineToCadHoleSeg(lineCoeffs, linePoints, candidateSegIdxs, mapHole2Lines);

	std::map<size_t, std::map<size_t, std::vector<Eigen::Vector3d>>> mapHole2Nodes;
	getIntersectionNodeOfLines(mapHole2Lines, 
			lineCoeffs, linePoints, mapHole2Nodes);

//debug files		
#ifdef VISUALIZATION_ENABLED
	{	
		std::string file_name = "origLine-" + name_ + ".pcd";
		savePcdFileOfLines(file_name, detectLinePoints);				
	}
	{	
		std::string file_name = "leftLine-" + name_ + ".pcd";
		savePcdFileOfLines(file_name, leftLinePoints);				
	}
	{	
		std::string file_name = "findLine-" + name_ + ".pcd";
		savePcdFileOfLines(file_name, linePoints);				
	}
#endif

	return mapHole2Nodes;
}

void CloudBorder::prepareData(
			const std::string &name, 
			const PointCloud::Ptr cloud,
			const Eigen::Vector4d &cloudPlane,
			const std::vector<std::vector<seg_pair_t>> &cadBorder,
			const Eigen::Vector4d &cadPlane,
			const std::vector<seg_pair_t> &cloudOuterSegs)
{
	name_ = name;
	cloud_ = cloud;
	cloudPlane_ = cloudPlane;
	cadBorder_ = cadBorder;
	cadPlane_ = cadPlane;
	cloudOuterSegs_ = cloudOuterSegs;

	//get cad outerSegs and cad holeSegs
	cadOuterSegs_ = cadBorder.front();
	cadHoleSegs_.clear();
	for (size_t k = 1; k < cadBorder.size(); k++)
	{
		const auto &vecSegs = cadBorder[k];
		cadHoleSegs_.insert(cadHoleSegs_.end(), vecSegs.begin(), vecSegs.end());
	}
	
	//get Neighbour cadSeg pairs
	cadSegNeighbours_.clear();
	for (size_t k = 1; k < cadBorder_.size(); k++)
	{
		std::vector<std::pair<size_t, size_t>> segPairs;
		const auto &vecSegs = cadBorder_[k];
		for (size_t j = 0; j < (vecSegs.size() - 1); j++)
		{
			segPairs.push_back(std::make_pair(j,j+1));
		}
		segPairs.push_back(std::make_pair(vecSegs.size() - 1, 0));
		cadSegNeighbours_[k] = segPairs;
	}

	PointCloud::Ptr cloud_filter(new pcl::PointCloud<pcl::PointXYZ>());
	projectionToPlane(cloudPlane_, cloud_, cloud_filter);

	PointCloud::Ptr cloud_sampling(new pcl::PointCloud<pcl::PointXYZ>);
	uniformSampling(0.01, cloud_filter, cloud_sampling);
	cloudSampling_ = cloud_sampling;
}

double CloudBorder::dist_to_seg(const Eigen::Vector3d &point, const seg_pair_t &seg)
{
	Eigen::Vector3d p = seg.first;
	Eigen::Vector3d n = (seg.second - seg.first).normalized();
	double dist = (point - p).cross(n).norm();	
	return dist;	
}

void CloudBorder::getOverlapOfOuterSegAndHoleSeg(
		std::vector<seg_pair_t> &candidateSegs)
{
	//find cloud outerSeg candidate which has overlap with holeSeg
	for (size_t j = 0; j < cadOuterSegs_.size(); j++)
	{
		auto &seg1 = cadOuterSegs_[j];
		bool find = false;
		for (auto &seg2 : cadHoleSegs_)
		{
			if (dist_to_seg(seg1.first, seg2) < 1e-5 
				&& dist_to_seg(seg1.second, seg2) < 1e-5)
			{
				find = true;
				break;
			}
		}
		if (find) candidateSegs.push_back(cloudOuterSegs_[j]);
	}
}

bool CloudBorder::detectLinesInCloudWall(
			std::vector<Eigen::VectorXf> &detectLineCoeffs,
			std::vector<PointCloud::Ptr> &detectLinePoints)
{
	//detect lines in cloud wall
	pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);
	pcl::Normal ptNormal(cloudPlane_(0), cloudPlane_(1), cloudPlane_(2));
	for (size_t i = 0; i < cloudSampling_->size(); i++)
	{
		normals->push_back(ptNormal);
	}
	std::vector<int> boundIndices;
	searchBoundaries(cloudSampling_, normals, boundIndices);
	auto boundPoints = geo::getSubSet(cloudSampling_, boundIndices, false);

	PointCloud::Ptr inputPoints(new pcl::PointCloud<pcl::PointXYZ>());
	pcl::copyPointCloud(*boundPoints, *inputPoints);

	float disthresh = 0.01;
	float connect_thresh = 0.3;
	while (inputPoints->size() > 10)
	{
		std::vector<int> indices;
		Eigen::VectorXf params;
		std::tie(indices, params) = geo::detectOneLineRansac(inputPoints, disthresh);
		if (indices.size() < 10) break;

		//check direction
		Eigen::Vector3f n = params.block<3, 1>(3, 0).normalized();
		bool isGoodDir = false;
		for (const auto &seg : cadHoleSegs_)
		{
			Eigen::Vector3f segDir 
				= (seg.second.cast<float>() - seg.first.cast<float>()).normalized();
			float dot = std::fabs(n.dot(segDir));
			if (dot > 0.996) 
			{
				isGoodDir = true;
				break;
			}
		}

		//get multi-section of inliers, remove too small section outliers
		auto inliers = geo::getSubSet(inputPoints, indices, false);
		if (isGoodDir)
		{
			PointCloud::Ptr goodInliers(new PointCloud());
			PointCloud::Ptr tmpInput(new PointCloud());
			pcl::copyPointCloud(*inliers, *tmpInput);
			while (tmpInput->size() >= 10)
			{
				std::vector<int> sectionIndex = clusterMainStructure(tmpInput, connect_thresh);
				if (sectionIndex.size() < 10) break;

				auto section = geo::getSubSet(tmpInput, sectionIndex, false);
				goodInliers->insert(goodInliers->end(), section->begin(), section->end());
				
				auto tmpLeft = geo::getSubSet(tmpInput, sectionIndex, true);
				tmpInput->swap(*tmpLeft);
			}
			if (goodInliers->size() >= 10)
			{
				detectLineCoeffs.push_back(params);
				detectLinePoints.push_back(goodInliers);  
			}				
		}
		
		auto leftPoints = geo::getSubSet(inputPoints, indices, true);
		inputPoints->swap(*leftPoints);
	}

#ifdef VISUALIZATION_ENABLED	
	{	
		std::string file_name = "boundPoints-" + name_ + ".pcd";
		savePcdFileOfBoundPoints(file_name, boundPoints);			
	}
#endif
	return true;
}

PointCloud::Ptr CloudBorder::getArea( 
		const seg_pair_t &seg, PointCloud::Ptr inputCloud, 
		const Eigen::Vector4d &cloudPlane, 
		const double calcLength, bool left)
{
	Eigen::Vector3d lineNormal = (seg.second - seg.first).normalized();
	Eigen::Vector3d calcDir = lineNormal.cross(cloudPlane.block<3,1>(0,0)); 
	if (true == left) 
	{
		calcDir = -calcDir;
	}

	Eigen::Vector3d pt2 = seg.second + calcLength * calcDir;
	Eigen::Vector3d pt3 = seg.first + calcLength * calcDir;
	std::vector<Eigen::Vector3d> filerPt;
	filerPt.push_back(seg.first);
	filerPt.push_back(seg.second);
	filerPt.push_back(pt2);
	filerPt.push_back(pt3);
	auto  pCloud = filerCloudByConvexHull(inputCloud, filerPt); 
	return pCloud;
}

int CloudBorder::judge_flag(const PointCloud::Ptr pCloudLeft, 
							const PointCloud::Ptr pCloudRight)
{
	int flag = 0;
	int sum = pCloudLeft->size() + pCloudRight->size();
	if (0 == sum) return flag;

	float rateLeft = float(pCloudLeft->size()) / float(sum);
	float rateRight = float(pCloudRight->size()) / float(sum);
	if (rateLeft > 0.8)
	{
		flag = -1;
	}
	else if (rateRight > 0.8)
	{
		flag = 1;
	}
	else
	{
		flag = 0;
	}
	return flag;	
}

int CloudBorder::judgeCloudInWhichSideOfLine(
			const std::string &fileName,
			const Eigen::VectorXf &line,
			const PointCloud::Ptr inliers)
{
	auto project_to_line = [](const Eigen::VectorXf &line, 
			const PointCloud::Ptr inliers)->std::map<float, size_t> {
		Eigen::Vector3f p = line.block<3,1>(0,0);
		Eigen::Vector3f n = line.block<3,1>(3,0).normalized();

		std::map<float, size_t> projectValue2Idx;
		for (size_t k = 0; k < inliers->size(); k++)
		{
			const auto &pt = inliers->points[k];
			Eigen::Vector3f point(pt.x, pt.y, pt.z);
			float projectValue = (point - p).dot(n);
			projectValue2Idx[projectValue] = k;
		}
		return projectValue2Idx;
	};

	std::map<float, size_t> projectValue2Idx = project_to_line(line, inliers);
	pcl::PointXYZ s = inliers->points[projectValue2Idx.begin()->second];
	pcl::PointXYZ e = inliers->points[projectValue2Idx.rbegin()->second];

	const double calcLength = 0.1;
	seg_pair_t tmpSeg(Eigen::Vector3d(s.x, s.y, s.z), Eigen::Vector3d(e.x, e.y, e.z));
	auto pCloudLeft = getArea(tmpSeg, cloudSampling_, cloudPlane_, calcLength, true);
	auto pCloudRight = getArea(tmpSeg, cloudSampling_, cloudPlane_, calcLength, false);
	int flag = judge_flag(pCloudLeft, pCloudRight);

	// Eigen::vector<Eigen::Vector3d> vec_tmp = ininterpolateSeg(Eigen::Vector3d(s.x, s.y, s.z), 
	// 				Eigen::Vector3d(e.x, e.y, e.z), 0.05);
	// int flag = 0;
	// for (size_t k = 0; k < vec_tmp.size() - 1; k++)
	// {
	// 	const auto &pt0 = vec_tmp[k];
	// 	const auto &pt1 = vec_tmp[k+1];
	// 	seg_pair_t tmpSeg(pt0, pt1);
	// 	auto pCloudLeft = getArea(tmpSeg, cloudSampling_, cloudPlane_, calcLength, true);
	// 	auto pCloudRight = getArea(tmpSeg, cloudSampling_, cloudPlane_, calcLength, false);
	// 	flag = judge_flag(pCloudLeft, pCloudRight);
	// }

#ifdef VISUALIZATION_ENABLED
	{			
		std::string file_name = "area-" + name_ + "-" + fileName 
				+ ".pcd";
		savePcdFileOfBothSideAreaOfLines(file_name, pCloudLeft, pCloudRight, tmpSeg);			
	}	
#endif
	return flag;
}

std::pair<bool, int> CloudBorder::inSameSideOfLine(
		const Eigen::Vector4d &cloudPlane, const Eigen::VectorXf &line, 
		const PointCloud::Ptr inliers, const float rate_Th)
{
	Eigen::Vector3f p = line.block<3, 1>(0, 0);
	Eigen::Vector3f n = line.block<3, 1>(3, 0).normalized();
	Eigen::Vector3f planeNormal(cloudPlane(0), cloudPlane(1), cloudPlane(2));
	planeNormal /= planeNormal.norm();
	
	int leftSum = 0;
	int rightSum = 0;					
	for (const auto &pt : inliers->points)
	{
		Eigen::Vector3f point(pt.x, pt.y, pt.z);
		Eigen::Vector3f cross = (point - p).cross(n);
		float dist = cross.norm();
		if (dist < 0.05) continue;

		float flag = cross.dot(planeNormal);
		int &side = (flag < 0) ? leftSum : rightSum;	
		side++;
	}
	int total = leftSum + rightSum;
	float leftRate = (float)leftSum / (float)total;
	float rightRate = (float)rightSum / (float)total;
	if (leftRate > rate_Th)
	{
		return std::make_pair(true, -1);
	}
	else if (rightRate > rate_Th)
	{
		return std::make_pair(true, 1);
	}
	else
	{
		return std::make_pair(false, 0);
	}
}

bool CloudBorder::tooCloseToLine(const Eigen::VectorXf &line, 
	const PointCloud::Ptr inliers, const float rate_Th) 
{
	Eigen::Vector3f p = line.block<3, 1>(0, 0);
	Eigen::Vector3f n = line.block<3, 1>(3, 0).normalized();	
	int sum = 0;
	for (const auto &pt : inliers->points)
	{
		Eigen::Vector3f point(pt.x, pt.y, pt.z);
		Eigen::Vector3f cross = (point - p).cross(n);
		float dist = cross.norm();
		if (dist < 0.05) 
		{
			sum++;
		}
	}	
	float rate = (float)sum / (float)inliers->size();
	if (rate > rate_Th)
	{
		return true;
	}
	return false;
}

bool CloudBorder::removeOuterSideLines(
		const std::vector<Eigen::VectorXf> &detectLineCoeffs,
		const std::vector<PointCloud::Ptr> &detectLinePoints,
		std::vector<Eigen::VectorXf> &leftLineCoeffs,
		std::vector<PointCloud::Ptr> &leftLinePoints)
{
	//check inside line, delete outerside line
	auto dist_to_line = [](const Eigen::Vector3f &point, const Eigen::VectorXf &line)->double {
		Eigen::Vector3f p = line.block<3, 1>(0, 0);
		Eigen::Vector3f n = line.block<3, 1>(3, 0).normalized();
		float dist = (point - p).cross(n).norm();	
		return dist;	
	};

	const float rate_Th = 0.95;
	const float distTh = 0.4;
	const float dotTh = 0.996;
	for (size_t i = 0; i < detectLineCoeffs.size(); i++)
	{
		Eigen::VectorXf line = detectLineCoeffs[i];
		bool isInsideLine = false;
		int lastSide = 0;
		for (size_t j = 0; j < detectLinePoints.size(); j++)
		{
			if (i == j) continue;
			PointCloud::Ptr inliers = detectLinePoints[j];
			if(true == tooCloseToLine(line, inliers, rate_Th)) continue;

			auto ret = inSameSideOfLine(cloudPlane_, line, inliers, rate_Th);
			if (false == ret.first) 
			{
				isInsideLine = true;
				break;
			}
			else if (0 != lastSide && lastSide != ret.second)
			{
				isInsideLine = true;
				break;
			}
			
			lastSide = ret.second;
		}

		if (false == isInsideLine)
		{
			isInsideLine = true;
			Eigen::Vector3f n = line.block<3, 1>(3, 0).normalized();
			float minDist = float(RAND_MAX);
			for (const auto &seg : cloudOuterSegs_)
			{
				Eigen::Vector3f segDir = (seg.second.cast<float>() - seg.first.cast<float>()).normalized();
				float dot = std::fabs(n.dot(segDir));
				if (dot < dotTh) continue;

				float dist1 = dist_to_line(seg.first.cast<float>(), line);
				float dist2 = dist_to_line(seg.second.cast<float>(), line);
				minDist = std::min(minDist, (dist1+dist2)/2.0f);
			}
			if (minDist < distTh)
			{
				isInsideLine = false;			
			}				
		}

		if (true == isInsideLine)
		{
			leftLineCoeffs.push_back(line);
			leftLinePoints.push_back(detectLinePoints[i]);
		}
	}

	return true;	
}

void CloudBorder::addCandidateSegs(const std::vector<seg_pair_t> &candidateSegs,
		std::vector<Eigen::VectorXf> &lineCoeffs,
		std::vector<PointCloud::Ptr> &linePoints,
		std::set<int> &candidateSegIdxs)
{
	for (const auto &seg : candidateSegs)
	{
		Eigen::VectorXf line(6);
		Eigen::Vector3d n = (seg.second - seg.first).normalized();
		line << seg.first(0), seg.first(1), seg.first(2), n(0), n(1), n(2);
		lineCoeffs.push_back(line);

		PointCloud::Ptr points(new pcl::PointCloud<pcl::PointXYZ>());
		auto vec_tmp = ininterpolateSeg(seg.first, seg.second, 0.001);
		for (const auto &p : vec_tmp)
		{
			points->push_back(pcl::PointXYZ(p(0), p(1), p(2)));
		}
		linePoints.push_back(points);
		candidateSegIdxs.insert(linePoints.size() - 1);
	}
} 

PointCloud::Ptr CloudBorder::pointsOverlapWithSeg(const PointCloud::Ptr inputCloud, 
							const seg_pair_t &seg)
{
	PointCloud::Ptr overlapPoints(new PointCloud());

	Eigen::Vector3d seg_vec = seg.second - seg.first;
	for (const auto &p : inputCloud->points)
	{
		Eigen::Vector3d point(p.x, p.y, p.z);
		double a = (point - seg.first).dot(seg_vec) 
						/ seg_vec.squaredNorm();
		if (a > 0.0 && a < 1.0)
		{
			overlapPoints->push_back(p);
		}
		else
		{
			double dist1 = (point - seg.first).norm();
			double dist2 = (point - seg.second).norm();
			if (dist1 < 0.4 || dist2 < 0.4)
			{
				overlapPoints->push_back(p);
			}
		}
	}
	return overlapPoints;
}

double CloudBorder::aveDistToSeg(const PointCloud::Ptr inputCloud, const seg_pair_t &seg) 
{
	double aveDist = 0.0;
	for (const auto &p : inputCloud->points)
	{
		Eigen::Vector3d point(p.x, p.y, p.z);
		double dist = dist_to_seg(point, seg);
		aveDist += dist;
	}
	aveDist /= double(inputCloud->size());
	return aveDist;
};

int CloudBorder::findNextSeg(const std::vector<seg_pair_t> &vecSegs,
						const std::vector<bool> &vecValid,
						const seg_pair_t &currtSeg,  seg_pair_t &nextSeg)
{
	for (size_t i = 0; i < vecSegs.size(); i++)
	{
		if (false == vecValid[i]) continue;

		const auto &seg = vecSegs[i];
		if ((currtSeg.second - seg.first).norm() < 1e-5)
		{
			nextSeg = seg;
			return i;
		}
		if ((currtSeg.second - seg.second).norm() < 1e-5)
		{
			nextSeg = std::make_pair(seg.second, seg.first);
			return i;
		}
	}

	return -1;
}

bool CloudBorder::anticlockwiseSortSegVec(
		const Eigen::Vector3d &planeNormal, const std::vector<seg_pair_t> &vecSegs,
		std::vector<seg_pair_t> &sortSegs, std::vector<int> &matchIdxs)
{
	if (vecSegs.empty()) return false;

	std::vector<seg_pair_t> connectSegs;
	std::vector<int> tmpMatchIdxs;
	std::vector<bool> vecValid(vecSegs.size(), true);

	connectSegs.push_back(vecSegs.front());
	tmpMatchIdxs.push_back(0);
	vecValid[0] =  false;
	while (connectSegs.size() < vecSegs.size())
	{
		const auto &currSeg = connectSegs.back();

		seg_pair_t nextSeg;
		int idx = findNextSeg(vecSegs, vecValid, currSeg, nextSeg);
		if (-1 == idx)
		{
			break;
		}
		connectSegs.push_back(nextSeg);
		tmpMatchIdxs.push_back(idx);
		vecValid[idx] = false;
	}

	const auto &firstSeg = connectSegs.front();
	Eigen::Vector3d calcDir = (firstSeg.second - firstSeg.first).cross(planeNormal); 
	Eigen::Vector3d pt = (firstSeg.first + firstSeg.second) / 2.0 + 0.02 * calcDir.normalized();
	Eigen::vector<seg_pair_t> segments(connectSegs.begin(), connectSegs.end());
	if (false == isPointInPolygon3D(pt, segments))
	{
		sortSegs = connectSegs;
		matchIdxs = tmpMatchIdxs;
	}
	else
	{
		sortSegs.push_back(std::make_pair(firstSeg.second, firstSeg.first));
		matchIdxs.clear();
		matchIdxs.push_back(tmpMatchIdxs.front());
		for (int j = connectSegs.size() - 1; j > 0; j--)
		{
			const auto &seg = connectSegs[j];
			sortSegs.push_back(std::make_pair(seg.second, seg.first));
			matchIdxs.push_back(tmpMatchIdxs[j]);
		}
		LOG(INFO) << "reverse segVec";
	}
	
	if (vecSegs.size() != sortSegs.size() || sortSegs.size() != matchIdxs.size())
	{
		LOG(WARNING) << "anticlockwiseSortSegVec Failed, vecSegs.size:" << vecSegs.size()
			<<" sortSegs.size:" << sortSegs.size()
			<< " matchIdxs.size:" << matchIdxs.size();
		return false;
	}
	else
	{
		return true;
	}
}

bool CloudBorder::matchLineToCadHoleSeg(
	const std::vector<Eigen::VectorXf> &lineCoeffs,
	const std::vector<PointCloud::Ptr> &linePoints,
	const std::set<int> &candidateSegIdxs,
	std::map<size_t, std::map<size_t, std::vector<size_t>>> &mapHole2Lines)
{
	for (size_t k = 1; k < cadBorder_.size(); k++)
	{
		std::map<size_t, std::vector<size_t>> mapSeg2Lines;
		const auto &origVecSegs = cadBorder_[k];
		std::vector<seg_pair_t> sortSegs; 
		std::vector<int> matchIdxs;
		if (false == anticlockwiseSortSegVec(cadPlane_.block<3,1>(0,0), 
											origVecSegs, sortSegs, matchIdxs))
		{
			continue;
		}

		for (size_t j = 0; j < sortSegs.size(); j++)
		{
			size_t origIdx = matchIdxs[j];
			const auto &seg = sortSegs[j];
			Eigen::Vector3f seg_vec = seg.second.cast<float>() - seg.first.cast<float>();
			Eigen::Vector3f segDir = seg_vec.normalized();
			//NOte: cloud is allways in right side of segDir,
			//rule: dir of segDir.cross(cadPlaneNormal) is right(flag = 1), the reverse dir is left(flag = -1)
			int segFlag = 1;

			std::vector<size_t> lineIdxs;
			for (size_t i = 0; i < lineCoeffs.size(); i++)
			{
				Eigen::VectorXf line = lineCoeffs[i];
				auto inliers = linePoints[i];
				Eigen::Vector3f n = line.block<3, 1>(3, 0).normalized();
				float dot = n.dot(segDir);
				if (std::fabs(dot) < 0.996) continue;
				
				double aveDist = aveDistToSeg(inliers, seg);
				if (aveDist > 1.0) continue;

				if (candidateSegIdxs.count(i))
				{
					lineIdxs.push_back(i);
					continue;
				}

				auto belongPoints = pointsOverlapWithSeg(inliers, seg);
				if (belongPoints->size() < 5) continue;	
				
				//correct line dir to be same with seg dir
				Eigen::VectorXf correctLine = line;
				if (dot < 0) correctLine.block<3,1>(3,0) = -line.block<3,1>(3,0);

				std::string fileName = "hole" + std::to_string(k) 
					+ "-sort" + std::to_string(j)
					+ "-seg" + std::to_string(origIdx) + "-line" + std::to_string(i);
				int line_flag = judgeCloudInWhichSideOfLine(fileName, correctLine, belongPoints);
				LOG(INFO) << fileName << " segFlag:" << segFlag
					<< " line_flag:" << line_flag << " " << (line_flag != segFlag ? "diff" : "same");
				if (0 != line_flag && (line_flag != segFlag)) continue;

				lineIdxs.push_back(i);
			}
			LOG(INFO) << "-----hole:" << k << " sort:" << j << " seg:" << origIdx 
						<< " lineIdxs size:" << lineIdxs.size();
			mapSeg2Lines[origIdx] = lineIdxs;

#ifdef VISUALIZATION_ENABLED
			{			
				std::string file_name = "belongLines-" + name_ + "-hole" + std::to_string(k) 
						+ "-sort" + std::to_string(j)
						+ "-seg" + std::to_string(origIdx) + ".pcd";
				savePcdFileOfBelongLines(file_name, linePoints, lineIdxs, seg);			
			}				
		}
#endif
		mapHole2Lines[k] = mapSeg2Lines;
	}

	return true;
}

void CloudBorder::getIntersectionNodeOfLines(
	const std::map<size_t, std::map<size_t, std::vector<size_t>>> &mapHole2Lines,
	const std::vector<Eigen::VectorXf> &lineCoeffs,
	const std::vector<PointCloud::Ptr> &linePoints,
	std::map<size_t, std::map<size_t, std::vector<Eigen::Vector3d>>> &mapHole2Nodes)
{
	//calc intersection node between two lines
	for (auto &it : cadSegNeighbours_)
	{
		size_t holeIdx = it.first;
		const auto findItor = mapHole2Lines.find(holeIdx);
		if (findItor == mapHole2Lines.end()) continue;
		std::map<size_t, std::vector<size_t>> mapSeg2Lines = findItor->second;

		std::map<size_t, std::vector<Eigen::Vector3d>> mapSeg2Nodes;
		for (auto &idxPair : it.second)
		{
			auto &segIdx1 = idxPair.first;
			auto &segIdx2 = idxPair.second;
			if (!mapSeg2Lines.count(segIdx1) || !mapSeg2Lines.count(segIdx2)) continue;

			std::vector<size_t>	lineIdxs1 = mapSeg2Lines[segIdx1];
			std::vector<size_t>	lineIdxs2 = mapSeg2Lines[segIdx2];

			std::vector<Eigen::Vector3d> vecNodes;
			for (size_t i = 0; i < lineIdxs1.size(); i++)
			{
				auto idx1 = lineIdxs1[i];
				auto tmp1 = lineCoeffs[idx1];
				Eigen::VectorXd line1(6);
				line1 << tmp1(0), tmp1(1), tmp1(2), tmp1(3), tmp1(4), tmp1(5);
				for (size_t j = 0; j < lineIdxs2.size(); j++)
				{
					auto idx2 = lineIdxs2[j];
					auto tmp2 = lineCoeffs[idx2];
					Eigen::VectorXd line2(6);
					line2 << tmp2(0), tmp2(1), tmp2(2), tmp2(3), tmp2(4), tmp2(5);
					Eigen::Vector3d interSectionPt;
					if (false == interSectionOfLineToLine(line1, line2, interSectionPt)) continue;
					auto vecPts1 = convertCloudToEigenVec(linePoints[idx1]);
					auto vecPts2 = convertCloudToEigenVec(linePoints[idx2]);

					auto ret1 = findNearestPt(vecPts1, interSectionPt);
					auto ret2 = findNearestPt(vecPts2, interSectionPt);
					if (ret1.first > 1.0 || ret2.first > 1.0) continue;

					vecNodes.push_back(interSectionPt);
				}
			}
			LOG(INFO) << "holeIdx:" << holeIdx 
				<< " segIdx1:" << segIdx1 << " segIdx2:" << segIdx2 
				<< " vecNodes:" << vecNodes.size();
			if (!vecNodes.empty()) 
			{
				if (!mapSeg2Nodes.count(segIdx1))
				{
					mapSeg2Nodes[segIdx1] = vecNodes;
				}
				if (!mapSeg2Nodes.count(segIdx2))
				{
					mapSeg2Nodes[segIdx2] = vecNodes;
				}
				auto &tmpVec1 = mapSeg2Nodes[segIdx1];
				tmpVec1.insert(tmpVec1.end(), vecNodes.begin(), vecNodes.end());

				auto &tmpVec2 = mapSeg2Nodes[segIdx2];
				tmpVec2.insert(tmpVec2.end(), vecNodes.begin(), vecNodes.end());					
			}			
		}
		if (!mapSeg2Nodes.empty()) mapHole2Nodes[holeIdx] = mapSeg2Nodes;
	}
}

void CloudBorder::savePcdFileOfLines(const std::string &file_name,
			const std::vector<PointCloud::Ptr> &linePoints)
{	
	std::default_random_engine e;
	std::uniform_real_distribution<double> random(0,1);
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr pCloud(new pcl::PointCloud<pcl::PointXYZRGB>());
	for (size_t i = 0; i < linePoints.size(); ++i)
	{
		auto cloud = linePoints[i];
		int r = int(random(e)*255);
		int g = int(random(e)*255);
		int b = int(random(e)*255);
		for (auto &p1 : cloud->points)
		{
			pcl::PointXYZRGB p2;
			p2.x = p1.x;
			p2.y = p1.y;
			p2.z = p1.z;
			p2.r = r;
			p2.g = g;
			p2.b = b;
			pCloud->push_back(p2);
		}	
	}
	pcl::io::savePCDFile(file_name, *pCloud);				
}

void CloudBorder::savePcdFileOfBoundPoints(const std::string &file_name,
			const PointCloud::Ptr boundPoints)
{	
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr pCloud(new pcl::PointCloud<pcl::PointXYZRGB>());
	for (auto &p1 : boundPoints->points)
	{
		pcl::PointXYZRGB p2;
		p2.x = p1.x;
		p2.y = p1.y;
		p2.z = p1.z;
		p2.r = 255;
		p2.g = 255;
		p2.b = 255;
		pCloud->push_back(p2);
	}	
	
	pcl::io::savePCDFile(file_name, *pCloud);				
}

void CloudBorder::savePcdFileOfBelongLines(const std::string &file_name,
		const std::vector<PointCloud::Ptr> &linePoints,
		const std::vector<size_t> &lineIdxs, const seg_pair_t &seg)
{	
	std::default_random_engine e;
	std::uniform_real_distribution<double> random(0,1);
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr pCloud(new pcl::PointCloud<pcl::PointXYZRGB>());
	for (auto i : lineIdxs)
	{
		auto cloud = linePoints[i];
		int r = int(random(e)*255);
		int g = int(random(e)*255);
		int b = int(random(e)*255);
		for (auto &p1 : cloud->points)
		{
			pcl::PointXYZRGB p2;
			p2.x = p1.x;
			p2.y = p1.y;
			p2.z = p1.z;
			p2.r = r;
			p2.g = g;
			p2.b = b;
			pCloud->push_back(p2);
		}	
	}
	//cad border
	Eigen::vector<Eigen::Vector3d> vecPoints 
		= ininterpolateSeg(seg.first, seg.second, 0.05);
	for (const auto &point : vecPoints)
	{
			pcl::PointXYZRGB p2;
			p2.x = point(0);
			p2.y = point(1);
			p2.z = point(2);
			p2.r = 255;
			p2.g = 255;
			p2.b = 255;
			pCloud->push_back(p2);
	}			
	pcl::io::savePCDFile(file_name, *pCloud);				
}

void CloudBorder::savePcdFileOfBothSideAreaOfLines(const std::string &file_name,
		const PointCloud::Ptr &pCloudLeft, const PointCloud::Ptr &pCloudRight,
		const seg_pair_t &seg)
{	
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr pCloud(new pcl::PointCloud<pcl::PointXYZRGB>());
	for (auto &p1 : pCloudLeft->points)
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
	for (auto &p1 : pCloudRight->points)
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

	Eigen::vector<Eigen::Vector3d> vecPoints 
		= ininterpolateSeg(seg.first, seg.second, 0.05);
	for (const auto &point : vecPoints)
	{
			pcl::PointXYZRGB p2;
			p2.x = point(0);
			p2.y = point(1);
			p2.z = point(2);
			p2.r = 255;
			p2.g = 255;
			p2.b = 255;
			pCloud->push_back(p2);
	}			
	pcl::io::savePCDFile(file_name, *pCloud);				
}

}

