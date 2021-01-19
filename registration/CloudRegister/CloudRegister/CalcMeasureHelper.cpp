#include "CalcMeasureHelper.h"

#include "funHelper.h"
#include <pcl/common/common.h>
#include <pcl/filters/crop_box.h>
#include <pcl/filters/passthrough.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/filters/crop_hull.h>
#include <pcl/surface/concave_hull.h>

namespace CloudReg
{

	void writePCDFile(const std::string& name, const std::vector<seg_pair_t>& segA, const std::vector<seg_pair_t>& segB)
	{
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr pCloudRGB(new pcl::PointCloud<pcl::PointXYZRGB>());
		for(auto& seg : segA)
		{
			auto vecPts = ininterpolateSeg(seg.first,seg.second,0.01);
			for(auto& pt : vecPts)
			{
				pcl::PointXYZRGB p2;
				p2.x = pt[0];	
				p2.y = pt[1];	
				p2.z = pt[2];	
				p2.r = 0;
				p2.g = 255;
				p2.b = 0;

				pCloudRGB->push_back(p2);
			}
		}
		
		for(auto& seg : segB)
		{
			auto vecPts = ininterpolateSeg(seg.first,seg.second,0.01);
			for(auto& pt : vecPts)
			{
				pcl::PointXYZRGB p2;
				p2.x = pt[0];	
				p2.y = pt[1];	
				p2.z = pt[2];	
				p2.r = 255;
				p2.g = 0;
				p2.b = 0;

				pCloudRGB->push_back(p2);
			}
		}
		pcl::io::savePCDFile(name, *pCloudRGB);
	}

	void writePCDFile(const std::string& name, const PointCloud::Ptr pCloud, std::vector<seg_pair_t>& border)
	{
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr pCloudRGB(new pcl::PointCloud<pcl::PointXYZRGB>());

		for(auto& pt : pCloud->points)
		{
			pcl::PointXYZRGB p2;
			p2.x = pt.x;	
			p2.y = pt.y;	
			p2.z = pt.z;	
			p2.r = 0;
			p2.g = 255;
			p2.b = 0;

			pCloudRGB->push_back(p2);
		}

		for(auto& seg : border)
		{
			auto vecPts = ininterpolateSeg(seg.first,seg.second,0.01);
			for(auto& pt : vecPts)
			{
				pcl::PointXYZRGB p2;
				p2.x = pt[0];	
				p2.y = pt[1];	
				p2.z = pt[2];	
				p2.r = 255;
				p2.g = 0;
				p2.b = 0;

				pCloudRGB->push_back(p2);
			}
		}
		pcl::io::savePCDFile(name, *pCloudRGB);

	}

	Eigen::Vector3d calcPerpendicular(const Eigen::Vector3d &pt,
			const Eigen::Vector3d &start, const Eigen::Vector3d &end)
	{
		Eigen::Vector3d retVal(0,0,0);
		Eigen::Vector3d segAB = start - end; 
		Eigen::Vector3d segPA = pt - start; 
		
		if(segAB.norm() < EPS_FLOAT_DOUBLE)
		{
			retVal = start;
			return retVal;
		}

		double u = segAB.dot(segPA)/segAB.squaredNorm();
		retVal = start + u * segAB;
		return retVal;
	}

	void groupDirection(const Eigen::Vector3d& horizenSeg, const std::vector<seg_pair_t>& border, 
			std::vector<seg_pair_t>& vecVertical, std::vector<seg_pair_t>& vecHorizen)
	{
		for(auto& seg : border)
		{
			auto curSeg = seg.first - seg.second;
			if (curSeg.norm() < EPS_FLOAT_DOUBLE) continue;

			double cos = horizenSeg.dot(curSeg)/(horizenSeg.norm() * curSeg.norm());
			// std::cout << "cos " <<fabs(cos) << std::endl;
			if (fabs(cos) < 0.1) 
				vecVertical.emplace_back(seg);
			else 
				vecHorizen.emplace_back(seg);
		}

		std::sort(vecHorizen.begin(), vecHorizen.end(), [&](const seg_pair_t& left, const seg_pair_t& right){
				return left.first[2] < right.first[2];});
	}

	PointCloud::Ptr filerCloudByConvexHull(pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud, 
											const std::vector<Eigen::Vector3d>& corners, const bool negative)
	{
        pcl::PointCloud<pcl::PointXYZ>::Ptr boundingbox_ptr (new pcl::PointCloud<pcl::PointXYZ>);
		for(const auto &corner : corners)
		{
			boundingbox_ptr->push_back(pcl::PointXYZ(corner[0], corner[1], corner[2]));
		}
        
        pcl::ConvexHull<pcl::PointXYZ> hull;                  
        hull.setInputCloud(boundingbox_ptr);                 
        hull.setDimension(2);                                 
        std::vector<pcl::Vertices> polygons;                 
        pcl::PointCloud<pcl::PointXYZ>::Ptr surface_hull (new pcl::PointCloud<pcl::PointXYZ>);
        hull.reconstruct(*surface_hull, polygons);           

        pcl::PointCloud<pcl::PointXYZ>::Ptr objects (new pcl::PointCloud<pcl::PointXYZ>);
        pcl::CropHull<pcl::PointXYZ> bb_filter;               
        bb_filter.setDim(2);
		bb_filter.setNegative(negative);
        bb_filter.setInputCloud(pCloud);                       
        bb_filter.setHullIndices(polygons);                   
        bb_filter.setHullCloud(surface_hull);                 
        bb_filter.filter(*objects);                           

		// LOG(INFO) << "inPut: " << pCloud->points.size() << " outPut: " << objects->points.size();
		return objects;
	}

	PointCloud::Ptr filerCloudByRange(pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud,
			const pcl::PointXYZ& min, const pcl::PointXYZ& max)
	{
		pcl::CropBox<pcl::PointXYZ> clipper;
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_final(new pcl::PointCloud<pcl::PointXYZ>());
		
		clipper.setMin(Eigen::Vector4f(min.x,min.y,min.z,1.0));
		clipper.setMax(Eigen::Vector4f(max.x,max.y,max.z,1.0));
        clipper.setNegative(false);
        clipper.setInputCloud(pCloud);
        clipper.filter(*cloud_final);
		LOG(INFO) << "inPut: " << pCloud->points.size() << " outPut: " << cloud_final->points.size();
		return cloud_final;
	}


	Eigen::vector<Eigen::Vector3d> getNearestPt(const Eigen::vector<Eigen::Vector3d>& vecCalcPt, 
			const PointCloud::Ptr pCloud, const double maxDist)
	{
		Eigen::vector<Eigen::Vector3d> vecCalcRaw(vecCalcPt.size(),Eigen::Vector3d(0,0,0));
		std::vector<double> vecDist(vecCalcPt.size(), 9999.0);
		for(auto& pt : pCloud->points)
		{
			Eigen::Vector3d ePt(pt.x, pt.y, pt.z);
			for(std::size_t i = 0; i < vecCalcPt.size(); i++)
			{
				double tmp = (ePt - vecCalcPt[i]).squaredNorm();
				if (tmp < maxDist && tmp < vecDist[i])
				{
					vecDist[i] = tmp;
					vecCalcRaw[i] = ePt;
				}
			}
		}

		for(std::size_t i = 0; i < vecCalcPt.size(); i++)
		{
			if (vecCalcRaw[i].norm() < EPS_FLOAT_DOUBLE)
			{
				vecCalcRaw[i] = vecCalcPt[i];
				LOG(WARNING) << "can't getNearestPt, use border instead";
			}
		}

		return vecCalcRaw;
	}

	void groupDirectionIndex(const Eigen::Vector3d& horizenSeg, const std::vector<seg_pair_t>& border, 
			std::vector<std::size_t>& vecVertical, std::vector<std::size_t>& vecHorizen)
	{
		for(std::size_t i = 0 ; i < border.size(); i++)
		{
			auto& seg = border[i];
			auto curSeg = seg.first - seg.second;
			double cos = horizenSeg.dot(curSeg)/(horizenSeg.norm() * curSeg.norm());
			//std::cout << "cos: " << std::setprecision(8) << fabs(cos) << std::endl;

			if (fabs(cos) < 0.1) 
				vecVertical.emplace_back(i);
			else 
				vecHorizen.emplace_back(i);
		}

	}

	bool isRootInSeg(const seg_pair_t& seg, const Eigen::Vector3d& p)
	{

		auto segAB = (seg.first - seg.second).normalized();
		auto segAP = (seg.first - p).normalized();
		double dotA = segAB.dot(segAP);

		auto segBA = seg.second - seg.first;
		auto segBP = seg.second - p;
		double dotB = segBA.dot(segBP);
		//LOG(INFO) << (std::setprecision(8)) << dotA << " " << dotB;
		if (dotA < -0.1 || dotB < -0.1)
			return false;
		
		return true;
	}

	std::tuple<std::size_t, std::size_t, int> getWallGrowAxisAndDir(const Eigen::Vector3d& sPt, const Eigen::Vector3d& ePt)
	{
		std::size_t index = std::fabs(sPt[0] - ePt[0]) > 0.1 ? 0 : 1;
		std::size_t indexOther = index == 0 ? 1 : 0;
		int dir = sPt[index] < ePt[index] ? 1 : -1;
		return std::make_tuple(index, indexOther, dir);

	}

	
	std::tuple<bool, Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d>
		calcOverlap(seg_pair_t& toSeg, seg_pair_t& calcSeg)
	{
		if ((toSeg.first - toSeg.second).dot(calcSeg.first - calcSeg.second) < 0)
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

		Eigen::Vector3d s1Pt(.0,.0,.0), e1Pt(.0,.0,.0);
		Eigen::Vector3d s2Pt(.0,.0,.0), e2Pt(.0,.0,.0);
		bool ret = false;
		if (findE1 && findS1)
		{
			ret = true;
			s1Pt = calcPerpendicular(calcSeg.first, toSeg.first, toSeg.second);
			e1Pt = calcPerpendicular(calcSeg.second, toSeg.first, toSeg.second);

			s2Pt = calcSeg.first;
			e2Pt = calcSeg.second;


		}
		else if (findE2 && findS2)
		{
			ret = true;
			s1Pt = toSeg.first;
			e1Pt = toSeg.second;

			s2Pt = calcPerpendicular(toSeg.first, calcSeg.first, calcSeg.second);
			e2Pt = calcPerpendicular(toSeg.second, calcSeg.first, calcSeg.second);

		}
		else if (findS1 && findE2)// 1 0 0 1
		{
			ret = true;
			s1Pt = calcPerpendicular(calcSeg.first, toSeg.first, toSeg.second);
			e1Pt = toSeg.second;

			s2Pt = calcSeg.first;
			e2Pt = calcPerpendicular(toSeg.second, calcSeg.first, calcSeg.second);

		}
		else if (findS2 && findE1)//0 1 1 0
		{
			ret = true;
			s1Pt = toSeg.first;
			e1Pt = calcPerpendicular(calcSeg.second, toSeg.first, toSeg.second);


			s2Pt = calcPerpendicular(toSeg.first, calcSeg.first, calcSeg.second);
			e2Pt = calcSeg.second;

		}

		return std::make_tuple(ret, s1Pt, e1Pt, s2Pt, e2Pt);
	}

	std::vector<seg_pair_t> calcBoxSegPair(std::vector<Eigen::Vector3d>& vecPt)
	{
		std::vector<seg_pair_t> vecSeg;
		for (std::size_t i = 1; i < 4; i++)
		{
			vecSeg.emplace_back(std::make_pair(vecPt[i - 1], vecPt[i]));
			if (i == 3)
				vecSeg.emplace_back(std::make_pair(vecPt[i], vecPt[0]));
		}

		for (std::size_t i = 5; i < 8; i++)
		{
			vecSeg.emplace_back(std::make_pair(vecPt[i - 1], vecPt[i]));
			if (i == 7)
				vecSeg.emplace_back(std::make_pair(vecPt[i], vecPt[4]));
		}
		for (std::size_t i = 0; i < 4; i++)
		{
			vecSeg.emplace_back(std::make_pair(vecPt[i], vecPt[i + 4]));
		}

		return vecSeg;
	}

	std::vector<Eigen::Vector3d> createRulerBox(seg_pair_t ruler, int thicknessDir, double thickness, double width)
	{
		if ((ruler.first - ruler.second).norm() < 1e-8)
		{
			LOG(ERROR) << "ruler length is 0";
			std::vector<Eigen::Vector3d> rulerPoints;
			return rulerPoints;
		}

		//thickness Direction
		Eigen::Vector3d thickn(0,0,0);
		thickn[thicknessDir] = 1;
		//ruler direction
		Eigen::Vector3d rulerab = ruler.second - ruler.first;
		Eigen::Vector3d rulern = rulerab.normalized();
		//width direction
		Eigen::Vector3d widthn = thickn.cross(rulern);

		//first 4 points
		Eigen::Vector3d pta = ruler.first + widthn * (width / 2);
		Eigen::Vector3d ptb = ruler.first - widthn * (width / 2);

		Eigen::Vector3d pt1 = pta + thickn * (thickness / 2);
		Eigen::Vector3d pt2 = pta - thickn * (thickness / 2);
		Eigen::Vector3d pt3 = ptb + thickn * (thickness / 2);
		Eigen::Vector3d pt4 = ptb - thickn * (thickness / 2);

		//last 4 points
		Eigen::Vector3d pt5 = pt1 + rulerab;
		Eigen::Vector3d pt6 = pt2 + rulerab;
		Eigen::Vector3d pt7 = pt3 + rulerab;
		Eigen::Vector3d pt8 = pt4 + rulerab;

		std::vector<Eigen::Vector3d> rulerPoints = {pt1, pt2, pt3, pt4, pt5, pt6, pt7, pt8};
		return rulerPoints;
	}

	std::vector<Eigen::Vector3d> getRulerCorners(const std::vector<Eigen::Vector3d>& rPoints)
	{
		std::vector<Eigen::Vector3d> corners;
		if (rPoints.empty())
		{
			LOG(ERROR) << "empty ruler points";
			return corners;
		}
		corners = {rPoints[0], rPoints[2], rPoints[4], rPoints[6]};
		return corners;
	}

	bool calIntersection(seg_pair_t line1, seg_pair_t line2, Eigen::Vector3d& intersec)
	{
		Eigen::Vector3d s1e1 = line1.second - line1.first;
		Eigen::Vector3d s2e2 = line2.second - line2.first;
		Eigen::Vector3d c12 = s1e1.cross(s2e2);
		if (c12.norm() <= 0.1)
		{
			LOG(INFO) << "Parallel line";
			return false;
		}

		Eigen::Vector3d s1s2 = line2.first - line1.first;
		Eigen::Vector3d cross1 = s1s2.cross(s1e1);
		Eigen::Vector3d cross2 = -s2e2.cross(s1e1);
		double scale = cross2.dot(cross1)/cross2.squaredNorm();
		if (scale < 0 || scale > 1)
		{
			LOG(INFO) << "There is no intersec";
			return false;
		}

		intersec = line2.first + scale*s2e2;
		return true;
	}

	bool calRuler3d(const std::vector<seg_pair_t>& wallBorder, const std::vector<seg_pair_t>& holeBorder, 
					 const seg_pair_t& rotateLine, const Eigen::Vector3d& P0,
					 const float& theta, seg_pair_t& ruler)
	{
		if (wallBorder.empty())
		{
			LOG(ERROR) << "empty input for ruler";
			return false;
		}
		
		if (P0 != rotateLine.first && P0 != rotateLine.second)
		{
			LOG(ERROR) << "P0 is not on the rotateLine for ruler";
			return false;
		}
		
		if (theta == 0)
		{
			ruler = rotateLine;
			return true;
		}

		//step0: Calculate normal vector
		seg_pair_t firstSeg = wallBorder[0];
		seg_pair_t bottomSeg = wallBorder.back();
		Eigen::Vector3d firstV = firstSeg.first - firstSeg.second;
		Eigen::Vector3d bottomV = bottomSeg.first - bottomSeg.second;
		Eigen::Vector3d nV = bottomV.cross(firstV);
		nV = nV.normalized();

		//step1: create R
		float radian = theta * geo::PI/180;
		Eigen::AngleAxisd t_V(radian, nV);
		Eigen::Matrix3d R = t_V.matrix();

		//step3: rotate
		Eigen::Vector3d rotateLineV = rotateLine.second - rotateLine.first;
		Eigen::Vector3d Q = R * rotateLineV + P0;
		seg_pair_t rotateQ = std::make_pair(P0, Q);

		//step4: Calculate the intersection point
		Eigen::Vector3d new3d;
		double rulerLength = 1000000.0;
		std::vector<seg_pair_t> allBorder = wallBorder;
		allBorder.insert(allBorder.end(), holeBorder.begin(), holeBorder.end());
		for (std::size_t i = 0; i < allBorder.size(); i++)
		{
			seg_pair_t line = allBorder[i];
			if (line == rotateLine)
				continue;

			Eigen::Vector3d s1e1 = line.second - line.first;
			Eigen::Vector3d s2e2 = rotateQ.second - rotateQ.first;
			Eigen::Vector3d c12 = s1e1.cross(s2e2);
			if (c12.norm() <= 0 + 1e-6 && ((line.first - P0).norm() < 1e-3 || (line.second - P0).norm() < 1e-3))
			{
				LOG(INFO) << "find the Vertical line"; // only for theta is 90
				ruler = line;
				LOG(INFO) <<"find the ruler end point is " << line.first[0] << " " <<line.first[1] 
					<< " " << line.first[2];
				LOG(INFO) <<"find the ruler end point is " << line.second[0] << " " <<line.second[1] 
					<< " " << line.second[2];
				return true;
			}
			
			Eigen::Vector3d intersec;
			if (calIntersection(rotateQ, line, intersec))
			{
				double length = (intersec - P0).norm();
				if (length <= 0 + 1e-6)
					continue;
				if (length < rulerLength)
				{
					rulerLength = length;
					new3d = intersec;
				}
			}
		}

		if (std::fabs(1000000.0 - rulerLength)< 1)
		{
			LOG(ERROR) << "can not find ruler, please check theta";
			return false;
		}
		LOG(INFO) <<"find the ruler end point is " << new3d[0] << " " <<new3d[1] << " " << new3d[2];
		ruler = std::make_pair(P0, new3d);
		return true;
	}

	//
	std::vector<seg_pair_t> calValidHoleVertical(
                        const std::vector<std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>>& holeBorders,
                        const std::pair<Eigen::Vector3d, Eigen::Vector3d>& horizen, int hAixs)
    {
        std::vector<seg_pair_t> validHoleVertical;
        const auto horizenSeg = horizen.first - horizen.second;
        for (size_t i = 0; i < holeBorders.size(); ++i)
        {
            auto &holeBorder = holeBorders[i];
            std::vector<seg_pair_t> vecHoleVertical, vecHoleHorizen;
            groupDirection(horizenSeg, holeBorder, vecHoleVertical, vecHoleHorizen);
            if (vecHoleHorizen.size() < 2 || vecHoleVertical.size() < 2)
            {
                LOG(ERROR) << "group Hole Direction Failed: " << vecHoleHorizen.size() << " -- " << vecHoleVertical.size();
			    continue;
            }
            std::sort(vecHoleVertical.begin(), vecHoleVertical.end(), [&](const seg_pair_t& left, const seg_pair_t& right){
				return left.first[hAixs] < right.first[hAixs];});
            validHoleVertical.emplace_back(vecHoleVertical[0]);
            validHoleVertical.emplace_back(vecHoleVertical.back());
        }

        return validHoleVertical;
    }

	int calWallHorizontalAxis(const seg_pair_t& seg)
    {
        double length = (seg.first - seg.second).norm();
		double len0 = std::fabs(length - std::fabs(seg.first[0] - seg.second[0]));
		double len1 = std::fabs(length - std::fabs(seg.first[1] - seg.second[1]));
		if (len0 < len1)
			return 0;
        return 1;
    }

	void cutOffRuler(seg_pair_t& ruler, double length)
    {
        double rLength = (ruler.first - ruler.second).norm();
        if (rLength <= length + 1e-4)  //200cm
            return;

        Eigen::Vector3d rulern = (ruler.second - ruler.first).normalized();
        Eigen::Vector3d rulerSecondEnd = ruler.first + length * rulern; 
        ruler.second = rulerSecondEnd;
    }

	std::vector<std::vector<Eigen::Vector3d>> getAllRulerBox(seg_pair_t ruler, int thicknessDir, 
									double thickness, double step, double boxLen, double boxWidth)
    {
        std::vector<std::vector<Eigen::Vector3d>> rulerBoxes;
        Eigen::Vector3d rulern = (ruler.second - ruler.first).normalized();
        double halfLen = boxLen / 2;
        auto vecPts = ininterpolateSeg(ruler.first,ruler.second,step);
        for (size_t i = 0; i < vecPts.size() - 1; ++i) //box length is 10mm
        {
            Eigen::Vector3d midPoint = vecPts[i];
            Eigen::Vector3d pt1 = midPoint - halfLen * rulern;
            Eigen::Vector3d pt2 = midPoint + halfLen * rulern;
            std::vector<Eigen::Vector3d> rulerB = createRulerBox(std::make_pair(pt1, pt2), 
                                                    thicknessDir, thickness, boxWidth); //
            rulerBoxes.emplace_back(rulerB);
        }

        return rulerBoxes;
    }

    calcMeassurment_t calFlatness(std::vector<seg_pair_t> rulers, int thicknessDir, Eigen::Vector4d plane, 
                                    pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud)
    {
        calcMeassurment_t measure;
        std::vector<std::vector<Eigen::Vector3d>> allBoxes;
		for(auto& ruler : rulers)
		{
			if ((ruler.first - ruler.second).norm() < 0.005)
				continue;
			auto boxes = getAllRulerBox(ruler, thicknessDir, 0., 0.005, 0.01, 0.025);
			allBoxes.insert(allBoxes.end(), boxes.begin(), boxes.end());
		}
        if (allBoxes.empty())
        {
            LOG(ERROR) << "empty boxes";
            return measure;
        }
        LOG(INFO) << "ruler get boxes num: " << allBoxes.size();

        std::vector<double> sumAll;
        for(size_t i = 0; i < allBoxes.size(); ++i)
        {
            auto box = allBoxes[i];
			std::vector<Eigen::Vector3d> corners = getRulerCorners(box);
			auto rangeCloud = filerCloudByConvexHull(pCloud, corners);

            if (rangeCloud->points.empty()) 
            {
                // LOG(ERROR) << "filerCloudByRange failed";
                continue;
            }
            double sum = 0;
            for (auto &p : rangeCloud->points)
                sum += std::fabs(pointToPLaneDist(plane, p));
            double avg = sum / rangeCloud->points.size();
            sumAll.emplace_back(avg);
        }

        if (sumAll.empty())
        {
            LOG(ERROR) << "empty sumAll, please check!";
            return measure;
        }
            
        double max =  *std::max_element(sumAll.begin(),sumAll.end());
        double min = *std::min_element(sumAll.begin(),sumAll.end());
        double difference = std::fabs(max - min);
        LOG(INFO) << "max avg: " << max << " min avg: " << min
                << " difference: " << difference;
        measure.value = difference;
        measure.rangeSeg.insert(measure.rangeSeg.end(), rulers.begin(), rulers.end());

        return measure;
    }

	bool judgeHoleBorder(const std::vector<std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>>& holeBorders,
                    std::pair<seg_pair_t, seg_pair_t> validWalls)
    {
        if (holeBorders.empty())
            return false;
       
        for(std::size_t i = 0; i < holeBorders.size(); i++)
        {
            auto &hole = holeBorders[i];
            auto iter1 = std::find(hole.begin(), hole.end(), validWalls.first);
            auto iter2 = std::find(hole.begin(), hole.end(), validWalls.second);
            if (iter1 != hole.end() && iter2 != hole.end())
                return true;
        }
        return false;
    }

}
