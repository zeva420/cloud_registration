#include "CalcRootFlatness.h"

#include "funHelper.h"
#include <pcl/common/common.h>
#include <pcl/filters/crop_box.h>
#include <pcl/filters/passthrough.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>

namespace CloudReg
{
    
    int calHorizontalAxisRoot(const seg_pair_t& seg)
    {
        double length = (seg.first - seg.second).norm();
        if (std::fabs(length - std::fabs(seg.first[0] - seg.second[0])) < 1e-4)
            return 1;
        else if (std::fabs(length - std::fabs(seg.first[1] - seg.second[1])) < 1e-4)
            return 0;
        return 2;
    }

    void cutOffRuler(seg_pair_t& ruler)
    {
        double length = (ruler.first - ruler.second).norm();
        if (length + 1e-4 <= 2 )  //200cm
            return;

        Eigen::Vector3d rulern = (ruler.second - ruler.first).normalized();
        Eigen::Vector3d rulerSecondEnd = ruler.first + 2*rulern;
        ruler.second = rulerSecondEnd;
    }

    std::vector<seg_pair_t> getDiagonalRuler(const std::pair<seg_pair_t, seg_pair_t>& zoomSeg)
    {
        seg_pair_t seg1 = zoomSeg.first;
        seg_pair_t seg2 = std::make_pair(seg1.second, zoomSeg.second.second);
        seg_pair_t seg3 = std::make_pair(zoomSeg.second.second, zoomSeg.second.first);
        seg_pair_t seg4 = std::make_pair(zoomSeg.second.first, seg1.first);
        std::vector<seg_pair_t> zoomBorder = {seg1, seg2, seg3, seg4}; //Clockwise direction

        seg_pair_t ruler1, ruler2;
        std::vector<seg_pair_t> holeBorder;
        calRuler3d(zoomBorder, holeBorder, seg2, seg2.first, 45, ruler1);
        Eigen::Vector3d rulern1 = (ruler1.second - ruler1.first).normalized();
        ruler1.first = ruler1.first + 0.05*rulern1; // Move down 50mm

        calRuler3d(zoomBorder, holeBorder, seg4, seg4.first, 45, ruler2);
        Eigen::Vector3d rulern2 = (ruler2.second - ruler2.first).normalized();
        ruler2.first = ruler2.first + 0.05*rulern2; // Move up 50mm

        cutOffRuler(ruler1);
        cutOffRuler(ruler2);
        
        std::vector<seg_pair_t> ruler = {ruler1, ruler2};
        return ruler;
    }

    std::vector<seg_pair_t> getMiddleRuler(const std::pair<seg_pair_t, seg_pair_t>& zoomSeg, 
                                            int vAixs, int hAixs)
    {
        std::vector<seg_pair_t> ruler;
        double length = (zoomSeg.first.second - zoomSeg.second.second).norm();
        double height = (zoomSeg.first.first - zoomSeg.first.second).norm();
        if (length + 1e-3 < 3)
            return ruler;
        LOG(INFO) << "Zoom length " << length << " height " << height;
        Eigen::Vector3d hypRuler0 = zoomSeg.first.first;
        hypRuler0[vAixs] = hypRuler0[vAixs] + height / 2;
        Eigen::Vector3d hypRuler1 = zoomSeg.second.first;
        hypRuler1[vAixs] = hypRuler1[vAixs] + height / 2;
        Eigen::Vector3d hypRuler = hypRuler1 - hypRuler0;
        Eigen::Vector3d rulern = hypRuler.normalized();
        if (length + 1e-3 >= 3  && length + 1e-4 < 6 ) //Measure once
        {
            double resDis = (length - 2) / 2;
            Eigen::Vector3d p1 = hypRuler0;
            p1[hAixs] = p1[hAixs] + resDis;
            Eigen::Vector3d p2 = p1 + 2*rulern;    //200 cm
            ruler.emplace_back(std::make_pair(p1, p2));
            return ruler;
        }

        double resDis = (length - 2) / 3;
        Eigen::Vector3d p1 = hypRuler0;
        p1[hAixs] = p1[hAixs] + resDis;
        Eigen::Vector3d p2 = p1 + 2*rulern;
        Eigen::Vector3d p3 = p2;
        p3[hAixs] = p3[hAixs] + resDis;
        Eigen::Vector3d p4 = p3 + 2*rulern;
        ruler.emplace_back(std::make_pair(p1, p2));
        ruler.emplace_back(std::make_pair(p3, p4));
        return ruler;
    }

    std::vector<std::vector<Eigen::Vector3d>> getAllRulerBox(seg_pair_t ruler, int thicknessDir, 
                                        double step, double boxLen, double boxWidth)
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
                                                    thicknessDir, 0., boxWidth); //thickness is tmp 0
            rulerBoxes.emplace_back(rulerB);
        }

        return rulerBoxes;
    }

    calcMeassurment_t calFlatness(seg_pair_t ruler, int thicknessDir, Eigen::Vector4d plane, 
                                    pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud)
    {
        calcMeassurment_t measure;
        std::vector<std::vector<Eigen::Vector3d>> allBoxes;
        allBoxes = getAllRulerBox(ruler, thicknessDir, 0.005, 0.01, 0.025);
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
            
            seg_pair_t seg1, seg2;
            seg1 = std::make_pair(box[0], box[2]); //pt1 pt3
            seg2 = std::make_pair(box[2], box[6]); //pt3 pt7
            std::vector<seg_pair_t> calcSeg = {seg1, seg2};

            pcl::PointXYZ min;
		    pcl::PointXYZ max;
            pcl::PointCloud<pcl::PointXYZ> cloud;
            cloud.width = calcSeg.size()*2;
            cloud.height = 1;
            cloud.is_dense = false;
            cloud.points.resize(cloud.width * cloud.height);
            
            for (size_t i = 0; i < calcSeg.size(); ++i)
            {
                auto& ptA = calcSeg[i].first;
                auto& ptB = calcSeg[i].second;

                cloud.points[i*2].x = ptA[0];
                cloud.points[i*2].y = ptA[1];
                cloud.points[i*2].z = ptA[2];
                
                cloud.points[i*2 + 1].x = ptB[0];
                cloud.points[i*2 + 1].y = ptB[1];
                cloud.points[i*2 + 1].z = ptB[2];
            }
            pcl::getMinMax3D(cloud, min, max);
            auto rangeCloud = filerCloudByRange(pCloud,min,max);
            if (rangeCloud->points.empty()) 
            {
                LOG(ERROR) << "filerCloudByRange failed";
                continue;
            }
            double sum = 0;
            for (auto &p : rangeCloud->points)
                sum += pointToPLaneDist(plane, p);
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
        measure.rangeSeg.emplace_back(ruler);

        return measure;
    }

    void calcRootFlatness(const std::vector<seg_pair_t>& rootBorder,
                          Eigen::Vector4d plane, pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud,
                          const double calcLengthTh)
    {
        // step0: cal Horizontal Axis
        const Eigen::Vector3d& baseSeg = rootBorder.front().first - rootBorder.front().second;
        int hAxis = calHorizontalAxisRoot(rootBorder.front());
        int vAxis = (hAxis == 0) ? 1 : 0;
        LOG(INFO) << "Horizontal Axis is " << hAxis;
        LOG(INFO) << "Horizontal vAxis is " << vAxis;

        //step1:  get vecCutBay
		std::vector<std::size_t> vecVerticalIndex;
		std::vector<std::size_t> vecHorizenIndex;
		groupDirectionIndex(baseSeg, rootBorder, vecVerticalIndex, vecHorizenIndex);
		
		std::vector<seg_pair_t> vecCutSeg;
		for(std::size_t i = 0; i< vecVerticalIndex.size(); i++)
		{
			seg_pair_t toSeg = rootBorder[vecVerticalIndex[i]];
			if ((toSeg.first - toSeg.second).norm() < calcLengthTh)
				continue;

			for(std::size_t j = i+1; j< vecVerticalIndex.size(); j++)
			{
				seg_pair_t calcSeg = rootBorder[vecVerticalIndex[j]];
				if ((calcSeg.first - calcSeg.second).norm() < calcLengthTh)
					continue;
				
				bool hasOverlap;
				Eigen::Vector3d s1Pt, e1Pt, s2Pt,e2Pt;
				std::tie(hasOverlap, s1Pt, e1Pt, s2Pt,e2Pt) = calcOverlap(toSeg,calcSeg);
				if (!hasOverlap) continue;

				if((s1Pt-e1Pt).norm() < calcLengthTh 
					|| (s2Pt-e2Pt).norm() < calcLengthTh)
				{
					LOG(INFO)<< "too short seg:" << vecHorizenIndex[i] <<" " << vecHorizenIndex[j];
					continue;
				}
                vecCutSeg.emplace_back(std::make_pair(s1Pt,s2Pt));
				vecCutSeg.emplace_back(std::make_pair(e1Pt,e2Pt));
            }
        }

        // step2: Sort the data of the zone division and calculate
        if (vecCutSeg.empty() || (vecCutSeg.size() % 2) != 0)
        {
            LOG(ERROR) << "root bay size is wrong " << vecCutSeg.size();
            return;
        }
       
        std::vector<calcMeassurment_t> allMeasure;
        std::vector<seg_pair_t> rulerAll;
        for (std::size_t i = 0; i < vecCutSeg.size() / 2; ++i)
        {
            auto HorizenSeg1 = vecCutSeg[2*i];
            auto HorizenSeg2 = vecCutSeg[2*i + 1];
            seg_pair_t firstSeg = (HorizenSeg1.first[hAxis] < HorizenSeg2.first[hAxis]) ?
                                    HorizenSeg1 : HorizenSeg2;
            seg_pair_t secondSeg = (firstSeg == HorizenSeg1) ? HorizenSeg2 : HorizenSeg1;

            seg_pair_t validSeg1, validSeg2;
            validSeg1.first = (firstSeg.first[vAxis] < firstSeg.second[vAxis]) ? firstSeg.first : firstSeg.second;
            validSeg1.second = (validSeg1.first == firstSeg.first) ? firstSeg.second : firstSeg.first;

            validSeg2.first = (secondSeg.first[vAxis] < secondSeg.second[vAxis]) ? secondSeg.first : secondSeg.second;
            validSeg2.second = (validSeg2.first == secondSeg.first) ? secondSeg.second : secondSeg.first;

            std::vector<seg_pair_t> rulers = getDiagonalRuler(std::make_pair(validSeg1, validSeg2));
            std::vector<seg_pair_t> middleRuler = getMiddleRuler(std::make_pair(validSeg1, validSeg2), vAxis, hAxis);
            if (!middleRuler.empty())
                rulers.insert(rulers.end(), middleRuler.begin(), middleRuler.end());
            rulerAll.insert(rulerAll.end(), rulers.begin(), rulers.end());
            LOG(INFO) << "zoom get rulers num: " << rulers.size();

            for (auto &ruler : rulers)
            {
                calcMeassurment_t measure;
                measure = calFlatness(ruler, 2, plane, pCloud);
                allMeasure.emplace_back(measure);
            }
        }

       //For visualization
        std::vector<seg_pair_t> vecRange;
        for(auto ruler : rulerAll)
        {
            std::vector<Eigen::Vector3d> rPoints =  createRulerBox(ruler, 2, 0.025, 0.025);
            for (size_t i = 0; i < rPoints.size(); ++i)
            {
                for(size_t j = 0; j < rPoints.size(); ++j)
                {
                    if (i == j) continue;
                    seg_pair_t pp = std::make_pair(rPoints[i], rPoints[j]);
                    vecRange.emplace_back(pp);
                }
            }
        }
        writePCDFile("testFla.pcd", pCloud, vecRange);
    }
}