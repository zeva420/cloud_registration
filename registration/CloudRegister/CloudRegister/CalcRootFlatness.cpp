#include "CalcRootFlatness.h"

#include "funHelper.h"
// #define VISUALIZATION_ENABLED
namespace CloudReg
{

    int calHorizontalAxisRoot(const seg_pair_t& seg)
    {
        double length = (seg.first - seg.second).norm();
        if (std::fabs(length - std::fabs(seg.first[0] - seg.second[0])) < 0.2)
            return 1;
        else if (std::fabs(length - std::fabs(seg.first[1] - seg.second[1])) < 0.2)
            return 0;
        return 2;
    }

    std::vector<seg_pair_t> getDiagonalRuler(const std::pair<seg_pair_t, seg_pair_t>& zoomSeg)
    {
        std::vector<seg_pair_t> ruler;
        seg_pair_t seg1 = zoomSeg.first;
        seg_pair_t seg2 = std::make_pair(seg1.second, zoomSeg.second.second);
        seg_pair_t seg3 = std::make_pair(zoomSeg.second.second, zoomSeg.second.first);
        seg_pair_t seg4 = std::make_pair(zoomSeg.second.first, seg1.first);
        std::vector<seg_pair_t> zoomBorder = {seg1, seg2, seg3, seg4}; //Clockwise direction

        seg_pair_t ruler1, ruler2;
        std::vector<seg_pair_t> holeBorder;
        if(!calRuler3d(zoomBorder, holeBorder, seg2, seg2.first, 45, ruler1))
            return ruler;
        Eigen::Vector3d rulern1 = (ruler1.second - ruler1.first).normalized();
        ruler1.first = ruler1.first + 0.05*rulern1; // Move down 50mm
        ruler1.second = ruler1.second - 0.05*rulern1; // Move up 50mm

        if(!calRuler3d(zoomBorder, holeBorder, seg4, seg4.first, 45, ruler2))
            return ruler;
        Eigen::Vector3d rulern2 = (ruler2.second - ruler2.first).normalized();
        ruler2.first = ruler2.first + 0.05*rulern2; // Move up 50mm
        ruler2.second = ruler2.second - 0.05*rulern2; // Move down 50mm

        cutOffRuler(ruler1, 2);
        cutOffRuler(ruler2, 2);
        
        ruler.emplace_back(ruler1);
        ruler.emplace_back(ruler2);
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
        double baseHeight = (zoomSeg.first.first[vAixs] < zoomSeg.first.second[vAixs]) ? 
                            zoomSeg.first.first[vAixs] : zoomSeg.first.second[vAixs];
        Eigen::Vector3d hypRuler0 = zoomSeg.first.first;
        hypRuler0[vAixs] = baseHeight + height / 2;
        Eigen::Vector3d hypRuler1 = zoomSeg.second.first;
        hypRuler1[vAixs] = baseHeight + height / 2;
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

        double resDis = (length - 4) / 3;
        Eigen::Vector3d p1 = hypRuler0;
        p1[hAixs] = p1[hAixs] + resDis;
        Eigen::Vector3d p2 = p1 + 2*rulern;   //200cm
        Eigen::Vector3d p3 = p2;
        p3[hAixs] = p3[hAixs] + resDis;
        Eigen::Vector3d p4 = p3 + 2*rulern;
        ruler.emplace_back(std::make_pair(p1, p2));
        ruler.emplace_back(std::make_pair(p3, p4));
        return ruler;
    }

    std::vector<std::tuple<std::vector<calcMeassurment_t>, std::vector<seg_pair_t>>>
    calRootFlatness(const std::vector<seg_pair_t>& rootBorder,
                          Eigen::Vector4d plane, pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud,
                          const double calcLengthTh)
    {
        // step0: cal Horizontal Axis
        std::vector<std::tuple<std::vector<calcMeassurment_t>, std::vector<seg_pair_t>>> returnValue;
        const Eigen::Vector3d& baseSeg = rootBorder.front().first - rootBorder.front().second;
        int hAxis = calHorizontalAxisRoot(rootBorder.front());
        int vAxis = (hAxis == 0) ? 1 : 0;
        LOG(INFO) << "Horizontal Axis is " << hAxis;
        LOG(INFO) << "Vertical Axis is " << vAxis;

        //step1:  get vecCutBay
        std::vector<seg_pair_t> vecCutSeg;
		std::vector<std::size_t> vecVerticalIndex;
		std::vector<std::size_t> vecHorizenIndex;
		groupDirectionIndex(baseSeg, rootBorder, vecVerticalIndex, vecHorizenIndex);
		
		for(std::size_t i = 0; i< vecVerticalIndex.size(); i++)
		{
			seg_pair_t toSeg = rootBorder[vecVerticalIndex[i]];
			if ((toSeg.first - toSeg.second).norm() < calcLengthTh)
				continue;

			for (int j = vecVerticalIndex.size() - 1; j >= i + 1; j--)
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
            return returnValue;
        }
        LOG(INFO) << "The root is divided into areas num: " << vecCutSeg.size() / 2;
       
        std::vector<seg_pair_t> vecRange;
        for (std::size_t i = 0; i < vecCutSeg.size() / 2; ++i)
        {
            std::vector<calcMeassurment_t> tmpMeasure;
            auto HorizenSeg1 = vecCutSeg[2*i];
            auto HorizenSeg2 = vecCutSeg[2*i + 1];
            std::vector<seg_pair_t> tmpSegs = {HorizenSeg1, HorizenSeg2};

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
            LOG(INFO) << "zoom get rulers num: " << rulers.size();

            for (auto &ruler : rulers)
            {
                calcMeassurment_t measure;
                std::vector<seg_pair_t> tmp = {ruler};
                measure = calFlatness(tmp, 2, plane, pCloud);
                
                std::vector<Eigen::Vector3d> rPoints =  createRulerBox(ruler, 2, 0.025, 0.025); //width 2.5cm
                std::vector<seg_pair_t> pair =  calcBoxSegPair(rPoints);
                measure.rangeSeg.insert(measure.rangeSeg.end(), pair.begin(), pair.end());
                vecRange.insert(vecRange.end(), pair.begin(), pair.end());
                tmpMeasure.emplace_back(measure);
            }
            returnValue.emplace_back(std::make_tuple(tmpMeasure, tmpSegs));
        }
#ifdef VISUALIZATION_ENABLED
        pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud_filtered(new pcl::PointCloud<pcl::PointXYZ>());
        uniformSampling(0.01, pCloud, pCloud_filtered);
        writePCDFile("RootFlatness.pcd", pCloud_filtered, vecRange);
#endif
        return returnValue;
    }
}