#include "CalcWallVerticality.h"

#include "funHelper.h"
#include <pcl/common/common.h>
#include <pcl/filters/crop_box.h>
#include <pcl/filters/passthrough.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>

namespace CloudReg
{
    Eigen::Vector4d calcWallPlane(std::vector<seg_pair_t> vecVertical)
	{
		pcl::PointCloud<pcl::PointXYZ> cloud;
		cloud.width = vecVertical.size()*2;
		cloud.height = 1;
		cloud.is_dense = false;
		cloud.points.resize(cloud.width * cloud.height);
		
		for (size_t i = 0; i < vecVertical.size(); ++i)
		{
			auto& ptA = vecVertical[i].first;
			auto& ptB = vecVertical[i].second;

			cloud.points[i*2].x = ptA[0];
			cloud.points[i*2].y = ptA[1];
			cloud.points[i*2].z = ptA[2];
			
			cloud.points[i*2 + 1].x = ptB[0];
			cloud.points[i*2 + 1].y = ptB[1];
			cloud.points[i*2 + 1].z = ptB[2];
		}
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_ptr(new pcl::PointCloud<pcl::PointXYZ>);
        cloud_ptr=cloud.makeShared();
        Eigen::Vector4d cadPlane =  calcPlaneParam(cloud_ptr);
        return cadPlane;
	}

    std::vector<Eigen::Vector3d> calBox(Eigen::Vector3d midPt, int type, int hAixs, Eigen::Vector3d rulern, calcMeassurment_t& measure)
    {
        double boxWidth = 0.025;
        double boxHight = 0.025;

        Eigen::Vector3d p1 = midPt;
        Eigen::Vector3d p2;
        if (type == 1)
            p2 = midPt + boxHight/2 * rulern;
        else
            p2 = midPt - boxHight/2 * rulern;

        int thicknessDir = (hAixs == 0) ? 1 : 0;
        std::vector<Eigen::Vector3d> rPoints =  createRulerBox(std::make_pair(p1, p2), thicknessDir, 0., boxWidth);
        // std::vector<seg_pair_t> pair =  calcBoxSegPair(rPoints);
        // // measure.rangeSeg.insert(measure.rangeSeg.end(), pair.begin(), pair.end());
        return getRulerCorners(rPoints);
    }

    void adjustHeight(const std::vector<seg_pair_t>& vecWallHorizen, double adjDis, 
                        Eigen::Vector3d pt1, Eigen::Vector3d& pt2, Eigen::Vector3d rulern)
    {
        for (size_t i = 1; i < vecWallHorizen.size(); ++i) //first one is the bottom
        {
            seg_pair_t horizen = vecWallHorizen[i];
            seg_pair_t pt12 = std::make_pair(pt1, pt2);
            Eigen::Vector3d intersec;
			if (calIntersection(pt12, horizen , intersec))
            {
                pt2[2] = horizen.second[2] - adjDis;
                LOG(INFO) << "adjust ruler height success"; 
                return;
            }
        }
        LOG(ERROR) << "adjust ruler height filed";
    }

    calcMeassurment_t calOneRulerData(const std::vector<seg_pair_t>& vecWallHorizen,
                                        double length,
                                        seg_pair_t Vertical1, seg_pair_t Vertical2, int type,
                                        const PointCloud::Ptr pCloud,
                                        Eigen::Vector4d cadPlane, int hAixs)
    {
        Eigen::Vector3d pt1 = (Vertical1.first[2] < Vertical1.second[2] + 1e-6) ? Vertical1.first : Vertical1.second;
        Eigen::Vector3d pt2 = (Vertical1.first == pt1) ? Vertical1.second : Vertical1.first;
        Eigen::Vector3d pt3 = (Vertical2.first[2] < Vertical2.second[2] + 1e-6) ? Vertical2.first : Vertical2.second;
        Eigen::Vector3d rulern = (pt2 - pt1).normalized();
        Eigen::Vector3d horizenn = (pt3 - pt1).normalized();
        if (type == 1)  //left ruler
        {
            pt1 = pt1 + 0.3 * horizenn;
            pt2 = pt2 + 0.3 * horizenn;
            pt2 = pt2 - 0.2 * rulern;
            adjustHeight(vecWallHorizen, 0.2, pt1, pt2, rulern); //For heterosexual walls
        }
        else if (type == 2)  //right ruler
        {
            pt1 = pt1 + 0.3 * horizenn;
            pt2 = pt2 + 0.3 * horizenn;
            pt1 = pt1 + 0.2 * rulern;
            adjustHeight(vecWallHorizen, 0., pt1, pt2, rulern);
        }
        else   //middle ruler
        {
            pt1 = pt1 + length / 2 * horizenn;
            pt2 = pt2 + length / 2 * horizenn;
            pt1 = pt1 + 0.1 * rulern;
            pt2 = pt2 - 0.1 * rulern;
            adjustHeight(vecWallHorizen, 0.1, pt1, pt2, rulern);
        }
        
       
        Eigen::Vector3d valid1, valid2;
        std::vector<double> sumAll;
        calcMeassurment_t item;
        item.rangeSeg.emplace_back(std::make_pair(pt1, pt2));
        auto vecRulerPts = ininterpolateSeg(pt1, pt2, 0.025);
        LOG(INFO) << "ruler get boxes num: " << vecRulerPts.size();

        for(size_t i = 0; i < vecRulerPts.size() - 1; ++i)
        {
            auto pt = vecRulerPts[i];
            auto corners = calBox(pt, type, hAixs, rulern, item);
            auto rangeCloud = filerCloudByConvexHull(pCloud, corners);
            if (rangeCloud->points.empty()) 
            {
                // LOG(ERROR) << "filerCloudByRange failed";
                continue;
            }
            double sum = 0;
            for (auto &p : rangeCloud->points)
                sum += std::fabs(pointToPLaneDist(cadPlane, p));
            double avg = sum / rangeCloud->points.size();
            sumAll.emplace_back(avg);
            valid1 = pt;
            break;
        }

        for(size_t i = vecRulerPts.size() - 1; i > 0 ; --i)
        {
            auto pt = vecRulerPts[i];
            auto corners = calBox(pt, type, hAixs, rulern, item);
            auto rangeCloud = filerCloudByConvexHull(pCloud, corners);
            if (rangeCloud->points.empty()) 
            {
                // LOG(ERROR) << "filerCloudByRange failed";
                continue;
            }
            double sum = 0;
            for (auto &p : rangeCloud->points)
                sum += std::fabs(pointToPLaneDist(cadPlane, p));
            double avg = sum / rangeCloud->points.size();
            sumAll.emplace_back(avg);
            valid2 = pt;
            break;
        }

        item.value = 0.;
        int minIndex = (hAixs == 0) ? 1 : 0;
        if (sumAll.size() < 2 || (valid1 - valid2).norm() < 0.005)
            LOG(ERROR) << "ruler has less than 2 endpoints";
        else
        {
            item.value = std::fabs(sumAll[0] - sumAll[1]);
            LOG(INFO) << "verticality avg is " << item.value;
            std::vector<Eigen::Vector3d> rPoints = createRulerBox(std::make_pair(valid1, valid2), minIndex, 0.025, 0.025);
            std::vector<seg_pair_t> pair =  calcBoxSegPair(rPoints);
            item.rangeSeg.insert(item.rangeSeg.end(), pair.begin(), pair.end());
        }
        return item;
    }

    std::map<int, std::vector<calcMeassurment_t>> calWallVerticality(const std::vector<seg_pair_t>& vecWallHorizen,
                                std::pair<seg_pair_t, seg_pair_t> validWalls,
                                const PointCloud::Ptr pCloud,
                                Eigen::Vector4d cadPlane, int hAixs)
    {
        std::map<int, std::vector<calcMeassurment_t>> vecRet;
        auto &Vertical1 = validWalls.first;
        auto &Vertical2 = validWalls.second;

        double length = std::fabs(Vertical1.first[hAixs] - Vertical2.first[hAixs]);
        LOG(INFO) << "wall length: " << length;
		if (length <= 0.1f + 1e-3)
		{
			LOG(INFO) << "too short length: " << length;
			return vecRet;
		}

        if (length > 0.1f && length  <= 0.6 + 1e-3)
        {
            auto item = calOneRulerData(vecWallHorizen, length, Vertical1, Vertical2, 3, pCloud, cadPlane, hAixs);
            if(!item.rangeSeg.empty())
                vecRet[3].emplace_back(item);
            return vecRet;
        }
        
        auto item1 = calOneRulerData(vecWallHorizen, length, Vertical1, Vertical2, 1, pCloud, cadPlane, hAixs);
        if(!item1.rangeSeg.empty())
            vecRet[1].emplace_back(item1);
        auto item2 = calOneRulerData(vecWallHorizen, length, Vertical2, Vertical1, 2, pCloud, cadPlane, hAixs);
        if(!item2.rangeSeg.empty())
            vecRet[2].emplace_back(item2);
        if (length >= 3 + 1e-6)
        {
            auto item = calOneRulerData(vecWallHorizen, length, Vertical1, Vertical2, 3, pCloud, cadPlane, hAixs);
            if(!item.rangeSeg.empty())
                vecRet[3].emplace_back(item);
        }
        return vecRet;
    }

    int calHorizontalAxis(const seg_pair_t& seg)
    {
        double length = (seg.first - seg.second).norm();
        if (std::fabs(length - std::fabs(seg.first[0] - seg.second[0])) < 1e-4)
            return 0;
        else if (std::fabs(length - std::fabs(seg.first[1] - seg.second[1])) < 1e-4)
            return 1;
        return 2;
    }

    std::tuple<std::vector<calcMeassurment_t>, std::vector<seg_pair_t>> calcVerticality(
			const std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>& wallBorder, 
            const std::vector<std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>>& holeBorders,
			pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud, Eigen::Vector4d plane, int wallIndex)
    {
        std::vector<calcMeassurment_t> returnMeasure;
        std::vector<seg_pair_t> returnSeg;
        // step0: cal Horizontal Axis
        auto horizen = wallBorder.back();
        int hAxis = calHorizontalAxis(horizen);
        LOG(INFO) <<"Wall " << wallIndex << " wallBorder size is " << wallBorder.size();
        LOG(INFO) << "Horizontal Axis is " << hAxis;

        // step1: Sort the data of the wall and hole
        std::vector<seg_pair_t> vecWallHorizen, vecWallVertical;	
		const auto horizenSeg = horizen.first - horizen.second;
        groupDirection(horizenSeg, wallBorder, vecWallVertical, vecWallHorizen);
        if (vecWallHorizen.size() < 2 || vecWallVertical.size() < 2)
		{
			LOG(ERROR) << "group Wall Direction Failed: " << vecWallHorizen.size() << " -- " << vecWallVertical.size();
			return std::make_tuple(returnMeasure, returnSeg);
		}
        LOG(INFO) << "plane " << plane[0] << " " << plane[1] << " " << plane[2] << " " << plane[3];
        LOG(INFO) << "group Wall Direction: " << vecWallHorizen.size() << " -- " << vecWallVertical.size();

        //step2: Count the number of separated walls
        auto validVertical = calValidHoleVertical(holeBorders, horizen, hAxis);
        LOG(INFO) << "valid hole Vertical num is " << validVertical.size();
        std::sort(validVertical.begin(), validVertical.end(), [&](const seg_pair_t& left, const seg_pair_t& right){
				return left.first[hAxis] < right.first[hAxis];});
        std::sort(vecWallVertical.begin(), vecWallVertical.end(), [&](const seg_pair_t& left, const seg_pair_t& right){
				return left.first[hAxis] < right.first[hAxis];});
        validVertical.insert(validVertical.begin(), vecWallVertical[0]); 
        validVertical.emplace_back(vecWallVertical.back());
        LOG(INFO) << "valid Wall num is " << validVertical.size() / 2;

        //step3: Calculate verticality
        std::map<int, std::vector<calcMeassurment_t>> allMeasure;
        for (std::size_t i = 0; i < validVertical.size() - 1; ++i)
        {
            auto &Vertical1 = validVertical[i];
            auto &Vertical2 = validVertical[i + 1];
            if (judgeHoleBorder(holeBorders, std::make_pair(Vertical1, Vertical2)))
                continue;
            
            auto measure = calWallVerticality(vecWallHorizen, std::make_pair(Vertical1, Vertical2), 
                                                pCloud, plane, hAxis);
            for (auto &map : measure)
                allMeasure[map.first].insert(allMeasure[map.first].end(), map.second.begin(), map.second.end());
        }

        LOG(INFO) << "Final output of wall verticality: ";
        std::vector<int> outOrder = {1, 2, 3};
        std::vector<seg_pair_t> vecRange;
        for(auto& type : outOrder)
        {
            if (!allMeasure.count(type))
                continue;
            auto measures = allMeasure[type];
            returnMeasure.insert(returnMeasure.end(), measures.begin(), measures.end());
            for(auto& item : measures)
            {
                vecRange.insert(vecRange.end(), item.rangeSeg.begin(), item.rangeSeg.end());
                LOG(INFO) << "Wall "<<wallIndex<<" verticality: " << item.value;
            }
        }
        writePCDFile("WallVerticality-" + std::to_string(wallIndex) + ".pcd", pCloud, vecRange);
        return std::make_tuple(returnMeasure, returnSeg);
    }

}