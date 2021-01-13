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

    std::vector<Eigen::Vector3d> calBox(Eigen::Vector3d midPt, int type, int hAixs)
    {
        double boxWidth = 0.025;
        double boxHight = 0.025;
       
        Eigen::Vector3d p1 = midPt;
        Eigen::Vector3d p2 = midPt;
        Eigen::Vector3d p3;
        Eigen::Vector3d p4;
        p1[hAixs] = midPt[hAixs] + boxWidth / 2;
        p2[hAixs] = midPt[hAixs] - boxWidth / 2;
        if (type == 1)
            p3 = Eigen::Vector3d(p2[0], p2[1], p2[2] + boxHight/2);
        else
            p3 = Eigen::Vector3d(p2[0], p2[1], p2[2] - boxHight/2);
        p4 = p3;
        p4[hAixs] = p3[hAixs] + boxWidth;
        
        std::vector<Eigen::Vector3d> corners = {p1, p2, p3, p4};
        return corners;
    }

    void adjustHeight(const std::vector<seg_pair_t>& vecWallHorizen, double adjDis, 
                        Eigen::Vector3d pt1, Eigen::Vector3d& pt2)
    {
        for (size_t i = 1; i < vecWallHorizen.size(); ++i) //first one is the bottom
        {
            seg_pair_t horizen = vecWallHorizen[i];
            seg_pair_t pt12 = std::make_pair(pt1, pt2);
            Eigen::Vector3d intersec;
			if (calIntersection(pt12, horizen , intersec))
            {
                pt2[2] = horizen.second[2] - adjDis;
                LOG(ERROR) << "adjust ruler height success"; 
                return;
            }
        }
        LOG(ERROR) << "adjust ruler height filed";
    }

    calcMeassurment_t calOneRulerData(const std::vector<seg_pair_t>& vecWallHorizen,
                                        double length,
                                        seg_pair_t Vertical1, int type,
                                        const PointCloud::Ptr pCloud,
                                        Eigen::Vector4d cadPlane, int hAixs)
    {
        Eigen::Vector3d pt1 = (Vertical1.first[2] < Vertical1.second[2] + 1e-6) ? Vertical1.first : Vertical1.second;
        Eigen::Vector3d pt2 = (Vertical1.first == pt1) ? Vertical1.second : Vertical1.first;
        if (type == 1)  //left ruler
        {
            pt1[hAixs] = pt1[hAixs] + 0.3;  //300mm
            pt2[hAixs] = pt2[hAixs] + 0.3;
            pt2[2] = pt2[2] - 0.2; //200mm
            adjustHeight(vecWallHorizen, 0.2, pt1, pt2); //For heterosexual walls
        }
        else if (type == 2)  //right ruler
        {
            pt1[hAixs] = pt1[hAixs] - 0.3;  //300mm
            pt1[2] = pt1[2] + 0.2;          //200mm
            pt2[hAixs] = pt2[hAixs] - 0.3;
            adjustHeight(vecWallHorizen, 0., pt1, pt2);
        }
        else   //middle ruler
        {
            pt1[hAixs] = pt1[hAixs] + length / 2;  //300mm
            pt1[2] = pt1[2] + 0.1;
            pt2[hAixs] = pt2[hAixs] + length / 2;
            pt2[2] = pt2[2] - 0.1;
            adjustHeight(vecWallHorizen, 0.1, pt1, pt2);
        }
        
        Eigen::vector<Eigen::Vector3d> vecCalcPt;
        vecCalcPt.emplace_back(pt1);
        vecCalcPt.emplace_back(pt2);
        
        std::vector<double> sumAll;
		for(auto& pt : vecCalcPt)
		{
            auto corners = calBox(pt, type, hAixs);
            auto rangeCloud = filerCloudByConvexHull(pCloud, corners);
            if (rangeCloud->points.empty()) 
            {
                LOG(ERROR) << "filerCloudByRange failed";
                continue;
            }
            double sum = 0;
            for (auto &p : rangeCloud->points)
                sum += pointToPLaneDist(cadPlane, p);
            double avg = sum / rangeCloud->points.size();
            sumAll.emplace_back(avg);
		}

        calcMeassurment_t item;
        item.value = 0.;
        if (sumAll.size() < 2)
            LOG(ERROR) << "ruler has less than 2 endpoints";
        else
            item.value = std::fabs(sumAll[0] - sumAll[1]);
        
        LOG(INFO) << "verticality avg is " << item.value;
        item.rangeSeg.emplace_back(std::make_pair(pt1, pt2));

        //just for Visualization
        int minIndex = (hAixs == 0) ? 1 : 0;
        std::vector<Eigen::Vector3d> rPoints = createRulerBox(std::make_pair(pt1, pt2), minIndex, 0.025, 0.025);
        std::vector<seg_pair_t> pair =  calcBoxSegPair(rPoints);
        item.rangeSeg.insert(item.rangeSeg.end(), pair.begin(), pair.end());

        return item;
    }

    std::vector<calcMeassurment_t> calWallVerticality(const std::vector<seg_pair_t>& vecWallHorizen,
                                std::pair<seg_pair_t, seg_pair_t> validWalls,
                                const PointCloud::Ptr pCloud,
                                Eigen::Vector4d cadPlane, int hAixs)
    {
        std::vector<calcMeassurment_t> vecRet;
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
            auto item = calOneRulerData(vecWallHorizen, length, Vertical1, 3, pCloud, cadPlane, hAixs);
            vecRet.emplace_back(item);
            return vecRet;
        }
        
        auto item1 = calOneRulerData(vecWallHorizen, length, Vertical1, 1, pCloud, cadPlane, hAixs);
        vecRet.emplace_back(item1);
        auto item2 = calOneRulerData(vecWallHorizen, length, Vertical2, 2, pCloud, cadPlane, hAixs);
        vecRet.emplace_back(item2);
        if (length >= 3 + 1e-6)
        {
            auto item = calOneRulerData(vecWallHorizen, length, Vertical1, 3, pCloud, cadPlane, hAixs);
            vecRet.emplace_back(item);
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

    void calcVerticality(const std::pair<Eigen::Vector3d, Eigen::Vector3d>& horizen,
			const std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>& wallBorder, 
            const std::vector<std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>>& holeBorders,
			pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud, int index)
    {
        // step0: cal Horizontal Axis
        int hAxis = calHorizontalAxis(horizen);
        LOG(INFO) << "Horizontal Axis is " << hAxis;

        // step1: Sort the data of the wall and hole
        std::vector<seg_pair_t> vecWallHorizen, vecWallVertical;	
		const auto horizenSeg = horizen.first - horizen.second;
        groupDirection(horizenSeg, wallBorder, vecWallVertical, vecWallHorizen);
        if (vecWallHorizen.size() < 2 || vecWallVertical.size() < 2)
		{
			LOG(ERROR) << "group Wall Direction Failed: " << vecWallHorizen.size() << " -- " << vecWallVertical.size();
			return;
		}

        Eigen::Vector4d cadPlane = calcWallPlane(vecWallVertical);
        LOG(INFO) << "cadPlane " << cadPlane[0] << " " << cadPlane[1] << " " << cadPlane[2] << " "
                    << cadPlane[3];

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
        std::vector<calcMeassurment_t> allMeasure;
        for (std::size_t i = 0; i < validVertical.size() - 1; ++i)
        {
            auto &Vertical1 = validVertical[i];
            auto &Vertical2 = validVertical[i + 1];
            if (judgeHoleBorder(holeBorders, std::make_pair(Vertical1, Vertical2)))
                continue;
            
            auto measure = calWallVerticality(vecWallHorizen, std::make_pair(Vertical1, Vertical2), 
                                                pCloud, cadPlane, hAxis);
            allMeasure.insert(allMeasure.end(), measure.begin(), measure.end());
        }

        std::vector<seg_pair_t> vecRange;
        for(auto& item : allMeasure)
            vecRange.insert(vecRange.end(), item.rangeSeg.begin(), item.rangeSeg.end());
        writePCDFile(std::to_string(index) +"-testWall.pcd", pCloud, vecRange);
    }

}