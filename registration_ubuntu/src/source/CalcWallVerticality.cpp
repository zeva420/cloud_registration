#include "CalcWallVerticality.h"

#include "funHelper.h"
#include <pcl/common/common.h>
#include <pcl/filters/crop_box.h>
#include <pcl/filters/passthrough.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>

namespace CloudReg
{
    std::vector<seg_pair_t> calValidHoleVertical(
                        std::vector<std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>>& holeBorders,
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

    void calBox(Eigen::Vector3d midPt, pcl::PointXYZ& min, pcl::PointXYZ& max, int type, int hAixs)
    {
        double boxWidth = 0.025;
        double boxHight = 0.025;
        std::vector<seg_pair_t> calcSeg;
        seg_pair_t seg1, seg2;
        
        seg1.first = midPt;
        seg1.second = midPt;
        seg1.first[hAixs] = midPt[hAixs] - boxWidth / 2;
        seg1.second[hAixs] = midPt[hAixs] + boxWidth / 2;
        calcSeg.emplace_back(seg1);
    
        seg2.first = seg1.first;
        if (type == 1)
            seg2.second = Eigen::Vector3d(seg1.first[0], seg1.first[1], seg1.first[2] + boxHight/2);
        else
            seg2.second = Eigen::Vector3d(seg1.first[0], seg1.first[1], seg1.first[2] - boxHight/2);
        calcSeg.emplace_back(seg2);

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
        // //// get the wall thickness?
        // pcl::PointXYZ min0;
		// pcl::PointXYZ max0;		
		// pcl::getMinMax3D(*pCloud, min0, max0);
        // double thickness = (hAixs == 0) ? std::fabs(max0.y - min0.y) : std::fabs(max0.x - min0.x);

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
        
        auto vecRawPt = getNearestPt(vecCalcPt, pCloud, 0.1*0.1);
        std::vector<double> sumAll;
		for(auto& pt : vecRawPt)
		{
			pcl::PointXYZ min;
		    pcl::PointXYZ max;
            calBox(pt, min, max, type, hAixs);
            auto rangeCloud = filerCloudByRange(pCloud,min,max);
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
        for (size_t i = 0; i < rPoints.size(); ++i)
        {
            for(size_t j = 0; j < rPoints.size(); ++j)
            {
                if (i == j) continue;
                seg_pair_t pp = std::make_pair(rPoints[i], rPoints[j]);
                item.rangeSeg.emplace_back(pp);
            }
        }

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

    bool judgeHoleBorder(std::vector<std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>>& holeBorders,
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
			std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>& wallBorder, 
            std::vector<std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>>& holeBorders,
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

    //test ruler
    /*
    void testRuler(std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>& wallBorder,
                    std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>& holeBorder,
                    pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud, int index)
    {
        if (wallBorder.empty())
        {
            LOG(ERROR) << "empty input for ruler";
        }
        std::cout << "pCloud " << pCloud->points.size() << std::endl;
        std::cout << "wallBorder " << wallBorder.size() << std::endl;
        std::cout << "holeBorder " << holeBorder.size() << std::endl;
        std::vector<double> vecSum = {0,0,0};
		for (std::size_t i = 0; i < wallBorder.size(); i++)
		{
			seg_pair_t line = wallBorder[i];
			vecSum[0] += std::fabs(line.first[0] - line.second[0]);
			vecSum[1] += std::fabs(line.first[1] - line.second[1]);
			vecSum[2] += std::fabs(line.first[2] - line.second[2]);
		}
		std::size_t minIndex = std::min_element(vecSum.begin(), vecSum.end()) - vecSum.begin();

        seg_pair_t baseLine = wallBorder[1];
        Eigen::Vector3d P0 = baseLine.first;
        float theta = 45;
        seg_pair_t ruler;
        
        if(calRuler3d(wallBorder, holeBorder,  baseLine, P0, theta, ruler))
        {
            std::vector<Eigen::Vector3d> rPoints =  createRulerBox(ruler, minIndex, 0.025, 0.025);
            std::vector<seg_pair_t> vecRange;
            vecRange.emplace_back(ruler);
            for (size_t i = 0; i < rPoints.size(); ++i)
            {
                for(size_t j = 0; j < rPoints.size(); ++j)
                {
                    if (i == j) continue;
                    seg_pair_t pp = std::make_pair(rPoints[i], rPoints[j]);
                    vecRange.emplace_back(pp);
                }
            }
            writePCDFile(std::to_string(index) +"testruler.pcd", pCloud, vecRange);
        }
    }*/
    
}