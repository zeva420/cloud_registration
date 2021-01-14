#include "CalcWallFlatness.h"

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
    bool checkAdjDis(std::vector<seg_pair_t> allBorder, int lhAxis, seg_pair_t hole1, seg_pair_t hole2,
                                 double baseHight, Eigen::Vector3d& resP)
    {
        double shortLen = 10000.;
        double longLen = 10000.;
        Eigen::Vector3d sd, ld;
        Eigen::Vector3d holep1 = hole1.first;
        holep1[2] = baseHight;
        Eigen::Vector3d holep2 = hole2.first;
        holep2[2] = baseHight;
       
        for(size_t i = 0; i < allBorder.size(); ++i)
        {
            double h = allBorder[i].first[lhAxis];
            double diff1 = std::fabs(h - holep1[lhAxis]);
            double diff2 = std::fabs(h - holep2[lhAxis]);
            if(h + 1e-6 < holep1[lhAxis] && diff1 < shortLen)
            {
                shortLen = diff1;
                sd = allBorder[i].first;
                sd[2] = baseHight;
            }
            else if(holep2[lhAxis] + 1e-6 < h && diff2 < longLen)
            {
                longLen = diff2;
                ld = allBorder[i].first;
                ld[2] = baseHight;
            }
        }

        if(longLen + 1e-6 < 10000. && longLen > 0.0125) //10cm + 2.5cm
        {
            Eigen::Vector3d longn = (ld - holep2).normalized();
            resP = holep2 + 0.01*longn;   //0.01
            return true;
        }

        if(shortLen + 1e-6 < 10000. && shortLen > 0.0125) //10cm + 2.5cm
        {
            Eigen::Vector3d shortn = (sd - holep1).normalized();
            resP = holep1 + 0.01*shortn;  //0.01
            return true;
        }
        return false;
    }
    bool adjustMeasurePt(Eigen::Vector3d& pt, std::vector<seg_pair_t> holes, 
                            std::vector<seg_pair_t> wallBorder, int lhAxis, double baseHight)
    {
        if (holes.empty())
            return true;

        std::vector<seg_pair_t> allBorder = wallBorder;
        allBorder.insert(allBorder.end(), holes.begin(), holes.end());

        for(size_t i = 0; i < holes.size() / 2; ++i)
        {
            double h1 = holes[2*i].first[lhAxis];
            double h2 = holes[2*i + 1].first[lhAxis];
            if(pt[lhAxis] > h1 && pt[lhAxis] < h2)
            {
                if (!checkAdjDis(allBorder, lhAxis, holes[2*i], holes[2*i + 1], baseHight, pt))
                    return false;
                else
                    return true;
            }
        }
        return true;
    }

    calcMeassurment_t calcSinglePt(Eigen::Vector3d pt, Eigen::Vector3d shortMeasureP, pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud, int lhAxis)
    {
        int thicknessDir = (lhAxis == 0) ? 1 : 0;
        double width = 0.025; //5cm
        double height = 0.025; //5cm

        Eigen::Vector3d p1 = pt;
        p1[2] = p1[2] + height / 2;
        Eigen::Vector3d p2 = pt;
        p2[2] = p2[2] - height / 2;

        std::vector<Eigen::Vector3d> rPoints =  createRulerBox(std::make_pair(p1, p2), thicknessDir, 0.025, width); //short wall
        std::vector<seg_pair_t> pair =  calcBoxSegPair(rPoints);
        std::vector<Eigen::Vector3d> corners = getRulerCorners(rPoints);
        auto rangeCloud = filerCloudByConvexHull(pCloud, corners);

        calcMeassurment_t measure;
        if (rangeCloud->points.empty()) 
        {
            LOG(ERROR) << "filerCloudByRange failed";
            return measure;
        }

        double sum = 0;
        for (auto &p : rangeCloud->points)
        {
            Eigen::Vector3d point(p.x, p.y, p.z);
            sum += std::fabs(point[thicknessDir] - shortMeasureP[thicknessDir]);
        }
        double avg = sum / rangeCloud->points.size();
        measure.value = avg;
        measure.rangeSeg.insert(measure.rangeSeg.end(), pair.begin(), pair.end());
        LOG(INFO) << "one measure value is " << avg;
        return measure;
    }

    void writePCDFileSq(const std::string& name, std::vector<PointCloud::Ptr> pClouds, std::vector<seg_pair_t>& border)
	{
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr pCloudRGB(new pcl::PointCloud<pcl::PointXYZRGB>());

		for(auto& pCloud : pClouds)
        {
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

    std::vector<calcMeassurment_t> calcWallSquareness(std::pair<size_t, size_t> wallPair, std::vector<PointCloud::Ptr> pClouds,
                                        const std::vector<vec_seg_pair_t>& vecWall, const std::vector<vec_seg_pair_t>& holeBorders,
                                        std::vector<std::size_t> lHolesIndex, std::vector<std::size_t> sHolesIndex)
    {
        vec_seg_pair_t shortWall = vecWall[wallPair.first];
        vec_seg_pair_t longWall = vecWall[wallPair.second];
        pcl::PointCloud<pcl::PointXYZ>::Ptr spCloud = pClouds[wallPair.first];
        pcl::PointCloud<pcl::PointXYZ>::Ptr lpCloud = pClouds[wallPair.second];
        auto seg1 = shortWall.back();
        auto seg2 = longWall.back();
        const auto shSeg = seg1.first - seg1.second;
        const auto lhSeg = seg2.first - seg2.second;
    
        Eigen::Vector3d lfirstV = longWall[0].first - longWall[0].second;
		Eigen::Vector3d lnV = lhSeg.cross(lfirstV);
		lnV = -1 * lnV.normalized(); //Normal vector towards the inside of the rectangle

        Eigen::Vector3d sfirstV = shortWall[0].first - shortWall[0].second;
		Eigen::Vector3d snV = shSeg.cross(sfirstV);
		snV = -1 * snV.normalized();

        Eigen::Vector3d connPt;
        if(seg1.first == seg2.first || seg1.first == seg2.second)
            connPt = seg1.first;
        else
            connPt = seg1.second;

        Eigen::Vector3d sother, lother;
        sother = (seg1.first == connPt) ? seg1.second : seg1.first;
        lother = (seg2.first == connPt) ? seg2.second : seg2.first;
        Eigen::Vector3d shortn = (sother - connPt).normalized();
        Eigen::Vector3d longn = (lother - connPt).normalized();

        Eigen::Vector3d shortMeasureP = connPt + 0.03 * shortn;  //30cm
        double lbaseHeight = seg2.first[2] + 0.03; //30cm
        double sbaseHeight = seg1.first[2] + 0.03; //30cm
        shortMeasureP[2] = sbaseHeight;

        std::vector<vec_seg_pair_t> lHoles, sHoles;
        std::vector<seg_pair_t> lHoleBorder, sHoleBorder;
        int lhAxis = calWallHorizontalAxis(seg2);
        int shAxis = calWallHorizontalAxis(seg1);
        std::vector<seg_pair_t> lvecWallHorizen, lvecWallVertical;
        std::vector<seg_pair_t> svecWallHorizen, svecWallVertical;
        groupDirection(lhSeg, longWall, lvecWallVertical, lvecWallHorizen);
        groupDirection(shSeg, shortWall, svecWallVertical, svecWallHorizen);
        for(auto index : lHolesIndex) 
            lHoles.emplace_back(holeBorders[index]);
        if (!lHoles.empty())
            lHoleBorder = calValidHoleVertical(lHoles, std::make_pair(seg2.first, seg2.second), lhAxis);
        for(auto index : sHolesIndex) 
            sHoles.emplace_back(holeBorders[index]);
        if (!sHoles.empty())
            sHoleBorder = calValidHoleVertical(sHoles, std::make_pair(seg1.first, seg1.second), shAxis);

        std::vector<calcMeassurment_t> allMeasure;
        std::vector<double> sumAll;
        if (!adjustMeasurePt(shortMeasureP, sHoleBorder, svecWallVertical, shAxis, sbaseHeight))
        {
            LOG(ERROR) << "can not find measure point in short wall";
            return allMeasure;
        }

        double lLen = (seg2.first - seg2.second).norm();
        Eigen::Vector3d p1 = connPt + 0.05 * longn;
        Eigen::Vector3d p2 = connPt + (lLen / 2) * longn;
        Eigen::Vector3d p3 = connPt + (lLen - 0.05) * longn;

        auto getResult = [&]( Eigen::Vector3d& p)
        { 
            p[2] = p[2] + lbaseHeight;
            if (adjustMeasurePt(p, lHoleBorder, lvecWallVertical, lhAxis, lbaseHeight))
            {
                calcMeassurment_t measure = calcSinglePt(p, shortMeasureP, lpCloud, lhAxis);
                if (!measure.rangeSeg.empty())
                {
                    Eigen::Vector3d po = p + (shortMeasureP[shAxis]-connPt[shAxis]) * lnV; //30cm
                    measure.rangeSeg.emplace_back(std::make_pair(p, po));
                    sumAll.emplace_back(measure.value);
                    allMeasure.emplace_back(measure);
                }
                    
            }
        };

        getResult(p1);
        getResult(p2);
        getResult(p3);

        if (!sumAll.empty())
        {
            double max =  *std::max_element(sumAll.begin(),sumAll.end());
            double min = *std::min_element(sumAll.begin(),sumAll.end());
            double difference = std::fabs(max - min);
            LOG(INFO) << "max avg: " << max << " min avg: " << min
                << " difference: " << difference;
            calcMeassurment_t measure;
            measure.value = difference;
            Eigen::Vector3d po = shortMeasureP + lLen * snV;
            measure.rangeSeg.emplace_back(std::make_pair(shortMeasureP, po));
            allMeasure.emplace_back(measure);
        }

        std::vector<PointCloud::Ptr> tClouds;
        std::vector<seg_pair_t> allSegs;
        for(auto& m : allMeasure)
            allSegs.insert(allSegs.end(), m.rangeSeg.begin(), m.rangeSeg.end());
        tClouds.emplace_back(spCloud);
        tClouds.emplace_back(lpCloud);
        std::cout << std::to_string(wallPair.first)+ "testSq.pcd" << std::endl;
        writePCDFileSq(std::to_string(wallPair.first)+ "testSq.pcd", tClouds , allSegs);

        return allMeasure;
    }

    void calcSquareness(const std::vector<vec_seg_pair_t>& vecWall, const std::vector<vec_seg_pair_t>& holeBorders,
			std::vector<PointCloud::Ptr> pClouds, std::map<std::size_t, std::vector<std::size_t>> holeMap, const double calcLengthTh)
    {
        if (vecWall.size() < 2)
        {
            LOG(ERROR) << "wall num is small " << vecWall.size();
            return;
        }
        LOG(INFO) << "input wall num is " << vecWall.size();
        
        //step1: Team up for the wall
        std::vector<std::pair<size_t, size_t>> wallPairs;
        for(size_t i = 0; i < vecWall.size() / 2; ++i)
        {
            auto seg1 = vecWall[2*i].back();
            auto seg2 = vecWall[2*i+1].back();
            double length1 = (seg1.first - seg1.second).norm();
            double length2 = (seg2.first - seg2.second).norm();

            if(length1 < calcLengthTh || length2 < calcLengthTh)
                continue;

            if (length1 < length2)
                wallPairs.emplace_back(std::make_pair(2*i, 2*i+1));
            else
                wallPairs.emplace_back(std::make_pair(2*i + 1, 2*i));
        }

        if(wallPairs.empty())
        {
            LOG(INFO) << "no valid wall pair";
            return;
        }
        LOG(INFO) << "valid wallPairs " << wallPairs.size();

        //step2: calc
        std::vector<calcMeassurment_t> allMeasure;
        for (size_t i = 0; i < wallPairs.size(); ++i)
        {
            auto wallPair = wallPairs[i];
            LOG(INFO) << "Handle wall " << wallPair.first << "---" << wallPair.second;
            std::vector<std::size_t> lHolesIndex, sHolesIndex;
            if (holeMap.count(wallPair.second))
                lHolesIndex = holeMap[wallPair.second];
            if (holeMap.count(wallPair.first))
                sHolesIndex = holeMap[wallPair.first];
            auto oneMeasure =  calcWallSquareness(wallPair, pClouds, vecWall, holeBorders, lHolesIndex, sHolesIndex);
            if (!oneMeasure.empty())
                allMeasure.insert(allMeasure.end(), oneMeasure.begin(), oneMeasure.end());
        }
    }
}