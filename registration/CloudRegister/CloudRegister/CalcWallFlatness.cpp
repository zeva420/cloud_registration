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
    std::map<int,std::vector<seg_pair_t>> calcSingleWall(const std::vector<seg_pair_t>& wallBorder, const std::vector<vec_seg_pair_t>& holeBorders,
                                    int hAxis, std::pair<seg_pair_t, seg_pair_t> wallData, int type)
    {
        auto seg1 = wallData.first;
        auto seg2 = wallData.second;
        std::map<int,std::vector<seg_pair_t>> rulerMap;

        std::vector<seg_pair_t> hborder;
        for(auto b : holeBorders)
            hborder.insert(hborder.end(), b.begin(), b.end());
        auto getDiagonalRuler = [&](int index, std::map<int,std::vector<seg_pair_t>>& rulerMap, int type)
        {
            seg_pair_t ruler;
            const auto& seg = wallBorder[index];
            
            if(calRuler3d(wallBorder, hborder, seg, seg.first, 45, ruler))
            {
                Eigen::Vector3d rulern = (ruler.second - ruler.first).normalized();
                ruler.first = ruler.first + 0.05*rulern; // Move 50mm
                ruler.second = ruler.second - 0.05*rulern; // Move 50mm
                cutOffRuler(ruler, 2);
                rulerMap[type].emplace_back(ruler);
            }
        };

        double length = std::fabs(seg1.first[hAxis] - seg2.first[hAxis]);
        if (type == 1 || type == 3) //Upper left ruler
        {
            int index = 1;
            getDiagonalRuler(index, rulerMap, 1);
        }

        if (type == 2 || type == 3) //bottom right ruler
        {
            int index = wallBorder.size()-1;
            getDiagonalRuler(index, rulerMap, 2);
        }

        if (length + 1e-4 >= 3)  //middle ruler
        {
            double len1 = (seg1.first - seg1.second).norm();
            double len2 = (seg2.first - seg2.second).norm();
            double height = (len1 > len2) ? len1 : len2;
            Eigen::Vector3d hypRuler0 = seg1.first;
            double baseHeight = (seg1.first[2] < seg1.second[2]) ? seg1.first[2] : seg1.second[2];
            hypRuler0[2] = baseHeight + height / 2;
            Eigen::Vector3d hypRuler1 = seg2.first;
            hypRuler1[2] = baseHeight + height / 2;
            Eigen::Vector3d hypRuler = hypRuler1 - hypRuler0;
            Eigen::Vector3d rulern = hypRuler.normalized();
            
            double resDis = (length - 2) / 2;
            Eigen::Vector3d p1 = hypRuler0;
            p1[hAxis] = p1[hAxis] + resDis;
            Eigen::Vector3d p2 = p1 + 2*rulern;    //200 cm
            rulerMap[3].emplace_back(std::make_pair(p1, p2));
        }

        return rulerMap;
    }

    std::map<int,std::vector<seg_pair_t>> calcWallRuler(const vec_seg_pair_t& wallBorder, const std::vector<vec_seg_pair_t>& holeBorders,
                                    int hAxis, std::vector<seg_pair_t> calValidVertical)
    {
        std::size_t wallNum = calValidVertical.size() / 2;
        LOG(INFO) << "valid Wall num is " << wallNum;
        std::map<int,std::vector<seg_pair_t>> rulers;
        bool symmetric = (wallBorder.back().first[hAxis] > wallBorder.back().second[hAxis]) ? 0 : 1;
        for (std::size_t i = 0; i < wallNum; ++i)
        {
            auto seg1 = calValidVertical[2*i];
            auto seg2 = calValidVertical[2*i + 1];
            double length = std::fabs(seg1.first[hAxis] - seg2.first[hAxis]);
            if (length + 1e-6 < 0.6)
            {
                LOG(INFO) << "wall length is small " << length;
                continue;
            }

            if (length >= 0.6)
            {
                int type = 3;
                if (wallNum > 1 && i == 0)
                    type = 1;
                else if (wallNum > 1 && i != 0)
                    type = 2;

                if (symmetric && type == 1)
                    type = 2;
                else if (symmetric && type == 2)
                    type = 1;
               
                auto ruler = calcSingleWall(wallBorder, holeBorders, hAxis, std::make_pair(seg1, seg2), type);
                for (auto &map : ruler)
                    rulers[map.first].insert(rulers[map.first].end(), map.second.begin(), map.second.end());
            }
        }
        return rulers;
    }

    std::vector<seg_pair_t> calHoleZoom(seg_pair_t seg1, seg_pair_t seg2, int hAxis)
    {
        Eigen::Vector3d p1 = (seg1.first[2] < seg1.second[2]) ? seg1.first : seg1.second;
        Eigen::Vector3d p2 = (p1 == seg1.first) ? seg1.second : seg1.first;
        Eigen::Vector3d p5 = (seg2.first[2] < seg2.second[2]) ? seg2.first : seg2.second;
        Eigen::Vector3d p6 = (p5 == seg2.first) ? seg2.second : seg2.first;
        Eigen::Vector3d p3 = p1;
        p3[hAxis] = (p1[hAxis] + p5[hAxis]) / 2;
        Eigen::Vector3d p4 = p2;
        p4[hAxis] = (p2[hAxis] + p6[hAxis]) / 2;

        seg_pair_t sega = std::make_pair(p1, p2);
        seg_pair_t segb = std::make_pair(p2, p4);
        seg_pair_t segc = std::make_pair(p4, p3);
        seg_pair_t segd = std::make_pair(p3, p1);

        seg_pair_t sege = std::make_pair(p3, p4);
        seg_pair_t segf = std::make_pair(p4, p6);
        seg_pair_t segg = std::make_pair(p6, p5);
        seg_pair_t segh = std::make_pair(p5, p3);

        std::vector<seg_pair_t> zoom = {sega, segb, segc, segd, sege, segf, segg, segh};
        return zoom;
    }

    double adjustHRulerBorder(seg_pair_t seg, std::vector<seg_pair_t> checkBorder, double resDis)
    {
        if (checkBorder.empty())
            return resDis;

        Eigen::Vector3d P0 = seg.second;
        Eigen::Vector3d new3d;
		double rulerLength = 1000000.0;
        for (std::size_t i = 0; i < checkBorder.size(); i++)
		{
			seg_pair_t line = checkBorder[i];
            Eigen::Vector3d intersec;
            if (calIntersection(seg, line, intersec))
            {
                double length = (intersec - P0).norm();
                if (length < rulerLength)
                {
                    rulerLength = length;
                    new3d = intersec;
                }
            }
        }
        if (std::fabs(1000000.0 - rulerLength) < 1)
		{
			LOG(ERROR) << "can not adjust hole ruler";
			return resDis;
		}

        double lengthr = rulerLength - 0.05; //50mm;
        if (lengthr <= 0)
            return 0;
        
        return (lengthr < resDis) ? lengthr : resDis;
    }

    std::vector<seg_pair_t> transfHoleRuler(seg_pair_t ruler, std::vector<seg_pair_t> checkBorder)
    {
        std::vector<seg_pair_t> rulersOut;
        
        double length = (ruler.first - ruler.second).norm();
        if (length + 1e-4 > 2)
        {
            LOG(ERROR) << "ruler length " << length << " can not transf";
            return rulersOut;
        }

        double resDis = (2 - length) / 2;
        seg_pair_t newRuler1, newRuler2;
        Eigen::Vector3d rulern = (ruler.second - ruler.first).normalized();
        double leftRes = adjustHRulerBorder(std::make_pair(ruler.second, ruler.first), checkBorder, resDis);
        Eigen::Vector3d p1 = ruler.first - leftRes * rulern;

        if ((resDis - leftRes) > 1e-6)
            resDis += (resDis - leftRes);
        double rightRes = adjustHRulerBorder(ruler, checkBorder, resDis);
        Eigen::Vector3d p2 = ruler.second + rightRes * rulern;
        
        newRuler1 = std::make_pair(p1, ruler.first);
        newRuler2 = std::make_pair(ruler.second, p2);
        rulersOut.emplace_back(newRuler1);
        rulersOut.emplace_back(newRuler2);
        
        return rulersOut;
    }

    std::map<int,std::vector<seg_pair_t>> calcHoleRuler(const seg_pair_t& horizenSeg, std::vector<seg_pair_t> checkBorder,
                            const std::vector<seg_pair_t>& holeBorder, int hAxis, int type)
    {
        std::map<int,std::vector<seg_pair_t>> rulers;
        std::vector<seg_pair_t> vecHoleHorizen, vecHoleVertical;
        bool symmetric = (horizenSeg.first[hAxis] > horizenSeg.second[hAxis]) ? 0 : 1;
        const auto hSeg = horizenSeg.first - horizenSeg.second;
        groupDirection(hSeg, holeBorder, vecHoleVertical, vecHoleHorizen);
        if (vecHoleHorizen.size() < 2 || vecHoleVertical.size() < 2)
		{
			LOG(ERROR) << "group hole Direction Failed: " << vecHoleHorizen.size() << " -- " << vecHoleVertical.size();
			return rulers;
		}
        std::sort(vecHoleVertical.begin(), vecHoleVertical.end(), [&](const seg_pair_t& left, const seg_pair_t& right){
				return left.first[hAxis] < right.first[hAxis];});
        
        std::vector<seg_pair_t> holeZoom = calHoleZoom(vecHoleVertical.front(), vecHoleVertical.back(), hAxis);
        if (type == 1 || type == 3)  //left ruler
        {
            seg_pair_t ruler1;
            std::vector<seg_pair_t> holeBorder;
            std::vector<seg_pair_t> zoomBorder = {holeZoom[0], holeZoom[1], holeZoom[2], holeZoom[3]};
            if(calRuler3d(zoomBorder, holeBorder, zoomBorder[2], zoomBorder[2].first, 45, ruler1))
            {
                Eigen::Vector3d rulern1 = (ruler1.second - ruler1.first).normalized();
                ruler1.first = ruler1.first - 0.03*rulern1; // Move up 30mm
                ruler1.second = ruler1.second + 0.03*rulern1; // Move down 30mm
                std::vector<seg_pair_t> rulersOut = transfHoleRuler(ruler1, checkBorder);
                int rulerNum = (symmetric == 1) ? 5 : 4;
                if (!rulersOut.empty())
                    rulers[rulerNum].insert(rulers[rulerNum].end(), rulersOut.begin(), rulersOut.end());
            }
        }

        if(type == 2 || type == 3) //right ruler
        {
            seg_pair_t ruler2;
            std::vector<seg_pair_t> holeBorder;
            std::vector<seg_pair_t> zoomBorder = {holeZoom[4], holeZoom[5], holeZoom[6], holeZoom[7]};
            if (calRuler3d(zoomBorder, holeBorder, zoomBorder[1], zoomBorder[1].first, 45, ruler2))
            {
                Eigen::Vector3d rulern2 = (ruler2.second - ruler2.first).normalized();
                ruler2.first = ruler2.first - 0.03*rulern2; // Move up 30mm
                ruler2.second = ruler2.second + 0.03*rulern2; // Move down 30mm
                std::vector<seg_pair_t> rulersOut = transfHoleRuler(ruler2, checkBorder);
                int rulerNum = (symmetric == 1) ? 4 : 5;
                if (!rulersOut.empty())
                    rulers[rulerNum].insert(rulers[rulerNum].end(), rulersOut.begin(), rulersOut.end());
            }
        }

        return rulers;
    }

    std::tuple<std::vector<calcMeassurment_t>, std::vector<seg_pair_t>>
    calWallFlatness(const vec_seg_pair_t& wallBorder, const std::vector<vec_seg_pair_t>& holeBorders,
			pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud, Eigen::Vector4d cadPlane, int wallIndex)
    {
        // step0: cal Horizontal and thickness Axis
        std::vector<calcMeassurment_t> allMeasure;
        std::vector<seg_pair_t> returnSeg;
        const std::pair<Eigen::Vector3d, Eigen::Vector3d>& horizenSeg = wallBorder.back();
        int hAxis = calWallHorizontalAxis(horizenSeg);
        int thicknessDir = (hAxis == 0) ? 1 : 0;
        LOG(INFO) << "Horizontal Axis is " << hAxis;
        LOG(INFO) << "thicknessDir Axis is " << thicknessDir;

        //step1: Separate walls according to doors and windows
        std::vector<seg_pair_t> vecWallHorizen, vecWallVertical;
        const auto hSeg = horizenSeg.first - horizenSeg.second;
        groupDirection(hSeg, wallBorder, vecWallVertical, vecWallHorizen);
        if (vecWallHorizen.size() < 2 || vecWallVertical.size() < 2)
		{
			LOG(ERROR) << "group Wall Direction Failed: " << vecWallHorizen.size() << " -- " << vecWallVertical.size();
			return std::make_tuple(allMeasure, returnSeg);
		}
        std::sort(vecWallVertical.begin(), vecWallVertical.end(), [&](const seg_pair_t& left, const seg_pair_t& right){
				return left.first[hAxis] < right.first[hAxis];});

        std::vector<seg_pair_t> calValidVertical, allVertical;
        auto validHoleVertical = calValidHoleVertical(holeBorders, horizenSeg, hAxis);
       
        if (!validHoleVertical.empty())
        {
            std::sort(validHoleVertical.begin(), validHoleVertical.end(), [&](const seg_pair_t& left, const seg_pair_t& right){
				return left.first[hAxis] < right.first[hAxis];});

            calValidVertical.emplace_back(validHoleVertical.front());
            calValidVertical.emplace_back(validHoleVertical.back());
            allVertical.insert(allVertical.end(), validHoleVertical.begin(), validHoleVertical.end());
        }
        calValidVertical.insert(calValidVertical.begin(), vecWallVertical[0]); //for wall
        calValidVertical.emplace_back(vecWallVertical.back());

        allVertical.insert(allVertical.begin(), vecWallVertical[0]);  //wall and hole
        allVertical.emplace_back(vecWallVertical.back());

        //step2: Handle walls
        std::map<int,std::vector<seg_pair_t>> wallRulers;
        wallRulers =  calcWallRuler(wallBorder, holeBorders, hAxis, calValidVertical);
        
        //step3: Handle holes
        std::map<int,std::vector<seg_pair_t>> holeRulers;
        if (!holeBorders.empty())
        {
            std::vector<double> lengths;
            for (std::size_t i = 0; i < allVertical.size() - 1; ++i)
            {
                const auto &Vertical1 = allVertical[i];
                const auto &Vertical2 = allVertical[i + 1];

                if (judgeHoleBorder(holeBorders, std::make_pair(Vertical1, Vertical2)))
                    continue;
                double length = std::fabs(Vertical1.first[hAxis] - Vertical2.first[hAxis]);
                lengths.emplace_back(length);
            }

            for(size_t i = 0; i < holeBorders.size(); ++i)
            {
                double leftLen = lengths[i];
                double rightLen = lengths[i + 1];
                if(leftLen < 0.4 + 1e-4 && rightLen < 0.4 + 1e-4) //The wall next to the hole is short
                    continue;

                int type = 3;
                if(leftLen >= 0.4 + 1e-4 && rightLen < 0.4 + 1e-4)
                    type = 1;
                if(rightLen >= 0.4 + 1e-4 && leftLen < 0.4 + 1e-4)
                    type = 2;

                std::vector<seg_pair_t> checkBorder;
                checkBorder = wallBorder;
                for(size_t j = 0; j < holeBorders.size(); ++j)
                {
                    if (i == j) continue;
                    checkBorder.insert(checkBorder.end(), holeBorders[j].begin(), holeBorders[j].end());
                }
               
                auto rulerMap = calcHoleRuler(horizenSeg, checkBorder, holeBorders[i], hAxis, type);
                for (auto &map : rulerMap)
                    holeRulers[map.first].insert(holeRulers[map.first].end(), map.second.begin(), map.second.end());
            }

        }

        std::vector<int> outOrder = {1, 2, 3};
        std::vector<seg_pair_t> vecRange;
        for(auto& type : outOrder)
        {
            if (!wallRulers.count(type))
                continue;
            auto rulers = wallRulers[type];
            for (auto& ruler : rulers)
            {
                std::vector<seg_pair_t> tmp = {ruler};
                calcMeassurment_t measure = calFlatness(tmp, thicknessDir, cadPlane, pCloud);
                if (!measure.rangeSeg.empty())
                {
                    std::vector<Eigen::Vector3d> rPoints =  createRulerBox(ruler, thicknessDir, 0.025, 0.025);
                    std::vector<seg_pair_t> pair =  calcBoxSegPair(rPoints);
                    vecRange.insert(vecRange.end(), pair.begin(), pair.end());
                    LOG(INFO) << "Wall "<<wallIndex<<" flatness: " << measure.value;
                    allMeasure.emplace_back(measure);
                }
            }
        }

        std::vector<int> outOrder1 = {4, 5};
        for(auto& type : outOrder1)
        {
            if (!holeRulers.count(type))
                continue;
            auto rulers = holeRulers[type];
            for(size_t i = 0; i < rulers.size() / 2; ++i)
            {
                auto ruler1 = rulers[2*i];
                auto ruler2 = rulers[2*i + 1];
                std::vector<seg_pair_t> tmp = {ruler1, ruler2};
                calcMeassurment_t measure = calFlatness(tmp, thicknessDir, cadPlane, pCloud);
                if (!measure.rangeSeg.empty())
                {
                    std::vector<Eigen::Vector3d> rPoints =  createRulerBox(ruler1, thicknessDir, 0.025, 0.025);
                    std::vector<seg_pair_t> pair =  calcBoxSegPair(rPoints);
                    vecRange.insert(vecRange.end(), pair.begin(), pair.end());
                    rPoints =  createRulerBox(ruler2, thicknessDir, 0.025, 0.025);
                    pair =  calcBoxSegPair(rPoints);
                    vecRange.insert(vecRange.end(), pair.begin(), pair.end());
                    LOG(INFO) << "Wall "<<wallIndex<<" flatness: " << measure.value;
                    allMeasure.emplace_back(measure);
                }
            }
        }
#ifdef VISUALIZATION_ENABLED
        pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud_filtered(new pcl::PointCloud<pcl::PointXYZ>());
        uniformSampling(0.01, pCloud, pCloud_filtered);
        writePCDFile("WallFlatness-" + std::to_string(wallIndex) + ".pcd", pCloud_filtered, vecRange);
#endif
        return std::make_tuple(allMeasure, returnSeg);
    }
}