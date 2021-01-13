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
    int calWallHorizontalAxis(const seg_pair_t& seg)
    {
        double length = (seg.first - seg.second).norm();
        if (std::fabs(length - std::fabs(seg.first[0] - seg.second[0])) < 1e-4)
            return 0;
        else if (std::fabs(length - std::fabs(seg.first[1] - seg.second[1])) < 1e-4)
            return 1;
        return 2;
    }

    std::vector<seg_pair_t> calcSingleWall(const std::vector<seg_pair_t>& wallBorder, const std::vector<vec_seg_pair_t>& holeBorders,
                                    int hAxis, std::pair<seg_pair_t, seg_pair_t> wallData, int type)
    {
        auto seg1 = wallData.first;
        auto seg2 = wallData.second;
        std::vector<seg_pair_t> rulers;

        std::vector<seg_pair_t> hborder;
        for(auto b : holeBorders)
            hborder.insert(hborder.end(), b.begin(), b.end());
        auto getDiagonalRuler = [&](int index, std::vector<seg_pair_t>& rulers)
        {
            seg_pair_t ruler;
            const auto& seg = wallBorder[index];
            
            if(calRuler3d(wallBorder, hborder, seg, seg.first, 45, ruler))
            {
                Eigen::Vector3d rulern = (ruler.second - ruler.first).normalized();
                ruler.first = ruler.first + 0.05*rulern; // Move 50mm
                ruler.second = ruler.second - 0.05*rulern; // Move 50mm
                cutOffRuler(ruler, 2);
                rulers.emplace_back(ruler);
            }
        };

        double length = std::fabs(seg1.first[hAxis] - seg2.first[hAxis]);
        if (type == 1 || type == 3) //Upper left ruler
        {
            int index = 1;
            getDiagonalRuler(index, rulers);
        }

        if (type == 2 || type == 3) //bottom right ruler
        {
            int index = wallBorder.size()-1;
            getDiagonalRuler(index, rulers);
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
            rulers.emplace_back(std::make_pair(p1, p2));
        }

        LOG(INFO) << "wall get rulers num: " << rulers.size();
        return rulers;
    }

    std::vector<seg_pair_t> calcWallRuler(const vec_seg_pair_t& wallBorder, const std::vector<vec_seg_pair_t>& holeBorders,
                                    int hAxis, std::vector<seg_pair_t> calValidVertical)
    {
        std::size_t wallNum = calValidVertical.size() / 2;
        LOG(INFO) << "valid Wall num is " << wallNum;
        std::vector<seg_pair_t> rulers;
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
                std::vector<seg_pair_t> ruler;
                int type = 3;
                if (wallNum > 1 && i == 0)
                    type = 1;
                else if (wallNum > 1 && i != 0)
                    type = 2;

                if (symmetric && type == 1)
                    type = 2;
                else if (symmetric && type == 2)
                    type = 1;
               
                ruler = calcSingleWall(wallBorder, holeBorders, hAxis, std::make_pair(seg1, seg2), type);
                rulers.insert(rulers.end(), ruler.begin(), ruler.end());
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

    std::vector<seg_pair_t> calcHoleRuler(const seg_pair_t& horizenSeg, 
                            const std::vector<seg_pair_t>& holeBorder, int hAxis, int type)
    {
        std::vector<seg_pair_t> rulers;
        std::vector<seg_pair_t> vecHoleHorizen, vecHoleVertical;
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
                rulers.emplace_back(ruler1);
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
                rulers.emplace_back(ruler2);
            }
        }

        return rulers;
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

    std::vector<seg_pair_t> transfHoleRuler(std::vector<seg_pair_t> rulersIn, std::vector<seg_pair_t> checkBorder)
    {
        std::vector<seg_pair_t> rulersOut;
        for(auto ruler : rulersIn)
        {
            double length = (ruler.first - ruler.second).norm();
            if (length + 1e-4 > 2)
            {
                LOG(ERROR) << "ruler length " << length << " can not transf";
                continue;
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
        }
        return rulersOut;
    }

    void calcWallFlatness(const vec_seg_pair_t& wallBorder, const std::vector<vec_seg_pair_t>& holeBorders,
			pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud, Eigen::Vector4d cadPlane, int index)
    {
        // step0: cal Horizontal and thickness Axis
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
			return;
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
        auto wallRulers =  calcWallRuler(wallBorder, holeBorders, hAxis, calValidVertical);
        
        //step3: Handle holes
        std::vector<seg_pair_t> allEmptRulers;
        std::vector<seg_pair_t> allHoleRulers;
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
               
                std::vector<seg_pair_t> ruler = calcHoleRuler(horizenSeg, holeBorders[i], hAxis, type);
                allEmptRulers.insert(allEmptRulers.end(), ruler.begin(), ruler.end());
                
                std::vector<seg_pair_t> checkBorder;
                checkBorder = wallBorder;
                for(size_t j = 0; j < holeBorders.size(); ++j)
                {
                    if (i == j) continue;
                    checkBorder.insert(checkBorder.end(), holeBorders[j].begin(), holeBorders[j].end());
                }
                auto realHoleRulers = transfHoleRuler(ruler, checkBorder);
                allHoleRulers.insert(allHoleRulers.end(), realHoleRulers.begin(), realHoleRulers.end());
            }

        }

        std::vector<calcMeassurment_t> allMeasure;
        for (auto &ruler : wallRulers)
        {
            calcMeassurment_t measure;
            measure = calFlatness(ruler, thicknessDir, cadPlane, pCloud);
            allMeasure.emplace_back(measure);
        }
        for (auto &ruler : allHoleRulers)
        {
            calcMeassurment_t measure;
            measure = calFlatness(ruler, thicknessDir, cadPlane, pCloud);
            allMeasure.emplace_back(measure);
        }

        std::vector<seg_pair_t> vecRange;
        for(auto& item : allMeasure)
        {
            for(auto& ruler : item.rangeSeg)
            {
                std::vector<Eigen::Vector3d> rPoints =  createRulerBox(ruler, 2, 0.025, 0.025);
                std::vector<seg_pair_t> pair =  calcBoxSegPair(rPoints);
                vecRange.insert(vecRange.end(), pair.begin(), pair.end());
            }
        }
        writePCDFile(std::to_string(index) +"-testWallFl.pcd", pCloud, vecRange);
    }
}