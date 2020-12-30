#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkRenderingOpenGL)
#include <pcl/visualization/pcl_visualizer.h>

#include "Segmentation.h"
#include "glog/logging.h"
#include "funHelper.h"

#include <pcl/common/transforms.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/features/normal_3d.h>
#include <pcl/segmentation/region_growing.h>

#include <random>


namespace CloudReg
{
bool Segmentation::run(pcl::PointCloud<pcl::PointXYZ>::Ptr origCloud,
			                CADModel &cadModel)
{   
    pcl::PointCloud<pcl::PointXYZ>::Ptr samplingCloud(new pcl::PointCloud<pcl::PointXYZ>());
    downSampling(0.05, origCloud, samplingCloud);

    Eigen::Matrix4d T(Eigen::Matrix4d::Identity());
    calibrateDirectionToAxisZ(samplingCloud, cadModel, T);

    pcl::PointCloud<pcl::PointXYZ>::Ptr transformed_cloud(new pcl::PointCloud<pcl::PointXYZ>());
    pcl::transformPointCloud(*origCloud, *transformed_cloud, T);
    origCloud->swap(*transformed_cloud);

    Eigen::Matrix4d bestT = Eigen::Matrix4d::Identity();
    calibrateCloudToCadPose(samplingCloud, cadModel, bestT);

    transformed_cloud->clear();
    pcl::transformPointCloud(*origCloud, *transformed_cloud, bestT);
    origCloud->swap(*transformed_cloud);

    doSegmentation(samplingCloud, cadModel);

    return true;
}

bool Segmentation::downSampling(double radius, pcl::PointCloud<pcl::PointXYZ>::Ptr cloud,
                        pcl::PointCloud<pcl::PointXYZ>::Ptr samplingCloud)
{
    LOG(INFO) << "********downSampling*******";
    uniformSampling(radius, cloud, samplingCloud);
    LOG(INFO) << "uniformSampling, radius:" << radius
        << " cloud size before:" << cloud->size()
        << " after: " << samplingCloud->size();
	return true;
}

bool Segmentation::statisticsForPointZ(float binSizeTh,
                pcl::PointCloud<pcl::PointXYZ>::Ptr cloud,
                std::vector<std::pair<int, std::vector<int>>> &zToNumVec)
{
    struct comp
    {
        bool operator()(const std::pair<int, std::vector<int>> l, const std::pair<int, std::vector<int>> r)
        {
            return (l.second.size() > r.second.size());
        }
    };
    std::map<int, std::vector<int>> zToNumMap;
    for (int i = 0; i < cloud->size(); i++)
    {
        auto &p = cloud->points[i];
        int z = std::floor(p.z / binSizeTh);
        zToNumMap[z].push_back(i);
    }
    if (zToNumMap.empty())
    {
        return false;
    }

    zToNumVec.clear();
    zToNumVec.insert(zToNumVec.end(), zToNumMap.begin(), zToNumMap.end());
    std::sort(zToNumVec.begin(), zToNumVec.end(), comp());
    return true;
}

bool Segmentation::calibrateDirectionToAxisZ(
            pcl::PointCloud<pcl::PointXYZ>::Ptr inputCloud, CADModel &cad, Eigen::Matrix4d &T)
{
    LOG(INFO) << "********calibrateDirectionToAxisZ*******";
    float binSizeTh = 0.2;
    std::vector<std::pair<int, std::vector<int>>> zToNumVec;
    statisticsForPointZ(binSizeTh, inputCloud, zToNumVec);

    auto it = zToNumVec.begin();
    auto &vecIdxs = it->second;
    auto subSet = geo::getSubSet(inputCloud, vecIdxs, false);
    LOG(INFO) << "zToNumVec:" << zToNumVec.size() << " z:" << it->first << " subSet:" << subSet->size();
    for (auto &it : zToNumVec)
    {
        LOG(INFO) << "z:" << it.first << " num:" << it.second.size();
    }
    Eigen::VectorXf coeff;
    std::vector<int> inlierIdxs;
    planeFitting(0.1, subSet, coeff, inlierIdxs);  
    Eigen::Vector4d plane(coeff(0), coeff(1), coeff(2), coeff(3));
    
    Eigen::Vector3d axisZ(0,0,1);
    double flag = (plane.block<3,1>(0,0).dot(axisZ) > 0) ? 1.0 : -1.0;  
    plane = flag * plane; 
    Eigen::Vector3d n = plane.block<3,1>(0,0); 
    double angle = std::acos(n.dot(axisZ));
    Eigen::AngleAxisd rotation_vector(angle, n.cross(axisZ));
    Eigen::Matrix3d R = rotation_vector.matrix();
    LOG(INFO) << "angle:" << (angle * 180.0 / 3.14) << " n:" << n(0) << "," << n(1) << "," << n(2) << " d:" << plane(3);

    T.block<3,3>(0,0) = R;

    pcl::PointCloud<pcl::PointXYZ>::Ptr calibratedCloud(new pcl::PointCloud<pcl::PointXYZ>());
    pcl::transformPointCloud(*inputCloud, *calibratedCloud, T);
    inputCloud->swap(*calibratedCloud);

    // remove some outliers by z
    auto get_z_range = [](const Eigen::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> &segments)
            ->std::pair<double, double> {
        double minZ = 2000;
        double maxZ = 0;
        for (auto &seg : segments)
        {
            minZ = std::min(minZ, seg.first(2));
            minZ = std::min(minZ, seg.second(2));
            maxZ = std::max(maxZ, seg.first(2));
            maxZ = std::max(maxZ, seg.second(2));
        }
        return std::make_pair(minZ, maxZ);
    };
    auto Bottons = cad.getTypedModelItems(ITEM_BOTTOM_E);
    auto Tops = cad.getTypedModelItems(ITEM_TOP_E);
    double minZ = 2000;
    double maxZ = 0;
    for (auto &it : Bottons)
    {
        auto highRange = get_z_range(it.segments_);
        minZ = std::min(minZ, highRange.first);
    }
    for (auto &it : Tops)
    {
        auto highRange = get_z_range(it.segments_);
        maxZ = std::max(maxZ, highRange.second);
    }

    double height = maxZ - minZ;
    LOG(INFO) << "minZ:" << minZ << " maxZ:" << maxZ << " height:" << height;
    binSizeTh = 0.1;
    zToNumVec.clear();
    statisticsForPointZ(binSizeTh, inputCloud, zToNumVec);
    std::vector<int> outliers;
    for (auto &it : zToNumVec)
    {
        auto &tmpVec = it.second;
        if (std::fabs(it.first) > (height/binSizeTh) && tmpVec.size() < 50)
        {
            outliers.insert(outliers.end(), tmpVec.begin(),tmpVec.end());
        }
    }
    
    auto leftCloud = geo::getSubSet(inputCloud, outliers, true);
    LOG(INFO) << "inputCloud:" << inputCloud->size() << " leftCloud:" << leftCloud->size();
    inputCloud->swap(*leftCloud);

    int number = Bottons.size() + Tops.size();
    if (zToNumVec.size() < number)
    {
        return false;
    }
    return true;
}

void Segmentation::estimateNormals(float radius, pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, 
                            pcl::PointCloud<pcl::Normal>::Ptr normals)
{
    pcl::search::KdTree<pcl::PointXYZ>::Ptr kdtree(new pcl::search::KdTree<pcl::PointXYZ>);
    pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> est;
    est.setInputCloud(cloud);   
    est.setRadiusSearch(radius);
    est.setSearchMethod(kdtree);  
    est.compute(*normals);
}

void Segmentation::regionGrowingSegmentation(
                            pcl::PointCloud<pcl::PointXYZ>::Ptr cloud,
                            pcl::PointCloud <pcl::Normal>::Ptr normals,
                            std::vector <pcl::PointIndices> &clusters)
{
    pcl::search::Search<pcl::PointXYZ>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZ>);
    pcl::RegionGrowing<pcl::PointXYZ, pcl::Normal> reg;
    reg.setMinClusterSize (700);
    // reg.setMaxClusterSize (20);
    reg.setSearchMethod (tree);
    reg.setNumberOfNeighbours (30);
    reg.setInputCloud (cloud);
    reg.setInputNormals (normals);
    reg.setSmoothnessThreshold (3.0 / 180.0 * M_PI);
    reg.setCurvatureThreshold (0.4);
    reg.extract (clusters);
}

pcl::PointCloud<pcl::Normal>::Ptr Segmentation::getSubSet(pcl::PointCloud<pcl::Normal>::Ptr cloud, 
                const std::vector<int>& indices, bool negative)
{
	pcl::PointCloud<pcl::Normal>::Ptr sub(new pcl::PointCloud<pcl::Normal>());

	pcl::PointIndices::Ptr pi(new pcl::PointIndices());
	pi->indices = indices;

	pcl::ExtractIndices<pcl::Normal> ei;
	ei.setInputCloud(cloud);
	ei.setIndices(pi);
	ei.setNegative(negative);
	ei.filter(*sub);

	return sub;
}

Eigen::Matrix<double, 2, 3> Segmentation::estimateTransform(
                const Eigen::Vector2d& s1, const Eigen::Vector2d& e1,
	            const Eigen::Vector2d& s2, const Eigen::Vector2d& e2) 
{
	Eigen::Matrix<double, 2, 3> T;
	Eigen::Vector2d n1 = (e1 - s1).normalized();
	Eigen::Vector2d n2 = (e2 - s2).normalized();
	double c = n1.dot(n2);
	double s = -(n1[0] * n2[1] - n1[1] * n2[0]);
	T(0, 0) = c;
	T(0, 1) = -s;
	T(1, 0) = s;
	T(1, 1) = c;

	Eigen::Vector2d cen1 = (s1 + e1) / 2.0;
	Eigen::Vector2d cen2 = (s2 + e2) / 2.0;

	T.block<2, 1>(0, 2) = cen1 - T.block<2, 2>(0, 0) * cen2;

	return T;
}

Eigen::Vector2d Segmentation::to_point_xy(const Eigen::Vector3d &point)
{
    return Eigen::Vector2d (point(0), point(1));
}

bool Segmentation::calibrateCloudByUpOrDown(CADModel &cad, 
                    pcl::PointCloud<pcl::PointXYZ>::Ptr ptSubSet1, 
                    pcl::PointCloud<pcl::Normal>::Ptr normalSubSet1, 
                    double &moveZ)
{
    //detect plane and their coeff
    pcl::PointCloud<pcl::PointXYZ>::Ptr copySet1(new pcl::PointCloud<pcl::PointXYZ>());
    pcl::copyPointCloud(*ptSubSet1, *copySet1);
    // pcl::PointCloud<pcl::Normal>::Ptr copyNormals1(new pcl::PointCloud<pcl::Normal>());
    // copyNormals1->insert(copyNormals1->end(), normalSubSet1->begin(), normalSubSet1->end());
    std::vector<Eigen::Vector4d> planeCoeffs;
    std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> planePoints;
    while (copySet1->size() > 200)
    {
        Eigen::Vector4d fitPlane;
        std::vector<int> inlierIndices;
        if (false == planeFitAndCheck(copySet1, fitPlane, inlierIndices)) break;

        auto inliers = geo::getSubSet(copySet1, inlierIndices, false);
        planeCoeffs.push_back(fitPlane);
        planePoints.push_back(inliers);
        auto tmpLeft = geo::getSubSet(copySet1, inlierIndices, true);
        copySet1->swap(*tmpLeft);
    }
    LOG(INFO) << "ptSubSet1:" << ptSubSet1->size() << " clusters:" << planePoints.size() << " copySet1:" << copySet1->size();

    //get planes' z value
    std::vector<double> vecPlaneZ;
    {
        pcl::PointCloud<pcl::PointXYZ>::Ptr tmpCloud(new pcl::PointCloud<pcl::PointXYZ>());
        Eigen::VectorXd line(6);
        line << 0, 0, 0, 0, 0, 1;
        for (auto &tmpPlane : planeCoeffs)
        {
            Eigen::Vector3d interSectionPt;
            if (false == interSectionOfLineToPlane(line, tmpPlane, interSectionPt)) continue;
            vecPlaneZ.push_back(interSectionPt(2));
        }
    }
    std::sort(vecPlaneZ.begin(), vecPlaneZ.end());

    //get cad planes' z value
    std::vector<double> vecCadZ;
    auto Bottoms = cad.getTypedModelItems(ITEM_BOTTOM_E);
    for (auto &it : Bottoms)
    {
        auto &p = it.points_.front();
        vecCadZ.push_back(p(2));
    }
    auto Tops = cad.getTypedModelItems(ITEM_TOP_E);
    for (auto &it : Tops)
    {
        auto &p = it.points_.front();
        vecCadZ.push_back(p(2));
    }
    // auto Beams = cad.getTypedModelItems(ITEM_BEAM_E);
    // for (auto &it : Beams)
    // {
    //     auto &p = it.points_.front();
    //     vecCadZ.push_back(p(2));
    // }
    std::sort(vecCadZ.begin(), vecCadZ.end());
    LOG(INFO) << "vecCadZ:" << vecCadZ.size() << " vecPlaneZ:" << vecPlaneZ.size();

    auto findMatchedZ = [](std::vector<double> &vec1, 
            std::vector<double> &vec2, double moveZ)->double {
        std::vector<double> vecMinDist;
        for (int i = 0; i < vec1.size(); i++)
        {
            auto &z1 = vec1[i];
            double minDist = double(RAND_MAX);
            for (int j = 0; j < vec2.size(); j++)
            {
                auto z2 = vec2[j] + moveZ;
                double dist = std::fabs(z1 - z2);
                minDist = std::min(minDist, dist);
            }

            vecMinDist.push_back(minDist);
        }
        double sum = std::accumulate(std::begin(vecMinDist), std::end(vecMinDist), 0.0);  
        double aveDist =  (vecMinDist.empty()) ? 2000 : (sum / vecMinDist.size());
        return aveDist;
    };

    //match and calc dist
    std::map<double, double> distToMoveZ;
    {
        for (int i = 0; i < vecCadZ.size(); i++)
        {
            auto &z1 = vecCadZ[i];
            for (int j = 0; j < vecPlaneZ.size(); j++)
            {
                auto &z2 = vecPlaneZ[j];
                double moveZ = z1 - z2;
                double aveDist = findMatchedZ(vecCadZ, vecPlaneZ, moveZ);
                distToMoveZ.insert(std::make_pair(aveDist, moveZ));
            }
        }
    }
    if (distToMoveZ.empty()) return false;
    if (distToMoveZ.begin()->first > 10) return false;
    LOG(INFO) << "moveZ:" << distToMoveZ.begin()->second << " dist:" << distToMoveZ.begin()->first;
    moveZ = distToMoveZ.begin()->second;

// #ifdef VISUALIZATION_ENABLED
    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
    T(2,3) = moveZ;

    std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> changedSubsets;
    for (int i = 0; i < planePoints.size(); i++)
    {
        auto tmp = planePoints[i];
        pcl::PointCloud<pcl::PointXYZ>::Ptr transformed_cloud(new pcl::PointCloud<pcl::PointXYZ>());
        pcl::transformPointCloud(*tmp, *transformed_cloud, T);
        changedSubsets.push_back(transformed_cloud);
    }

    {
        std::default_random_engine e;
        std::uniform_real_distribution<double> random(0,1);
        pcl::PointCloud<pcl::PointXYZRGB>::Ptr pCloudRgb(new pcl::PointCloud<pcl::PointXYZRGB>());
        for (int i = 0; i < changedSubsets.size(); i++)
        {
            auto tmp = changedSubsets[i];
            int r = int(random(e)*255);
            int g = int(random(e)*255);
            int b = int(random(e)*255);
            for (auto &p1 : tmp->points)
            {
                pcl::PointXYZRGB p2;
                p2.x = p1.x;
                p2.y = p1.y;
                p2.z = p1.z;
                p2.r = r;
                p2.g = g;
                p2.b = b;
                pCloudRgb->push_back(p2);
            }
        }
        std::string file_name = "changed-Subset1.pcd";
        pcl::io::savePCDFile(file_name, *pCloudRgb);	
    }
// #endif

    return true;
}

bool Segmentation::calibrateCloudByWall(CADModel &cad,
                        pcl::PointCloud<pcl::PointXYZ>::Ptr ptSubSet2, 
                        pcl::PointCloud<pcl::Normal>::Ptr normalSubSet2, 
                        Eigen::Matrix<double, 2, 3> &T)
{
    //detect planes and their coeff
    pcl::PointCloud<pcl::PointXYZ>::Ptr copySet2(new pcl::PointCloud<pcl::PointXYZ>());
    pcl::copyPointCloud(*ptSubSet2, *copySet2);
    // pcl::PointCloud<pcl::Normal>::Ptr copyNormals2(new pcl::PointCloud<pcl::Normal>());
    // copyNormals2->insert(copyNormals2->end(), normalSubSet2->begin(), normalSubSet2->end());
    std::vector<Eigen::Vector4d> wallCoeffs;
    std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> wallPoints;
    while (copySet2->size() > 200)
    {
        Eigen::Vector4d fitPlane;
        std::vector<int> inlierIndices;
        if (false == planeFitAndCheck(copySet2, fitPlane, inlierIndices)) break;

        auto inliers = geo::getSubSet(copySet2, inlierIndices, false);
        wallCoeffs.push_back(fitPlane);
        wallPoints.push_back(inliers);
        auto tmpLeft = geo::getSubSet(copySet2, inlierIndices, true);
        copySet2->swap(*tmpLeft);
    }
    LOG(INFO) << "ptSubSet2:" << ptSubSet2->size() << " clusters:" << wallPoints.size() << " copySet2:" << copySet2->size();
    
    //project to XOY
    std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> wallProjectPts;
    for (auto wall : wallPoints)
    {   
        pcl::PointCloud<pcl::PointXYZ>::Ptr projectPts(new pcl::PointCloud<pcl::PointXYZ>());
        for (auto &p : wall->points)
        {
            projectPts->push_back(pcl::PointXYZ(p.x, p.y, 0.0));
        }
        wallProjectPts.push_back(projectPts);
    }
    
    //detect lines in XOY
    float disthresh = 0.03;
    float connect_thresh = 0.1;
    float radius = 0.1;
    std::vector<int> parentIdxVec;
    std::vector<Eigen::VectorXf> lineCoeffs;
    std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> linePoints;
    std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> lineSegs;
    for (int i = 0; i < wallProjectPts.size(); i++)
    {   
        auto currProjectPts = wallProjectPts[i];
        std::vector<int> indices;
        Eigen::VectorXf params;
        std::tie(indices, params) = geo::detectOneLineRansac(currProjectPts, disthresh);
        auto inliers = geo::getSubSet(currProjectPts, indices, false);

        bool find = false;
        pcl::PointCloud<pcl::PointXYZ>::Ptr inlierPoints(new pcl::PointCloud<pcl::PointXYZ>());
        pcl::copyPointCloud(*inliers, *inlierPoints);
        while (inlierPoints->size() > 10)
        {   
            std::vector<int> segIndex = clusterMainStructure(inlierPoints, connect_thresh);
            auto segInliers = geo::getSubSet(inlierPoints, segIndex, false);

            std::pair<Eigen::Vector3d, Eigen::Vector3d> segment;
            if (true == detectLineEndPoints(segInliers, params, radius, segment))
            {
                find = true;
                lineSegs.push_back(segment);
            }
            auto tmpLeft = geo::getSubSet(inlierPoints, segIndex, true);
            inlierPoints->swap(*tmpLeft);			
        }
        if (find == true)
        {
            lineCoeffs.push_back(params);
            linePoints.push_back(inliers);
            parentIdxVec.push_back(i);            
        }
    }
    LOG(INFO) << "lineCoeffs size:" << lineCoeffs.size();

    auto pt_MinDist_to_cloud = [](Eigen::Vector3d &point, 
            pcl::PointCloud<pcl::PointXYZ>::Ptr cloud)->double {
        double minDist = double(RAND_MAX);
        for (auto &p : cloud->points)
        {
            Eigen::Vector3d pt(p.x, p.y, p.z);
            double dist = (point - pt).norm();
            minDist = std::min(minDist, dist);
        }
        return minDist;
    };

    //calc intersection node between two lines in XOY
    std::vector<Eigen::Vector3d> vecNodes;
    for (int i = 0; i < lineCoeffs.size(); i++)
    {
        auto tmp1 = lineCoeffs[i];
        Eigen::VectorXd line1(6);
        line1 << tmp1(0), tmp1(1), tmp1(2), tmp1(3), tmp1(4), tmp1(5);
        for (int j = i+1; j < lineCoeffs.size(); j++)
        {
            auto tmp2 = lineCoeffs[j];
            Eigen::VectorXd line2(6);
            line2 << tmp2(0), tmp2(1), tmp2(2), tmp2(3), tmp2(4), tmp2(5);
            Eigen::Vector3d interSectionPt;
            if (false == interSectionOfLineToLine(line1, line2, interSectionPt)) continue;

            double minDist1 = pt_MinDist_to_cloud(interSectionPt, linePoints[i]);
            double minDist2 = pt_MinDist_to_cloud(interSectionPt, linePoints[j]);
            if (minDist1 > 1.5 || minDist2 > 1.5) continue;

            vecNodes.push_back(interSectionPt);
        }
    }

    //get cad nodes in XOY
    std::vector<Eigen::Vector3d> vecCadNodes;
    auto Walls = cad.getTypedModelItems(ITEM_WALL_E);
    for (auto &it : Walls)
    {
        for (auto &seg : it.segments_)
        {
            if (seg.first(2) < 1e-5 && seg.second(2) < 1e-5)
            {
                vecCadNodes.push_back(seg.first);
                break;
            }
        }
    }
    LOG(INFO) << "vecCadNodes:" << vecCadNodes.size() << " vecNodes:" << vecNodes.size();

    auto get_IdxPairs = [](int vecSize)->std::vector<std::pair<int, int>> {
        std::vector<std::pair<int, int>> idxPairs;
        for (int i = 0; i < vecSize; i++)
        {
            for (int j = i+1; j < vecSize; j++)
            {
                idxPairs.push_back(std::make_pair(i,j));
            }
        }
        return idxPairs;
    };

    auto nodeMatching = [&](const std::vector<Eigen::Vector3d> &nodesVec1, 
            const std::vector<Eigen::Vector3d> &nodesVec2, 
            const Eigen::Matrix<double, 2, 3> &T)->std::map<std::pair<int, int>, double> {
        
        std::map<std::pair<int, int>, double> idxPairToDist;
        for (int i = 0; i < nodesVec1.size(); i++)
        {
            Eigen::Vector2d point = to_point_xy(nodesVec1[i]);
            double minDist = 20000;
            int bestIdx = -1;
            for (int j = 0; j < nodesVec2.size(); j++)
            {
                Eigen::Vector2d p = to_point_xy(nodesVec2[j]);
                Eigen::Vector2d q = T.block<2, 2>(0, 0) * p + T.block<2, 1>(0, 2);
                double dist = (point - q).norm();

                if (dist < minDist)
                {
                    minDist = dist;
                    bestIdx = j;
                }
            }
            if (-1 == bestIdx) continue;

            auto idxPair = std::make_pair(i,bestIdx);
            idxPairToDist.insert(std::make_pair(idxPair, minDist));
        }
        return idxPairToDist;
    };

    auto ave_dist = [](std::map<std::pair<int, int>, double> &idxPairToDist)->double {
        if (idxPairToDist.empty()) return double(RAND_MAX);

        double aveDist = 0.0;
        for (auto &it : idxPairToDist)
        {
            aveDist += it.second;
        }
        aveDist /= idxPairToDist.size();
        return aveDist;
    };

    //match nodes and calc dist
    std::map<double, Eigen::Matrix<double, 2, 3>> distToTrans;
    std::vector<std::pair<int, int>> idxPairsVec1 = get_IdxPairs(vecCadNodes.size());
    std::vector<std::pair<int, int>> idxPairsVec2 = get_IdxPairs(vecNodes.size());
    for (int i = 0; i < idxPairsVec1.size(); i++)
    {
        auto &pair1 = idxPairsVec1[i];
        Eigen::Vector2d s1 = to_point_xy(vecCadNodes[pair1.first]);
        Eigen::Vector2d e1 = to_point_xy(vecCadNodes[pair1.second]);
        double len1 = (s1 - e1).norm();
        for (int j = 0; j < idxPairsVec2.size(); j++)
        {
            auto &pair2 = idxPairsVec2[j];
            Eigen::Vector2d s2 = to_point_xy(vecNodes[pair2.first]);
            Eigen::Vector2d e2 = to_point_xy(vecNodes[pair2.second]);
            double len2 = (s2 - e2).norm();
            if (std::fabs(len1-len2) > 0.2) continue;
              
            Eigen::Matrix<double, 2, 3> Ts = estimateTransform(s1, e1, s2, e2); 
            Eigen::Matrix<double, 2, 3> Te = estimateTransform(e1, s1, s2, e2);

            auto idxPairToDist1 = nodeMatching(vecCadNodes, vecNodes, Ts);
            auto idxPairToDist2 = nodeMatching(vecCadNodes, vecNodes, Te);

            distToTrans.insert(std::make_pair(ave_dist(idxPairToDist1), Ts));
            distToTrans.insert(std::make_pair(ave_dist(idxPairToDist2), Te));
        }
    }
    if (distToTrans.empty())
    {
        return false;
    }
    auto bestItor = distToTrans.begin();
    LOG(INFO) << "bestMatch dist:" << bestItor->first;
    if (bestItor->first > 10) return false;

    T = bestItor->second;

// #ifdef VISUALIZATION_ENABLED
    Eigen::Matrix4d tmpT = Eigen::Matrix4d::Identity();
	tmpT.block<2, 2>(0, 0) = T.block<2, 2>(0, 0);
	tmpT.block<2, 1>(0, 3) = T.block<2, 1>(0, 2);
    pcl::PointCloud<pcl::PointXYZ>::Ptr nodeCloud(new pcl::PointCloud<pcl::PointXYZ>());
    for (auto &p : vecNodes)
    {
        nodeCloud->push_back(pcl::PointXYZ(p(0), p(1), p(2)));
    }
    {
        pcl::PointCloud<pcl::PointXYZ>::Ptr transformed_cloud(new pcl::PointCloud<pcl::PointXYZ>());
        pcl::transformPointCloud(*nodeCloud, *transformed_cloud, tmpT);
        nodeCloud->swap(*transformed_cloud);
    }

    std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> changedSubsets;
    for (int i = 0; i < wallPoints.size(); i++)
    {
        auto tmp = wallPoints[i];
        pcl::PointCloud<pcl::PointXYZ>::Ptr transformed_cloud(new pcl::PointCloud<pcl::PointXYZ>());
        pcl::transformPointCloud(*tmp, *transformed_cloud, tmpT);
        changedSubsets.push_back(transformed_cloud);
    }

    {
        auto pCloudRgb = to_rgp_cloud(nodeCloud, 255, 0, 0);
        std::string file_name = "nodeCloud.pcd";
        pcl::io::savePCDFile(file_name, *pCloudRgb);	
    }

    {
        std::default_random_engine e;
        std::uniform_real_distribution<double> random(0,1);
        pcl::PointCloud<pcl::PointXYZRGB>::Ptr pCloudRgb(new pcl::PointCloud<pcl::PointXYZRGB>());
        for (int i = 0; i < changedSubsets.size(); i++)
        {
            auto tmp = changedSubsets[i];
            int r = int(random(e)*255);
            int g = int(random(e)*255);
            int b = int(random(e)*255);
            auto pCurrCloud = to_rgp_cloud(tmp, r, g, b);
            pCloudRgb->insert(pCloudRgb->end(), pCurrCloud->begin(), pCurrCloud->end());
        }
        std::string file_name = "changed-Subset2.pcd";
        pcl::io::savePCDFile(file_name, *pCloudRgb);	
    }
// #endif

    return true;
}

bool Segmentation::calibrateCloudToCadPose(
                        pcl::PointCloud<pcl::PointXYZ>::Ptr inputCloud,
                        CADModel &cad, Eigen::Matrix4d &bestT)
{
    LOG(INFO) << "********calibrateCloudToCadPose*******";
    pcl::PointCloud<pcl::Normal>::Ptr cloudNormals(new pcl::PointCloud<pcl::Normal>());
    estimateNormals(0.1, inputCloud, cloudNormals);

    std::vector<int> upOrDown;
    std::vector<int> walls;
    Eigen::Vector3d axisZ(0,0,1);
    for (int i = 0; i < cloudNormals->size(); i++)
    {
        auto &p = cloudNormals->points[i];
        Eigen::Vector3d n(p.normal_x, p.normal_y, p.normal_z);
        double dot = std::fabs(n.dot(axisZ));
        if (dot > 0.9) upOrDown.push_back(i);
        if (dot < 0.2) walls.push_back(i);
    }

    auto ptSubSet1 = geo::getSubSet(inputCloud, upOrDown, false);
    auto ptSubSet2 = geo::getSubSet(inputCloud, walls, false);
    auto normalSubSet1 = getSubSet(cloudNormals, upOrDown, false);
    auto normalSubSet2 = getSubSet(cloudNormals, walls, false);

    double moveZ = 0.0;
    if (false == calibrateCloudByUpOrDown(cad, ptSubSet1, normalSubSet1, moveZ)) return false;

    Eigen::Matrix<double, 2, 3> T = Eigen::Matrix<double, 2, 3>::Identity();
    if (false == calibrateCloudByWall(cad, ptSubSet2, normalSubSet2, T)) return false;

    bestT = Eigen::Matrix4d::Identity();
	bestT.block<2, 2>(0, 0) = T.block<2, 2>(0, 0);
	bestT.block<2, 1>(0, 3) = T.block<2, 1>(0, 2);
    bestT(2, 3) = moveZ;

    pcl::PointCloud<pcl::PointXYZ>::Ptr transformed_cloud(new pcl::PointCloud<pcl::PointXYZ>());
    pcl::transformPointCloud(*inputCloud, *transformed_cloud, bestT);
    inputCloud->swap(*transformed_cloud);

// #ifdef VISUALIZATION_ENABLED
    {
        std::default_random_engine e;
        std::uniform_real_distribution<double> random(0,1);
        int r = int(random(e)*255);
        int g = int(random(e)*255);
        int b = int(random(e)*255);
        auto pCloudRgb = to_rgp_cloud(inputCloud, r, g, b);
        std::string file_name = "changed-inputCloud.pcd";
        pcl::io::savePCDFile(file_name, *pCloudRgb);	
    }
// #endif
}

pcl::PointCloud<pcl::PointXYZ>::Ptr Segmentation::segmentateSingleModelItem(
        pcl::PointCloud<pcl::PointXYZ>::Ptr inputCloud, const ModelItem &item,
        const std::string &name)
{
    pcl::PointCloud<pcl::PointXYZ>::Ptr cadPoints(new pcl::PointCloud<pcl::PointXYZ>());
    for (auto &p : item.points_)
    {
        cadPoints->push_back(pcl::PointXYZ(p(0), p(1), p(2)));
    }
    LOG(INFO) << "cadPoints:" << cadPoints->size();
    Eigen::Vector4d cadPlane;
    {
        Eigen::VectorXf coeff;
        std::vector<int> inlierIdxs;
        planeFitting(0.003, cadPoints, coeff, inlierIdxs);
        cadPlane = Eigen::Vector4d(coeff(0), coeff(1), coeff(2), coeff(3));
    }
    auto findMaxAxis = [](const Eigen::Vector3d &n)->int
    {
        int maxIdx = -1;
        if (n.norm() > 1e-5)
        {
            std::vector<double> valueVec;
            valueVec.push_back(std::fabs(n(0)));
            valueVec.push_back(std::fabs(n(1)));
            valueVec.push_back(std::fabs(n(2)));
            auto maxItor = std::max_element(valueVec.begin(), valueVec.end());
            maxIdx = std::distance(valueVec.begin(), maxItor);
        }
        return maxIdx;
    };
    int maxAxis = findMaxAxis(cadPlane.block<3,1>(0,0));
    if (-1 == maxAxis) return nullptr;
    LOG(INFO) << "maxAxis:" << maxAxis;

    auto to_2D_point = [](int deletCoordinate, const Eigen::Vector3d &p)->Eigen::Vector2d {
        Eigen::Vector2d pt2D;
        int j = 0;
        for (int i = 0; i < p.size(); i++)
        {
            if (i == deletCoordinate) continue;
            pt2D(j) = p(i);
            j++;
        }
        return pt2D;
    };
    Eigen::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>> segments2D;
    for (auto &seg : item.segments_)
    {
        auto pair2D = std::make_pair(to_2D_point(maxAxis, seg.first), to_2D_point(maxAxis, seg.second));
        segments2D.push_back(pair2D);
    }

    auto getNearPointsToPlane = [](pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, 
        const Eigen::Vector4d &plane, double distTh)->std::vector<int> {
        Eigen::Vector3d n = plane.block<3, 1>(0, 0);
        double d = plane(3);

        std::vector<int> indexVec;
        for (int i = 0; i < cloud->size(); i++)
        {
            auto &p = cloud->points[i];
            Eigen::Vector3d point(p.x, p.y, p.z);
            double dist = std::fabs(n.dot(point) + d);
            if (dist < distTh)
            {
                indexVec.push_back(i);
            }
        }
        return indexVec;
    };
    std::vector<int> indexVec = getNearPointsToPlane(inputCloud, cadPlane, 0.25);
    auto nearPoints = geo::getSubSet(inputCloud, indexVec, false);
    LOG(INFO) << "nearPoints:" << nearPoints->size();
// #ifdef VISUALIZATION_ENABLED
    {
        std::default_random_engine e;
        std::uniform_real_distribution<double> random(0,1);
        int r = int(random(e)*255);
        int g = int(random(e)*255);
        int b = int(random(e)*255);
        auto pCloudRgb = to_rgp_cloud(nearPoints, r, g, b, 0.01);
        std::string file_name = name + "-nearPoints.pcd";
        pcl::io::savePCDFile(file_name, *pCloudRgb);	
    }
// #endif
    Eigen::Vector4d fitPlane;
    pcl::PointCloud<pcl::PointXYZ>::Ptr fitPoints(new pcl::PointCloud<pcl::PointXYZ>());
    while (nearPoints->size() > 200)
    {
        std::vector<int> inlierIndices;
        if (false == planeFitAndCheck(nearPoints, fitPlane, inlierIndices)) break;

        LOG(INFO) << "inlierIndices:" << inlierIndices.size();
        fitPoints = geo::getSubSet(nearPoints, inlierIndices, false);
        if (std::fabs(cadPlane.block<3,1>(0,0).dot(fitPlane.block<3,1>(0,0))) > 0.9)
        {
            break;
        }
        auto tmpLeft = geo::getSubSet(nearPoints, inlierIndices, true);
        nearPoints->swap(*tmpLeft);
    }
    LOG(INFO) << "cadPlane:" << cadPlane(0) << " " << cadPlane(1) << " " << cadPlane(2) << " " << cadPlane(3);
    LOG(INFO) << "fitPlane:" << fitPlane(0) << " " << fitPlane(1) << " " << fitPlane(2) << " " << fitPlane(3);
    LOG(INFO) << "fitPoints:" << fitPoints->size();
    if (fitPoints->empty()) return false;

    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_filter(new pcl::PointCloud<pcl::PointXYZ>());
    projectionToPlane(cadPlane, fitPoints, cloud_filter);

    std::vector<int> indices;
    for (int i = 0; i < cloud_filter->size(); i++)
    {
        auto &p = cloud_filter->points[i];
        Eigen::Vector3d point(p.x, p.y, p.z);
        Eigen::Vector2d point2D = to_2D_point(maxAxis, point);
        if (isPointInPolygon2D(point2D, segments2D))
        {
            indices.push_back(i);
        }
    }
    LOG(INFO) << "indices:" << indices.size();
    auto segmentedCloud = geo::getSubSet(fitPoints, indices, false);
// #ifdef VISUALIZATION_ENABLED
    {
        std::default_random_engine e;
        std::uniform_real_distribution<double> random(0,1);
        int r = int(random(e)*255);
        int g = int(random(e)*255);
        int b = int(random(e)*255);
        auto pCloudRgb = to_rgp_cloud(fitPoints, r, g, b, 0.01);
        std::string file_name = name + "-fitPoints.pcd";
        pcl::io::savePCDFile(file_name, *pCloudRgb);	
    }
// #endif

    return segmentedCloud;
}

bool Segmentation::doSegmentation(
        pcl::PointCloud<pcl::PointXYZ>::Ptr inputCloud, CADModel &cad)
{
    std::default_random_engine e;
    std::uniform_real_distribution<double> random(0,1);

    LOG(INFO) << "********doSegmentation*******";
    LOG(INFO) << "inputCloud:" << inputCloud->size();
    std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> segmentWallPoints;
    auto Walls = cad.getTypedModelItems(ITEM_WALL_E);
    LOG(INFO) << "Walls:" << Walls.size();
    for (int k = 0; k < Walls.size(); k++)
    {
        auto &item = Walls[k];
        LOG(INFO) << "=========== segmentate wall:" << k;
        auto segmentPoints = segmentateSingleModelItem(inputCloud, item, "wall-" + std::to_string(k));
        if (nullptr == segmentPoints || segmentPoints->empty()) continue;

        segmentWallPoints.push_back(segmentPoints);       
        {
            int r = int(random(e)*255);
            int g = int(random(e)*255);
            int b = int(random(e)*255);
            auto pCloudRgb = to_rgp_cloud(segmentPoints, r, g, b, 0.01);
            std::string file_name = "wall-" + std::to_string(k) + "-segmentPoints.pcd";
            pcl::io::savePCDFile(file_name, *pCloudRgb);	
        } 
    } 

    std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> segmentTopPoints;
    auto Tops = cad.getTypedModelItems(ITEM_TOP_E);
    LOG(INFO) << "Tops:" << Tops.size();
    for (int k = 0; k < Tops.size(); k++)
    {
        auto &item = Tops[k];
        LOG(INFO) << "=========== segmentate top:" << k;
        auto segmentPoints = segmentateSingleModelItem(inputCloud, item, "top-" + std::to_string(k));
        if (nullptr == segmentPoints || segmentPoints->empty()) continue;

        segmentTopPoints.push_back(segmentPoints);
        {
            int r = int(random(e)*255);
            int g = int(random(e)*255);
            int b = int(random(e)*255);
            auto pCloudRgb = to_rgp_cloud(segmentPoints, r, g, b, 0.01);
            std::string file_name = "top-" + std::to_string(k) + "-segmentPoints.pcd";
            pcl::io::savePCDFile(file_name, *pCloudRgb);	
        } 
    }   

    std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> segmentBottomPoints;
    auto Bottoms = cad.getTypedModelItems(ITEM_BOTTOM_E);
    LOG(INFO) << "Bottoms:" << Bottoms.size();
    for (int k = 0; k < Bottoms.size(); k++)
    {
        auto &item = Bottoms[k];
        LOG(INFO) << "=========== segmentate bottom:" << k;
        auto segmentPoints = segmentateSingleModelItem(inputCloud, item, "bottom-" + std::to_string(k));
        if (nullptr == segmentPoints || segmentPoints->empty()) continue;

        segmentBottomPoints.push_back(segmentPoints);
        {
            int r = int(random(e)*255);
            int g = int(random(e)*255);
            int b = int(random(e)*255);
            auto pCloudRgb = to_rgp_cloud(segmentPoints, r, g, b, 0.01);
            std::string file_name = "bottom-" + std::to_string(k) + "-segmentPoints.pcd";
            pcl::io::savePCDFile(file_name, *pCloudRgb);	
        } 
    }
}


bool Segmentation::planeFitAndCheck(pcl::PointCloud<pcl::PointXYZ>::Ptr nearPoints,
                        Eigen::Vector4d &fitPlane, std::vector<int> &inlierIndices)
{
    pcl::PointCloud<pcl::PointXYZ>::Ptr sampleCloud(new pcl::PointCloud<pcl::PointXYZ>());
    uniformSampling(0.03, nearPoints, sampleCloud);

    pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>());
    estimateNormals(0.06, sampleCloud, normals);
    LOG(INFO) << "normals:" << normals->size();

    Eigen::VectorXf coeff;
    std::vector<int> indices;
    planeFitting(0.05, nearPoints, coeff, indices);
    fitPlane = Eigen::Vector4d(coeff(0), coeff(1), coeff(2), coeff(3));
    auto inliers = geo::getSubSet(nearPoints, indices, false);
    LOG(INFO) << "inliers:" << inliers->size();

    auto findNearestPtIdx = [](pcl::PointXYZ &point, 
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud)->int {
        double minDist = double(RAND_MAX);
        int index = -1;
        for (int i = 0; i < cloud->size(); i++)
        {
            auto &p = cloud->points[i];
            double xDist = std::fabs(point.x - p.x);
            if (xDist > 0.03) continue;
            double yDist = std::fabs(point.y - p.y);
            if (yDist > 0.03) continue;
            double zDist = std::fabs(point.z - p.z);
            if (zDist > 0.03) continue;

            double dist = xDist*xDist + yDist*yDist + zDist*zDist;
            if (dist < minDist)
            {
                minDist = dist;
                index = i;
            }
        }
        return index;
    };

    inlierIndices.clear();
    for (int i = 0; i < inliers->size(); i++)
    {
        auto &p = inliers->points[i];
        int idx = findNearestPtIdx(p, sampleCloud);
        if (-1 == idx) continue;
        auto &np = normals->points[idx];
        Eigen::Vector3d n(np.normal_x, np.normal_y, np.normal_z);
        if (std::fabs(fitPlane.block<3,1>(0,0).dot(n)) < 0.7) continue;
        inlierIndices.push_back(indices[i]);
    }
    if (inlierIndices.size() < 10) return false;

    return true;
}

pcl::PointCloud<pcl::PointXYZRGB>::Ptr Segmentation::to_rgp_cloud(
                            const pcl::PointCloud<pcl::PointXYZ>::Ptr cloud,
                            int r, int g, int b, double sampleRadius)
{
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr pCloudRgb(new pcl::PointCloud<pcl::PointXYZRGB>());
    pcl::PointCloud<pcl::PointXYZ>::Ptr tmpCloud(new pcl::PointCloud<pcl::PointXYZ>());
    if (sampleRadius > 1e-5)
    {
        uniformSampling(sampleRadius, cloud, tmpCloud);
    }
    else
    {
        tmpCloud = cloud;
    }
    for (const auto &p1 : tmpCloud->points)
    {
        pcl::PointXYZRGB p2;
        p2.x = p1.x;
        p2.y = p1.y;
        p2.z = p1.z;
        p2.r = r;
        p2.g = g;
        p2.b = b;
        pCloudRgb->push_back(p2);
    } 
    return pCloudRgb;
}

} //namespace CloudReg
