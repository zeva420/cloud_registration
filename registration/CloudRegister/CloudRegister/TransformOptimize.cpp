#include "TransformOptimize.h"
#include "glog/logging.h"
#include "funHelper.h"

#include <pcl/common/transforms.h>

//#include <pcl/keypoints/uniform_sampling.hpp>
#include <pcl/filters/uniform_sampling.h>

#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/sample_consensus/ransac.h>
#include <pcl/sample_consensus/sac_model_plane.h>
#include <pcl/registration/transformation_estimation_svd.h>



#include <random>


namespace CloudReg
{
bool TransformOptimize::run(
                const std::map<ModelItemType, std::vector<PointCloud::Ptr>> &mapCloudItem,
                const CADModel &cadModel, const Eigen::Vector3d &center)
{
    for (auto &it : mapCloudItem)
    {
        if (ITEM_HOLE_E == it.first) continue;
        LOG(INFO) << "input " <<  toModelItemName(it.first) << " vecSize:" << it.second.size();
        Eigen::vector<PointsAndPlane> vecItems;
        for (auto &cloud : it.second)
        {
            PointsAndPlane item;
            item.cloudPtr_ = cloud;
            vecItems.push_back(item);
        }
        type2CloudItems_[it.first] = vecItems;
    }

    //get model plane coeff
    getModelPlaneCoeff(cadModel, center);

    matchCloudToMode();

    downSampling();

    //optimize with sampling Cloud
    Eigen::Matrix4d transform;
    optimize(transform);

    //transform with input Cloud
    transformCloud(transform);

    //get plane coeff with input Cloud
    Eigen::Vector3d newCenter = transform.block<3,3>(0,0) * center + transform.block<3,1>(0,3);
    getCloudPlaneCoeff(newCenter);

    optRets_.mapClouds_.clear();
    fillResult(transform, optRets_);

    //view Dist with sampling Cloud
    viewModelAndChangedCloud();

    return optRets_.valid() ? true : false;
}

bool TransformOptimize::downSampling()
{
    LOG(INFO) << "********downSampling*******";
	double radius = 0.05;
    
    for (auto &it : type2CloudItems_)
    {
        auto &vecItems = it.second;
		Eigen::vector<PointsAndPlane> vecSamplingItems;
        for (auto &item : vecItems)
        {
            PointCloud::Ptr cloud_sampling(new pcl::PointCloud<pcl::PointXYZ>);
            uniformSampling(radius, item.cloudPtr_, cloud_sampling);

            PointsAndPlane samplingItem;
            samplingItem.cloudPtr_ = cloud_sampling;
            vecSamplingItems.push_back(samplingItem);
            LOG(INFO) << "uniformSampling, radius:" << radius << ", "
                << toModelItemName(it.first) << " cloud size before:" << item.cloudPtr_->size()
                << " after: " << samplingItem.cloudPtr_->size();
        }
        type2SamplingItems_[it.first] = vecSamplingItems;
    }

	return true;
}

bool TransformOptimize::getModelPlaneCoeff(const CADModel &cadModel,
                                            const Eigen::Vector3d &center)
{
    LOG(INFO) << "********getModelPlaneCoeff*******";
    std::vector<ModelItem> walls3d = cadModel.getTypedModelItems(ITEM_WALL_E);
    std::vector<ModelItem> bottoms3d = cadModel.getTypedModelItems(ITEM_BOTTOM_E);
    std::vector<ModelItem> tops3d = cadModel.getTypedModelItems(ITEM_TOP_E);
    std::vector<ModelItem> beams3d = cadModel.getTypedModelItems(ITEM_BEAM_E);

    auto getModelItemPoints = [](std::vector<ModelItem> &modelItems, const Eigen::Vector3d &center)
        ->Eigen::vector<PointsAndPlane>
    {
		Eigen::vector<PointsAndPlane> vecItems;
        for (auto &it : modelItems)
        {
            PointCloud::Ptr cloud(new PointCloud());
            cloud->resize(it.points_.size());
            std::transform(it.points_.begin(), it.points_.end(), cloud->begin(),
                [](const Eigen::Vector3d& v) { return pcl::PointXYZ(v(0), v(1), v(2)); });

            PointsAndPlane item;
            item.cloudPtr_ = cloud;
            vecItems.push_back(item);
        }
        return vecItems;     
    };

    auto wallItems = getModelItemPoints(walls3d, center);
    auto bottomItems = getModelItemPoints(bottoms3d, center);
    auto topItems = getModelItemPoints(tops3d, center);
    auto beamItems = getModelItemPoints(beams3d, center);
    type2ModelItems_[ITEM_WALL_E] = wallItems;
    type2ModelItems_[ITEM_BOTTOM_E] = bottomItems;
    type2ModelItems_[ITEM_TOP_E] = topItems;
    type2ModelItems_[ITEM_BEAM_E] = beamItems;

    for (auto &it : type2ModelItems_)
    {
        LOG(INFO) << "for " << toModelItemName(it.first);
        auto &vecModelItems = it.second;
        for (auto &item : vecModelItems)
        {
            Eigen::Vector4d plane;    
            calcPlaneCoeff(item.cloudPtr_, center, plane);
            item.plane_ = plane;  
        }
    }
	return true;
}

bool TransformOptimize::getCloudPlaneCoeff(const Eigen::Vector3d &center)
{
    LOG(INFO) << "********getCloudPlaneCoeff*******";
    for (auto &it : type2CloudItems_)
    {
        LOG(INFO) << "for " << toModelItemName(it.first);
        auto &vecCloudItems = it.second;
        for (auto &item : vecCloudItems)
        {
            Eigen::Vector4d plane;    
            calcPlaneCoeff(item.cloudPtr_, center, plane);
            item.plane_ = plane;  
        }
    }

	return true;
}

bool TransformOptimize::calcPlaneCoeff(PointCloud::Ptr inputCloud, 
                    const Eigen::Vector3d &center, Eigen::Vector4d &planeCoeff)
{
    Eigen::VectorXf coeff;
    std::vector<int> inlierIdxs;
    planeFitting(0.003, inputCloud, coeff, inlierIdxs);
    // auto inliers = geo::getSubSet(inputCloud, inlierIdxs, false);
    // Eigen::Vector4d plane = calcPlaneParam(inliers);

    Eigen::Vector4d plane(coeff(0), coeff(1), coeff(2), coeff(3));
    Eigen::Vector3d n = plane.block<3,1>(0,0);
    double d = plane(3);
    double flag = ((n.dot(center) + d) > 0) ? 1.0 : -1.0;  
    plane = flag * plane;     
    LOG(INFO) << "plane coeff " << plane[0] << " " << plane[1] 
        << " " << plane[2] << " " << plane[3] 
        << ", inlierRate:" << 100.0 * double(inlierIdxs.size()) / double(inputCloud->size()) << "%";
    planeCoeff = plane;

	return true;
}

std::pair<double,double> TransformOptimize::calcCloudToPLaneAveDist(Eigen::Vector4d &plane,
                                PointCloud::Ptr cloud, bool bMedian)
{
    double aveDist = 0.0;
	double medianDist = 0.0;
	std::vector<double> vecDist;
    for (auto &p : cloud->points)
    {
        double dist = std::abs(pointToPLaneDist(plane, p));
        aveDist += dist;
		vecDist.emplace_back(dist);
    }
    if (0 < cloud->size())
    {
        aveDist /= double(cloud->size());


		if (bMedian && vecDist.size() > 1)
		{
			std::sort(vecDist.begin(), vecDist.end());
			if (vecDist.size() % 2 == 1)
			{
				medianDist = vecDist[vecDist.size() / 2 - 1];
			}
			else
			{
				medianDist = 0.5 *(vecDist[vecDist.size() / 2 - 1] + vecDist[vecDist.size() / 2]);

			}
		}

    }

    return std::make_pair(aveDist, medianDist);
}

bool TransformOptimize::matchCloudToMode()
{
    LOG(INFO) << "********matchCloudToMode*******";
    for (auto &it1 : type2ModelItems_)
    {
        auto type = it1.first;
		Eigen::vector<PointsAndPlane> &vecModelItems = it1.second;
        auto it2 = type2CloudItems_.find(type);
        if (it2 == type2CloudItems_.end()) continue; 
		Eigen::vector<PointsAndPlane> &vecCloudItems = it2->second;
		if (vecModelItems.size() != vecCloudItems.size())
		{
			LOG(WARNING) << "the vecModelItems size mismatch";
			continue;
		}

        std::map<int, std::map<double, int>> model2CloudDists;
        for (int i = 0; i < vecModelItems.size(); i++)
        {
            auto &plane = vecModelItems[i].plane_;
            Eigen::Vector3d n = plane.block<3,1>(0,0);
            double d  = plane(3);
            for (int j = 0; j < vecCloudItems.size(); j++)
            {
                auto aveDist = calcCloudToPLaneAveDist(plane, vecCloudItems[j].cloudPtr_,false);
                model2CloudDists[i][aveDist.first] = j;
            }  
        }

		Eigen::vector<PointsAndPlane> matchedCloudItems;
        matchedCloudItems.resize(vecCloudItems.size());
        for (auto &it : model2CloudDists)
        {
            int i = it.first;
            int j = it.second.begin()->second;
            double dist = it.second.begin()->first;
            matchedCloudItems[i] = vecCloudItems[j];
            LOG(INFO) << toModelItemName(type) << " cloud to model plane aveDist:" 
                << dist << " match idxPair: " << i << "-" << j;
        }
        it2->second = matchedCloudItems;
    }

    return true;
}

bool TransformOptimize::optimize(Eigen::Matrix4d &transform)
{
    LOG(INFO) << "********optimize*******";
    clear();
    initSolver();

    //add vertex
    LOG(INFO) << "********addSE3Vertex*******";
    addSE3Vertex(&transform, g2o::SE3Quat(), false, false);

    //add edges
    LOG(INFO) << "********addWallPointToModelPlaneEdges*******";
    addWallPointToModelPlaneEdges(transform);

    LOG(INFO) << "********optData*******";
    optData(10, false, true);

    LOG(INFO) << "********getSE3Transfor*******";
    getSE3Transfor(transform);

	return true;
}

bool TransformOptimize::addWallPointToModelPlaneEdges(Eigen::Matrix4d &transform)
{
	const VertexID_G2O_t id = getVertexID(&transform);
	if (INVAILD_G2O_VERTEXID == id)
	{
		LOG(INFO) << "can't get transform id";
		return false;
	}

    for (auto &it : type2ModelItems_)
    {
        if (ITEM_BOTTOM_E == it.first) continue;

        auto &vecModelItems = it.second;
        auto &vecSamplingItems = type2SamplingItems_[it.first];
        if (vecSamplingItems.size() != vecModelItems.size()) continue;

        double weight = 10000.0;        		
        if (ITEM_TOP_E == it.first) weight = 700000.0;

        for (int i = 0; i < vecModelItems.size(); i++)
        {
            Eigen::Vector4d plane = vecModelItems[i].plane_;
            auto cloud = vecSamplingItems[i].cloudPtr_;
            weight /= double(cloud->size());
            for (auto &v : cloud->points)
            {
                Eigen::Vector3d pt(v.x, v.y, v.z);
                const Eigen::Matrix<double, 1, 1> information
                    =  weight* Eigen::Matrix<double, 1, 1>::Identity();
                g2o::EdgePtToPlaneDist* e = new g2o::EdgePtToPlaneDist(pt, plane);
                e->setVertex(0, optimizer_.vertex(id));
                e->setMeasurement(0.0);
                e->information() = information;
                addEdge("EdgePtToPlaneDist", e);
            } 
        }
    }

    return true;
}

bool TransformOptimize::getSE3Transfor(Eigen::Matrix4d &transform)
{
	const VertexID_G2O_t id = getVertexID(&transform);
	if (INVAILD_G2O_VERTEXID == id)
	{
		LOG(INFO) << "can't get transform id";
		return false;
	}
	g2o::VertexSE3Expmap *pVertexSE3 = static_cast<g2o::VertexSE3Expmap *>(optimizer_.vertex(id));
	g2o::SE3Quat estimate = pVertexSE3->estimate();
	transform = estimate.to_homogeneous_matrix();

	return true;
}

bool TransformOptimize::transformCloud(Eigen::Matrix4d &finalT)
{
    LOG(INFO) << "********transformCloud*******";
    // std::stringstream ss;
    // ss << std::fixed << std::setprecision(8);
    // ss << "\n------final R------\n";
    // ss << finalT (0,0) << " " << finalT (0,1) << " " << finalT (0,2) << "\n";
    // ss << finalT (1,0) << " " << finalT (1,1) << " " << finalT (1,2) << "\n";
    // ss << finalT (2,0) << " " << finalT (2,1) << " " << finalT (2,2) << "\n";

    // ss << "------final t------\n";
    // ss << finalT (0,3) << " " << finalT (1,3) << " " << finalT (2,3) << "\n";
    // LOG(INFO) << ss.str();
    LOG(INFO) << "------final T------\n" << finalT;

	auto transfor = [&](Eigen::Matrix4d& finalT, PointCloud::Ptr cloud) {
		
		PointCloud::Ptr transformed_cloud(new pcl::PointCloud<pcl::PointXYZ>());
		pcl::transformPointCloud(*cloud, *transformed_cloud, finalT);
		cloud->swap(*transformed_cloud);
	};

	Eigen::Matrix4d newT = Eigen::Matrix4d::Identity();
    for (auto &it : type2ModelItems_)
    {
        auto &vecModelItems = it.second;
        auto &vecCloudItems = type2CloudItems_[it.first];
        if (vecModelItems.size() != vecCloudItems.size()) continue;

        for (int i = 0; i < vecModelItems.size(); i++)
        {              
            transfor(finalT, vecCloudItems[i].cloudPtr_);
            auto distError = calcCloudToPLaneAveDist(vecModelItems[i].plane_, vecCloudItems[i].cloudPtr_,true);
            LOG(INFO) << "first: " << toModelItemName(it.first) << " cloud to model plane, aveDist:" 
                << distError.first << " medianDist:" << distError.second;
            
            if (ITEM_BOTTOM_E == it.first)
            {
                newT(2, 3) = -distError.second;
            }
        }

        
    }
	for (auto &it : type2ModelItems_)
	{
		auto &vecModelItems = it.second;
		auto &vecCloudItems = type2CloudItems_[it.first];
		if (vecModelItems.size() != vecCloudItems.size()) continue;

		for (int i = 0; i < vecModelItems.size(); i++)
		{
			transfor(newT, vecCloudItems[i].cloudPtr_);
			//auto distError = calcCloudToPLaneAveDist(vecModelItems[i].plane_, vecCloudItems[i].cloudPtr_.true);
			//LOG(INFO) << "second: " << toModelItemName(it.first) << " cloud to model plane, aveDist:"
			//	<< distError.first << " medianDist:" << distError.second;
		}
	}
    finalT *= newT;

	return true;
}

bool TransformOptimize::fillResult(Eigen::Matrix4d &finalT, optCloudRets &optRets)
{
    for (auto &it : type2ModelItems_)
    {
        auto &vecModelItems = it.second;
        auto &vecCloudItems = type2CloudItems_[it.first];
        if (vecModelItems.size() != vecCloudItems.size()) continue;

		Eigen::vector<OptPlane> vecOptPlane;
        for (int i = 0; i < vecModelItems.size(); i++)
        {
            OptPlane piece;
            piece.cloud_ = vecCloudItems[i].cloudPtr_;
            piece.cloudPlane_ = vecCloudItems[i].plane_;
            piece.cadPlane_ = vecModelItems[i].plane_;
            vecOptPlane.push_back(piece);
        }
        optRets.mapClouds_[it.first] = vecOptPlane;
    }
    optRets.T_ = finalT;

	return true;
}



void TransformOptimize::projectCloudToXOYPlane(Eigen::Vector3d &startPt,
                PointCloud::Ptr model,
                Eigen::Matrix4f &T)
{
    pcl::PointXYZ firstPt = model->front();
    pcl::PointXYZ secondPt = model->points[1];
    pcl::PointXYZ lastPt = model->back();
    PointCloud::Ptr source(new pcl::PointCloud<pcl::PointXYZ>());
    source->push_back(firstPt);
    source->push_back(secondPt);
    source->push_back(lastPt);

    double l12 = std::sqrt(std::pow(secondPt.x-firstPt.x,2) 
            + std::pow(secondPt.y-firstPt.y,2) + std::pow(secondPt.z-firstPt.z,2));
    double l13 = std::sqrt(std::pow(lastPt.x-firstPt.x,2) 
            + std::pow(lastPt.y-firstPt.y,2) + std::pow(lastPt.z-firstPt.z,2));

    pcl::PointXYZ A= pcl::PointXYZ(startPt(0),startPt(1),startPt(2));
    pcl::PointXYZ B = pcl::PointXYZ(startPt(0)+l12,startPt(1),startPt(2));
    pcl::PointXYZ C = pcl::PointXYZ(startPt(0),startPt(1)+l13,startPt(2));
    PointCloud::Ptr traget(new pcl::PointCloud<pcl::PointXYZ>());
    traget->push_back(A);
    traget->push_back(B);
    traget->push_back(C);

    pcl::registration::TransformationEstimationSVD<pcl::PointXYZ, pcl::PointXYZ> TESVD;  
    // pcl::registration::TransformationEstimationSVD<pcl::PointXYZ, pcl::PointXYZ>::Matrix4 T;  
    TESVD.estimateRigidTransformation(*source, *traget, T); 

    std::stringstream ss;
    ss << std::fixed << std::setprecision(5);
    ss << "startPt:" << startPt(0) << "," << startPt(1) << "," << startPt(2) << "\n";
    ss << "project transform R ==============\n";
    ss << T (0,0) << " " << T (0,1) << " " << T (0,2) << "\n";
    ss << T (1,0) << " " << T (1,1) << " " << T (1,2) << "\n";
    ss << T (2,0) << " " << T (2,1) << " " << T (2,2) << "\n";

    ss << "project transform t ==============\n";
    ss << T (0,3) << " " << T (1,3) << " " << T (2,3) << "\n";
    // LOG(INFO) << ss.str();

    double maxLen = 0.0;
    for (int i = 1; i < model->size(); i++)
    {
        auto p = model->points[i];
        double l = std::sqrt(std::pow(p.x-firstPt.x,2) 
            + std::pow(p.y-firstPt.y,2) + std::pow(p.z-firstPt.z,2));
        maxLen = std::max(maxLen, l);
    }
    startPt(0) += maxLen;
}

bool TransformOptimize::viewModelAndChangedCloud()
{
    LOG(INFO) << "********viewModelAndChangedCloud*******";

#ifdef VISUALIZATION_ENABLED
    auto& tmpVisualMap = type2SamplingItems_;
    for (auto &it : type2ModelItems_)
    {
        auto &vecModelItems = it.second;
        auto &vecCloudItems = tmpVisualMap[it.first];
        if (vecModelItems.size() != vecCloudItems.size()) continue;

        for (int i = 0; i < vecModelItems.size(); i++)
        {
            Eigen::Vector4d &model_plane = vecModelItems[i].plane_;
            Eigen::Vector4d &cloud_plane = vecCloudItems[i].plane_;
            Eigen::Vector3d n1 = model_plane.block<3,1>(0,0);
            Eigen::Vector3d n2 = cloud_plane.block<3,1>(0,0);
            double d1 = model_plane(3);
            double d2 = cloud_plane(3);
            double dot = n1.dot(n2);
            double angle = std::acos(std::fabs(dot)) * 180.0 / 3.1415;
            LOG(INFO) << "i:" << i << " dot:" << dot << " angle:" << angle 
                << " d1:" << d1 << " d2:" << d2 << " diff:" << std::fabs(d1-d2);
        }
    }

    

    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud3D_dist2Model(new pcl::PointCloud<pcl::PointXYZRGB>);
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud2D_dist2Model(new pcl::PointCloud<pcl::PointXYZRGB>);
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud3D_dist2Cloud(new pcl::PointCloud<pcl::PointXYZRGB>);
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud2D_dist2Cloud(new pcl::PointCloud<pcl::PointXYZRGB>);
  
    Eigen::Vector3d startPt(0,0,0);
    for (auto &it : type2ModelItems_)
    {
        auto &vecModelItems = it.second;
        auto &vecCloudItems = tmpVisualMap[it.first];
        if (vecModelItems.size() != vecCloudItems.size()) continue;

        for (int i = 0; i < vecModelItems.size(); i++)
        {
            auto model = vecModelItems[i].cloudPtr_;
            auto cloud = vecCloudItems[i].cloudPtr_;
            Eigen::Vector4d model_plane = vecModelItems[i].plane_;
            Eigen::Vector4d cloud_plane = vecCloudItems[i].plane_;
            Eigen::Matrix4f T;
            projectCloudToXOYPlane(startPt, model, T);
            
            {
                pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_rgb(new pcl::PointCloud<pcl::PointXYZRGB>);
                for (auto &p : cloud->points)
                {
                    double dist = pointToPLaneDist(model_plane, p);
                    pcl::PointXYZRGB p_rgb = getColorPtByDist(p, dist);
                    cloud_rgb->push_back(p_rgb);
                }
                cloud3D_dist2Model->insert(cloud3D_dist2Model->end(), cloud_rgb->begin(), cloud_rgb->end());

                pcl::PointCloud<pcl::PointXYZRGB>::Ptr transforCloud(new pcl::PointCloud<pcl::PointXYZRGB>());
                pcl::transformPointCloud(*cloud_rgb, *transforCloud, T);
                cloud2D_dist2Model->insert(cloud2D_dist2Model->end(), transforCloud->begin(), transforCloud->end());
            }

            {
                pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_rgb(new pcl::PointCloud<pcl::PointXYZRGB>);
                for (auto &p : cloud->points)
                {
                    double dist = pointToPLaneDist(cloud_plane, p);
                    pcl::PointXYZRGB p_rgb = getColorPtByDist(p, dist);
                    cloud_rgb->push_back(p_rgb);
                }
                cloud3D_dist2Cloud->insert(cloud3D_dist2Cloud->end(), cloud_rgb->begin(), cloud_rgb->end());

                std::string label = "finalT-cloud" + std::to_string(i);

                pcl::PointCloud<pcl::PointXYZRGB>::Ptr transforCloud(new pcl::PointCloud<pcl::PointXYZRGB>());
                pcl::transformPointCloud(*cloud_rgb, *transforCloud, T);
                cloud2D_dist2Cloud->insert(cloud2D_dist2Cloud->end(), transforCloud->begin(), transforCloud->end());
            }
        }        

    }
    savePCDFile<pcl::PointXYZRGB>("dist-to-modelPlane-3D.pcd", *cloud3D_dist2Model);
    savePCDFile<pcl::PointXYZRGB>("dist-to-modelPlane-2D.pcd", *cloud2D_dist2Model);
    savePCDFile<pcl::PointXYZRGB>("dist-to-cloudPlane-3D.pcd", *cloud3D_dist2Cloud);
    savePCDFile<pcl::PointXYZRGB>("dist-to-cloudPlane-2D.pcd", *cloud2D_dist2Cloud);

    for (auto &it : type2ModelItems_)
    {
        auto &vecModelItems = it.second;
        auto &vecCloudItems = tmpVisualMap[it.first];
        if (vecModelItems.size() != vecCloudItems.size()) continue;

        for (int i = 0; i < vecModelItems.size(); i++)
        {
            auto model = vecModelItems[i].cloudPtr_;
            auto cloud = vecCloudItems[i].cloudPtr_;
            
            std::vector<Eigen::Vector3d> modelPoints;
            for (int j = 0; j < model->size()-1 ; j++)
            {
                pcl::PointXYZ &p1 = model->points[j];
                pcl::PointXYZ &p2 = model->points[j+1];

                auto vec_tmp = ininterpolateSeg(Eigen::Vector3d(p1.x, p1.y, p1.z), 
                    Eigen::Vector3d(p2.x, p2.y, p2.z), 0.05);
                modelPoints.insert(modelPoints.end(), vec_tmp.begin(), vec_tmp.end());
            }
            pcl::PointXYZ &p1 = model->back();
            pcl::PointXYZ &p2 = model->front();
            auto vec_tmp = ininterpolateSeg(Eigen::Vector3d(p1.x, p1.y, p1.z), 
                Eigen::Vector3d(p2.x, p2.y, p2.z), 0.05);
            modelPoints.insert(modelPoints.end(), vec_tmp.begin(), vec_tmp.end());

            pcl::PointCloud<pcl::PointXYZRGB>::Ptr optCloud(new pcl::PointCloud<pcl::PointXYZRGB>);
            for (auto &p : modelPoints)
            {
                pcl::PointXYZRGB p_rgb;
                p_rgb.x = p(0);
                p_rgb.y = p(1);
                p_rgb.z = p(2);
                p_rgb.r = 255;
                p_rgb.g = 0;
                p_rgb.b = 0;            
                optCloud->push_back(p_rgb);
            }
            for (auto &p : cloud->points)
            {
                pcl::PointXYZRGB p_rgb;
                p_rgb.x = p.x;
                p_rgb.y = p.y;
                p_rgb.z = p.z;
                p_rgb.r = 0;
                p_rgb.g = 255;
                p_rgb.b = 0; 
                optCloud->push_back(p_rgb);
            }
            savePCDFile<pcl::PointXYZRGB>("optimized-" + toModelItemName(it.first) 
                                + "-" + std::to_string(i)  + ".pcd", *optCloud); 
        }
    }

 
#endif
	return true;
}

} //namespace CloudReg
