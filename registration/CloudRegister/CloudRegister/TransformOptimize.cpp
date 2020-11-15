#include "TransformOptimize.h"
#include "glog/logging.h"

#include <pcl/common/transforms.h>

#include <pcl/keypoints/impl/uniform_sampling.hpp>

#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/sample_consensus/ransac.h>
#include <pcl/sample_consensus/sac_model_plane.h>

#include <pcl/visualization/pcl_visualizer.h>

namespace CloudReg
{
bool TransformOptimize::run(std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &vecCloudPtr,
			                CADModel &cadModel)
{
    downSampling(vecCloudPtr);

    std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> model_vec;
    convertToPclCloud(cadModel, model_vec);

    //get model plane coeff
    std::vector<Eigen::Vector4d> modelPlanes;
    getModelPlaneCoeff(model_vec, modelPlanes);

    matchCloudToMode(modelPlanes, vecCloudPtr);

    //optimize
    Eigen::Matrix4d transform;
    optimize(vecCloudPtr, modelPlanes, transform);

    transformCloud(modelPlanes, vecCloudPtr, transform);

    viewModelAndChangedCloud(model_vec, vecCloudPtr);

    return true;
}

void TransformOptimize::uniformSampling(double radius, 
                    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, 
                    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_filtered)
{
    int cloudSize = cloud->size();
    
    // Uniform sampling object.
    pcl::UniformSampling<pcl::PointXYZ> filter;
    filter.setInputCloud(cloud);
    filter.setRadiusSearch(radius);
    // We need an additional object to store the indices of surviving points.
    pcl::PointCloud<int> keypointIndices;

    filter.compute(keypointIndices);
    pcl::copyPointCloud(*cloud, keypointIndices.points, *cloud_filtered);

    LOG(INFO) << "uniformSampling, radius:" << radius 
        << " cloud size before:" << cloudSize 
        << " after: "<< cloud_filtered->size();
}

bool TransformOptimize::downSampling(std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &vecCloudPtr)
{
    LOG(INFO) << "********downSampling*******";
    for (auto cloud : vecCloudPtr)
    {
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_sampling(new pcl::PointCloud<pcl::PointXYZ>);
        uniformSampling(0.01, cloud, cloud_sampling);
        cloud->swap(*cloud_sampling);
    }
}

bool TransformOptimize::convertToPclCloud(CADModel &cadModel, 
                        std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &model_vec)
{
    const auto& bottom3d = cadModel.getTypedModelItems(ITEM_BOTTOM_E).front();
    const auto& top3d = cadModel.getTypedModelItems(ITEM_BOTTOM_E).back();
    std::vector<ModelItem> walls3d = cadModel.getTypedModelItems(ITEM_WALL_E);

    for (auto &item : walls3d)
    {
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
        cloud->resize(item.points_.size());
        std::transform(item.points_.begin(), item.points_.end(), cloud->begin(),
            [](const Eigen::Vector3d& v) { return pcl::PointXYZ(v(0), v(1), v(2)); });
        model_vec.push_back(cloud);
    }
    {
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
        cloud->resize(bottom3d.points_.size());
        std::transform(bottom3d.points_.begin(), bottom3d.points_.end(), cloud->begin(),
            [](const Eigen::Vector3d& v) { return pcl::PointXYZ(v(0), v(1), v(2)); });
        model_vec.push_back(cloud);
    }
    {
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
        cloud->resize(top3d.points_.size());
        std::transform(top3d.points_.begin(), top3d.points_.end(), cloud->begin(),
            [](const Eigen::Vector3d& v) { return pcl::PointXYZ(v(0), v(1), v(2)); });
        model_vec.push_back(cloud);
    }
    return true;
}

void TransformOptimize::planeFitting(double distTh, 
                pcl::PointCloud<pcl::PointXYZ>::Ptr ground, 
                Eigen::VectorXf &coeff, std::vector<int> &inlierIdxs)
{
    pcl::SampleConsensusModelPlane<pcl::PointXYZ>::Ptr model(
                        new pcl::SampleConsensusModelPlane<pcl::PointXYZ>(ground));
    pcl::RandomSampleConsensus<pcl::PointXYZ> ransac(model);
    ransac.setDistanceThreshold(distTh);
    ransac.computeModel();
    ransac.getInliers(inlierIdxs);
    
    //ax+by_cz_d=0，coeff: a,b,c,d
    ransac.getModelCoefficients(coeff);
    LOG(INFO) << "plane coeff " << coeff[0] << " " <<coeff[1] 
        << " " << coeff[2] << " " << coeff[3];
}

bool TransformOptimize::getModelPlaneCoeff(
                    std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &model_vec, 
                    std::vector<Eigen::Vector4d> &modelPlanes)
{
    LOG(INFO) << "********get model plane coeff*******";
    int index = -1;
    for (auto cloud : model_vec)
    {
        index++;
        pcl::PointCloud<pcl::PointXYZ>::Ptr ground_inliers(new pcl::PointCloud<pcl::PointXYZ>);
        Eigen::VectorXf coeff;
        std::vector<int> inlierIdxs;
        planeFitting(0.1, cloud, coeff, inlierIdxs);
        Eigen::Vector4d plane(coeff(0), coeff(1), coeff(2), coeff(3));
        modelPlanes.push_back(plane);  
    }
}

double TransformOptimize::calcCloudToPLaneAveDist(Eigen::Vector4d &plane,
                                pcl::PointCloud<pcl::PointXYZ>::Ptr cloud)
{
    Eigen::Vector3d n = plane.block<3,1>(0,0);
    double d  = plane(3);

    double aveDist = 0.0;
    for (auto &p : cloud->points)
    {
        Eigen::Vector3d point(p.x, p.y, p.z);
        double dist = std::abs(n.dot(point) + d) / n.norm();
        aveDist += dist;
    }
    if (0 < cloud->size())
    {
        aveDist /= double(cloud->size());
    }

    return aveDist;
}

bool TransformOptimize::matchCloudToMode(
                    std::vector<Eigen::Vector4d> &modePlanes,
                    std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &vecCloudPtr)
{
    LOG(INFO) << "********matchCloudToMode*******";
    if (modePlanes.size() != vecCloudPtr.size())
    {
        LOG(INFO) << "size is not eqaul!!";
        return false;
    }

    std::map<int, std::map<double, int>> model2CloudDists;
    for (int i = 0; i < modePlanes.size(); i++)
    {
        auto &plane = modePlanes[i];
        Eigen::Vector3d n = plane.block<3,1>(0,0);
        double d  = plane(3);
        for (int j = 0; j < vecCloudPtr.size(); j++)
        {
            double aveDist = calcCloudToPLaneAveDist(plane, vecCloudPtr[j]);
            model2CloudDists[i][aveDist] = j;
        }
    }

    std::stringstream ss;
    	ss << std::fixed << std::setprecision(4);
    ss << "====print model2CloudDists info====\n";
    for (auto &it : model2CloudDists)
    {
        ss << "model " << it.first << " match to:\n";
        for (auto &item : it.second)
        {
            ss << "   cloud " << item.second << " dist:" << item.first << "\n";
        }
    }
    // LOG(INFO) << ss.str();

    std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> tmpCloudPtr;
    tmpCloudPtr.resize(vecCloudPtr.size());
    for (auto &it : model2CloudDists)
    {
        int i = it.first;
        int j = it.second.begin()->second;
        tmpCloudPtr[i] = vecCloudPtr[j];
        LOG(INFO) << "cloud to plane aveDist:" << it.second.begin()->first;
    }
    vecCloudPtr.clear();
    vecCloudPtr.insert(vecCloudPtr.end(), tmpCloudPtr.begin(), tmpCloudPtr.end());
    return true;
}

bool TransformOptimize::optimize(std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &vecCloudPtr,
                                std::vector<Eigen::Vector4d> &modelPlanes,
                                Eigen::Matrix4d &transform)
{
    LOG(INFO) << "********optimize*******";
    clear();
    initSolver();

    //add vertex
    LOG(INFO) << "********addSE3Vertex*******";
    addSE3Vertex(&transform, g2o::SE3Quat(), false, false);

    //add edges
    LOG(INFO) << "********addWallPointToModelPlaneEdges*******";
    addWallPointToModelPlaneEdges(vecCloudPtr, modelPlanes, transform);

    LOG(INFO) << "********optData*******";
    optData(10, false, true);

    LOG(INFO) << "********getSE3Transfor*******";
    getSE3Transfor(transform);
}

bool TransformOptimize::addWallPointToModelPlaneEdges(
                        std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &vecCloudPtr,
                        std::vector<Eigen::Vector4d> &modelPlanes,
                        Eigen::Matrix4d &transform)
{
	const VertexID_G2O_t id = getVertexID(&transform);
	if (INVAILD_G2O_VERTEXID == id)
	{
		LOG(INFO) << "can't get transform id";
		return false;
	}

    for (int i = 0; i < vecCloudPtr.size(); i++)
    {
        Eigen::Vector4d plane = modelPlanes[i];
        auto cloud = vecCloudPtr[i];

        for (auto &v : cloud->points)
        {
            Eigen::Vector3d pt(v.x, v.y, v.z);
            const Eigen::Matrix<double, 1, 1> information
                =  Eigen::Matrix<double, 1, 1>::Identity();
            g2o::EdgePtToPlaneDist* e = new g2o::EdgePtToPlaneDist(pt, plane);
            e->setVertex(0, optimizer_.vertex(id));
            e->setMeasurement(0.0);
            e->information() = information;
            addEdge("EdgePtToPlaneDist", e);
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

bool TransformOptimize::transformCloud(
                    std::vector<Eigen::Vector4d> &modePlanes,
                    std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &vecCloudPtr,
                    Eigen::Matrix4d &finalT)
{
    LOG(INFO) << "********transformCloud*******";
    std::stringstream ss;
    ss << std::fixed << std::setprecision(5);
    ss << "------final R------\n";
    ss << finalT (0,0) << " " << finalT (0,1) << " " << finalT (0,2) << "\n";
    ss << finalT (1,0) << " " << finalT (1,1) << " " << finalT (1,2) << "\n";
    ss << finalT (2,0) << " " << finalT (2,1) << " " << finalT (2,2) << "\n";

    ss << "------final t------\n";
    ss << finalT (0,3) << " " << finalT (1,3) << " " << finalT (2,3) << "\n";
    LOG(INFO) << ss.str();

    for (int i = 0; i < vecCloudPtr.size(); i++)
    {
        auto &cloud = vecCloudPtr[i];
        pcl::PointCloud<pcl::PointXYZ>::Ptr transformed_cloud(new pcl::PointCloud<pcl::PointXYZ>());
        pcl::transformPointCloud(*cloud, *transformed_cloud, finalT);    
        cloud->swap(*transformed_cloud); 

        Eigen::Vector4d plane = modePlanes[i];
        double aveDist = calcCloudToPLaneAveDist(plane, cloud);
        LOG(INFO) << "cloud to plane, aveDist:" << aveDist;

        std::string fileName = "optimized-cloud-" + convertToSimpleIDStr(i)  + ".pcd";
        pcl::io::savePCDFile(fileName, *cloud);   
    }
}

bool TransformOptimize::viewModelAndChangedCloud(
                std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &model_vec,
                std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &vecCloudPtr)
{
    LOG(INFO) << "********viewModelAndChangedCloud*******";
    //viewer window init
    pcl::visualization::PCLVisualizer viewer("demo");
    // int v1(0);
    // viewer.createViewPort(0.0, 0.0, 1.0, 1.0, v1);
    // float bckgr_gray_level = 0.0;  // Black
    // float txt_gray_lvl = 1.0 - bckgr_gray_level;
    // viewer.setBackgroundColor(bckgr_gray_level, bckgr_gray_level, bckgr_gray_level, v1);     
    // viewer.setSize(1280, 1024);  // Visualiser window size
    std::default_random_engine e;
    std::uniform_real_distribution<double> random(0,1);

    int index = -1;
    for (auto cloud : model_vec)
    {
        index++;
        for (int i = 0; i < cloud->size()-1 ; i++)
        {
            pcl::PointXYZ &p1 = cloud->points[i];
            pcl::PointXYZ &p2 = cloud->points[i+1];
            std::string lineName = "wall" + convertToSimpleIDStr(index) 
                + "-line" + convertToSimpleIDStr(i) + "-" + convertToSimpleIDStr(i+1);
            // viewer.addLine(p1, p2, 255, 255, 255, lineName, v1);
            viewer.addLine(p1, p2, 255, 255, 255, lineName);
        }
        pcl::PointXYZ &p1 = cloud->back();
        pcl::PointXYZ &p2 = cloud->front();
        std::string lineName = "wall" + convertToSimpleIDStr(index) + "-line-b-f";
        // viewer.addLine(p1, p2, 255, 255, 255, lineName, v1);   
        viewer.addLine(p1, p2, 255, 255, 255, lineName);       
    }

    for (auto cloud : vecCloudPtr)
    {
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_sampling(new pcl::PointCloud<pcl::PointXYZ>);
        uniformSampling(0.07, cloud, cloud_sampling);
        cloud->swap(*cloud_sampling);
    }

    index = -1;
    for (auto cloud : vecCloudPtr)
    {
        index++;
        std::string label = "finalT-cloud" + convertToSimpleIDStr(index);
        pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> cloud_color(cloud, 0, 255, 0);
        // viewer.addPointCloud(cloud, cloud_color, label, v1);  
        viewer.addPointCloud(cloud, cloud_color, label);
    }

    while (!viewer.wasStopped())
    {
        viewer.spinOnce();
    }
}

} //namespace CloudReg
