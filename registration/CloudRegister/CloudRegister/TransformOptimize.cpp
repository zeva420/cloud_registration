#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkRenderingOpenGL)
#include <pcl/visualization/pcl_visualizer.h>

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
bool TransformOptimize::run(std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &vecCloudPtr,
			                CADModel &cadModel)
{
    std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> model_vec;
    convertToPclCloud(cadModel, model_vec);
    if (vecCloudPtr.size() != model_vec.size())
    {
        LOG(INFO) << "vecCloudPtr size:" << vecCloudPtr.size() 
            << " != model_vec size:" << model_vec.size();
        return false;
    }

    //get model plane coeff
    Eigen::vector<Eigen::Vector4d> modelPlanes;
    getModelPlaneCoeff(model_vec, modelPlanes);

    matchCloudToMode(modelPlanes, vecCloudPtr);

    std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> vecSamplingCloud;
    downSampling(vecCloudPtr, vecSamplingCloud);

    //optimize with sampling Cloud
    Eigen::Matrix4d transform;
    optimize(vecSamplingCloud, modelPlanes, transform);

    //transform with input Cloud
    transformCloud(modelPlanes, vecCloudPtr, transform);

    //get plane coeff with input Cloud
    Eigen::vector<Eigen::Vector4d> cloudPlanes;
    getModelPlaneCoeff(vecCloudPtr, cloudPlanes);

    optRets_.clear();
    fillResult(modelPlanes, cloudPlanes, model_vec, vecCloudPtr, optRets_);

    //view Dist with sampling Cloud
    viewModelAndChangedCloud(modelPlanes, cloudPlanes, model_vec, vecSamplingCloud);

    return optRets_.empty() ? false : true;
}

bool TransformOptimize::downSampling(std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &vecCloudPtr,
                        std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &vecSamplingCloud)
{
    LOG(INFO) << "********downSampling*******";
	double radius = 0.01;
    for (auto cloud : vecCloudPtr)
    {
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_sampling(new pcl::PointCloud<pcl::PointXYZ>);
        uniformSampling(radius, cloud, cloud_sampling);
        vecSamplingCloud.push_back(cloud_sampling);
		LOG(INFO) << "uniformSampling, radius:" << radius
			<< " cloud size before:" << cloud->size()
			<< " after: " << cloud_sampling->size();
    }

	return true;
}

bool TransformOptimize::convertToPclCloud(CADModel &cadModel, 
                        std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &model_vec)
{
    std::vector<ModelItem> walls3d = cadModel.getTypedModelItems(ITEM_WALL_E);
    std::vector<ModelItem> bottoms3d = cadModel.getTypedModelItems(ITEM_BOTTOM_E);
    std::vector<ModelItem> tops3d = cadModel.getTypedModelItems(ITEM_TOP_E);

    for (auto &item : walls3d)
    {
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
        cloud->resize(item.points_.size());
        std::transform(item.points_.begin(), item.points_.end(), cloud->begin(),
            [](const Eigen::Vector3d& v) { return pcl::PointXYZ(v(0), v(1), v(2)); });
        model_vec.push_back(cloud);
    }
    for (auto &item : bottoms3d)
    {
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
        cloud->resize(item.points_.size());
        std::transform(item.points_.begin(), item.points_.end(), cloud->begin(),
            [](const Eigen::Vector3d& v) { return pcl::PointXYZ(v(0), v(1), v(2)); });
        model_vec.push_back(cloud);
    }
    for (auto &item : tops3d)
    {
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
        cloud->resize(item.points_.size());
        std::transform(item.points_.begin(), item.points_.end(), cloud->begin(),
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
        << " " << coeff[2] << " " << coeff[3] 
        << ", inlierRate:" << 100.0 * double(inlierIdxs.size()) / double(ground->size()) << "%";
}

bool TransformOptimize::getModelPlaneCoeff(
                    std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &model_vec, 
                    Eigen::vector<Eigen::Vector4d> &modelPlanes)
{
    LOG(INFO) << "********get model plane coeff*******";
    int index = -1;
    for (auto cloud : model_vec)
    {
        index++;
        pcl::PointCloud<pcl::PointXYZ>::Ptr ground_inliers(new pcl::PointCloud<pcl::PointXYZ>);
        Eigen::VectorXf coeff;
        std::vector<int> inlierIdxs;
        planeFitting(0.05, cloud, coeff, inlierIdxs);
        Eigen::Vector4d plane(coeff(0), coeff(1), coeff(2), coeff(3));
        modelPlanes.push_back(plane);  
    }

	return true;
}

double TransformOptimize::pointToPLaneDist(Eigen::Vector4d &plane,
                                            pcl::PointXYZ &p)
{
    Eigen::Vector3d n = plane.block<3,1>(0,0);
    double d  = plane(3);
    Eigen::Vector3d point(p.x, p.y, p.z);
    double dist = std::abs(n.dot(point) + d) / n.norm();    
    return dist;
}

double TransformOptimize::calcCloudToPLaneAveDist(Eigen::Vector4d &plane,
                                pcl::PointCloud<pcl::PointXYZ>::Ptr cloud)
{
    double aveDist = 0.0;
    for (auto &p : cloud->points)
    {
        double dist = pointToPLaneDist(plane, p);
        aveDist += std::abs(dist);
    }
    if (0 < cloud->size())
    {
        aveDist /= double(cloud->size());
    }

    return aveDist;
}

bool TransformOptimize::matchCloudToMode(
                    Eigen::vector<Eigen::Vector4d> &modePlanes,
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
        LOG(INFO) << "cloud to model plane aveDist:" << it.second.begin()->first;
    }
    vecCloudPtr.clear();
    vecCloudPtr.insert(vecCloudPtr.end(), tmpCloudPtr.begin(), tmpCloudPtr.end());

    return true;
}

bool TransformOptimize::optimize(std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &vecCloudPtr,
                                Eigen::vector<Eigen::Vector4d> &modelPlanes,
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
    optData(15, false, true);

    LOG(INFO) << "********getSE3Transfor*******";
    getSE3Transfor(transform);

	return true;
}

bool TransformOptimize::addWallPointToModelPlaneEdges(
                        std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &vecCloudPtr,
                        Eigen::vector<Eigen::Vector4d> &modelPlanes,
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

        double weight = 10000.0 / double(cloud->size());
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
                    Eigen::vector<Eigen::Vector4d> &modePlanes,
                    std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &vecCloudPtr,
                    Eigen::Matrix4d &finalT)
{
    LOG(INFO) << "********transformCloud*******";
    std::stringstream ss;
    ss << std::fixed << std::setprecision(5);
    ss << "\n------final R------\n";
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
        LOG(INFO) << "cloud to model plane, aveDist:" << aveDist;

    }

	return true;
}

bool TransformOptimize::fillResult(
                Eigen::vector<Eigen::Vector4d> &modePlanes,
                Eigen::vector<Eigen::Vector4d> &cloudPlanes,
                std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &model_vec,
                std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &vecCloudPtr, 
                optCloudRets &optRets)
{
    OptCloud walls;
    walls.type_ = WALL_E;
    for (int i = 0; i < vecCloudPtr.size() - 2; i++)
    {
        auto cloud = vecCloudPtr[i];
        walls.vecCloud_.push_back(cloud);
        walls.vecCloudPlane_.push_back(cloudPlanes[i]);
        walls.vecCadPlane_.push_back(modePlanes[i]);

    }
    optRets[WALL_E] = walls;

    if (vecCloudPtr.size() > 2)
    {
        OptCloud bottom;
        bottom.type_ = BOTTOM_E;
        int index = vecCloudPtr.size() - 2;
        auto cloud  = vecCloudPtr[index];
        bottom.vecCloud_.push_back(cloud);
        bottom.vecCloudPlane_.push_back(cloudPlanes[index]);
        bottom.vecCadPlane_.push_back(modePlanes[index]);   
        optRets[BOTTOM_E] = bottom; 
    }

    if (vecCloudPtr.size() > 2)
    {
        OptCloud top;
        top.type_ = TOP_E;
        int index = vecCloudPtr.size() - 1;
        auto cloud  = vecCloudPtr[index];
        top.vecCloud_.push_back(cloud);
        top.vecCloudPlane_.push_back(cloudPlanes[index]);
        top.vecCadPlane_.push_back(modePlanes[index]);   
        optRets[TOP_E] = top; 
    }

	return true;
}



void TransformOptimize::projectCloudToXOYPlane(Eigen::Vector3d &startPt,
                pcl::PointCloud<pcl::PointXYZ>::Ptr model,
                Eigen::Matrix4f &T)
{
    pcl::PointXYZ firstPt = model->front();
    pcl::PointXYZ secondPt = model->points[1];
    pcl::PointXYZ lastPt = model->back();
    pcl::PointCloud<pcl::PointXYZ>::Ptr source(new pcl::PointCloud<pcl::PointXYZ>());
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
    pcl::PointCloud<pcl::PointXYZ>::Ptr traget(new pcl::PointCloud<pcl::PointXYZ>());
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

bool TransformOptimize::viewModelAndChangedCloud(
                Eigen::vector<Eigen::Vector4d> &modePlanes,
                Eigen::vector<Eigen::Vector4d> &cloudPlanes,
                std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &model_vec,
                std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &vecCloudPtr)
{
    LOG(INFO) << "********viewModelAndChangedCloud*******";
#if 1
    // pcl::visualization::PCLVisualizer viewer("demo");
    // for (int i = 0; i < model_vec.size(); i++)
    // {
    //     auto cloud = model_vec[i];
    //     for (int j = 0; j < cloud->size(); j++)
    //     {
    //         auto &p = cloud->points[j];
    //         viewer.addSphere(p, 0.3f, std::to_string(i) + std::to_string(j));
    //     }
    // }

    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud3D_dist2Model(new pcl::PointCloud<pcl::PointXYZRGB>);
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud2D_dist2Model(new pcl::PointCloud<pcl::PointXYZRGB>);
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud3D_dist2Cloud(new pcl::PointCloud<pcl::PointXYZRGB>);
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud2D_dist2Cloud(new pcl::PointCloud<pcl::PointXYZRGB>);
    // viewer.addText(
    //     "[-inf, -0.01]:blue, [-0.01, 0]:cyan, [0, 0.01]:green, [0.01, 0.02]:yellow, [0.02, 0.03]:pink, [0.03, +inf]:red", 
    //     0, 0.8);
    Eigen::Vector3d startPt(0,0,0);
    for (int i = 0; i < vecCloudPtr.size(); i++)
    {
        auto cloud = vecCloudPtr[i];
        auto model = model_vec[i];
        Eigen::Vector4d model_plane = modePlanes[i];
        Eigen::Vector4d cloud_plane = cloudPlanes[i];

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

            pcl::PointCloud<pcl::PointXYZRGB>::Ptr transformed_cloud (new pcl::PointCloud<pcl::PointXYZRGB>());
            pcl::transformPointCloud(*cloud_rgb, *transformed_cloud, T);
            cloud2D_dist2Model->insert(cloud2D_dist2Model->end(), transformed_cloud->begin(), transformed_cloud->end());
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

            std::string label = "finalT-cloud" + convertToSimpleIDStr(i);
            // viewer.addPointCloud(cloud_rgb, label);

            pcl::PointCloud<pcl::PointXYZRGB>::Ptr transformed_cloud (new pcl::PointCloud<pcl::PointXYZRGB>());
            pcl::transformPointCloud(*cloud_rgb, *transformed_cloud, T);
            cloud2D_dist2Cloud->insert(cloud2D_dist2Cloud->end(), transformed_cloud->begin(), transformed_cloud->end());
        
        }
    }
    savePCDFile<pcl::PointXYZRGB>("dist-to-modelPlane-3D.pcd", *cloud3D_dist2Model);
    savePCDFile<pcl::PointXYZRGB>("dist-to-modelPlane-2D.pcd", *cloud2D_dist2Model);
    savePCDFile<pcl::PointXYZRGB>("dist-to-cloudPlane-3D.pcd", *cloud3D_dist2Cloud);
    savePCDFile<pcl::PointXYZRGB>("dist-to-cloudPlane-2D.pcd", *cloud2D_dist2Cloud);

    for (int i = 0; i < model_vec.size(); i++)
    {
        auto modelCloud = model_vec[i];
        std::vector<Eigen::Vector3d> modelPoints;
        for (int j = 0; j < modelCloud->size()-1 ; j++)
        {
            pcl::PointXYZ &p1 = modelCloud->points[j];
            pcl::PointXYZ &p2 = modelCloud->points[j+1];

            auto vec_tmp = ininterpolateSeg(Eigen::Vector3d(p1.x, p1.y, p1.z), 
                Eigen::Vector3d(p2.x, p2.y, p2.z), 0.05);
            modelPoints.insert(modelPoints.end(), vec_tmp.begin(), vec_tmp.end());
        }
        pcl::PointXYZ &p1 = modelCloud->back();
        pcl::PointXYZ &p2 = modelCloud->front();
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
        for (auto &p : vecCloudPtr[i]->points)
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
        savePCDFile<pcl::PointXYZRGB>("optimized-cloud-" + std::to_string(i)  + ".pcd", *optCloud); 
    }

    // while (!viewer.wasStopped())
    // {
    //     viewer.spinOnce();
    // }
#endif
	return true;
}

} //namespace CloudReg
