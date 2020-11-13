
#include "BaseType.h"
#include "TransformOptimize.h"


namespace cloudReg
{

void TransformOptimize::run(std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &walls_vec,
                                std::vector<Eigen::Vector4d> &modelPlanes,
                                Eigen::Matrix4d &transform)
{
    clear();
    initSolver();

    //add vertex
    addSE3Vertex(&transform, g2o::SE3Quat(), false, false);

    //add edges
    addWallPointToModelPlaneEdges(walls_vec, modelPlanes, transform);

    optData(10, false, true);

    getSE3Transfor(transform);
}

bool TransformOptimize::addWallPointToModelPlaneEdges(
                        std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &walls_vec,
                        std::vector<Eigen::Vector4d> &modelPlanes,
                        Eigen::Matrix4d &transform)
{
	const VertexID_G2O_t id = getVertexID(&transform);
	if (INVAILD_G2O_VERTEXID == id)
	{
		std::cerr << "can't get transform id" << std::endl;
		return false;
	}

    for (int i = 0; i < walls_vec.size(); i++)
    {
        Eigen::Vector4d plane = modelPlanes[i];
        auto cloud = walls_vec[i];

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
		std::cerr << "can't get transform id" << std::endl;
		return false;
	}
	g2o::VertexSE3Expmap *pVertexSE3 = static_cast<g2o::VertexSE3Expmap *>(optimizer_.vertex(id));
	g2o::SE3Quat estimate = pVertexSE3->estimate();
	transform = estimate.to_homogeneous_matrix();

	return true;
}

} //namespace CloudReg
