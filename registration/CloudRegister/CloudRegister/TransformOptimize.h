#ifndef TRANSFORM_OPTIMIZE_H
#define TRANSFORM_OPTIMIZE_H

#include "BaseType.h"
#include "glog/logging.h"
#include "CADModel.h"

//edge
#include "EdgeConstraint.h"

//pcl
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>

/*
#include "g2o/core/factory.h"
#include "g2o/stuff/macros.h"
#include "g2o/types/slam3d/types_slam3d.h"
#include "g2o/core/sparse_optimizer.h"
#include "g2o/core/optimization_algorithm_factory.h"
#include "g2o/core/robust_kernel.h"
#include "g2o/core/robust_kernel_factory.h"
*/

#include "g2o/core/block_solver.h"
#include "g2o/solvers/eigen/linear_solver_eigen.h"
#include "g2o/core/optimization_algorithm_levenberg.h"
#include "g2o/core/robust_kernel_impl.h"

namespace CloudReg
{

typedef int32_t VertexID_G2O_t;

#define INVAILD_G2O_VERTEXID 0

class TransformOptimize
{
public:
    TransformOptimize(const std::string& name, const std::string& logStr)
    	: name_(name)
		, logStr_(logStr)
		, optimizer_() 
	{
	}

    ~TransformOptimize()
    {
		clear();
    }

    bool run(std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &vecCloudPtr,
			CADModel &cadModel);

private:

	void uniformSampling(double radius, 
                    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, 
                    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_filtered);

	bool downSampling(std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &vecCloudPtr);

	bool convertToPclCloud(CADModel &cadModel, 
					std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &model_vec);

	void planeFitting(double distTh, 
					pcl::PointCloud<pcl::PointXYZ>::Ptr ground, 
					Eigen::VectorXf &coeff, std::vector<int> &inlierIdxs);

	bool getModelPlaneCoeff(
					std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &model_vec,
					std::vector<Eigen::Vector4d> &modelPlanes);

	double calcCloudToPLaneAveDist(Eigen::Vector4d &plane,
                                pcl::PointCloud<pcl::PointXYZ>::Ptr cloud);

	bool matchCloudToMode(
                    std::vector<Eigen::Vector4d> &modePlanes,
                    std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &vecCloudPtr);

	bool optimize(std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &vecCloudPtr,
					std::vector<Eigen::Vector4d> &modelPlanes,
					Eigen::Matrix4d &transform);

    bool addWallPointToModelPlaneEdges(
                    std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &vecCloudPtr,
                    std::vector<Eigen::Vector4d> &modelPlanes,
					Eigen::Matrix4d &transform);

    bool getSE3Transfor(Eigen::Matrix4d &transform);

	bool transformCloud(
					std::vector<Eigen::Vector4d> &modePlanes,
					std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &vecCloudPtr,
					Eigen::Matrix4d &finalT);

	bool viewModelAndChangedCloud(
					std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &model_vec,
					std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &vecCloudPtr);

	std::string convertToSimpleIDStr(uint64_t id)
	{
		std::stringstream ss;
		ss << std::setfill('0') << std::setw(4) << std::to_string((id % 10000));
		return ss.str();
	}

private:

	void clear()
	{
        optimizer_.clear();
		mapValue2Id_.clear();
		mapEdge_.clear();	
	}

	void initSolver()
	{
#if 0
		g2o::BlockSolverX::LinearSolverType * linearSolver = 
                        new g2o::LinearSolverEigen<g2o::BlockSolverX::PoseMatrixType>();
        g2o::BlockSolverX * blockSolver = new g2o::BlockSolverX(linearSolver);
        g2o::OptimizationAlgorithmLevenberg* algorithm 
				= new g2o::OptimizationAlgorithmLevenberg(blockSolver);
#endif

#if 1
		auto linearSolver = g2o::make_unique<g2o::LinearSolverEigen<g2o::BlockSolverX::PoseMatrixType>>();
		auto blockSolver = g2o::make_unique<g2o::BlockSolverX>(std::move(linearSolver));
		g2o::OptimizationAlgorithmLevenberg* algorithm 
				= new g2o::OptimizationAlgorithmLevenberg(std::move(blockSolver));
#endif
		optimizer_.setAlgorithm(algorithm);
    }

	//be carefull ,the id must begin with 1
    const VertexID_G2O_t getFreeID() const
    {
    	return static_cast<const VertexID_G2O_t>(mapValue2Id_.size() + 1);
    }

	template<class T>
	const VertexID_G2O_t createG2OID(const T* pValue)
	{
		if (nullptr == pValue)
			return INVAILD_G2O_VERTEXID;

		uint64_t key = (uint64_t)pValue;

		auto iter = mapValue2Id_.find(key);
		if (iter != mapValue2Id_.end())
			return INVAILD_G2O_VERTEXID;

		VertexID_G2O_t id = getFreeID();
		mapValue2Id_[key] = id;
		return id;

	}

	template<class T>
	const VertexID_G2O_t getVertexID(const T* pValue) const
	{
		if (nullptr == pValue)
			return INVAILD_G2O_VERTEXID;

		uint64_t key = (uint64_t)pValue;
		auto iter = mapValue2Id_.find(key);

		return (iter != mapValue2Id_.end() ? iter->second : INVAILD_G2O_VERTEXID);
	}

	template<class T>
	VertexID_G2O_t addSE3Vertex(const T* pValue, const g2o::SE3Quat& estimate, bool bMarginalized, bool bFix)
	{
		const VertexID_G2O_t id = createG2OID(pValue);
		if (INVAILD_G2O_VERTEXID == id)
		{
			LOG(INFO) << "createG2OID Failed: " << pValue;
			return id;
		}

		g2o::VertexSE3Expmap *pVertex = new g2o::VertexSE3Expmap();
    	pVertex->setEstimate(estimate);
    	pVertex->setMarginalized(bMarginalized);
    	pVertex->setId(id);
    	pVertex->setFixed(bFix);
    	optimizer_.addVertex(pVertex);
		return id;
		
	}

	template<class T>
	void addEdge(const std::string& name, T* e)
	{
		optimizer_.addEdge(e);
		mapEdge_[name].insert(e);
	}

    template<class T>
    double computeSumError(std::set<T*> &edges)
    {
        double sumErr = 0.0;
        for (auto &e : edges)
        {
            e->computeError();
            double chi2 = e->chi2();
            sumErr += chi2;
        }
        return sumErr;
    }

	double statisticsError(const std::string& prex = "", bool bPrint = false)
	{
    	std::stringstream ss;
    	ss << std::fixed << std::setprecision(4)
        	<< "########### statistics opt "<< prex  << " " << name_ << " ##########\n";

		double errSum = 0.0f;
		for(auto& value : mapEdge_)
		{
			std::string edgeName = value.first;
			std::size_t edgeNum = value.second.size();
			double err = computeSumError(value.second);

			errSum += err;
			ss << "------ " << edgeName << " : " << err << "/" << edgeNum <<"/" <<(double)err/edgeNum << "\n"; 
		}

		ss << "vertices = " << optimizer_.vertices().size() << " edges = " << optimizer_.edges().size()
			<< " sumErr = " << errSum <<"\n";
		ss << "########################################################################################" << "\n";

		if(bPrint)
			LOG(INFO) << ss.str();

		return errSum;
	}

    int optData(uint64_t maxIterations, bool bVerbose, bool bStatistics = false)
    {
        if (optimizer_.vertices().empty() || optimizer_.edges().empty())
        {
            return false;
        }

		if (bStatistics) statisticsError("before", bStatistics);

        optimizer_.setVerbose(bVerbose);
        optimizer_.initializeOptimization();
        int optTime = optimizer_.optimize(maxIterations);

		if (bStatistics) statisticsError("after", bStatistics);

        return optTime;
    }

private:
    std::string name_;
	std::string logStr_;
    
	g2o::SparseOptimizer optimizer_;

    //added Vertex
	std::map<uint64_t, VertexID_G2O_t>  mapValue2Id_;

	std::map<std::string, std::set<g2o::OptimizableGraph::Edge*>> mapEdge_;

};

} //namespace CloudReg

#endif // TRANSFORM_OPTIMIZE_H
