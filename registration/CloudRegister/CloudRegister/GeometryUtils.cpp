#include "GeometryUtils.h"

#include <pcl/common/transforms.h>
#include <pcl/filters/passthrough.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/filters/radius_outlier_removal.h>
#ifdef UBUNTU_SWITCH
#include <pcl/keypoints/impl/uniform_sampling.hpp>
#else
#include <pcl/filters/uniform_sampling.h>
#endif
#include <pcl/kdtree/kdtree.h>
#include <pcl/segmentation/extract_clusters.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/ransac.h>
#include <pcl/sample_consensus/sac_model_line.h>
#include <pcl/sample_consensus/sac_model_plane.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/filters/conditional_removal.h>
#include "funHelper.h"
#include "glog/logging.h"

// trans2d
#define RBLOCK(T) T.block<2, 2>(0, 0)
#define tBLOCK(T) T.block<2, 1>(0, 2)

namespace trans2d {

Matrix2x3f buildTransform(float theta, const Eigen::Vector2f& t) {
	float c{ std::cos(theta) }, s{ std::sin(theta) };

	Matrix2x3f T;
	RBLOCK(T) << c, -s, s, c;
	tBLOCK(T) = t;

	return T;
}

Matrix2x3f inverse(const Matrix2x3f& T) {
	Matrix2x3f invT;
	RBLOCK(invT) = RBLOCK(T).transpose();
	tBLOCK(invT) = -RBLOCK(invT) * tBLOCK(T);
	return invT;
}

Eigen::Vector2f transform(const Eigen::Vector2f& p, const Matrix2x3f& T) {
	return RBLOCK(T) * p + tBLOCK(T);
}

// return T where T* s2 = s1, T* e2 = e1
Matrix2x3f estimateTransform(const Eigen::Vector2f& s1, const Eigen::Vector2f& e1,
	const Eigen::Vector2f& s2, const Eigen::Vector2f& e2) {
	Matrix2x3f T;

	Eigen::Vector2f n1 = (e1 - s1).normalized();
	Eigen::Vector2f n2 = (e2 - s2).normalized();
	float c = n1.dot(n2);
	float s = -(n1[0] * n2[1] - n1[1] * n2[0]);
	T(0, 0) = c;
	T(0, 1) = -s;
	T(1, 0) = s;
	T(1, 1) = c;

	Eigen::Vector2f cen1 = (s1 + e1) / 2.0f;
	Eigen::Vector2f cen2 = (s2 + e2) / 2.0f;

	tBLOCK(T) = cen1 - RBLOCK(T) * cen2;

	return T;
}

// simply shift to 3d
Eigen::Matrix4f asTransform3d(const Matrix2x3f& T2d) {
	Eigen::Matrix4f T3d = Eigen::Matrix4f::Identity();
	T3d.block<2, 2>(0, 0) = RBLOCK(T2d);
	T3d.block<2, 1>(0, 3) = tBLOCK(T2d);
	return T3d;
}

#undef RBLOCK
#undef tBLOCK

}

namespace geo {

float distance_to_segment_2d(const Eigen::Vector2f& p, const Eigen::Vector2f& s, const Eigen::Vector2f& e) {
	Eigen::Vector2f sp = p - s;
	Eigen::Vector2f se = e - s;
	float se2 = se.squaredNorm();
	if (se2 < 1e-8) return sp.norm();

	float t = sp.dot(se) / se2;
	if (t < 0.f) return sp.norm();
	if (t > 1.f) return (p - e).norm();

	return std::fabs(sp(0) * se(1) - sp(1) * se(0)) / std::sqrt(se2);
}

std::vector<std::size_t> sort_points_counter_clockwise(const Eigen::vector<Eigen::Vector2f>& points) {
#ifdef UBUNTU_SWITCH
	auto cen_theta = ll::mapf([](const Eigen::Vector2f& p) { return (float)std::atan2(double(p[1]), double(p[0])); }, points);
#else
	auto cen_theta = ll::mapf([](const Eigen::Vector2f& p) { return std::atan2f(p[1], p[0]); }, points);
#endif

	auto indices = ll::range(cen_theta.size());
	std::sort(indices.begin(), indices.end(), [&cen_theta](std::size_t i, std::size_t j) { return cen_theta[i] < cen_theta[j]; });

	return indices;
}

SegmentIntersectInfo zIntersectSegment(const Eigen::Vector3d& a, const Eigen::Vector3d& b, double z) {
	SegmentIntersectInfo sii;

	double dz = b(2, 0) - a(2, 0);
	if (std::fabs(dz) < 1e-6) return sii;

	sii.lambda_ = (z - a(2, 0)) / dz;
	if (sii.valid()) sii.point_ = a * (1 - sii.lambda_) + b * sii.lambda_;

	return sii;
}

Eigen::vector<Eigen::Vector3d> intersectContour(const Eigen::Vector3d& s, const Eigen::Vector3d& e,
	const Eigen::vector<Eigen::Vector3d>& points) {
	Eigen::Vector3d se = e - s;

	Eigen::vector<Eigen::Vector3d> intersections;
	for (std::size_t i = 0; i < points.size(); ++i) {
		Eigen::Vector3d s2 = points[i];
		Eigen::Vector3d e2 = points[(i + 1) % points.size()];
		Eigen::Vector3d s2e2 = e2 - s2;

		Eigen::Vector3d c1 = (s2 - s).cross(se);
		Eigen::Vector3d c2 = s2e2.cross(se);
		double t = -c2.transpose().dot(c1) / c2.squaredNorm();

		if (t > 0. && t < 1.) intersections.emplace_back(s2 + s2e2 * t);
	}

	return intersections;
}

bool isInShape2D(const Eigen::Vector3d& point, const Eigen::vector<Eigen::Vector3d>& outline, const std::vector<Eigen::vector<Eigen::Vector3d> >& holes) {
	Eigen::Vector3d minp = Eigen::Vector3d::Ones() * 100.; //! 
	Eigen::Vector3d maxp = Eigen::Vector3d::Ones() * -100.;
	for (const Eigen::Vector3d& p : outline) {
		for (int i = 0; i < 3; ++i) {
			if (p(i) < minp(i)) minp(i) = p(i);
			if (p(i) > maxp(i)) maxp(i) = p(i);
		}
	}

	// this is not very readable, may refactor someday. or, just write more duplicate codes...
	int dim1{ 0 }, dim2{ 1 }, dimfix{ 2 };
	Eigen::Vector3d dp = maxp - minp;
	if (dp(0) < 1e-6) std::swap(dimfix, dim1);
	else if (dp(1) < 1e-6) std::swap(dimfix, dim2);

	//todo: should use margin check instead.
	for (int i = 0; i < 3; ++i)
		if (i != dimfix) 
			if (point[i] > maxp[i] || point[i] < minp[i]) return false;

	// we check dim1 & dim2
	Eigen::Vector3d s = point;
	Eigen::Vector3d e = Eigen::Vector3d::Zero();
	e[dim1] = 1.; // just need to pass the point
	e += point;

	Eigen::vector<Eigen::Vector3d> intersections = intersectContour(s, e, outline);
	if (intersections.empty()) return false;

	for (const auto& hole : holes) {
		Eigen::vector<Eigen::Vector3d> ins = intersectContour(s, e, hole);
		intersections.insert(intersections.end(), ins.begin(), ins.end());
	}

	if (!intersections.empty() && intersections.size() % 2 != 0) {
		LOG(INFO) << "intersections issue, cannot judge, intersection size: " << intersections.size();
		return false;
	}

	Eigen::Vector3d se = e - s;
	std::size_t left = std::count_if(intersections.begin(), intersections.end(), [&s, &se](const Eigen::Vector3d& ip) {
		return (ip - s).dot(se) < 0.;
	});

	return left % 2 == 1;
}

// cloud

PointCloud::Ptr passThrough(PointCloud::Ptr cloud, const std::string& field, float low, float high) {
	pcl::PassThrough<Point> filter;
	filter.setInputCloud(cloud);
	filter.setFilterFieldName(field);
	filter.setFilterLimits(low, high);

	PointCloud::Ptr filtered(new PointCloud());
	filter.filter(*filtered);
	return filtered;
}

PointCloud::Ptr downsampleUniformly(PointCloud::Ptr cloud, float radius) {
	PointCloud::Ptr filtered(new PointCloud());

	pcl::UniformSampling<Point> filter;
	filter.setInputCloud(cloud);
	filter.setRadiusSearch(radius);
#ifdef UBUNTU_SWITCH
	pcl::PointCloud<int> keypointIndices;
	filter.compute(keypointIndices);
	pcl::copyPointCloud(*cloud, keypointIndices.points, *filtered);
#else
	filter.filter(*filtered);
#endif

	return filtered;
}

PointCloud::Ptr getSubSet(PointCloud::Ptr cloud, const std::vector<int>& indices, bool negative) {
	PointCloud::Ptr sub(new PointCloud());

	pcl::PointIndices::Ptr pi(new pcl::PointIndices());
	pi->indices = indices;

	pcl::ExtractIndices<Point> ei;
	ei.setInputCloud(cloud);
	ei.setIndices(pi);
	ei.setNegative(negative);
	ei.filter(*sub);

	return sub;
}

PointCloud::Ptr clusterMainStructure(PointCloud::Ptr cloud, float distance) {
	pcl::search::KdTree<Point>::Ptr tree(new pcl::search::KdTree<Point>);
	tree->setInputCloud(cloud);

	std::vector<pcl::PointIndices> clusters;
	pcl::EuclideanClusterExtraction<Point> ece;
	ece.setClusterTolerance(distance);
	ece.setSearchMethod(tree);
	ece.setInputCloud(cloud);
	ece.extract(clusters);

	// cluster is already sorted, just use the first one.
	PointCloud::Ptr main(new PointCloud());
	pcl::copyPointCloud(*cloud, clusters.front().indices, *main);

	return main;
}

std::pair<std::vector<int>, Eigen::VectorXf> detectOneLineRansac(PointCloud::Ptr cloud, float disthresh) {
	pcl::SampleConsensusModelLine<Point>::Ptr lineModel(new pcl::SampleConsensusModelLine<Point>(cloud));
	pcl::RandomSampleConsensus<Point> ransac(lineModel);
	ransac.setDistanceThreshold(disthresh);
	ransac.computeModel();

	std::pair<std::vector<int>, Eigen::VectorXf> re;
	ransac.getInliers(re.first);
	ransac.getModelCoefficients(re.second);

	return re;
}

std::pair<PointCloud::Ptr, Eigen::Vector4f> refinePlanePattern(PointCloud::Ptr cloud, double disthresh) {
	pcl::SampleConsensusModelPlane<Point>::Ptr planeModel(new pcl::SampleConsensusModelPlane<Point>(cloud));
	pcl::RandomSampleConsensus<Point> ransac(planeModel);
	ransac.setDistanceThreshold(disthresh);
	ransac.computeModel();

	std::vector<int> inliers;
	Eigen::VectorXf params;
	ransac.getInliers(inliers);
	ransac.getModelCoefficients(params);

	auto plane = geo::getSubSet(cloud, inliers);
	Eigen::Vector4f abcd = CloudReg::calcPlaneParam(plane).cast<float>();

	// filter again.
	plane = geo::filterPoints(cloud, [&abcd, disthresh](const Point& p) {
		double dis = abcd(0) * p.x + abcd(1) * p.y + abcd(2) * p.z + abcd(3);
		return std::fabs(dis) <= disthresh;
	});

	return std::make_pair(plane, abcd);
}

PointCloud::Ptr transfromPointCloud(PointCloud::Ptr cloud, const Eigen::Matrix4f& T) {
	PointCloud::Ptr trans(new PointCloud());
	pcl::transformPointCloud<Point>(*cloud, *trans, T);
	return trans;
}

PointCloud::Ptr mapPoints(PointCloud::Ptr cloud, std::function<Point (const Point & p)> mapfunc) {
	PointCloud::Ptr newcloud(new PointCloud());
	newcloud->reserve(cloud->size());
	for (const auto& p : *cloud) newcloud->push_back(mapfunc(p));
	newcloud->width = newcloud->size();
	newcloud->height = 1;

	return newcloud;
}

class FunctorCondition : public pcl::ConditionBase<Point> {
public:
	using FilterFunc = std::function<bool (const Point&)>;
	FunctorCondition(FilterFunc filterfunc) :
		pcl::ConditionBase<Point>(), func_(filterfunc) {}

	bool evaluate(const Point& point) const override {
		return func_(point);
	}

private:
	FilterFunc func_;
};

PointCloud::Ptr filterPoints(PointCloud::Ptr cloud, std::function<bool (const Point&)> evafunc, bool negative) {
	pcl::ConditionalRemoval<Point> filter;
	{
		FunctorCondition::Ptr cond(new FunctorCondition(evafunc));
		filter.setCondition(cond);
	}
	filter.setInputCloud(cloud);

	PointCloud::Ptr filtered_cloud(new PointCloud());
	filter.filter(*filtered_cloud);
	filter.setKeepOrganized(!negative);

	return filtered_cloud;
}

}
