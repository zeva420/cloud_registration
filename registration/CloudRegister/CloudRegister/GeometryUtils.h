#pragma once
// geometry utilities

#include <pcl/common/common.h>
#include <pcl/point_types.h>

#include "ll.hpp"
#include "BaseType.h"

inline Point operator+(const Point& p, const Point& q) { return Point(p.x + q.x, p.y + q.y, p.z + q.z); }
inline Point operator-(const Point& p, const Point& q) { return Point(p.x - q.x, p.y - q.y, p.z - q.z); }
inline Point operator*(const Point& p, float a) { return Point(p.x * a, p.y * a, p.z * a); }
inline Point operator/(const Point& p, float a) { return p * (1.0f / a); }

// 2d transformation
namespace trans2d {
using Matrix2x3f = Eigen::Matrix<float, 2, 3>;

extern Matrix2x3f buildTransform(float theta, const Eigen::Vector2f& t);

extern Matrix2x3f inverse(const Matrix2x3f& T);

extern Eigen::Vector2f transform(const Eigen::Vector2f& p, const Matrix2x3f& T);

// return T where T* s2 = s1, T* e2 = e1
extern Matrix2x3f estimateTransform(const Eigen::Vector2f& s1, const Eigen::Vector2f& e1,
	const Eigen::Vector2f& s2, const Eigen::Vector2f& e2);

// T* src = dst, #src == #dst
extern Matrix2x3f estimateTransform(const Eigen::vector<Eigen::Vector2f>& src, const Eigen::vector<Eigen::Vector2f>& dst);

// shift to 3d
extern Eigen::Matrix4f asTransform3d(const Matrix2x3f& T2d);
};

namespace geo {
constexpr float PI = 3.14159265f;
constexpr float TAU = PI * 2.f;

/* point operations */
inline double random(double max = 1.0) { return std::rand() / static_cast<double>(RAND_MAX)* max; }

inline float dot(const Point& p, const Point& q) { return p.x * q.x + p.y * q.y + p.z * q.z; }

inline Point cross(const Point& p, const Point& q) { return Point(p.y * q.z - p.z * q.y, p.z * q.x - p.x * q.z, p.x * q.y - p.y * q.x); }

inline float length_squared(const Point& p) { return dot(p, p); }

inline float length(const Point& p) { return std::sqrt(length_squared(p)); }

inline Point P_(const Eigen::Vector3f& v) { return Point(v[0], v[1], v[2]); }
inline Eigen::Vector3f V_(const Point& p) { Eigen::Vector3f v; v << p.x, p.y, p.z; return v; }

// 2d
inline float length2d_squared(const Point& p) { return p.x * p.x + p.y * p.y; }

inline float length2d(const Point& p) { return std::sqrt(length2d_squared(p)); }

inline float cross2d(const Point& p, const Point& q) { return p.x * q.y - p.y * q.x; }

extern float distance_to_segment_2d(const Eigen::Vector2f& p, const Eigen::Vector2f& s, const Eigen::Vector2f& e);

extern std::vector<std::size_t> sort_points_counter_clockwise(const Eigen::vector<Eigen::Vector2f>& points);

// simple geometric
struct SegmentIntersectInfo {
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	bool valid() const { return lambda_ >= 0. && lambda_ <= 1.; }

	double lambda_{ -1. };
	Eigen::Vector3d point_;
};
extern SegmentIntersectInfo zIntersectSegment(const Eigen::Vector3d& a, const Eigen::Vector3d& b, double z);

// see CADModel::InterpolateShape
// warn, no error check.
extern bool isInShape2D(const Eigen::Vector3d& point, const Eigen::vector<Eigen::Vector3d>& outline, const std::vector<Eigen::vector<Eigen::Vector3d> >& holes);

/* cloud operations */

extern PointCloud::Ptr passThrough(PointCloud::Ptr cloud, const std::string& field, float low, float high);

extern PointCloud::Ptr downsampleUniformly(PointCloud::Ptr cloud, float radius);

// set negative to TRUE to inverse grab.
extern PointCloud::Ptr getSubSet(PointCloud::Ptr cloud, const std::vector<int>& indices, bool negative = false);

extern PointCloud::Ptr clusterMainStructure(PointCloud::Ptr cloud, float distance);

//note: this method does NOT ensure to get the 'main line'.
//return: indices of inliers, [p, n] of the line (dim = 6).
extern std::pair<std::vector<int>, Eigen::VectorXf> detectOneLineRansac(PointCloud::Ptr cloud, float disthresh);

// note the cloud should looks like a plane.
extern std::pair<PointCloud::Ptr, Eigen::Vector4f> refinePlanePattern(PointCloud::Ptr cloud, double disthresh);

extern PointCloud::Ptr transfromPointCloud(PointCloud::Ptr cloud, const Eigen::Matrix4f& T);

extern PointCloud::Ptr mapPoints(PointCloud::Ptr cloud, std::function<Point (const Point & p)> mapfunc);

extern PointCloud::Ptr filterPoints(PointCloud::Ptr cloud, std::function<bool (const Point&)> evafunc, bool negative = false);
};

