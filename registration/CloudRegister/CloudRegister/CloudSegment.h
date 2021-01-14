#pragma once

#include "BaseType.h"
#include "CADModel.h"
#include "GeometryUtils.h"

namespace CloudReg {

class CloudSegment {
public:
	CloudSegment(PointCloud::Ptr orgcloud, const CADModel& model) :
		orgCloud_(orgcloud), cadModel_(model) {}

	struct PlaneCloud {
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

		PointCloud::Ptr cloud_{ nullptr };
		Eigen::Vector4f abcd_;

		Eigen::Vector3f n() const { return abcd_.block<3, 1>(0, 0); }
	};

	struct SegmentResult {
		explicit SegmentResult(const Eigen::Matrix4f& T): T_(T) {} // in case i forget to assign T_..., todo: fix this.

		std::map<ModelItemType, std::vector<PlaneCloud> > clouds_;

		Eigen::Matrix4f T_; // T* org = cur, note we already transformed after segment.

		bool valid() const { return !clouds_.empty(); }

		std::vector<PlaneCloud> allPlanes() const {
			std::vector<PlaneCloud> all;
			for (const auto& pr : clouds_)
				all.insert(all.end(), pr.second.begin(), pr.second.end());
			return all;
		}

		Eigen::Vector3f originalCenter() const { return -T_.block<3, 3>(0, 0).transpose()* T_.block<3, 1>(0, 3); }

		std::string to_string() const;
	};

	SegmentResult run();

private:
	struct Segment {
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

		Eigen::Vector3f s_, e_;
		inline Eigen::Vector3f mid() const { return (s_ + e_) * 0.5f; }
		inline Eigen::Vector3f dir() const { return (e_ - s_).normalized(); }
		inline float len() const { return (e_ - s_).norm(); }

		// 2d logic
		void beCounterClockwise(const Eigen::Vector3f& cen);
	};

	PointCloud::Ptr orgCloud_;
	const CADModel& cadModel_;

	// cad aabb box from model
	double cadx1_{ METER_100 }, cadx2_{ -METER_100 };
	double cady1_{ METER_100 }, cady2_{ -METER_100 };
	double cadz1_{ METER_100 }, cadz2_{ -METER_100 };

	double zfloor_, zroof_; // org cloud

	PointCloud::Ptr sparsedCloud_{ nullptr }; // most time we work on this.

	Eigen::Matrix4f T_; // T* org = cur, note we already transformed after segment.

	// constants
	static constexpr double METER_100 = 100.;
	static constexpr float DOWNSAMPLE_SIZE = 0.01f;

	// main processes
	bool calibrateDirectionToAxisZ();
	void recordModelBoundingBox();
	bool alignCloudToCADModel();
	SegmentResult segmentByCADModel();
	void refineSegmentResult(SegmentResult& sr) const;

	PointCloud::Ptr sliceMainBody();
	double detectRoofHeight(PointCloud::Ptr cloud) const;

	std::vector<PlaneCloud> detectPlanes(PointCloud::Ptr cloud,
		double disthresh, std::size_t inlier_count_thresh, std::size_t countthresh = 10000) const;
	std::vector<PlaneCloud> detectRegionPlanes(PointCloud::Ptr cloud, double anglediff, double curvediff, std::size_t min_points) const;
	void detectPlanesRecursively(PointCloud::Ptr cloud, std::vector<PlaneCloud>& planes,
		double disthresh, std::size_t inlier_count_thresh, std::size_t countthresh) const;

	trans2d::Matrix2x3f chooseTransformByHoles(const Eigen::vector<trans2d::Matrix2x3f>& Ts, const std::vector<PlaneCloud>& walls) const;

	// get SORTED intersection points on cad
	Eigen::vector<Eigen::Vector3d> intersectCADModelOnZ(const CADModel& cadModel, float z) const;

	// line: p, n, 
	// note that the cloud no need to be an "actual segment", either the detected segment cloud.
	std::vector<Segment> splitSegments(PointCloud::Ptr cloud, const Eigen::Vector3f p, const Eigen::Vector3f& n,
		float connect_thresh, float min_seg_length) const;

	// 2d, segments & blueprint shall be sorted.
	Eigen::vector<trans2d::Matrix2x3f> computeSegmentAlignCandidates(const std::vector<Segment>& segments,
		const std::vector<Segment>& blueprint, float disthresh) const;

	SegmentResult segmentCloudByCADModel(PointCloud::Ptr cloud) const; // after the cloud was aligned properly.

	inline PointCloud::Ptr sparsedCloud();

	void removeFarPoints(PointCloud::Ptr inputCloud, const CADModel& cad);

	bool statisticsForPointZ(float binSizeTh, PointCloud::Ptr cloud,
		std::vector<std::pair<int, std::vector<int>>>& zToNumVec);

	// debug
	void _show_result(const SegmentResult& sr) const;
};

}
