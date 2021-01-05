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
		std::vector<PlaneCloud> walls_;
		std::vector<PlaneCloud> beams_;
		PlaneCloud roof_, floor_;

		bool valid() const {
			return !walls_.empty() && !roof_.cloud_ && !floor_.cloud_;
		}

		//! careful 
		std::vector<PlaneCloud> allPlanes() const {
			std::vector<PlaneCloud> all;
			all.reserve(walls_.size() + beams_.size() + 2);
			all.insert(all.end(), walls_.begin(), walls_.end());
			all.insert(all.end(), beams_.begin(), beams_.end());
			all.push_back(roof_);
			all.push_back(floor_);

			return all;
		}
	};

	SegmentResult run();

private:
	struct Segment {
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

		Eigen::Vector3f s_, e_;
		inline Eigen::Vector3f mid() const { return (s_ + e_) * 0.5f; }
		inline Eigen::Vector3f dir() const { return (e_ - s_).normalized(); }
		inline float len() const { return (e_ - s_).norm(); }
	};

	PointCloud::Ptr orgCloud_;
	const CADModel& cadModel_;

	// cad aabb box from model
	double cadx1_{ METER_100 }, cadx2_{ -METER_100 };
	double cady1_{ METER_100 }, cady2_{ -METER_100 };
	double cadz1_{ METER_100 }, cadz2_{ -METER_100 };

	double zfloor_, zroof_; // org cloud

	PointCloud::Ptr sparsedCloud_{ nullptr }; // most time we work on this.

	// constants
	static constexpr double METER_100 = 100.;
	static constexpr float DOWNSAMPLE_SIZE = 0.01f;

	// main processes
	void recordModelBoundingBox();
	bool alignCloudToCADModel();
	SegmentResult segmentByCADModel();

	PointCloud::Ptr sliceMainBody();
	double detectRoofHeight(PointCloud::Ptr cloud) const;

	std::vector<PlaneCloud> detectPlanes(PointCloud::Ptr cloud,
		double disthresh, std::size_t inlier_count_thresh, std::size_t countthresh = 10000) const;
	std::vector<PlaneCloud> detectRegionPlanes(PointCloud::Ptr cloud, double anglediff, double curvediff, std::size_t min_points) const;
	void detectPlanesRecursively(PointCloud::Ptr cloud, std::vector<PlaneCloud>& planes,
		double disthresh, std::size_t inlier_count_thresh, std::size_t countthresh) const;

	// line: p, n, 
	// note that the cloud no need to be an "actual segment", either the detected segment cloud.
	std::vector<Segment> splitSegments(PointCloud::Ptr cloud, const Eigen::Vector3f p, const Eigen::Vector3f& n,
		float connect_thresh, float min_seg_length) const;

	// 2d, segments & blueprint shall be sorted.
	Eigen::vector<trans2d::Matrix2x3f> computeSegmentAlignCandidates(const std::vector<Segment>& segments,
		const Eigen::vector<Eigen::Vector2f>& blueprint, float disthresh) const;

	SegmentResult segmentCloudByCADModel(PointCloud::Ptr cloud) const; // after the cloud was aligned properly.

	inline PointCloud::Ptr sparsedCloud();
};

}
