#pragma once

#include "BaseType.h"
#include "CADModel.h"

namespace CloudReg {

class CoarseMatching {
	struct Segment {
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

		Eigen::Vector3f s_, e_;

		inline Eigen::Vector3f mid() const { return (s_ + e_) * 0.5f; }
		inline Eigen::Vector3f dir() const { return (e_ - s_).normalized(); }
	};

public:
	CoarseMatching();

	bool run(const std::vector<PointCloud::Ptr>& allPieces, const CADModel& cadModel);

private:
	static constexpr double PLANE_REFINE_DISTANCE = 0.02;
	static constexpr double LINE_DETECT_DIS_THRESH = 0.02;
	static constexpr double LINE_DETECT_CONNECT_THRESH = 0.05;

	static constexpr float OUTLINE_LINK_DIS_THRESH = 0.15f;
	static constexpr float OUTLINE_LINK_DIS_THRESH_SQUARED = OUTLINE_LINK_DIS_THRESH * OUTLINE_LINK_DIS_THRESH;

	static constexpr float BLUEPRINT_MATCH_DIS_THRESH = 0.25f;

	std::pair<PointCloud::Ptr, Eigen::Vector4f> refinePlanePattern(PointCloud::Ptr cloud, double disthresh) const;

	// we already get the plane parameters, use as initialization
	Segment detectSegementXoY(PointCloud::Ptr cloud, const Eigen::Vector4f& plane_param, double disthresh, double connect_thresh) const;

	Eigen::Matrix4f computeOutlineTransformation(const Eigen::vector<Eigen::Vector2f>& blueprint,
		const Eigen::vector<Eigen::Vector2f>& outline, float disthresh);
};

}
