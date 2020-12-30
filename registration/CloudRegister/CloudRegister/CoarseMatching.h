#pragma once

#include "BaseType.h"
#include "CADModel.h"
#include "GeometryUtils.h"

namespace CloudReg {

class CoarseMatching {
	struct Segment {
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

		Eigen::Vector3f s_, e_;

		inline Eigen::Vector3f mid() const { return (s_ + e_) * 0.5f; }
		inline Eigen::Vector3f dir() const { return (e_ - s_).normalized(); }
		inline float len() const { return (e_ - s_).norm(); }
	};

public:
	CoarseMatching();

	struct MatchResult {
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

		bool isValid() const { return !walls_.empty() && walls_.size() == matchIndices_.size(); }

		std::vector<PointCloud::Ptr> getAllPieces() const;

		std::string toString() const;

		// refined & transformed clouds
		std::vector<PointCloud::Ptr> walls_;
		PointCloud::Ptr floor_, roof_;

		// the index of wall in CADModel which the wall matched to
		// #matchIndices_ = #walls_ 
		// walls_[i] -> cadModel.getTypedModelItems(ITEM_WALL_E)[matchIndices_[i]]
		std::vector<std::size_t> matchIndices_;

		Eigen::Matrix4f T_; //! note that the clouds are already transformed by T_;
	};

	MatchResult run(const std::vector<PointCloud::Ptr>& allPieces, const CADModel& cadModel);

	std::vector<PointCloud::Ptr> testSegment(PointCloud::Ptr cloud, const CADModel& cadModel);

private:
	static constexpr double PLANE_REFINE_DISTANCE = 0.05;
	static constexpr double LINE_DETECT_DIS_THRESH = 0.02;
	static constexpr double LINE_DETECT_CONNECT_THRESH = 0.05;

	static constexpr float OUTLINE_LINK_DIS_THRESH = 0.15f;
	static constexpr float OUTLINE_LINK_DIS_THRESH_SQUARED = OUTLINE_LINK_DIS_THRESH * OUTLINE_LINK_DIS_THRESH;

	static constexpr float BLUEPRINT_MATCH_DIS_THRESH = 0.25f;

	std::pair<PointCloud::Ptr, Eigen::Vector4f> refinePlanePattern(PointCloud::Ptr cloud, double disthresh) const;

	// we already get the plane parameters, use as initialization, 
	Segment detectMainSegementXoY(PointCloud::Ptr cloud, const Eigen::Vector4f& plane_param, double disthresh, double connect_thresh) const;

	Eigen::vector<trans2d::Matrix2x3f> computeOutlineTransformCandidates(const Eigen::vector<Eigen::Vector2f>& blueprint,
		const Eigen::vector<Eigen::Vector2f>& outline, float disthresh);

	// wallparams: [Vector6f(p, n)] 2d line params
	trans2d::Matrix2x3f chooseTransformByHoles(const CADModel& cadModel, const std::vector<PointCloud::Ptr>& walls,
		const std::vector<Eigen::VectorXf>& wallparams, const Eigen::vector<trans2d::Matrix2x3f>& Ts, const float floorZ) const;

	std::vector<Segment> detectSegmentsXOY(PointCloud::Ptr cloud, float linethresh, float connect_thresh, float min_seg_length) const;
	// line_params, v6: [p, n]
	std::vector<Segment> splitSegments(PointCloud::Ptr cloud, Eigen::VectorXf line_params, float connect_thresh, float min_seg_length) const;

	// get SORTED intersection points on cad
	Eigen::vector<Eigen::Vector3d> intersectCADModelOnZ(const CADModel& cadModel, float z) const;
};

}
