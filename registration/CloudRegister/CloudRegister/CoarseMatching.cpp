#include "CoarseMatching.h"
#include "glog/logging.h"

#include "GeometryUtils.h"

#include <pcl/sample_consensus/model_types.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/ransac.h>
#include <pcl/sample_consensus/sac_model_line.h>
#include <pcl/sample_consensus/sac_model_plane.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/ModelCoefficients.h>
#include <pcl/visualization/pcl_visualizer.h>

namespace CloudReg {

CoarseMatching::CoarseMatching() {}

bool CoarseMatching::run(std::vector<PointCloud::Ptr>& allPieces, const CADModel& cadModel) {

	const auto& holes = cadModel.getTypedModelItems(ITEM_HOLE_E);
	for (const auto& hole : holes) LOG(INFO) << hole.toString();

	// return true;

	constexpr float HOR_CHECK_EPS = 0.01f;

	LOG(INFO) << "coarse match " << allPieces.size() << " walls.";

	//1. basically, we use walls to do coarse matching.
	std::vector<std::pair<PointCloud::Ptr, Segment>> wallPieces;
	wallPieces.reserve(allPieces.size()); // may -2
	PointCloud::Ptr floor, roof; // still we need floor and roof to get the proper Y
	float zFloor, zRoof;

	for (auto pr : ll::enumerate(allPieces)) {
		auto cloud = *pr.iter;

		PointCloud::Ptr plane;
		Eigen::Vector4f params;
		std::tie(plane, params) = refinePlanePattern(cloud, PLANE_REFINE_DISTANCE);

		float nz = params[2];
		LOG(INFO) << pr.index << "] nz: " << nz;
		if (std::fabs(nz) < HOR_CHECK_EPS) {
			// we got a wall piece
			auto segment = detectSegementXoY(plane, params, LINE_DETECT_DIS_THRESH, LINE_DETECT_CONNECT_THRESH);
			wallPieces.emplace_back(plane, segment);
		} else {
			float z = -params[3] / nz;
			if (z > 0.f) { roof = plane; zRoof = z; } else { floor = plane; zFloor = z; }
		}
	}

	// LOG(INFO) << wallPieces.size() << " walls found, floor -> roof: "<< ;
	LOG(INFO) << ll::unsafe_format("%d walls found, floor %.3f -> roof %.3f, total height: %3.f.",
		wallPieces.size(), zFloor, zRoof, (zRoof - zFloor));
	LOG_IF(INFO, wallPieces.size() % 2 != 0) << "note that wall size is odd, we may missed some.";

	//2. setup outline
	Eigen::vector<Eigen::Vector2f> cens(wallPieces.size());
	std::transform(wallPieces.begin(), wallPieces.end(), cens.begin(),
		[](const std::pair<PointCloud::Ptr, Segment>& pr) {
		return static_cast<Eigen::Vector2f>(pr.second.mid().block<2, 1>(0, 0));
	});
	auto sorted_indices = geo::sort_points_counter_clockwise(cens);

	Eigen::vector<Eigen::Vector2f> outline;
	{
		Eigen::vector<Eigen::Vector3f> lineends;
		for (auto i : sorted_indices) {
			const auto& seg = wallPieces[i].second;
			if (seg.s_.cross(seg.e_)[2] > 0.f) {
				lineends.push_back(seg.s_);
				lineends.push_back(seg.e_);
			} else {
				lineends.push_back(seg.e_);
				lineends.push_back(seg.s_);
			}
		}

		auto uiter = std::unique(lineends.begin(), lineends.end(), [](const Eigen::Vector3f& v1, const Eigen::Vector3f& v2) {
			return (v1 - v2).squaredNorm() < OUTLINE_LINK_DIS_THRESH_SQUARED;
		});

		outline.resize(std::distance(lineends.begin(), uiter));
		std::transform(lineends.begin(), uiter, outline.begin(), [](const Eigen::Vector3f& v) { return static_cast<Eigen::Vector2f>(v.block<2, 1>(0, 0));  });
		if ((outline.back() - outline.front()).squaredNorm() < OUTLINE_LINK_DIS_THRESH_SQUARED) outline.pop_back();
	}
	LOG(INFO) << "outline/wall size: " << outline.size() << " / " << wallPieces.size();

	//3. match
	const auto& blueprint3d = cadModel.getTypedModelItems(ITEM_BOTTOM_E).front();
	LOG(INFO) << blueprint3d.points_.size() << " points in blueprint.";

	Eigen::vector<Eigen::Vector2f> blueprint(blueprint3d.points_.size());
	std::transform(blueprint3d.points_.begin(), blueprint3d.points_.end(), blueprint.begin(),
		[](const Eigen::Vector3d& v) { return Eigen::Vector2f(v(0, 0), v(1, 0)); });

	Eigen::Matrix4f T = computeOutlineTransformation(blueprint, outline, BLUEPRINT_MATCH_DIS_THRESH);
	T(2, 3) = -zFloor; // shift to zero
	LOG(INFO) << "T: \n" << T << "\n";

	//4. transform point cloud
	allPieces.clear();
	for (auto& pr : wallPieces) allPieces.push_back(geo::transfromPointCloud(pr.first, T));
	allPieces.push_back(geo::transfromPointCloud(floor, T));
	allPieces.push_back(geo::transfromPointCloud(roof, T));

	// debug write.
	for (auto pr : ll::enumerate(allPieces))
		pcl::io::savePCDFile(ll::unsafe_format("coarse-wall-%d.pcd", pr.index), *(*pr.iter), true);

	// simple visualization
	pcl::visualization::PCLVisualizer viewer;

	for (auto pr : ll::enumerate(allPieces)) {
		auto cloud = *pr.iter;

		double r{ geo::random() }, g{ geo::random() }, b{ geo::random() };

		pcl::visualization::PointCloudColorHandlerCustom<Point> color(cloud, r * 255., g * 255., b * 255.);
		viewer.addPointCloud(cloud, color, "piece" + std::to_string(pr.index));
	}

	for (auto pr : ll::enumerate(outline)) {
		const Eigen::Vector2f& v = *pr.iter;
		Point p(v[0], v[1], 0.f);

		viewer.addSphere(p, 0.1f, "outline" + std::to_string(pr.index));
		viewer.addText3D(std::to_string(pr.index), Point(p.x, p.y, 0.2f), 0.3f);
	}

	{
		auto cad = cadModel.genTestCloud();
		pcl::visualization::PointCloudColorHandlerCustom<Point> color(cad, 255., 0., 0.);
		viewer.addPointCloud(cad, color, "cad");
	}

	while (!viewer.wasStopped()) {
		viewer.spinOnce(33);
		std::this_thread::sleep_for(std::chrono::microseconds(33));
	}

	return true;
}

std::pair<PointCloud::Ptr, Eigen::Vector4f> CoarseMatching::refinePlanePattern(PointCloud::Ptr cloud, double disthresh) const {
	pcl::SampleConsensusModelPlane<Point>::Ptr planeModel(new pcl::SampleConsensusModelPlane<Point>(cloud));
	pcl::RandomSampleConsensus<Point> ransac(planeModel);
	ransac.setDistanceThreshold(disthresh);
	ransac.computeModel();

	std::vector<int> inliers;
	Eigen::VectorXf params;
	ransac.getInliers(inliers);
	ransac.getModelCoefficients(params);

	auto plane = geo::getSubSet(cloud, inliers);
	Eigen::Vector4f abcd = params.block<4, 1>(0, 0);

	LOG(INFO) << "refined plane: " << plane->size() * 100.f / static_cast<float>(cloud->size()) << "% inliers. params: " << params.transpose() << "\n";

	return std::make_pair(plane, abcd);
#if 0
	pcl::ModelCoefficients coeffs;
	pcl::PointIndices inliers;
	pcl::SACSegmentation<Point> seg;
	seg.setOptimizeCoefficients(true);
	seg.setModelType(pcl::SACMODEL_PLANE);
	seg.setMethodType(pcl::SAC_RANSAC);
	seg.setDistanceThreshold(disthresh);
	LL_LOG_F("probability: %.2f, max iterations %d.\n", seg.getProbability(), seg.getMaxIterations());
	seg.setProbability(0.2f);
	seg.setMaxIterations(100);

	seg.setInputCloud(cloud);

	float ratio = 0.f;
	for (int i = 0; i < 1 && ratio < 0.3f; ++i) {
		seg.segment(inliers, coeffs);

		ratio = inliers.indices.size() / static_cast<float>(cloud->size());
		LL_LOG_F("%d: %.2f inliers.\n", i, ratio * 100.f);
	}

	auto plane = geo::getSubSet(cloud, inliers.indices);
	Eigen::Vector4f params;
	for (int i = 0; i < 4; ++i) params(i, 0) = coeffs.values[i];

	LL_LOG_F("refined plane: %.2f.\n", plane->size() * 100.f / static_cast<float>(cloud->size()));
	LL_LOG << params.transpose() << "\n";

	return std::make_pair(plane, params);
#endif
}

CoarseMatching::Segment CoarseMatching::detectSegementXoY(PointCloud::Ptr cloud, const Eigen::Vector4f& plane_param,
	double disthresh, double connect_thresh) const {

	if (std::abs(plane_param[2]) > 0.1f) {
		LOG(WARNING) << "initial plane n.z: " << plane_param[2] << ", may not a proper call.";
	}

	double ratio = std::min(0.1, 20000. / cloud->size()); // need sparse, todo: if we now the height, we can simply do passthrough
	PointCloud::Ptr sparsed = geo::filterPoints(cloud, [ratio](const Point& p) { return geo::random() < ratio; });
	PointCloud::Ptr cloud2d = geo::mapPoints(sparsed, [](const Point& p) { return Point(p.x, p.y, 0.f); });

	LOG(INFO) << "filter|map to 2d: " << cloud->size() << " -> " << cloud2d->size();

	std::vector<int> indices;
	Eigen::VectorXf params;
	std::tie(indices, params) = geo::detectOneLineRansac(cloud2d, disthresh);
	LOG_IF(INFO, std::fabs(params[5]) > 1e-6) << "detected line is not parallel to XoY!";

	auto inliers = geo::getSubSet(cloud2d, indices);
	auto body = geo::clusterMainStructure(inliers, connect_thresh);
	LOG(INFO) << "get segment cloud: " << cloud2d->size() << "->" << inliers->size() << " -> " << body->size();

	// we now get ends
	Point p = geo::P_(params.block<3, 1>(0, 0));
	Point n = geo::P_(params.block<3, 1>(3, 0));

	auto projlens = ll::mapf([p, n](const Point& q) {return geo::dot(q - p, n); }, body->points);
	auto pr = std::minmax_element(projlens.begin(), projlens.end());
	std::size_t i = std::distance(projlens.begin(), pr.first);
	std::size_t j = std::distance(projlens.begin(), pr.second);

	Segment segment;
	segment.s_ = geo::V_(body->points[i]);
	segment.e_ = geo::V_(body->points[j]);

	return segment;
}

Eigen::Matrix4f CoarseMatching::computeOutlineTransformation(const Eigen::vector<Eigen::Vector2f>& blueprint,
	const Eigen::vector<Eigen::Vector2f>& outline, float disthresh) {
	// try to align by each wall
	struct can {
		std::size_t i, j;
		trans2d::Matrix2x3f T;
		float maxmindis;

		std::vector<std::size_t> nearest_indices; // defines which indices were matched: #nearest_indices = #outline

		std::string to_string() const {
			std::stringstream ss;
			ss << "(" << i << ", " << j << ") " << maxmindis << " [" << ll::string_join(nearest_indices, ", ") << "]";
			// << "\n"  << T;
			return ss.str();
		}
	};

	const float disthresh_squared = disthresh * disthresh;

	std::vector<can> candidates;
	for (std::size_t i = 0; i < blueprint.size(); ++i) {
		Eigen::Vector2f s1 = blueprint[i];
		Eigen::Vector2f e1 = blueprint[(i + 1) % blueprint.size()];
		float len = (e1 - s1).norm();

		// check each segment
		for (std::size_t j = 0; j < outline.size(); ++j) {
			Eigen::Vector2f s2 = outline[j];
			Eigen::Vector2f e2 = outline[(j + 1) % outline.size()];

			// check wall lens
			//todo: check wall lens 2
			if (std::fabs((e2 - s2).norm() - len) > 2.0f * disthresh)
				continue;

			trans2d::Matrix2x3f T = trans2d::estimateTransform(s1, e1, s2, e2);

			// now check max-min distance
			// [(index, distance^ 2), ...] todo: may check one to one
			auto get_nearest_point_in_blueprint = [&blueprint](const Eigen::Vector2f& p)->std::pair<std::size_t, float> {
				auto prpr = ll::min_by([&p](const Eigen::Vector2f& q) { return (p - q).squaredNorm(); }, blueprint);
				return std::make_pair(std::distance(blueprint.begin(), prpr.first), prpr.second);
			};

			auto matchids = ll::mapf([&](const Eigen::Vector2f& p) {
				return get_nearest_point_in_blueprint(trans2d::transform(p, T));
			}, outline);

			float maxmindis = ll::max_by(&ll::get_value<std::size_t, float>, matchids).second;

			if (maxmindis > disthresh_squared) continue;

			can c;
			c.i = i;
			c.j = j;
			c.T = T;
			c.maxmindis = std::sqrt(maxmindis);
			c.nearest_indices = ll::mapf(&ll::get_key<std::size_t, float>, matchids);
			candidates.push_back(c);
		}
	}

	if (candidates.empty()) {
		LOG(INFO) << "no candidates found, match failed!\n";
		return Eigen::Matrix4f::Identity();
	}

	LOG(INFO) << candidates.size() << " candidates found:\n";
	for (auto pr : ll::enumerate(candidates))
		LOG(INFO) << "[" << pr.index << "] " << pr.iter->to_string() << "\n";

	// remove duplicate match pattern
	{
		auto is_same = [](const can& ci, const can& cj) {
			return ci.nearest_indices.size() == cj.nearest_indices.size() &&
				std::equal(ci.nearest_indices.begin(), ci.nearest_indices.end(), cj.nearest_indices.begin());
		};

		// std::make_unique need a sort, simple brute-force. this can be optimized
		std::vector<can> unique_candidates;
		for (const auto& ci : candidates) {
			auto iter = std::find_if(unique_candidates.begin(), unique_candidates.end(), [&ci, &is_same](const can& cj) { return is_same(ci, cj); });
			if (iter == unique_candidates.end()) {
				unique_candidates.push_back(ci);
			} else if (iter->maxmindis > ci.maxmindis) {
				iter->i = ci.i;
				iter->j = ci.j;
				iter->T = ci.T;
				iter->maxmindis = ci.maxmindis;
			}
		}

		candidates.swap(unique_candidates);
	}

	LOG(INFO) << candidates.size() << " UNIQUE candidates found:\n";
	for (auto pr : ll::enumerate(candidates))
		LOG(INFO) << "[" << pr.index << "] " << pr.iter->to_string() << "\n";

	//todo: may remove duplicated transform
	auto iter = std::min_element(candidates.begin(), candidates.end(), [](const can& ci, const can& cj) { return ci.maxmindis < cj.maxmindis; });

	LOG(INFO) << "choosing: " << std::distance(candidates.begin(), iter) << "\n";

	return trans2d::asTransform3d(iter->T);
}

}
