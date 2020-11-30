#include "CoarseMatching.h"
#include "glog/logging.h"

#include <pcl/sample_consensus/model_types.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/ransac.h>
#include <pcl/sample_consensus/sac_model_line.h>
#include <pcl/sample_consensus/sac_model_plane.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/ModelCoefficients.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/filters/radius_outlier_removal.h>

namespace CloudReg {

///

std::vector<PointCloud::Ptr> CoarseMatching::MatchResult::getAllPieces() const {
	auto pieces = walls_;
	pieces.push_back(floor_);
	pieces.push_back(roof_);
	return pieces;
}

std::string CoarseMatching::MatchResult::toString() const {
	std::stringstream ss;
	ss << "MatchResult: \n"
		<< walls_.size() << " walls, match index: [" << ll::string_join(matchIndices_, ", ") << "]\n"
		<< "T: \n" << T_;

	return ss.str();
}

///

CoarseMatching::CoarseMatching() {}

CoarseMatching::MatchResult CoarseMatching::run(const std::vector<PointCloud::Ptr>& allPieces, const CADModel& cadModel) {
	constexpr float HOR_CHECK_EPS = 0.01f;

	LOG(INFO) << "coarse match " << allPieces.size() << " walls.";

#if 0
	{
		pcl::visualization::PCLVisualizer viewer;

		auto frame = cadModel.genTestFrameCloud();
		auto frag = cadModel.genTestFragCloud(0.01);

		pcl::visualization::PointCloudColorHandlerCustom<Point> red(frame, 255., 0., 0.);
		viewer.addPointCloud(frame, red, "frame");
		viewer.addPointCloud(frag, "fragment");

		while (!viewer.wasStopped()) {
			viewer.spinOnce(33);
		}
	}
#endif

	//0. downsample
	const auto sparsedPieces = ll::mapf([](PointCloud::Ptr cloud) { return geo::downsampleUniformly(cloud, 0.01f); }, allPieces);
	// debug write.
	for (auto pr : ll::enumerate(sparsedPieces))
		pcl::io::savePCDFile(ll::unsafe_format("coarse-wall-%d.pcd", pr.index), *(*pr.iter), true);

	//1. basically, we use walls to do coarse matching.
	std::vector<std::pair<PointCloud::Ptr, Segment>> wallPieces;
	wallPieces.reserve(sparsedPieces.size()); // may -2
	PointCloud::Ptr floor, roof; // still we need floor and roof to get the proper Y
	float zFloor, zRoof;

	for (auto pr : ll::enumerate(sparsedPieces)) {
		auto cloud = *pr.iter;

		PointCloud::Ptr plane;
		Eigen::Vector4f params;
		std::tie(plane, params) = refinePlanePattern(cloud, PLANE_REFINE_DISTANCE);

		float nz = params[2];
		LOG(INFO) << pr.index << "] nz: " << nz;
		if (std::fabs(nz) < HOR_CHECK_EPS) {
			// we got a wall piece
			auto segment = detectMainSegementXoY(plane, params, LINE_DETECT_DIS_THRESH, LINE_DETECT_CONNECT_THRESH);
			// sort counter clock wise
			if (segment.s_.cross(segment.e_)[2] < 0.f) {
				Eigen::Vector3f s = segment.s_;
				segment.s_ = segment.e_;
				segment.e_ = s;
			}
			wallPieces.emplace_back(plane, segment);
		} else {
			float z = -params[3] / nz;
			if (z > 0.f) { roof = plane; zRoof = z; } else { floor = plane; zFloor = z; }
		}
	}

	// LOG(INFO) << wallPieces.size() << " walls found, floor -> roof: "<< ;
	LOG(INFO) << ll::unsafe_format("%d walls found, floor %.3f -> roof %.3f, total height: %.3f.",
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
			lineends.push_back(seg.s_);
			lineends.push_back(seg.e_);
		}

		// debug
		{
			std::stringstream ss;
			for (const auto& p : lineends) ss << p.transpose() << ", ";
			LOG(INFO) << "sorted segments: " << ss.str();
		}

		auto uiter = std::unique(lineends.begin(), lineends.end(), [](const Eigen::Vector3f& v1, const Eigen::Vector3f& v2) {
			return (v1 - v2).squaredNorm() < OUTLINE_LINK_DIS_THRESH_SQUARED;
		});

		outline.resize(std::distance(lineends.begin(), uiter));
		std::transform(lineends.begin(), uiter, outline.begin(), [](const Eigen::Vector3f& v) { return static_cast<Eigen::Vector2f>(v.block<2, 1>(0, 0));  });
		if ((outline.back() - outline.front()).squaredNorm() < OUTLINE_LINK_DIS_THRESH_SQUARED) outline.pop_back();

		// debug
		std::stringstream ss;
		for (std::size_t i = 0; i < outline.size(); ++i) {
			const auto& cur = outline[i];
			const auto& nxt = outline[(i + 1) % outline.size()];

			ss << "(" << cur.transpose() << ") -- " << (nxt - cur).squaredNorm() << " --> ";
		}
		LOG(INFO) << "outline chain (threshold: " << OUTLINE_LINK_DIS_THRESH_SQUARED << "): " << ss.str();
	}
	LOG(INFO) << "outline/wall size: " << outline.size() << " / " << wallPieces.size();

	//3. match
	const auto& blueprint3d = cadModel.getTypedModelItems(ITEM_BOTTOM_E).front();
	LOG(INFO) << blueprint3d.points_.size() << " points in blueprint.";

	Eigen::vector<Eigen::Vector2f> blueprint(blueprint3d.points_.size());
	std::transform(blueprint3d.points_.begin(), blueprint3d.points_.end(), blueprint.begin(),
		[](const Eigen::Vector3d& v) { return Eigen::Vector2f(v(0, 0), v(1, 0)); });

	Eigen::Matrix4f T = Eigen::Matrix4f::Identity();
	{
		const auto& Ts = computeOutlineTransformCandidates(blueprint, outline, BLUEPRINT_MATCH_DIS_THRESH);
		if (Ts.size() == 1) T = trans2d::asTransform3d(Ts.front());
		else {
			std::vector<PointCloud::Ptr> walls;
			std::vector<Eigen::VectorXf> wallparams;
			{
				walls.reserve(wallPieces.size());
				wallparams.reserve(wallPieces.size());

				for (const auto& pr : wallPieces) {
					walls.push_back(pr.first);
					Eigen::VectorXf params(6, 1);
					params.block<3, 1>(0, 0) = pr.second.s_;
					params.block<3, 1>(3, 0) = pr.second.dir();
					wallparams.push_back(params);
				}
			}
			T = trans2d::asTransform3d(chooseTransformByHoles(cadModel, walls, wallparams, Ts, zFloor));
		}
	}
	T(2, 3) = -zFloor; // shift to zero
	LOG(INFO) << "T: \n" << T << "\n";

	//4. transform point cloud, prepare result
	CoarseMatching::MatchResult result;
	{
		LL_ASSERT(blueprint3d.segments_.size() == blueprint3d.points_.size() && "havnt build segment?");

		const auto& bpSegments = blueprint3d.segments_;
		// search by ends distance
		auto get_nearst_wall_index_in_blueprint = [&bpSegments](const Segment& seg, const Eigen::Matrix4f& T)-> std::size_t {
			Eigen::Vector3f s = T.block<3, 3>(0, 0) * seg.s_ + T.block<3, 1>(0, 3);
			Eigen::Vector3f e = T.block<3, 3>(0, 0) * seg.e_ + T.block<3, 1>(0, 3);

			auto pr = ll::min_by([&s, &e](const std::pair<Eigen::Vector3d, Eigen::Vector3d>& bpseg) {
				Eigen::Vector3d sd(s(0, 0), s(1, 0), s(2, 0)), ed(e(0, 0), e(1, 0), e(2, 0));
				return (sd - bpseg.first).squaredNorm() + (ed - bpseg.second).squaredNorm();
			}, bpSegments);

			return std::distance(bpSegments.cbegin(), pr.first);
		};

		result.T_ = T;

		result.walls_.reserve(wallPieces.size());
		result.matchIndices_.reserve(wallPieces.size());
		for (auto& pr : wallPieces) {
			result.walls_.push_back(geo::transfromPointCloud(pr.first, T));
			result.matchIndices_.push_back(get_nearst_wall_index_in_blueprint(pr.second, T));
		}

		result.floor_ = geo::transfromPointCloud(floor, T);
		result.roof_ = geo::transfromPointCloud(roof, T);
	}

	LOG(INFO) << result.toString();

	// simple visualization
	pcl::visualization::PCLVisualizer viewer;

	for (auto pr : ll::enumerate(result.walls_)) {
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
		auto cad = cadModel.genTestFrameCloud();
		pcl::visualization::PointCloudColorHandlerCustom<Point> color(cad, 255., 0., 0.);
		viewer.addPointCloud(cad, color, "cad");
	}

	while (!viewer.wasStopped()) {
		viewer.spinOnce(33);
	}

	return result;
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

CoarseMatching::Segment CoarseMatching::detectMainSegementXoY(PointCloud::Ptr cloud, const Eigen::Vector4f& plane_param,
	double disthresh, double connect_thresh) const {

	if (std::abs(plane_param[2]) > 0.1f) {
		LOG(WARNING) << "initial plane n.z: " << plane_param[2] << ", may not a proper call.";
	}

	PointCloud::Ptr cloud2d = geo::mapPoints(cloud, [](const Point& p) { return Point(p.x, p.y, 0.f); });

	LOG(INFO) << "filter|map to 2d: " << cloud->size() << " -> " << cloud2d->size();

	std::vector<int> indices;
	Eigen::VectorXf params;
	std::tie(indices, params) = geo::detectOneLineRansac(cloud2d, disthresh);
	LOG_IF(INFO, std::fabs(params[5]) > 1e-6) << "detected line is not parallel to XoY!";

	auto inliers = geo::getSubSet(cloud2d, indices);
	// auto body = geo::clusterMainStructure(inliers, connect_thresh);
	PointCloud::Ptr body(new PointCloud());
	{
		pcl::RadiusOutlierRemoval<Point> ror;
		ror.setInputCloud(inliers);
		ror.setRadiusSearch(0.05);
		ror.setMinNeighborsInRadius(10);
		ror.filter(*body);
	}
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

Eigen::vector<trans2d::Matrix2x3f> CoarseMatching::computeOutlineTransformCandidates(const Eigen::vector<Eigen::Vector2f>& blueprint,
	const Eigen::vector<Eigen::Vector2f>& outline, float disthresh) {
	// try to align by each wall
	struct can {
		std::size_t i, j;
		trans2d::Matrix2x3f T;
		float maxmindis;
		float avgdis;

		float a() const { return std::atan2f(T(1, 0), T(0, 0)); }

		std::string to_string() const {
			std::stringstream ss;
			ss << "(" << i << ", " << j << ") " << maxmindis << ", " << avgdis
				<< " [" << (a() * 180.f / geo::PI) << ", (" << t().transpose() << ")]";
			return ss.str();
		}

	private:
		Eigen::Vector2f t() const { return -T.block<2, 2>(0, 0).transpose() * T.block<2, 1>(0, 2); }
	};

	auto min_distance_to_blueprint = [&blueprint](const Eigen::Vector2f& p) {
		float mindis = std::numeric_limits<float>::max();
		for (std::size_t i = 0; i < blueprint.size(); ++i) {
			std::size_t j = (i + 1) % blueprint.size();
			float dis = geo::distance_to_segment_2d(p, blueprint[i], blueprint[j]);
			if (dis < mindis) mindis = dis;
		}

		return mindis;
	};

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
			auto mindis = ll::mapf([&T, &min_distance_to_blueprint](const Eigen::Vector2f& p) {
				Eigen::Vector2f q = T.block<2, 2>(0, 0) * p + T.block<2, 1>(0, 2);
				return min_distance_to_blueprint(q);
			}, outline);

			can c;
			c.i = i;
			c.j = j;
			c.T = T;
			c.maxmindis = *(std::max_element(mindis.begin(), mindis.end()));
			c.avgdis = ll::sum(mindis) / mindis.size();
			candidates.push_back(c);
		}
	}

	// sort
	std::sort(candidates.begin(), candidates.end(), [](const can& c1, const can& c2) {return c1.maxmindis < c2.maxmindis; });

	{
		std::stringstream ss;
		for (std::size_t i = 0; i < candidates.size() && i < 20; ++i) ss << i << "| " << candidates[i].to_string() << "\n";
		LOG(INFO) << "all candidates (threshold: " << disthresh << ") \n" << ss.str();
	}

	if (candidates.front().maxmindis > disthresh) {
		LOG(INFO) << "no candidates found, match failed!\n";
		auto I = trans2d::buildTransform(0.f, Eigen::Vector2f::Zero());
		return Eigen::vector<trans2d::Matrix2x3f>(1, I);
	}

	// remove duplicate match pattern by transformation
	{
		constexpr float ANGLE_THRESH = 5.f / geo::PI;
		constexpr float TRANSLATE_THRESH = 0.1f;

		auto is_same = [ANGLE_THRESH, TRANSLATE_THRESH](const can& ci, const can& cj) {
			float da = ci.a() - cj.a();
			if (da > geo::PI) da -= geo::TAU;
			else if (da < -geo::PI) da += geo::TAU;

			if (std::fabs(da) > ANGLE_THRESH) return false;

			// t dis
			Eigen::Matrix2f Ri = ci.T.block<2, 2>(0, 0);
			Eigen::Vector2f ti = ci.T.block<2, 1>(0, 2);
			Eigen::Vector2f tj = cj.T.block<2, 1>(0, 2);

			float dt = (Ri.transpose() * (tj - ti)).norm();
			return dt < TRANSLATE_THRESH;
		};

		std::cout << "===================================================\n";
		for (const auto& can : candidates) {
			is_same(candidates.front(), can);
		}

		// simple brute-force. this can be optimized
		// keep in mind that the candidates are already sorted.
		std::vector<can> unique_candidates;
		for (const auto& ci : candidates) {
			if (ci.maxmindis > disthresh) break;

			auto iter = std::find_if(unique_candidates.begin(), unique_candidates.end(), [&ci, &is_same](const can& cj) { return is_same(ci, cj); });
			if (iter == unique_candidates.end()) {
				unique_candidates.push_back(ci);
			} else if (iter->avgdis > ci.avgdis) {
				// copy.
				iter->i = ci.i;
				iter->j = ci.j;
				iter->T = ci.T;
				iter->maxmindis = ci.maxmindis;
				iter->avgdis = ci.avgdis;
			}
		}

		candidates.swap(unique_candidates);
	}

	LOG(INFO) << candidates.size() << " UNIQUE candidates found:\n";
	for (auto pr : ll::enumerate(candidates))
		LOG(INFO) << pr.index << "| " << pr.iter->to_string();

	Eigen::vector<trans2d::Matrix2x3f> Ts;
	Ts.reserve(candidates.size());
	for (const auto& c : candidates) Ts.emplace_back(c.T);

	return Ts;
}

trans2d::Matrix2x3f CoarseMatching::chooseTransformByHoles(const CADModel& cadModel,
	const std::vector<PointCloud::Ptr>& walls, const std::vector<Eigen::VectorXf>& wallparams,
	const Eigen::vector<trans2d::Matrix2x3f>& Ts, const float floorZ) const {

	constexpr float HALF_SLICE_THICKNESS = 0.02f;

	// todo: we may check order
	// T* src = dst, return #{matched srcs}
	auto get_matched_num = [](const Eigen::vector<Eigen::Vector2f>& src,
		const Eigen::vector<Eigen::Vector2f>& dst, const trans2d::Matrix2x3f& T) -> std::size_t {
		Eigen::vector<Eigen::Vector2f> transed_src(src.size());
		std::transform(src.begin(), src.end(), transed_src.begin(), [&T](const Eigen::Vector2f& v) { return trans2d::transform(v, T); });

		return std::count_if(transed_src.begin(), transed_src.end(), [&dst](Eigen::Vector2f& sp) {
			return std::any_of(dst.begin(), dst.end(), [&sp](const Eigen::Vector2f& dp) { return (sp - dp).norm() < 0.1f; });
		});
	};

	auto get_consistent_transformations_for = [&](const Eigen::vector<trans2d::Matrix2x3f>& checkTs, float z)-> Eigen::vector<trans2d::Matrix2x3f> {
		Eigen::vector<Eigen::Vector3d> cadHoleEnds = intersectCADModelOnZ(cadModel, z);
		if (cadHoleEnds.empty()) {
			LOG(WARNING) << "bad function call: no intersection points found at: " << z;
			return Ts;
		}

		// slice and detect ends
		Eigen::vector<Eigen::Vector3f> wallEnds;
		for (auto pr : ll::enumerate(walls)) {
			auto wall = *pr.iter;
			auto params = wallparams[pr.index];

			auto line = geo::passThrough(wall, "z", floorZ + z - HALF_SLICE_THICKNESS, floorZ + z + HALF_SLICE_THICKNESS);

			auto subsegs = splitSegments(line, params, 0.1f, 0.1f);

			for (const auto& seg : subsegs) {
				wallEnds.emplace_back(seg.s_);
				wallEnds.emplace_back(seg.e_);
			}
		}

		Eigen::vector<Eigen::Vector2f> hole2d, wall2d;
		{
			hole2d.reserve(cadHoleEnds.size());
			wall2d.reserve(wallEnds.size());
			for (const auto& v : cadHoleEnds) hole2d.emplace_back(v(0), v(1));
			for (const auto& v : wallEnds) wall2d.emplace_back(v(0), v(1));
		}

		// now check if matched
		Eigen::vector<trans2d::Matrix2x3f> re;
		for (const auto& T : checkTs) {
			// succeed if only each cad hole have candidate
			std::size_t cnt = get_matched_num(hole2d, wall2d, trans2d::inverse(T));

			if (cnt == cadHoleEnds.size()) re.emplace_back(T);
		}

		return re;
	};

	// choose height
	if (!cadModel.containModels(ITEM_HOLE_E)) {
		LOG(WARNING) << "no holes in cad model! choose the first one (with least max min distance).";
		return Ts.front();
	}

	const auto& holes = cadModel.getTypedModelItems(ITEM_HOLE_E);
	auto hranges = ll::mapf([](const ModelItem& mi) { return mi.highRange_; }, holes);

	const double low = ll::min_by(ll::get_key<double, double>, hranges).second;
	const double high = ll::max_by(ll::get_value<double, double>, hranges).second;
	LOG(INFO) << "hole high range: " << low << " -> " << high << " (" << (high - low) << ")";

	// this can be very inefficient for now.
	// todo: rethink & optimize this process, use grid for now.
	constexpr double GRID_SIZE = HALF_SLICE_THICKNESS * 2.5;
	const int N = static_cast<int>((high - low) / GRID_SIZE) + 1;
	std::vector<int> cell_counts(N, 0);
	for (const auto& hr : hranges) {
		int ilow = (hr.first - low) / GRID_SIZE;
		int ihigh = (hr.second - low) / GRID_SIZE;

		for (int i = ilow; i <= ihigh; ++i) cell_counts[i]++;
	}
	LOG(INFO) << "cell counts: " << ll::string_join(cell_counts, ", ");

	// only avoid jumps for now...
	Eigen::vector<trans2d::Matrix2x3f> leftTs = Ts;
	for (std::size_t i = 1; i < cell_counts.size() - 1; ++i) {
		if (cell_counts[i - 1] == cell_counts[i] && cell_counts[i] == cell_counts[i + 1] && cell_counts[i] != 0) {
			float mid = low + i * GRID_SIZE + GRID_SIZE / 2.f;

			auto lTs = get_consistent_transformations_for(leftTs, mid);
			if (lTs.empty()) {
				LOG(WARNING) << "no consistent transformation at: " << mid;
			} else if (lTs.size() == 1) {
				LOG(INFO) << "unique consistence check succeed at: " << mid;
				return lTs.front();
			} else {
				LOG(INFO) << lTs.size() << " consistent transformation found at: " << mid;
				leftTs.swap(lTs);
			}
		}
	}

	// we failed to find any
	LOG(WARNING) << "failed to remove ambiguity, left candidates: " << leftTs.size() << "/" << Ts.size() << ". choose the first one (with least max min distance).";
	return leftTs.front();
}

std::vector<CoarseMatching::Segment> CoarseMatching::detectSegmentsXOY(PointCloud::Ptr cloud,
	float linethresh, float connect_thresh, float min_seg_length) const {
	ll::TimeCounter tc([](int ms) { LOG(INFO) << ms << " ms passed in " << __FUNCTION__; });

	auto cloud2d = geo::mapPoints(cloud, [](const Point& p) { return Point(p.x, p.y, 0.f); });

	auto rest = cloud2d;
	LOG(INFO) << "left points: " << rest->size();

	std::vector<Segment> segments;

	std::vector<int> inliers;
	Eigen::VectorXf params;
	while (rest->size() > 1000) {
		std::tie(inliers, params) = geo::detectOneLineRansac(rest, linethresh);
		LOG(INFO) << "line found: " << params.transpose();

		auto line = geo::getSubSet(rest, inliers);
		auto subsegs = splitSegments(line, params, connect_thresh, min_seg_length);
		segments.insert(segments.end(), subsegs.begin(), subsegs.end());

		rest = geo::getSubSet(rest, inliers, true);

		LOG(INFO) << segments.size() << " segments, left points: " << rest->size();
	}

	return segments;
}

std::vector<CoarseMatching::Segment> CoarseMatching::splitSegments(PointCloud::Ptr cloud, Eigen::VectorXf line_params,
	float connect_thresh, float min_seg_length) const {
	Point p = geo::P_(line_params.block<3, 1>(0, 0));
	Point n = geo::P_(line_params.block<3, 1>(3, 0));

	auto projlens = ll::mapf([p, n](const Point& q) {return geo::dot(q - p, n); }, cloud->points);

	// a simple 1d cluster, todo: optimize this
	std::vector<std::unordered_set<std::size_t>> sets;
	for (std::size_t i = 0; i < projlens.size(); ++i) {
		// find connected sets
		std::vector<std::size_t> connected_sets;
		for (auto pr : ll::enumerate(sets)) {
			auto& ids = *pr.iter;
			if (std::any_of(ids.begin(), ids.end(), [&](std::size_t id) {
				return std::fabs(projlens[i] - projlens[id]) < connect_thresh; }))
				connected_sets.push_back(pr.index);
		}

		if (connected_sets.empty()) {
			// find an empty set to 
			std::unordered_set<std::size_t> ids;
			ids.insert(i);
			sets.push_back(ids);
		} else {
			// join them
			auto& joined = sets[connected_sets.front()];
			joined.insert(i);
			for (std::size_t j = 1; j < connected_sets.size(); ++j) {
				auto& tobe_joined = sets[connected_sets[j]];
				joined.insert(tobe_joined.begin(), tobe_joined.end());
				tobe_joined.clear(); // clear
			}
		}
	}

	auto cnt = std::count_if(sets.begin(), sets.end(), [](const std::unordered_set<std::size_t>& c) { return c.size() > 1; });

	std::vector<Segment> segments;
	segments.reserve(cnt);
	for (const auto& ids : sets) {
		if (ids.size() < 2) continue;

		auto pr = std::minmax_element(ids.begin(), ids.end(), [&projlens](std::size_t i, std::size_t j) {
			return projlens[i] < projlens[j];
		});

		Segment s;
		s.s_ = geo::V_(cloud->points[*pr.first]);
		s.e_ = geo::V_(cloud->points[*pr.second]);

		if ((s.s_ - s.e_).norm() < min_seg_length) continue;

		segments.push_back(s);
	}

	LOG(INFO) << segments.size() << " / " << cnt << " segment(s) found.";

	return segments;
}

Eigen::vector<Eigen::Vector3d> CoarseMatching::intersectCADModelOnZ(const CADModel& cadModel, float z) const {
	const auto& roof = cadModel.getTypedModelItems(ITEM_TOP_E).front();
	const auto& floor = cadModel.getTypedModelItems(ITEM_BOTTOM_E).front();

	const auto& bottom_outline = floor.points_;

	const float floorZ = floor.points_.front()(2, 0);
	const float roofZ = roof.points_.front()(2, 0);

	if (z <= floorZ || z >= roofZ || !cadModel.containModels(ITEM_HOLE_E)) return bottom_outline;

	// we now try intersect each hole
	// note that the hole can be non-convex, so we only need intersection points
	Eigen::vector<Eigen::Vector3d> points;

	const auto& holes = cadModel.getTypedModelItems(ITEM_HOLE_E); // we already checked this in above

	for (const auto& hole : holes) {
		const auto& hr = hole.highRange_;
		if (z <= hr.first || z >= hr.second) continue;

		Eigen::vector<Eigen::Vector3d> ips;
		for (const auto& seg : hole.segments_) {
			auto sii = geo::zIntersectSegment(seg.first, seg.second, z);
			if (sii.valid()) { //todo: use margin check
				ips.emplace_back(sii.point_);
			}
		}
		LOG_IF(WARNING, ips.size() % 2 != 0) << "intersection points is odd! may intersected with ends.";

		// todo: may just insert refer to its parent index, but we just sort here then
		points.insert(points.end(), ips.begin(), ips.end());
	}

	LOG(INFO) << points.size() << " intersection points found on z: " << z;

	if (!points.empty()) {
		// remove duplicated points
		// if points are too close, then this may NOT be an detectable point
		constexpr double CHECK_THRESH = 0.01;
		std::vector<Eigen::vector<Eigen::Vector3d>> sets; // a simple cluster
		for (const auto& p : points) {
			auto iter = std::find_if(sets.begin(), sets.end(), [&p, CHECK_THRESH](const Eigen::vector<Eigen::Vector3d>& set) {
				return std::any_of(set.begin(), set.end(), [&p, CHECK_THRESH](const Eigen::Vector3d& q) { return (p - q).norm() < CHECK_THRESH; });
			});
			if (iter == sets.end()) {
				sets.emplace_back(1, p);
			} else {
				iter->emplace_back(p);
			}
		}

		points.clear();
		for (const auto& s : sets) if (s.size() == 1) points.emplace_back(s.front());

		LOG(INFO) << points.size() << " valid points found.";
	}

	return points;
}

}
