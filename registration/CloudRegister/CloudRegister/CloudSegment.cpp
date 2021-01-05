#include "CloudSegment.h"

#include "glog/logging.h"

#include <pcl/sample_consensus/model_types.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/ransac.h>
#include <pcl/sample_consensus/sac_model_line.h>
#include <pcl/sample_consensus/sac_model_plane.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/segmentation/region_growing.h>
#include <pcl/ModelCoefficients.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/filters/radius_outlier_removal.h>
#include <pcl/features/normal_3d.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/search/search.h>

#include "GeometryUtils.h"
#include "funHelper.h"

namespace CloudReg {

CloudSegment::SegmentResult CloudSegment::run() {
	recordModelBoundingBox();

	if (!alignCloudToCADModel()) {
		LOG(ERROR) << "failed to align cloud to cadModel.";
	}

	return segmentByCADModel();
}

void CloudSegment::recordModelBoundingBox() {
	const auto& floor = cadModel_.getTypedModelItems(ITEM_BOTTOM_E).front();
	for (const auto& p : floor.points_) {
		if (p(0) < cadx1_) cadx1_ = p(0);
		if (p(0) > cadx2_) cadx2_ = p(0);
		if (p(1) < cady1_) cady1_ = p(1);
		if (p(1) > cady2_) cady2_ = p(1);
	}

	cadz1_ = floor.points_.front()(2, 0);
	const auto& roof = cadModel_.getTypedModelItems(ITEM_TOP_E).front();
	cadz2_ = roof.points_.front()(2, 0);

	LOG(INFO) << ll::unsafe_format("cad model AABB bounding box: %.3f -> %.3f, %.3f -> %.3f, %.3f -> %.3f.",
		cadx1_, cadx2_, cady1_, cady2_, cadz1_, cadz2_);
}

bool CloudSegment::alignCloudToCADModel() {
	//1. found vertical planes
	auto cloud = sliceMainBody(); // zroof, zfloor record here.

	//note: 10000 may relate to downsample in sliceMainBody
	// auto planes = detectPlanes(cloud, 0.02, 10000); 
	auto planes = detectRegionPlanes(cloud, 3. / 180. * geo::PI, 1., 10000);

	auto wallCandidates = ll::filter([](const PlaneCloud& pc) { return std::fabs(pc.n()(2)) < 0.02f; }, planes);
	LOG(INFO) << ll::unsafe_format("%d / %d wall candidates detected: ", wallCandidates.size(), planes.size());
	{
		std::stringstream ss;
		for (const auto& pc : wallCandidates)
			ss << pc.n().transpose() << "\n";
		LOG(INFO) << "candidates normals: " << ss.str();
	}

	//2. map to segment 
	auto wallSegments = ll::mapf([&](const PlaneCloud& pc) {
		const Eigen::Vector4f& params = pc.abcd_;
		Eigen::Vector3f n = Eigen::Vector3f(-params(1), params(0), 0.f).normalized();
		Eigen::Vector3f p;
		if (std::fabs(params(0)) < 1e-6) p = Eigen::Vector3f(0.f, -params(3) / params(1), 0.f);
		else p = Eigen::Vector3f(-params(3) / params(0), 0.f, 0.f);

		return splitSegments(pc.cloud_, p, n, 0.1f, 0.1f);
	}, wallCandidates);

	// join them all, 
	// todo: reconsider this, like, if the segments is too much, may try take some of them.
	std::vector<Segment> allSegments;
	{
		std::size_t cnt = ll::sum_by([](const auto& v) { return v.size(); }, wallSegments);
		allSegments.reserve(cnt);
		for (const auto& segs : wallSegments) allSegments.insert(allSegments.end(), segs.begin(), segs.end());

		for (auto& s : allSegments) {
			// we dont need z info, set to -2 just for visualization.
			s.s_(2, 0) = -2.;
			s.e_(2, 0) = -2.;
		}

		// to counter clockwise.	
		for (auto& seg : allSegments) {
			double ps = atan2(seg.s_(1), seg.s_(0));
			double pe = atan2(seg.e_(1), seg.e_(0));
			if (pe < ps) std::swap(seg.s_, seg.e_);
		}

		std::stringstream ss;
		for (auto pr : ll::enumerate(allSegments)) ss << ll::unsafe_format("%02d: %.3f\n", pr.index, pr.iter->len());
		LOG(INFO) << allSegments.size() << " segments detected, length: \n" << ss.str();
	}

	//3. match
	// use cad floor, since the beams may affect roof's shape.
	const auto& blueprint3d = cadModel_.getTypedModelItems(ITEM_BOTTOM_E).front();
	Eigen::vector<Eigen::Vector2f> blueprint(blueprint3d.points_.size());
	std::transform(blueprint3d.points_.rbegin(), blueprint3d.points_.rend(), blueprint.begin(),
		[](const Eigen::Vector3d& v) { return Eigen::Vector2f(v(0, 0), v(1, 0)); });
	//^^^ note the above "rbegin rend", by this we reverse the points to counterclockwise.

	constexpr double MATCH_THRESH = 0.25;
	Eigen::vector<trans2d::Matrix2x3f> Ts = computeSegmentAlignCandidates(allSegments, blueprint, MATCH_THRESH);
	if (Ts.empty()) {
		LOG(ERROR) << " no candidate transform found.";
		return false;
	}

	// use the first one(minimal error) for now
	Eigen::Matrix4f T = trans2d::asTransform3d(trans2d::inverse(Ts.front()));
	T(2, 3) = -zfloor_; // shift up

	LOG(INFO) << "cloud trans found: \n" << T;

	//! note we do transform here
	sparsedCloud_ = geo::transfromPointCloud(sparsedCloud(), T);
	orgCloud_ = geo::transfromPointCloud(orgCloud_, T);

	return true;
}

CloudSegment::SegmentResult CloudSegment::segmentByCADModel() {
	// bad thing is, due to efficiency, we cannot just segment the original cloud. 
	auto sr = segmentCloudByCADModel(sparsedCloud());

	// bounding box filter
	float x1{ METER_100 }, x2{ -METER_100 }, y1{ METER_100 }, y2{ -METER_100 }, z1{ METER_100 }, z2{ -METER_100 };
	{
		auto update_aabb = [&](PointCloud::Ptr cloud) {
			if (!cloud) return;

			Point p1, p2;
			pcl::getMinMax3D(*cloud, p1, p2);

			x1 = std::min(x1, p1.x);
			y1 = std::min(y1, p1.y);
			z1 = std::min(z1, p1.z);
			x2 = std::max(x2, p2.x);
			y2 = std::max(y2, p2.y);
			z2 = std::max(z2, p2.z);
		};

		for (const auto& pc : sr.walls_) update_aabb(pc.cloud_);
		for (const auto& pc : sr.beams_) update_aabb(pc.cloud_);
		update_aabb(sr.roof_.cloud_);
		update_aabb(sr.floor_.cloud_);

		LOG(INFO) << ll::unsafe_format("bounding box: %.3f -> %.3f, %.3f -> %.3f, %.3f -> %.3f.",
			x1, x2, y1, y2, z1, z2);

		// margin
		x1 -= DOWNSAMPLE_SIZE;
		x2 += DOWNSAMPLE_SIZE;
		y1 -= DOWNSAMPLE_SIZE;
		y2 += DOWNSAMPLE_SIZE;
		z1 -= DOWNSAMPLE_SIZE;
		z2 += DOWNSAMPLE_SIZE;
	}

	// we now classify original cloud
	const auto& allplanes = sr.allPlanes();
	std::vector<std::vector<int>> clusters(allplanes.size());
	std::vector<pcl::search::KdTree<Point>> trees(allplanes.size());
	for (std::size_t i = 0; i < allplanes.size(); ++i) {
		auto cloud = allplanes[i].cloud_;
		if (cloud && !cloud->empty()) trees[i].setInputCloud(cloud);
	}

	pcl::Indices searchIndices;
	std::vector<float> searchDis;
	constexpr double MAX_DIS_SQUARED = 0.015 * 0.015;

	for (std::size_t idx = 0; idx < orgCloud_->size(); ++idx) {
		const auto& p = orgCloud_->points[idx];
		if (p.x< x1 || p.x> x2 || p.y< y1 || p.y> y2 || p.z< z1 || p.z> z2) continue;

		Eigen::Vector4f phom(p.x, p.y, p.z, 1.);

		auto ins = ll::filter([&](std::size_t i) {
			const auto& pc = allplanes[i];
			if (!pc.cloud_ || pc.cloud_->empty()) return false;

			float dis = phom.transpose().dot(pc.abcd_);
			return std::fabs(dis) < 0.02f;
		}, ll::range(allplanes.size()));

		if (ins.empty()) continue;
		if (ins.size() == 1) clusters[ins.front()].emplace_back(idx);
		else {
			// then we have to judge by distance to cloud....
			auto pr = ll::min_by([&](std::size_t i) {
				int cnt = trees[i].nearestKSearch(p, 1, searchIndices, searchDis);
				return cnt >= 1 ? searchDis.front() : 100.f;
			}, ins);

			if (pr.second < MAX_DIS_SQUARED) {
				std::size_t i = *pr.first;
				clusters[i].emplace_back(idx);
			}
		}
	}

	{
		std::stringstream ss;
		for (std::size_t i = 0; i < allplanes.size(); ++i)
			ss << allplanes[i].cloud_->size() << " -> " << clusters[i].size() << ", ";
		LOG(INFO) << "re-cluster done: " << ss.str();
	}

	// out
	int offset = 0;
	for (std::size_t i = 0; i < sr.walls_.size(); ++i) {
		if (clusters[i + offset].size() > 0)
			sr.walls_[i].cloud_ = geo::getSubSet(orgCloud_, clusters[i + offset]);
	}
	offset += sr.walls_.size();
	for (std::size_t i = 0; i < sr.beams_.size(); ++i) {
		if (clusters[i + offset].size() > 0)
			sr.beams_[i].cloud_ = geo::getSubSet(orgCloud_, clusters[i + offset]);
	}
	offset += sr.beams_.size();
	if (!clusters[offset].empty()) sr.roof_.cloud_ = geo::getSubSet(orgCloud_, clusters[offset]);
	offset += 1;
	if (!clusters[offset].empty()) sr.floor_.cloud_ = geo::getSubSet(orgCloud_, clusters[offset]);

#ifdef VISUALIZATION_ENABLED
	pcl::visualization::PCLVisualizer viewer;

	auto add_cloud = [&viewer](const std::string& name, PointCloud::Ptr cloud, double r, double g, double b) {
		pcl::visualization::PointCloudColorHandlerCustom<Point> color(cloud, r, g, b);
		viewer.addPointCloud(cloud, color, name);
	};
	auto add_cloud_rc = [&viewer, &add_cloud](const std::string& name, PointCloud::Ptr cloud) {
		double r{ geo::random() }, g{ geo::random() }, b{ geo::random() };
		add_cloud(name, cloud, r * 255., g * 255., b * 255.);
	};

	// add_cloud("cloud", cloud, 100., 100., 100.);

	for (const auto& pr : ll::enumerate(sr.walls_))
		if (pr.iter->cloud_)
			add_cloud_rc("wall" + std::to_string(pr.index), pr.iter->cloud_);
	for (const auto& pr : ll::enumerate(sr.beams_))
		if (pr.iter->cloud_)
			add_cloud_rc("beam" + std::to_string(pr.index), pr.iter->cloud_);
	add_cloud_rc("roof", sr.roof_.cloud_);
	add_cloud_rc("floor", sr.floor_.cloud_);

	add_cloud("cad", cadModel_.genTestFrameCloud(), 255., 0., 0.);

	while (!viewer.wasStopped()) {
		viewer.spinOnce(33);
	}
#endif

	return sr;
}

PointCloud::Ptr CloudSegment::sliceMainBody() {
	constexpr double MACHINE_RADIUS = 0.8f;
	constexpr double ZCUT_LOW = 0.5; // always avoid cloud under 0.5m, relative to floor
	constexpr double ZCUT_MARGIN = 0.1;

	auto cloud = sparsedCloud();

	//1. simple prune by bounding box
	{
		const double dw = (cadx2_ - cadx1_ - MACHINE_RADIUS);
		const double dh = (cady2_ - cady1_ - MACHINE_RADIUS);
		const double BOUNDING_BOX_MAX_HALF_WIDTH = std::sqrt(dw * dw + dh * dh);

		cloud = geo::passThrough(cloud, "x", -BOUNDING_BOX_MAX_HALF_WIDTH, BOUNDING_BOX_MAX_HALF_WIDTH);
		cloud = geo::passThrough(cloud, "y", -BOUNDING_BOX_MAX_HALF_WIDTH, BOUNDING_BOX_MAX_HALF_WIDTH);
	}

	//2. we now cut height range
	const double wallHeight = cadz2_ - cadz1_;
	const double zRoof = detectRoofHeight(cloud);
	const double zFloor = zRoof - wallHeight;

	zfloor_ = zFloor;
	zroof_ = zRoof;

	double zBeamLow{ cadz2_ }; // refer to cad model
	if (cadModel_.containModels(ITEM_BEAM_E)) {
		auto beams = cadModel_.getTypedModelItems(ITEM_BEAM_E);
		for (const auto& beam : beams) {
			const auto& hr = beam.highRange_;
			zBeamLow = std::min(zBeamLow, hr.first);
		}

		LOG(INFO) << ll::unsafe_format("beam lowest z: %.3f, (roof: %.3f)", zBeamLow, cadz2_);
	}

	double zhigh = zFloor + zBeamLow; // avoid beams
	double zlow = zFloor + ZCUT_LOW;

	LL_ASSERT((zhigh - zlow) > (ZCUT_MARGIN * 3) && "no range available.");

	LOG(INFO) << ll::unsafe_format("slice main body: %.3f -> %.3f, with margin: %.3f", zlow, zhigh, ZCUT_MARGIN);

	// todo: we may find a "closed" range, which may faster our process, but that's now the general case

	return geo::passThrough(cloud, "z", zlow + ZCUT_MARGIN, zhigh - ZCUT_MARGIN);
}

double CloudSegment::detectRoofHeight(PointCloud::Ptr cloud) const {
	const float THRESH_COSINE = std::cos(5. / 180. * geo::PI);

	auto subRoof = geo::filterPoints(cloud, [THRESH_COSINE](const Point& p) {
		return p.z > 0. && p.z / geo::length(p) > THRESH_COSINE;
	});
	LL_ASSERT(!subRoof->empty() && "failed to estimate roof height.");

	//todo: may ensure calibrate (0, 0, 1) here
	//todo: may remove outliers

	double z = ll::sum_by([](const Point& p) { return p.z; }, subRoof->points) / static_cast<double>(subRoof->size());
	LOG(INFO) << ll::unsafe_format("roof height estimated with %d points: %.3f", subRoof->size(), z);

	return z;
}

std::vector<CloudSegment::PlaneCloud> CloudSegment::detectPlanes(PointCloud::Ptr cloud,
	double disthresh, std::size_t inlier_count_thresh, std::size_t countthresh) const {
	std::vector<CloudSegment::PlaneCloud> planes;
	detectPlanesRecursively(cloud, planes, disthresh, inlier_count_thresh, countthresh);
	return planes;
}

void CloudSegment::detectPlanesRecursively(PointCloud::Ptr cloud, std::vector<PlaneCloud>& planes,
	double disthresh, std::size_t inlier_count_thresh, std::size_t countthresh) const {
	//todo: we may try speed up this process

	if (cloud->size() < inlier_count_thresh || planes.size() >= countthresh) return;

	pcl::SampleConsensusModelPlane<Point>::Ptr planeModel(new pcl::SampleConsensusModelPlane<Point>(cloud));
	pcl::RandomSampleConsensus<Point> ransac(planeModel);
	ransac.setDistanceThreshold(disthresh);
	ransac.computeModel();

	std::vector<int> inliers;
	Eigen::VectorXf params;
	ransac.getInliers(inliers);
	ransac.getModelCoefficients(params);

	if (inliers.size() < inlier_count_thresh) return;

	PlaneCloud pc;
	pc.cloud_ = geo::getSubSet(cloud, inliers);
	pc.abcd_ = params.block<4, 1>(0, 0);
	planes.push_back(pc);

	auto left = geo::getSubSet(cloud, inliers, true);

	LOG(INFO) << ll::unsafe_format("found a plane: %d + %d", inliers.size(), left->size());
	detectPlanesRecursively(left, planes, disthresh, inlier_count_thresh, countthresh);
}

std::vector<CloudSegment::PlaneCloud> CloudSegment::detectRegionPlanes(PointCloud::Ptr cloud, double anglediff, double curvediff, std::size_t min_points) const {
	pcl::search::KdTree<Point>::Ptr tree(new pcl::search::KdTree<Point>());
	tree->setInputCloud(cloud);

	// normals
	pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>());
	{
		pcl::NormalEstimation<Point, pcl::Normal> n;
		n.setInputCloud(cloud);
		n.setSearchMethod(tree);
		n.setKSearch(10);
		n.compute(*normals);
	}

	std::vector<pcl::PointIndices> clusters;
	{
		pcl::RegionGrowing<Point, pcl::Normal> reg;
		reg.setMinClusterSize(min_points);
		reg.setSearchMethod(tree);
		reg.setNumberOfNeighbours(30);
		reg.setInputCloud(cloud);
		reg.setInputNormals(normals);
		reg.setSmoothnessThreshold(anglediff);
		reg.setCurvatureThreshold(curvediff);

		reg.extract(clusters);

		LOG(INFO) << clusters.size() << " clusters found.";
	}

	auto planes = ll::mapf([&](const pcl::PointIndices& pi) {
		PlaneCloud pc;
		auto piece = geo::getSubSet(cloud, pi.indices);
		std::tie(pc.cloud_, pc.abcd_) = geo::refinePlanePattern(piece, 0.05f);
		return pc;
	}, clusters);

#if 0 // show segmentation
	pcl::visualization::PCLVisualizer viewer;

	auto add_cloud = [&viewer](const std::string& name, PointCloud::Ptr cloud, double r, double g, double b) {
		pcl::visualization::PointCloudColorHandlerCustom<Point> color(cloud, r, g, b);
		viewer.addPointCloud(cloud, color, name);
	};
	auto add_cloud_rc = [&viewer, &add_cloud](const std::string& name, PointCloud::Ptr cloud) {
		double r{ geo::random() }, g{ geo::random() }, b{ geo::random() };
		add_cloud(name, cloud, r * 255., g * 255., b * 255.);
	};

	// add_cloud("cloud", cloud, 0., 0., 150.);
	// viewer.addPointCloudNormals<Point, pcl::Normal>(cloud, normals, 10, 0.2f, "normals");
	for (auto pr : ll::enumerate(clusters)) {
		auto piece = geo::getSubSet(cloud, pr.iter->indices);
		add_cloud_rc("piece_" + std::to_string(pr.index), piece);
	}

	while (!viewer.wasStopped()) {
		viewer.spinOnce(33);
	}
#endif

	return planes;
}

std::vector<CloudSegment::Segment> CloudSegment::splitSegments(PointCloud::Ptr cloud, const Eigen::Vector3f vp, const Eigen::Vector3f& vn,
	float connect_thresh, float min_seg_length) const {
	const auto p = geo::P_(vp);
	const auto n = geo::P_(vn);

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

		std::size_t sid = *pr.first;
		std::size_t eid = *pr.second;
		if (std::fabs(projlens[eid] - projlens[sid]) < min_seg_length) continue;

		Segment s;
		s.s_ = geo::V_(cloud->points[sid]);
		s.e_ = geo::V_(cloud->points[eid]);

		segments.push_back(s);
	}

	LOG(INFO) << segments.size() << " / " << cnt << " segment(s) found.";

	return segments;
}

Eigen::vector<trans2d::Matrix2x3f> CloudSegment::computeSegmentAlignCandidates(const std::vector<Segment>& segments,
	const Eigen::vector<Eigen::Vector2f>& blueprint, float disthresh) const {
	struct can {
		trans2d::Matrix2x3f T_;
		float error_;

		float a() const { return std::atan2f(T_(1, 0), T_(0, 0)); }

		std::string to_string() const {
			std::stringstream ss;
			ss << error_ << " [" << (a() * 180.f / geo::PI) << ", (" << t().transpose() << ")]";
			return ss.str();
		}

	private:
		Eigen::Vector2f t() const { return -T_.block<2, 2>(0, 0).transpose() * T_.block<2, 1>(0, 2); }
	};

	auto minimum_distance_to_wall_segments = [&segments](const Eigen::Vector2f& p)-> float {
		auto pr = ll::min_by([&p](const Segment& seg) {
			Eigen::Vector2f s(seg.s_(0), seg.s_(1));
			Eigen::Vector2f e(seg.e_(0), seg.e_(1));
			return geo::distance_to_segment_2d(p, s, e);
		}, segments);
		return pr.second;
	};

	// sum of point distance to segment
	auto match_error = [&](const trans2d::Matrix2x3f& T) {
		auto mindis = ll::mapf([&](const Eigen::Vector2f& p) {
			Eigen::Vector2f q = trans2d::transform(p, T);
			return minimum_distance_to_wall_segments(q);
		}, blueprint);
		return ll::sum(mindis);
	};

	std::vector<can> cans;
	for (const auto& seg : segments) {
		Eigen::Vector2f s1(seg.s_(0), seg.s_(1));
		Eigen::Vector2f e1(seg.e_(0), seg.e_(1));

		float len = (e1 - s1).norm();

		for (std::size_t i = 0; i < blueprint.size(); ++i) {
			Eigen::Vector2f s2 = blueprint[i];
			Eigen::Vector2f e2 = blueprint[(i + 1) % blueprint.size()];

			if (std::fabs((e2 - s2).norm() - len) > 2. * disthresh) continue;

			trans2d::Matrix2x3f T = trans2d::estimateTransform(s1, e1, s2, e2);
			float err = match_error(T);

			can c;
			c.T_ = T;
			c.error_ = err;

			cans.push_back(c);
		}
	}

	// sort and output.
	std::sort(cans.begin(), cans.end(), [](const auto& c1, const auto& c2) { return c1.error_ < c2.error_; });

	// filter & unique.
	{
		constexpr float ANGLE_THRESH = 5.f / geo::PI;
		constexpr float TRANSLATE_THRESH = 0.2f;
		constexpr float MAX_AVG_DISTANCE = 1.; // this is very loose, i think

		auto is_same = [ANGLE_THRESH, TRANSLATE_THRESH](const can& ci, const can& cj) {
			float da = ci.a() - cj.a();
			if (da > geo::PI) da -= geo::TAU;
			else if (da < -geo::PI) da += geo::TAU;

			if (std::fabs(da) > ANGLE_THRESH) return false;

			// t dis
			Eigen::Matrix2f Ri = ci.T_.block<2, 2>(0, 0);
			Eigen::Vector2f ti = ci.T_.block<2, 1>(0, 2);
			Eigen::Vector2f tj = cj.T_.block<2, 1>(0, 2);

			float dt = (Ri.transpose() * (tj - ti)).norm();
			return dt < TRANSLATE_THRESH;
		};

		std::vector<can> unique_cans;
		for (const auto& ci : cans) {
			if (ci.error_ > MAX_AVG_DISTANCE* blueprint.size()) continue;

			auto iter = std::find_if(unique_cans.begin(), unique_cans.end(), [&ci, &is_same](const can& cj) { return is_same(ci, cj); });
			if (iter == unique_cans.end())
				unique_cans.push_back(ci);
		}
		cans.swap(unique_cans);

		LOG(INFO) << ll::unsafe_format("%d / %d unique candidates found.", cans.size(), unique_cans.size());
	}

	{
		std::stringstream ss;
		for (auto pr : ll::enumerate(cans)) ss << pr.index << ":\t" << pr.iter->to_string() << "\n";
		LOG(INFO) << "match errors: \n" << ss.str();
	}

	Eigen::vector<trans2d::Matrix2x3f> Ts;
	Ts.reserve(cans.size());
	for (const auto& c : cans) Ts.emplace_back(c.T_);

	return Ts;
}

CloudSegment::SegmentResult CloudSegment::segmentCloudByCADModel(PointCloud::Ptr thecloud) const {
	// we now process the ORIGINAL cloud
	constexpr double SLICE_HALF_THICKNESS = 0.1f;

	// simple pass through
	auto cloud = geo::filterPoints(thecloud, [&](const Point& p) {
		return (p.x > cadx1_ - SLICE_HALF_THICKNESS && p.x < cadx2_ + SLICE_HALF_THICKNESS) &&
			(p.y > cady1_ - SLICE_HALF_THICKNESS && p.y < cady2_ + SLICE_HALF_THICKNESS);
		// && (p.z > cadz1_ - SLICE_HALF_THICKNESS && p.z < cadz2_ + SLICE_HALF_THICKNESS);
	});

	// a, b, c should defines a plane
	auto get_slice = [&](const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c)-> PlaneCloud {
		Eigen::Vector3d n = (b - a).cross(c - a).normalized();
		auto slice = geo::filterPoints(cloud, [&a, &n, SLICE_HALF_THICKNESS](const Point& p) {
			Eigen::Vector3d ap = Eigen::Vector3d(p.x, p.y, p.z) - a;
			double dis = std::fabs(ap.dot(n));
			return dis <= SLICE_HALF_THICKNESS;
		});

		return detectRegionPlanes(slice, 5. / 180. * geo::PI, 1., slice->size() / 2).front();
	};

	SegmentResult sr;
	// each wall
	const auto& walls = cadModel_.getTypedModelItems(ITEM_WALL_E);
	sr.walls_.reserve(walls.size());
	for (const auto& wi : walls) {
		PlaneCloud pc = get_slice(wi.points_[0], wi.points_[1], wi.points_[2]);
		sr.walls_.emplace_back(pc);
	}

	// floor & roof
	const auto& roof = cadModel_.getTypedModelItems(ITEM_TOP_E).front();
	sr.roof_ = get_slice(roof.points_[0], roof.points_[1], roof.points_[2]);

	// floor may not aligned properly.
	// todo: we can use the pre-extracted planes
	{
		constexpr double MIN_WALL_HEIGHT = 1.;
		constexpr double MAX_WALL_HEIGHT = 5.;
		constexpr double HIST_STEP = 0.1;
		constexpr std::size_t HIST_COUNT = static_cast<std::size_t>((MAX_WALL_HEIGHT - MIN_WALL_HEIGHT) / HIST_STEP);

		// hist 0:0.1:5
		std::array<std::size_t, HIST_COUNT> hists;
		hists.fill(0);
		for (const auto& p : (*cloud)) {
			int i = (cadz2_ - MIN_WALL_HEIGHT - p.z) / HIST_STEP;
			if (i < 0 || i >= HIST_COUNT) continue;
			hists[i] += 1;
		}

		auto iter = std::max_element(hists.begin(), hists.end());
		std::size_t i = std::distance(hists.begin(), iter);
		double floorz1 = cadz2_ - MIN_WALL_HEIGHT - (i + 1) * HIST_STEP;
		double floorz2 = cadz2_ - MIN_WALL_HEIGHT - (i - 1) * HIST_STEP;
		LOG(INFO) << ll::unsafe_format("slice floor %.3f -> %.3f.", floorz1, floorz2);

		auto slice = geo::passThrough(cloud, "z", floorz1, floorz2);
		sr.floor_ = detectRegionPlanes(slice, 5. / 180. * geo::PI, 1., slice->size() / 2).front();
	}

	//todo: beams

#if 0
	pcl::visualization::PCLVisualizer viewer;

	auto add_cloud = [&viewer](const std::string& name, PointCloud::Ptr cloud, double r, double g, double b) {
		pcl::visualization::PointCloudColorHandlerCustom<Point> color(cloud, r, g, b);
		viewer.addPointCloud(cloud, color, name);
	};
	auto add_cloud_rc = [&viewer, &add_cloud](const std::string& name, PointCloud::Ptr cloud) {
		double r{ geo::random() }, g{ geo::random() }, b{ geo::random() };
		add_cloud(name, cloud, r * 255., g * 255., b * 255.);
	};

	add_cloud("cloud", cloud, 100., 100., 100.);

	for (const auto& pr : ll::enumerate(sr.walls_))
		add_cloud_rc("wall" + std::to_string(pr.index), pr.iter->cloud_);
	add_cloud_rc("roof", sr.roof_.cloud_);
	add_cloud_rc("floor", sr.floor_.cloud_);

	add_cloud("cad", cadModel_.genTestFrameCloud(), 255., 0., 0.);

	while (!viewer.wasStopped()) {
		viewer.spinOnce(33);
	}
#endif

	return sr;
}

PointCloud::Ptr CloudSegment::sparsedCloud() {
	if (!sparsedCloud_) sparsedCloud_ = geo::downsampleUniformly(orgCloud_, DOWNSAMPLE_SIZE);
	return sparsedCloud_;
}

}