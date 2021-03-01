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
#include <pcl/filters/crop_box.h>
#include <pcl/filters/radius_outlier_removal.h>
#include <pcl/features/normal_3d.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/search/search.h>
#include <pcl/filters/passthrough.h>
#include <pcl/filters/extract_indices.h>

#include "GeometryUtils.h"
#include "funHelper.h"

#include "SimpleViewer.h"

namespace CloudReg {

void CloudSegment::Segment::beCounterClockwise(const Eigen::Vector3f& cen) {
	double ps = atan2(s_(1) - cen(1), s_(0) - cen(0));
	double pe = atan2(e_(1) - cen(1), e_(0) - cen(0));
	if (pe < ps) std::swap(s_, e_);
}

std::string CloudSegment::SegmentResult::to_string() const {
	std::stringstream ss;
	ss << "SegmentResult: { \nplanes: ";
	for (int i = 0; i < ITEM_MAX_E; ++i) {
		ModelItemType mit = static_cast<ModelItemType>(i);

		auto iter = clouds_.find(mit);
		if (iter == clouds_.end()) ss << toModelItemName(mit) << ": 0 / 0";
		else {
			const auto& pcs = iter->second;
			std::size_t valid = std::count_if(pcs.begin(), pcs.end(), [](const PlaneCloud& p) { return p.cloud_; });
			ss << toModelItemName(mit) << ": " << valid << " / " << pcs.size();
		}

		ss << ", ";
	}
	ss << "\nT: \n" << T_ << "\n}";

	return ss.str();
}

CloudSegment::SegmentResult CloudSegment::run() {
	if (!calibrateDirectionToAxisZ()) {
		LOG(ERROR) << "failed to calibrate direction to axis-Z.";
	}

#ifdef VISUALIZATION_ENABLED
	if (1) {
		SimpleViewer viewer;
		viewer.addCloud(sparsedCloud());
		viewer.show();
	}
#endif

	recordModelBoundingBox();

	if (!alignCloudToCADModel()) {
		LOG(ERROR) << "failed to align cloud to cadModel.";
	}

	auto sr = segmentByCADModel();
	return sr;
}

bool CloudSegment::calibrateDirectionToAxisZ() {
	LOG(INFO) << "********calibrateDirectionToAxisZ*******";
	removeFarPoints(orgCloud_, cadModel_);

	pcl::PointCloud<pcl::PointXYZ>::Ptr filteredCloud(new pcl::PointCloud<pcl::PointXYZ>());
	pcl::PassThrough<pcl::PointXYZ> pass;
	pass.setInputCloud(orgCloud_);
	pass.setFilterFieldName("z");
	pass.setFilterLimits(0.5, 1.5);
	pass.filter(*filteredCloud);

	PointCloud::Ptr samplingCloud(new PointCloud());
	uniformSampling(0.1, filteredCloud, samplingCloud);

	float binSizeTh = 0.5;
	std::vector<std::pair<int, std::vector<int>>> zToNumVec;
	if (false == statisticsForPointZ(binSizeTh, samplingCloud, zToNumVec)) {
		LOG(WARNING) << "statistics For Point.Z Failed.";
		return false;
	}

	auto& vecIdxs = zToNumVec.begin()->second;
	auto subSet = geo::getSubSet(samplingCloud, vecIdxs, false);
	LOG(INFO) << "zToNumVec:" << zToNumVec.size() << " firstZ:"
		<< zToNumVec.begin()->first << " subSet:" << subSet->size();

	for (auto pr : ll::enumerate(zToNumVec)) {
		auto it = *pr.iter;
		LOG(INFO) << "z:" << it.first << " num:" << it.second.size();
		if (3 <= pr.index) break;
	}

	Eigen::VectorXf coeff;
	std::vector<int> inlierIdxs;
	planeFitting(0.1, subSet, coeff, inlierIdxs);
	Eigen::Vector4f plane(coeff(0), coeff(1), coeff(2), coeff(3));

	Eigen::Vector3f axisZ(0, 0, 1);
	float flag = (plane.block<3, 1>(0, 0).dot(axisZ) > 0) ? 1.0 : -1.0;
	plane = flag * plane;
	Eigen::Vector3f n = plane.block<3, 1>(0, 0);
	float angle = std::acos(n.dot(axisZ));
	Eigen::AngleAxisf rotation_vector(angle, n.cross(axisZ));
	Eigen::Matrix3f R = rotation_vector.matrix();
	LOG(INFO) << "angle:" << (angle * 180.0 / 3.14) << " n:"
		<< n(0) << "," << n(1) << "," << n(2) << " d:" << plane(3);
	LOG(INFO) << "calibrated to AxisZ, R: \n" << R;

	Eigen::Matrix4f T(Eigen::Matrix4f::Identity());
	T.block<3, 3>(0, 0) = R;

	//transform orig cloud
	orgCloud_ = geo::transfromPointCloud(orgCloud_, T);
#if 0
	pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud_filtered(new pcl::PointCloud<pcl::PointXYZ>());
	uniformSampling(0.01, orgCloud_, pCloud_filtered);
	pcl::io::savePCDFile("origin.pcd", *pCloud_filtered);
#endif
	return true;
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
	auto planes = detectRegionPlanes(cloud, 3. / 180. * geo::PI, 1., 500);

	//auto nzs = ll::mapf([](const PlaneCloud& pc) { return std::fabs(pc.n()(2)); }, planes);
	//std::sort(nzs.begin(), nzs.end());
	//LOG(INFO) << "nzs: " << ll::string_join(nzs, ", ");

	auto wallCandidates = ll::filter([](const PlaneCloud& pc) { return std::fabs(pc.n()(2)) < 0.1f; }, planes);
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
		for (auto& seg : allSegments) seg.beCounterClockwise(Eigen::Vector3f::Zero());

		std::stringstream ss;
		for (auto pr : ll::enumerate(allSegments)) ss << ll::unsafe_format("%02d: %.3f\n", pr.index, pr.iter->len());
		LOG(INFO) << allSegments.size() << " segments detected, length: \n" << ss.str();
	}

	//3. match
	// we need get segments on blueprint first.
	std::vector<Segment> blueprint;
	{
		Eigen::Vector3f floorCenter = Eigen::Vector3f::Zero();

		// floor contours
		{
			const auto& blueprint3d = cadModel_.getTypedModelItems(ITEM_BOTTOM_E).front();
			const Eigen::vector<Eigen::Vector3d>& points = blueprint3d.points_;
			for (std::size_t i = 0; i < points.size(); ++i) {
				Segment s;
				s.s_ = points[(i + 1) % points.size()].cast<float>();
				s.e_ = points[i].cast<float>();
				blueprint.push_back(s);

				floorCenter += points[i].cast<float>();
			}

			floorCenter /= static_cast<float>(points.size());
			LOG(INFO) << "floor center: " << floorCenter.transpose();
		}

		// beams
		if (cadModel_.containModels(ITEM_BEAM_E)) {
			const auto& beams = cadModel_.getTypedModelItems(ITEM_BEAM_E);
			LOG(INFO) << ll::unsafe_format("%d beams in cad(means %d actually).", beams.size(), beams.size() / 2);

			for (std::size_t bi = 0; bi < beams.size(); bi += 2) {
				const Eigen::vector<Eigen::Vector3d>& points = beams[bi].points_;
				// get proj
				Eigen::Vector3d a, b;
				bool valid = false;
				for (std::size_t i = 0; i < points.size(); ++i) {
					const Eigen::Vector3d& s = points[i];
					const Eigen::Vector3d& e = points[(i + 1) % points.size()];

					if ((e - s).block<2, 1>(0, 0).squaredNorm() < 1e-4) continue;

					if (!valid) {
						a = s;
						b = e;
						valid = true;
					} else {
						// merge
						double projlen = (s - a).dot(b - a);
						if (projlen < 0) a = s;
						else if (projlen > (b - a).block<2, 1>(0, 0).squaredNorm()) b = s;

						projlen = (e - a).dot(b - a);
						if (projlen < 0) a = e;
						else if (projlen > (b - a).block<2, 1>(0, 0).squaredNorm()) b = e;
					}
				}

				Segment s;
				s.s_ = a.cast<float>();
				s.e_ = b.cast<float>();
				s.s_(2, 0) = 0.;
				s.e_(2, 0) = 0.;
				s.beCounterClockwise(floorCenter);

				LOG(INFO) << "found a beam segment: " << s.s_.transpose() << " -> " << s.e_.transpose() << " with length: " << s.len();

				blueprint.push_back(s);
			}
		}

		LOG(INFO) << blueprint.size() << " segments found from blueprint.";
	}

	constexpr double MATCH_THRESH = 0.15;
	Eigen::vector<trans2d::Matrix2x3f> Ts = computeSegmentAlignCandidates(allSegments, blueprint, MATCH_THRESH);
	if (Ts.empty()) {
		LOG(ERROR) << " no candidate transform found.";
		return false;
	}

	trans2d::Matrix2x3f T2d = Ts.size() == 1 ? Ts.front() : chooseTransformByHoles(Ts, wallCandidates);

	Eigen::Matrix4f T = trans2d::asTransform3d(trans2d::inverse(T2d));
	T(2, 3) = -zfloor_; // shift up

	LOG(INFO) << "cloud trans found: \n" << T;

	//! note we do transform here
	sparsedCloud_ = geo::transfromPointCloud(sparsedCloud(), T);
	orgCloud_ = geo::transfromPointCloud(orgCloud_, T);
	// and we assign T_ here, 
	T_ = T;

#if VISUALIZATION_ENABLED	
	// debug show
	if (0) {
		SimpleViewer viewer;
		for (const auto& pc : wallCandidates) viewer.addCloud(pc.cloud_);

		for (const auto& seg : allSegments)
			viewer.addSegment(geo::P_(seg.s_), geo::P_(seg.e_));

		for (const auto& seg : blueprint)
			viewer.addSegment(geo::P_(seg.s_), geo::P_(seg.e_));

		viewer.addCloud(cadModel_.genTestFrameCloud(), 255., 0., 0.);

		viewer.show();
	}
#endif
	return true;
}

CloudSegment::SegmentResult CloudSegment::segmentByCADModel() {
	// bad thing is, due to efficiency, we cannot just segment the original cloud. 
	auto sr = segmentCloudByCADModel(sparsedCloud());

	// bounding box filter
	float x1{ METER_100 }, x2{ -METER_100 }, y1{ METER_100 }, y2{ -METER_100 }, z1{ METER_100 }, z2{ -METER_100 };
	{
		for (const auto& pr : sr.clouds_)
			for (const auto& pc : pr.second) {
				if (!pc.cloud_) continue;

				Point p1, p2;
				pcl::getMinMax3D(*pc.cloud_, p1, p2);

				x1 = std::min(x1, p1.x);
				y1 = std::min(y1, p1.y);
				z1 = std::min(z1, p1.z);
				x2 = std::max(x2, p2.x);
				y2 = std::max(y2, p2.y);
				z2 = std::max(z2, p2.z);
			}

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

	LOG(INFO) << ll::unsafe_format("would cluster %d clusters.", clusters.size());

	std::vector<int> searchIndices;
	std::vector<float> searchDis;
	constexpr double MAX_DIS_SQUARED = DOWNSAMPLE_SIZE * DOWNSAMPLE_SIZE * 2.25;
	constexpr float ON_PLANE_CHECK_THRESH = 0.05f;


	for (std::size_t idx = 0; idx < orgCloud_->size(); ++idx) {
		const auto& p = orgCloud_->points[idx];
		if (p.x< x1 || p.x> x2 || p.y< y1 || p.y> y2 || p.z< z1 || p.z> z2) continue;

		Eigen::Vector4f phom(p.x, p.y, p.z, 1.);

#if 0
		auto ins = ll::filter([&](std::size_t i) {
			const auto& pc = allplanes[i];
			if (!pc.cloud_ || pc.cloud_->empty()) return false;

			float dis = phom.transpose().dot(pc.abcd_);
			return std::fabs(dis) < ON_PLANE_CHECK_THRESH;
		}, ll::range(allplanes.size()));

		if (ins.empty()) continue;
		if (ins.size() == 1) {
			// we still need a dis check
			std::size_t i = ins.front();
			int cnt = trees[i].nearestKSearch(p, 1, searchIndices, searchDis);
			if (cnt > 0 && searchDis.front() < MAX_DIS_SQUARED)
				clusters[ins.front()].emplace_back(idx);
		} else {
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
#else
		// we use the one with minimal distance
		auto indices = ll::range(allplanes.size());
		auto planedis = ll::mapf([&](std::size_t i)->float {
			const auto& pc = allplanes[i];
			if (!pc.cloud_ || pc.cloud_->empty()) return false;

			float dis = phom.transpose().dot(pc.abcd_);
			return std::fabs(dis);
		}, indices);

		auto ins = ll::filter([&](std::size_t i) {
			if (planedis[i] > ON_PLANE_CHECK_THRESH) return false;
			int cnt = trees[i].nearestKSearch(p, 1, searchIndices, searchDis);
			return cnt > 0 && searchDis.front() < MAX_DIS_SQUARED;
		}, indices);

		if (ins.empty()) continue;

		std::size_t i = *(ll::min_by([&planedis](std::size_t i) { return planedis[i]; }, ins).first);

		clusters[i].emplace_back(idx);
#endif
	}

	{
		std::stringstream ss;
		for (std::size_t i = 0; i < allplanes.size(); ++i) {
			std::size_t size = allplanes[i].cloud_ ? allplanes[i].cloud_->size() : 0;
			ss << size << " -> " << clusters[i].size() << ", ";
		}
		LOG(INFO) << "re-cluster done: " << ss.str();
	}

	// out
	int offset = 0;
	for (auto& pr : sr.clouds_) {
		for (std::size_t i = 0; i < pr.second.size(); ++i) {
			pr.second[i].cloud_ = geo::getSubSet(orgCloud_, clusters[offset++]);
		}
	}

	_show_result(sr);

	return sr;
}

void CloudSegment::refineSegmentResult(SegmentResult& sr) const {
	//todo: refactor this.
	// use wall segment to slice floor & roof
	auto itw = sr.clouds_.find(ITEM_WALL_E);

	auto is_all_wall_found = [&]() {
		if (itw == sr.clouds_.end()) return false;
		const auto& walls = itw->second;
		return walls.size() > 2 &&
			std::all_of(walls.begin(), walls.end(), [](const PlaneCloud& pc) { return pc.cloud_; });
	};

	if (!is_all_wall_found()) {
		LOG(WARNING) << "not all of the walls was segmented, skip the refine process.";
		return;
	}

	// get the contour
	Eigen::vector<Eigen::Vector2f> points;
	{
#define _A(n) abcd##n(0)
#define _B(n) abcd##n(1)
#define _D(n) abcd##n(3)
		auto intersect_plane_on_xoy = [](const Eigen::Vector4f& abcd1, const Eigen::Vector4f& abcd2)->Eigen::Vector2f {
			// no error check
			Eigen::Vector2f p;
			p(0) = (_B(1) * _D(2) - _B(2) * _D(1)) / (_B(2) * _A(1) - _B(1) * _A(2));
			p(1) = (_A(1) * _D(2) - _A(2) * _D(1)) / (_A(2) * _B(1) - _A(1) * _B(2));
			return p;
		};
#undef  _A
#undef  _B
#undef  _D

		const auto& walls = itw->second;
		points.reserve(walls.size());
		for (std::size_t i = 0; i < walls.size(); ++i) {
			Eigen::Vector2f p = intersect_plane_on_xoy(walls[i].abcd_, walls[(i + 1) % walls.size()].abcd_);
			points.emplace_back(p);
		}
	}

	auto is_in_contour = [&](const Point& p) {
		float y = p.y;
		std::vector<float> xs;
		for (std::size_t i = 0; i < points.size(); ++i) {
			const Eigen::Vector2f& s = points[i];
			const Eigen::Vector2f& e = points[(i + 1) % points.size()];
			if (std::fabs(e(1) - s(1)) < 1e-6) continue;

			float t = (y - s(1)) / (e(1) - s(1));
			if (t > 0 && t < 1) {
				float x = (1 - t) * s(0) + t * e(0);
				xs.push_back(x);
			}
		}

		std::sort(xs.begin(), xs.end());
		for (std::size_t i = 0; i < xs.size(); ++i) {
			if (xs[i] > p.x) return i % 2 == 1;
		}

		return false;
	};

	auto itf = sr.clouds_.find(ITEM_BOTTOM_E);
	auto itr = sr.clouds_.find(ITEM_TOP_E);
	if (itf != sr.clouds_.end() && !itf->second.empty()) {
		PlaneCloud& pc = itf->second.front();
		if (pc.cloud_) {
			pc.cloud_ = geo::filterPoints(pc.cloud_, is_in_contour);
			LOG(INFO) << "floor refined: " << pc.cloud_->size() << " points left.";
		}
	}
	if (itr != sr.clouds_.end() && !itr->second.empty()) {
		PlaneCloud& pc = itr->second.front();
		if (pc.cloud_) {
			pc.cloud_ = geo::filterPoints(pc.cloud_, is_in_contour);
			LOG(INFO) << "roof refined: " << pc.cloud_->size() << " points left.";
		}
	}
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
	if (false && cadModel_.containModels(ITEM_BEAM_E)) {
		//^^^^ false: just let beams in the matching.
		auto beams = cadModel_.getTypedModelItems(ITEM_BEAM_E);
		for (const auto& beam : beams) {
			const auto& hr = beam.highRange_;
			if (hr.first < 1e-4) {} else zBeamLow = std::min(zBeamLow, hr.first);
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

std::vector<CloudSegment::PlaneCloud> CloudSegment::detectRegionPlanes(
	PointCloud::Ptr cloud, double anglediff, double curvediff, std::size_t min_points) const {
	return detectRegionPlanes(cloud, nullptr, anglediff, curvediff, min_points);
}

std::vector<CloudSegment::PlaneCloud> CloudSegment::detectRegionPlanes(
	PointCloud::Ptr cloud, pcl::PointCloud<pcl::Normal>::Ptr normals, double anglediff, double curvediff, std::size_t min_points) const {
	// normals
	if (!normals) normals = computeNormals(cloud);

	std::vector<pcl::PointIndices> clusters;
	{
		pcl::search::KdTree<Point>::Ptr tree(new pcl::search::KdTree<Point>());
		tree->setInputCloud(cloud);

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
	SimpleViewer viewer;

	// add_cloud("cloud", cloud, 0., 0., 150.);
	// viewer.addPointCloudNormals<Point, pcl::Normal>(cloud, normals, 10, 0.2f, "normals");
	for (auto pr : ll::enumerate(clusters)) {
		auto piece = geo::getSubSet(cloud, pr.iter->indices);
		viewer.addCloud("piece_" + std::to_string(pr.index), piece);
	}

	viewer.show();
#endif

	return planes;
}

trans2d::Matrix2x3f CloudSegment::chooseTransformByHoles(const Eigen::vector<trans2d::Matrix2x3f>& Ts, const std::vector<PlaneCloud>& walls) const {
	constexpr float HALF_SLICE_THICKNESS = 0.02f;
	constexpr float MATCH_DIS_CHECK = 0.1f;

	// T*src -> dst
	auto get_matched_num = [MATCH_DIS_CHECK](const Eigen::vector<Eigen::Vector2f>& src,
		const Eigen::vector<Eigen::Vector2f>& dst, const trans2d::Matrix2x3f& T) {
		Eigen::vector<Eigen::Vector2f> transed_src(src.size());
		std::transform(src.begin(), src.end(), transed_src.begin(), [&T](const Eigen::Vector2f& v) { return trans2d::transform(v, T); });

		return std::count_if(transed_src.begin(), transed_src.end(), [MATCH_DIS_CHECK, &dst](Eigen::Vector2f& sp) {
			return std::any_of(dst.begin(), dst.end(), [MATCH_DIS_CHECK, &sp](const Eigen::Vector2f& dp) { return (sp - dp).norm() < MATCH_DIS_CHECK; });
		});
	};

	// for the case when we cannot remove ambiguity
	trans2d::Matrix2x3f canT;
	std::size_t cancnt{ 0 };

	auto get_consistent_transformations_for = [&](const Eigen::vector<trans2d::Matrix2x3f>& checkTs, float z) {
		Eigen::vector<Eigen::Vector3d> cadHoleEnds = intersectCADModelOnZ(cadModel_, z);
		if (cadHoleEnds.empty()) {
			LOG(WARNING) << "bad function call: no intersection points found at: " << z;
			return checkTs;
		}

		// slice and detect ends
		Eigen::vector<Eigen::Vector3f> wallEnds;
		for (const auto& pc : walls) {
			auto line = geo::passThrough(pc.cloud_, "z", zfloor_ + z - HALF_SLICE_THICKNESS, zfloor_ + z + HALF_SLICE_THICKNESS);
			// we simply ignore z here.
			const Eigen::Vector4f& params = pc.abcd_;
			Eigen::Vector3f n = Eigen::Vector3f(-params(1), params(0), 0.f).normalized();
			Eigen::Vector3f p;
			if (std::fabs(params(0)) < 1e-6) p = Eigen::Vector3f(0.f, -params(3) / params(1), 0.f);
			else p = Eigen::Vector3f(-params(3) / params(0), 0.f, 0.f);

			auto subsegs = splitSegments(line, p, n, 0.1f, 0.1f);
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
			std::size_t cnt = get_matched_num(hole2d, wall2d, T);

			if (cnt == cadHoleEnds.size()) re.emplace_back(T);

			//SimpleViewer viewer;
			//for (const Eigen::Vector2f& p : hole2d) {
			//	Eigen::Vector2f tp = trans2d::transform(p, T);
			//	auto ps = Point(p(0), p(1), 0.);
			//	auto pd = Point(tp(0), tp(1), 0.);

			//	viewer.addPoint(ps, 255., 255., 0.);
			//	viewer.addPoint(pd, 255., 0., 0.);
			//	viewer.addLine(ps, pd, 255., 255., 255.);
			//}
			//for (const Eigen::Vector2f& p : wall2d)
			//	viewer.addPoint(Point(p(0), p(1), 0.), 0., 255., 0.);

			//viewer.show();
			if (cnt > cancnt) {
				cancnt = cnt;
				canT = T;
			}
		}

		return re;
	};

	// choose height
	if (!cadModel_.containModels(ITEM_HOLE_E)) {
		LOG(WARNING) << "no holes in cad model! choose the first one (with least max min distance).";
		return Ts.front();
	}

	const auto& holes = cadModel_.getTypedModelItems(ITEM_HOLE_E);
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
	// for (std::size_t i = 1; i < cell_counts.size() - 1; ++i) {
	for (std::size_t i = cell_counts.size() - 2; i >= 1; --i) {
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
	LOG(WARNING) << "failed to remove ambiguity, left candidates: " << leftTs.size() << "/" << Ts.size();

	if (cancnt > 0) {
		LOG(INFO) << "use the one with max match count: " << cancnt;
		return canT;
	} else {
		LOG(INFO) << "use the first one (with minimal matching error).";
		return leftTs.front();
	}
}

// get SORTED intersection points on cad
Eigen::vector<Eigen::Vector3d> CloudSegment::intersectCADModelOnZ(const CADModel& cadModel, float z) const {
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
	const std::vector<Segment>& blueprint, float disthresh) const {
	struct can {
		trans2d::Matrix2x3f T_;
		float error_;

#ifdef UBUNTU_SWITCH
		float a() const { return std::atan2(T_(1, 0), T_(0, 0)); }
#else
		float a() const { return std::atan2f(T_(1, 0), T_(0, 0)); }
#endif

		std::string to_string() const {
			std::stringstream ss;
			ss << error_ << " [" << (a() * 180.f / geo::PI) << ", (" << t().transpose() << ")]";
			return ss.str();
		}

	private:
		Eigen::Vector2f t() const { return -T_.block<2, 2>(0, 0).transpose() * T_.block<2, 1>(0, 2); }
	};

	Eigen::vector<Eigen::Vector2f> blueprintpoints;
	blueprintpoints.reserve(blueprint.size() * 2);
	for (const auto& s : blueprint) {
		blueprintpoints.emplace_back(s.s_.block<2, 1>(0, 0));
		blueprintpoints.emplace_back(s.e_.block<2, 1>(0, 0));
	}

	auto minimum_distance_to_wall_segments = [&segments](const Eigen::Vector2f& p)-> float {
		auto pr = ll::min_by([&p](const Segment& seg) {
			Eigen::Vector2f s(seg.s_(0), seg.s_(1));
			Eigen::Vector2f e(seg.e_(0), seg.e_(1));
			return geo::distance_to_segment_2d(p, s, e);
		}, segments);
		return pr.second;
	};

	// max-min point distance
	auto match_error = [&](const trans2d::Matrix2x3f& T)->float {
		auto mindis = ll::mapf([&](const Eigen::Vector2f& p) {
			Eigen::Vector2f q = trans2d::transform(p, T);
			return minimum_distance_to_wall_segments(q);
		}, blueprintpoints);
		return *std::max_element(mindis.begin(), mindis.end());
	};

	std::vector<can> cans;
	for (const auto& seg : segments) {
		Eigen::Vector2f s1 = seg.s_.block<2, 1>(0, 0);
		Eigen::Vector2f e1 = seg.e_.block<2, 1>(0, 0);

		float len = (e1 - s1).norm();

		for (const auto& bseg : blueprint) {
			Eigen::Vector2f s2 = bseg.s_.block<2, 1>(0, 0);
			Eigen::Vector2f e2 = bseg.e_.block<2, 1>(0, 0);

			if (std::fabs((e2 - s2).norm() - len) > 2. * disthresh) continue;

			trans2d::Matrix2x3f T = trans2d::estimateTransform(s1, e1, s2, e2);
			float err = match_error(T);

			if (err > disthresh) continue;

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
	constexpr double SLICE_HALF_THICKNESS = 0.1f;
	constexpr double NORMAL_CHECK = 0.866; // 30
	constexpr double GROWTH_ANGLE = 10. / 180. * geo::PI;

	// simple boxcrop
	auto cloud = geo::filterPoints(thecloud, [&](const Point& p) {
		return (p.x > cadx1_ - SLICE_HALF_THICKNESS && p.x < cadx2_ + SLICE_HALF_THICKNESS) &&
			(p.y > cady1_ - SLICE_HALF_THICKNESS && p.y < cady2_ + SLICE_HALF_THICKNESS);
		// && (p.z > cadz1_ - SLICE_HALF_THICKNESS && p.z < cadz2_ + SLICE_HALF_THICKNESS);
	});

	// normals
	auto normals = computeNormals(cloud);

	auto get_slice = [&](const Eigen::vector<Eigen::Vector3d>& points,
		const std::vector<Eigen::vector<Eigen::Vector3d>>& holes) {
		LL_ASSERT(points.size() > 2 && "?");

		// crop box, 
		// we may try use crop hull, however, the shape can not be convex.
		PointCloud::Ptr slice(new PointCloud());
		NormalCloud::Ptr slicenormals(new NormalCloud());
		{
			std::vector<int> nearidx;
			{
				Eigen::Vector4f minp = Eigen::Vector4f::Ones() * 100.;
				Eigen::Vector4f maxp = Eigen::Vector4f::Ones() * -100.;
				for (const Eigen::Vector3d& p : points) {
					for (int i = 0; i < 3; ++i) {
						if (p(i) < minp(i)) minp(i) = p(i);
						if (p(i) > maxp(i)) maxp(i) = p(i);
					}
				}
				Eigen::Vector4f d = Eigen::Vector4f::Ones() * SLICE_HALF_THICKNESS;
				minp += -d;
				maxp += d;
				minp(3) = 1.;
				maxp(3) = 1.;


				pcl::CropBox<Point> cb;
				cb.setMin(minp);
				cb.setMax(maxp);
				cb.setInputCloud(cloud);
				// cb.filter(*slice);
				cb.filter(nearidx);
			}

			// we then check their normal
			std::vector<int> indices;
			Eigen::Vector3d n = (points[1] - points[0]).cross(points[2] - points[1]).normalized();
			for (auto i : nearidx) {
				const pcl::Normal& pn = normals->points[i];
				Eigen::Vector3d pnv(pn.normal_x, pn.normal_y, pn.normal_z);

				if (std::fabs(pnv.dot(n)) > NORMAL_CHECK) indices.push_back(i);
			}

			slice = geo::getSubSet(cloud, indices);
			{
				pcl::ExtractIndices<pcl::Normal> ei;
				ei.setInputCloud(normals);
				pcl::PointIndices::Ptr pi(new pcl::PointIndices());
				pi->indices = indices;
				ei.setIndices(pi);
				ei.setNegative(false);
				ei.filter(*slicenormals);
			}

			LOG(INFO) << slice->size() << " points sliced.";
		}

		//todo: we can optimize the normal calculation in the below process
		auto planes = detectRegionPlanes(slice, slicenormals, GROWTH_ANGLE, 1., slice->size() / 4);
		if (planes.empty()) return PlaneCloud();

		if (0) {
			SimpleViewer viewer;
			viewer.addCloud(slice);
			for (const auto& p : planes) viewer.addCloud(p.cloud_);

			viewer.show();
		}

		if (planes.size() == 1)  return planes.front();

		std::stringstream ss;
		for (auto pr : ll::enumerate(planes))
			ss << pr.index << "] " << pr.iter->cloud_->size() << ", \t" << pr.iter->abcd_.transpose() << "\n";
		LOG(INFO) << "more than one planes were detected: \n" << ss.str();

		if (0) {
			SimpleViewer viewer;
			auto add_outline = [&viewer](const Eigen::vector<Eigen::Vector3d>& outline) {
				for (std::size_t i = 0; i < outline.size(); ++i) {
					auto s = geo::P_(outline[i].cast<float>());
					auto e = geo::P_(outline[(i + 1) % outline.size()].cast<float>());
					viewer.addSegment(s, e, 255., 255., 255.);
				}
			};

			for (const auto& plane : planes) viewer.addCloud(plane.cloud_);
			add_outline(points);
			for (const auto& hole : holes) add_outline(hole);

			viewer.show();
		}

		// we now check if we should merge them into one.
		std::sort(planes.begin(), planes.end(), [](const PlaneCloud& pc1, const PlaneCloud& pc2) {
			return pc1.cloud_->size() > pc2.cloud_->size();
		});
		auto base = planes.front();

		for (std::size_t i = 1; i < planes.size(); ++i) {
			const auto& plane = planes[i];
			double p = plane.abcd_.block<3, 1>(0, 0).cross(base.abcd_.block<3, 1>(0, 0)).norm();
			if (p > 0.02) continue;

			// parallel, then we check each points
			std::size_t inners = std::count_if(plane.cloud_->points.begin(), plane.cloud_->points.end(),
				[&points, &holes](const Point& p) {
				Eigen::Vector3d vp(p.x, p.y, p.z);
				return geo::isInShape2D(vp, points, holes);
			});

			float ratio = static_cast<float>(inners) / plane.cloud_->size();
			LOG(INFO) << ll::unsafe_format("inner ratio: %d / %d = %.2f %%", inners, plane.cloud_->size(), ratio*100.f);
			if (ratio > 0.7f) {
				base.cloud_->points.insert(base.cloud_->points.end(), plane.cloud_->points.begin(), plane.cloud_->points.end());
				LOG(INFO) << "a piece joined: " << plane.cloud_->size();
			}
		}

		base.cloud_->width = 1;
		base.cloud_->height = base.cloud_->size();

		return base;
	};

	// for roof only, copy from refineSegmentResult.
	auto is_in_contour = [](const Eigen::vector<Eigen::Vector3d>& points, const Point& p) {
		float y = p.y;
		std::vector<float> xs;
		for (std::size_t i = 0; i < points.size(); ++i) {
			const Eigen::Vector3d& s = points[i];
			const Eigen::Vector3d& e = points[(i + 1) % points.size()];
			if (std::fabs(e(1) - s(1)) < 1e-6) continue;

			float t = (y - s(1)) / (e(1) - s(1));
			if (t > 0 && t < 1) {
				float x = (1 - t) * s(0) + t * e(0);
				xs.push_back(x);
			}
		}

		std::sort(xs.begin(), xs.end());
		for (std::size_t i = 0; i < xs.size(); ++i) {
			if (xs[i] > p.x) return i % 2 == 1;
		}

		return false;
	};

	// todo: copy from get_slice, refactor these later.
	auto get_roof_slice = [&](const Eigen::vector<Eigen::Vector3d>& points) {
		LL_ASSERT(points.size() > 2 && "?");

		// crop box, 
		// we may try use crop hull, however, the shape can not be convex.
		PointCloud::Ptr slice(new PointCloud());
		NormalCloud::Ptr slicenormals(new NormalCloud());
		{
			std::vector<int> nearidx;
			{
				Eigen::Vector4f minp = Eigen::Vector4f::Ones() * 100.;
				Eigen::Vector4f maxp = Eigen::Vector4f::Ones() * -100.;
				for (const Eigen::Vector3d& p : points) {
					for (int i = 0; i < 3; ++i) {
						if (p(i) < minp(i)) minp(i) = p(i);
						if (p(i) > maxp(i)) maxp(i) = p(i);
					}
				}
				Eigen::Vector4f d = Eigen::Vector4f::Ones() * SLICE_HALF_THICKNESS;
				minp += -d;
				maxp += d;
				minp(3) = 1.;
				maxp(3) = 1.;


				pcl::CropBox<Point> cb;
				cb.setMin(minp);
				cb.setMax(maxp);
				cb.setInputCloud(cloud);
				// cb.filter(*slice);
				cb.filter(nearidx);
			}

			// we then check their normal
			std::vector<int> indices;
			Eigen::Vector3d n = (points[1] - points[0]).cross(points[2] - points[1]).normalized();
			for (auto i : nearidx) {
				const pcl::Normal& pn = normals->points[i];
				Eigen::Vector3d pnv(pn.normal_x, pn.normal_y, pn.normal_z);

				if (std::fabs(pnv.dot(n)) > NORMAL_CHECK) indices.push_back(i);
			}

			slice = geo::getSubSet(cloud, indices);
			{
				pcl::ExtractIndices<pcl::Normal> ei;
				ei.setInputCloud(normals);
				pcl::PointIndices::Ptr pi(new pcl::PointIndices());
				pi->indices = indices;
				ei.setIndices(pi);
				ei.setNegative(false);
				ei.filter(*slicenormals);
			}

			LOG(INFO) << slice->size() << " points sliced.";
		}

		//todo: we can optimize the normal calculation in the below process
		auto planes = detectRegionPlanes(slice, slicenormals, GROWTH_ANGLE, 1., 30);
		if (planes.empty()) return PlaneCloud();

		//if (1) {
		//	SimpleViewer viewer;
		//	// viewer.addCloud(slice);
		//	for (const auto& p : planes) viewer.addCloud(p.cloud_);
		//	// viewer.addCloud(mainplane.cloud_);

		//	viewer.show();
		//}

		auto iter = std::max_element(planes.begin(), planes.end(), [](const PlaneCloud& pc1, const PlaneCloud& pc2) {
			return pc1.cloud_->size() < pc2.cloud_->size();
		});
		PlaneCloud& mainplane = *iter;

		auto is_on_mainplane = [&mainplane](const Point& p) {
			float dis = mainplane.abcd_(0) * p.x + mainplane.abcd_(1) * p.y +
				mainplane.abcd_(2) * p.z + mainplane.abcd_(3);
			return std::fabs(dis) < 0.02f;
		};

		std::size_t maincnt = mainplane.cloud_->size();

		for (auto it = planes.begin(); it != planes.end(); ++it) {
			if (it == iter) continue;

			for (const auto& p : it->cloud_->points)
				if (is_on_mainplane(p) && is_in_contour(points, p)) mainplane.cloud_->points.push_back(p);
		}
		mainplane.cloud_->height = mainplane.cloud_->points.size();
		mainplane.cloud_->width = 1;

		LOG(INFO) << ll::unsafe_format("roof recollect: %d -> %d.", maincnt, mainplane.cloud_->points.size());

		return mainplane;
	};

	SegmentResult sr(T_);

	std::vector<Eigen::vector<Eigen::Vector3d>> noholes;
	const auto& walls = cadModel_.getTypedModelItems(ITEM_WALL_E);
	std::vector<std::vector<ModelItem>> wallHoles(walls.size());
	if (cadModel_.containModels(ITEM_HOLE_E)) {
		const auto& holes = cadModel_.getTypedModelItems(ITEM_HOLE_E);
		for (const auto& hole : holes) {
			if (hole.parentIndex_ >= 0 && hole.parentIndex_ < walls.size()) wallHoles[hole.parentIndex_].push_back(hole);
			else LOG(WARNING) << "invalid hole parent index: " << hole.parentIndex_;
		}
	}
	for (auto pr : ll::zip(walls.begin(), walls.end(), wallHoles.begin(), wallHoles.end())) {
		const auto& wall = *pr.first;
		const auto& holes = *pr.second;

		auto slice = get_slice(wall.points_, ll::mapf([](const ModelItem& mi) { return mi.points_; }, holes));
		sr.clouds_[ITEM_WALL_E].emplace_back(slice);
	}


	for (int i = 0; i < ITEM_MAX_E; ++i) {
		ModelItemType mit = static_cast<ModelItemType>(i);
		if (mit == ITEM_HOLE_E || mit == ITEM_WALL_E) continue; // we dont care about holes. and, we alreay procced the walls.

		if (mit == ITEM_BOTTOM_E) {
			// floor need handle carefully.
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
			double floorz1 = cadz2_ - MIN_WALL_HEIGHT - (i + 2) * HIST_STEP;
			double floorz2 = cadz2_ - MIN_WALL_HEIGHT - (i - 1) * HIST_STEP;
			LOG(INFO) << ll::unsafe_format("slice floor %.3f -> %.3f.", floorz1, floorz2);

			PointCloud::Ptr slice(new PointCloud());
			pcl::PointCloud<pcl::Normal>::Ptr slicenormals(new pcl::PointCloud<pcl::Normal>());
			{
				std::vector<int> indices;
				{
					for (int i = 0; i < cloud->size(); ++i) {
						const auto& p = cloud->points[i];
						if (p.z > floorz1&& p.z < floorz2) {
							const auto& n = normals->points[i];
							if (std::fabs(n.normal_z) > NORMAL_CHECK) indices.push_back(i);
						}
					}
				}

				slice = geo::getSubSet(cloud, indices);
				{
					pcl::ExtractIndices<pcl::Normal> ei;
					ei.setInputCloud(normals);
					pcl::PointIndices::Ptr pi(new pcl::PointIndices());
					pi->indices = indices;
					ei.setIndices(pi);
					ei.setNegative(false);
					ei.filter(*slicenormals);
				}
			}

			// auto slice = geo::passThrough(cloud, "z", floorz1, floorz2);
			//todo: we can still use box crop here.
			auto planes = detectRegionPlanes(slice, slicenormals, GROWTH_ANGLE, 1., slice->size() / 4);
			auto pc = *std::max_element(planes.begin(), planes.end(), [](const PlaneCloud& pc1, const PlaneCloud& pc2) {
				return pc1.cloud_->size() > pc2.cloud_->size();
			});

			sr.clouds_[mit].emplace_back(pc);
		}
		//else if (mit == ITEM_TOP_E) {
		//	for (const auto& mi : cadModel_.getTypedModelItems(mit))
		//		sr.clouds_[mit].emplace_back(get_roof_slice(mi.points_));
		//} 
		else {
			for (const auto& mi : cadModel_.getTypedModelItems(mit))
				sr.clouds_[mit].emplace_back(get_slice(mi.points_, noholes));
		}
	}

	LOG(INFO) << "segmented: " << sr.to_string();

	//refineSegmentResult(sr);
#ifdef VISUALIZATION_ENABLED
	_show_result(sr);
#endif

	return sr;
}

PointCloud::Ptr CloudSegment::sparsedCloud() {
	if (!sparsedCloud_) sparsedCloud_ = geo::downsampleUniformly(orgCloud_, DOWNSAMPLE_SIZE);
	return sparsedCloud_;
}

void CloudSegment::removeFarPoints(PointCloud::Ptr inputCloud, const CADModel& cad) {
	auto Bottons = cad.getTypedModelItems(ITEM_BOTTOM_E);
	auto Walls = cad.getTypedModelItems(ITEM_WALL_E);

	double width = 0;
	Eigen::Vector3d firstPt = Bottons.front().points_.front();
	for (auto& p : Bottons.front().points_) {
		double len = (firstPt - p).norm();
		width = std::max(width, len);
	}
	double height = std::fabs(Walls.front().highRange_.second - Walls.front().highRange_.first);

	double distSquareTh = width * width + height * height;
	if (distSquareTh < 1.0) distSquareTh = 20. * 20.;

	PointCloud::Ptr leftCloud(new PointCloud());
	for (auto& p : inputCloud->points) {
		double dist = p.x * p.x + p.y * p.y + p.z * p.z;
		if (dist < distSquareTh) leftCloud->push_back(p);
	}
	LOG(INFO) << "removeFarPoints, inputCloud:" << inputCloud->size() << " leftCloud:"
		<< leftCloud->size() << " distTh:" << std::sqrt(distSquareTh)
		<< " width:" << width << " height:" << height;
	inputCloud->swap(*leftCloud);
	return;
}

bool CloudSegment::statisticsForPointZ(float binSizeTh, PointCloud::Ptr cloud,
	std::vector<std::pair<int, std::vector<int>>>& zToNumVec) {
	struct comp {
		bool operator()(const std::pair<int, std::vector<int>> l, const std::pair<int, std::vector<int>> r) {
			return (l.second.size() > r.second.size());
		}
	};
	std::map<int, std::vector<int>> zToNumMap;
	for (int i = 0; i < cloud->size(); i++) {
		auto& p = cloud->points[i];
		int z = std::floor(p.z / binSizeTh);
		zToNumMap[z].push_back(i);
	}
	if (zToNumMap.empty()) {
		return false;
	}

	zToNumVec.clear();
	zToNumVec.insert(zToNumVec.end(), zToNumMap.begin(), zToNumMap.end());
	std::sort(zToNumVec.begin(), zToNumVec.end(), comp());
	return true;
}

CloudSegment::NormalCloud::Ptr CloudSegment::computeNormals(PointCloud::Ptr cloud) const {
	pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>());

	pcl::search::KdTree<Point>::Ptr tree(new pcl::search::KdTree<Point>());
	tree->setInputCloud(cloud);

	pcl::NormalEstimation<Point, pcl::Normal> n;
	n.setInputCloud(cloud);
	n.setSearchMethod(tree);
	n.setKSearch(10);
	n.compute(*normals);

	return normals;
}

void CloudSegment::_show_result(const SegmentResult& sr) const {
#if VISUALIZATION_ENABLED	
	SimpleViewer viewer;

	// viewer.addCloud(sparsedCloud_, 128., 128., 128.);
	for (const auto& pr : sr.clouds_)
		for (const auto& pc : pr.second)
			if (pc.cloud_)
				viewer.addCloud(pc.cloud_);

	viewer.addCloud(cadModel_.genTestFrameCloud(), 255., 0., 0.);

	viewer.show();
#endif
}

}