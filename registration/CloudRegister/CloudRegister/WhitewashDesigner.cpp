#include "WhitewashDesigner.h"
#include "ll.hpp"
#include "LPWrapper.h"
#include "glog/logging.h"

namespace CloudReg {

std::string to_string(const Salient& s) {
	return ll::unsafe_format("Salient: {height: %.3f, area: %.3f, [%.3f, %.3f, %.3f] -> [%.3f, %.3f, %.3f]}",
		s.height_, s.area_, s.bbp1_.x, s.bbp1_.y, s.bbp1_.z, s.bbp2_.x, s.bbp2_.y, s.bbp2_.z);
}

void WhitewashDesigner::setupWalls(std::vector<Wall>&& walls) {
	walls_.swap(walls);
}

void WhitewashDesigner::addConstraint(const WallConstraint& wc) {
	LL_ASSERT(wc.i_ >= 0 && wc.i_ < walls_.size() && wc.j_ >= 0 && wc.j_ < walls_.size());
	// LL_ASSERT(&& "wrong order.");
	if (walls_[wc.i_].pos_ > walls_[wc.j_].pos_) {
		LOG(INFO) << " WallConstraint should be reversed.";
	} else constraints_.emplace_back(wc);
}

bool WhitewashDesigner::solve() {
	LL_ASSERT(!walls_.empty() && !constraints_.empty());

	// 
	const int N = walls_.size();
	const int COLS = N * 3;
	const int ROWS = walls_.size() * 2 + constraints_.size() * 2;

	Eigen::MatrixXd A(ROWS, COLS);
	Eigen::VectorXd b(ROWS);
	Eigen::VectorXd c(COLS), d1(COLS), d2(COLS);

	A.fill(0.);
	// [-1, 1] should be enough
	d1.fill(-1.);
	d2.fill(1.);

	// each wall
	for (int i = 0; i < N; ++i) {
		const auto& wall = walls_[i];
		int r = i * 2;
		int ii = i * 3;

		A(r, ii) = -1.;
		A(r, ii + 1) = -1.;
		A(r, ii + 2) = -1.;
		b(r) = -(wall.minSalientPaintThickness_ + wall.maxSalientHeight_);

		A(r + 1, ii) = -1.;
		A(r + 1, ii + 2) = -1.;
		b(r + 1) = -wall.minWallPaintThickness_;

		d1[ii + 1] = 0.;
		d2[ii + 1] = wall.maxSalientHeight_;
		d1[ii + 2] = 0.;

		c[ii] = wall.length_;
		c[ii + 1] = wall.length_ * 4.;
		c[ii + 2] = wall.length_ * 8.;
	}

	// each wall pair
	for (std::size_t i = 0; i < constraints_.size(); ++i) {
		const auto& wc = constraints_[i];
		double pij = walls_[wc.j_].pos_ - walls_[wc.i_].pos_;

		int r = N * 2 + i * 2;

		A(r, wc.i_ * 3) = 1.;
		A(r, wc.j_ * 3) = 1.;
		b(r) = pij - wc.expectedDistance_ - wc.lowBound_;

		A(r + 1, wc.i_ * 3) = -1.;
		A(r + 1, wc.j_ * 3) = -1.;
		b(r + 1) = -(pij - wc.expectedDistance_ - wc.highBound_);
	}

	//	Eigen::IOFormat fmt(4, 0, ", ", "\n", "[", "]");
	//#define P(a) #a<< ": "<< a.transpose().format(fmt) <<"\n"
	//	LOG(INFO) << "problem constructed: \na:\n" << A << "\n"
	//		<< P(b) << P(c) << P(d1) << P(d2);
	//#undef P


	auto result = LPWrapper::Solve(c, A, b, d1, d2);
	LOG(INFO) << result.to_string();

	if (!result.isOptimal()) {
		LOG(WARNING) << "failed to find optimal.";
		return false;
	}

	// fill result.
	// guidelines_.reserve(N);
	for (std::size_t i = 0; i < N; ++i) {
		float x = result.x_(i * 3);
		float y = result.x_(i * 3 + 1);
		float z = result.x_(i * 3 + 2);

		if (z > 0.f) y = walls_[i].maxSalientHeight_;

		auto& wall = walls_[i];
		wall.paintThickness_ = x + z;
		wall.salientChipping_ = y;
		wall.wallChipping_ = z;
	}

	// log result.
	{
		std::stringstream ss;
		for (auto pr : ll::enumerate(walls_))
			ss << ll::unsafe_format("[wall %d] chip wall %.3f, chip salient %.3f, paint wall %.3f.\n",
			pr.index, pr.iter->wallChipping_, pr.iter->salientChipping_, pr.iter->paintThickness_);

		ss << "\n";

		auto fnl_pos = [&](std::size_t i, bool left) {
			return walls_[i].pos_ + (-walls_[i].wallChipping_ + walls_[i].paintThickness_) * (left ? 1 : -1);
		};
		for (const auto& wc : constraints_) {
			float pij = (fnl_pos(wc.j_, false) - fnl_pos(wc.i_, true));
			ss << ll::unsafe_format("[pair %d, %d] expect %.3f, final %.3f. diff %.3f in [%.3f, %.3f].\n",
				wc.i_, wc.j_, wc.expectedDistance_, pij, (pij - wc.expectedDistance_), wc.lowBound_, wc.highBound_);
		}

		LOG(INFO) << "so the whitewash paint: \n" << ss.str();
	}

	return true;
}
}

#include "GeometryUtils.h"
#include "SimpleViewer.h"
#include <pcl/kdtree/kdtree.h>
#include <pcl/segmentation/extract_clusters.h>

namespace CloudReg {
void WhitewashDesigner::inputData(const CloudItem& floor, const vecItems_t& walls) {
	LL_ASSERT(walls.size() == floor.cadBorder_.front().size() && "wall size mismatch.");

	walls_.resize(walls.size());
	for (std::size_t i = 0; i < walls.size(); ++i) {
		const auto& wall = walls[i];
		auto& w = walls_[i];

		w.salients_ = detectSalients(wall);
		if (!w.salients_.empty()) {
			for (const auto& s : w.salients_) {
				if (s.area_ > params_.minSalientArea_)
					w.maxSalientHeight_ = std::max<double>(w.maxSalientHeight_, s.height_);
			}
			w.maxSalientHeight_ = std::min(w.maxSalientHeight_, params_.maxSalientHeight_); // just ignore
		}

		double err = { 0. };
		const Eigen::Vector4d& params = wall.cloudPlane_;
		const double a{ params(0) }, b{ params(1) }, c{ params(2) }, d{ params(3) };
		if (std::fabs(a) < 0.1f) {
			w.pos_ = -d / b; // y
			err = std::sqrt(a * a + c * c);
		} else if (std::fabs(b) < 0.1f) {
			w.pos_ = -d / a;
			err = std::sqrt(b * b + c * c);
		} else LOG(ERROR) << "?";

		LOG(INFO) << ll::unsafe_format("direction error of wall %d: %.3f", i, err);
	}

	// length set refer to CAD
	const auto& segs = floor.cadBorder_.front();
	for (std::size_t i = 0; i < segs.size(); ++i) {
		walls_[i].length_ = (segs[i].second - segs[i].first).norm();
	}
	LOG(INFO) << ll::unsafe_format("%d walls setup.", walls_.size());

	// now detect constraints and set estimate position
	getWallConstraintPair(floor.cloudBorder_.front(), floor.cadBorder_.front());

	// debug
	{
		std::stringstream ss;
		ss << ll::unsafe_format("%d walls, %d constrains: \n", walls_.size(), constraints_.size());
		for (auto pr : ll::enumerate(walls_))
			ss << ll::unsafe_format("[wall %d] pos: %.6f, salient: %.6f, length: %.6f.\n",
			pr.index, pr.iter->pos_, pr.iter->maxSalientHeight_, pr.iter->length_);
		ss << "-----------------------------------------------\n";
		for (const auto& wc : constraints_)
			ss << ll::unsafe_format("[pair %d, %d] distance: %.6f, [%.6f, %.6f].\n", wc.i_, wc.j_, wc.expectedDistance_, wc.lowBound_, wc.highBound_);
		ss << "===============================================\n";

		LOG(INFO) << "WhitewashDesigner input state: \n" << ss.str();
	}
}

std::vector<Salient> WhitewashDesigner::detectSalients(const CloudItem& wall) const {
	constexpr float HUMPED_THRESH = 0.002f;
	constexpr float DOWNSAMPLE = 0.01f;
	constexpr float CUSTER_DISTANCE = 0.03f;
	constexpr int MIN_HUMPED_POINTS = 100;

	//todo: this may need to return with sign.
	auto dis_to_plane = [&](const Point& p) {
		Eigen::Vector4d vp(p.x, p.y, p.z, 1.);
		float dis = std::fabs(vp.dot(wall.cloudPlane_));
		return dis;
	};

	auto sparsed = geo::downsampleUniformly(wall.pCloud_, DOWNSAMPLE);

	const float den = sparsed->size() / getSolidArea(wall);
	LOG(INFO) << ll::unsafe_format("wall density: %.3f (count in each m^2), num thresh: %d.", den, MIN_HUMPED_POINTS);

	auto hump = geo::filterPoints(sparsed, [&](const Point& p) { return dis_to_plane(p) > HUMPED_THRESH; });

	std::vector<Salient> salients;
	std::vector<PointCloud::Ptr> salientClouds;
	{
		pcl::search::KdTree<Point>::Ptr tree(new pcl::search::KdTree<Point>);
		tree->setInputCloud(hump);

		std::vector<pcl::PointIndices> clusters;
		pcl::EuclideanClusterExtraction<Point> ece;
		ece.setClusterTolerance(0.1f);
		ece.setSearchMethod(tree);
		ece.setInputCloud(hump);
		ece.extract(clusters);

		if (!clusters.empty()) {
			for (const auto& indices : clusters) {
				const std::size_t CNT = indices.indices.size();
				if (CNT < MIN_HUMPED_POINTS) continue;

				Eigen::Vector4f v1, v2;
				pcl::getMinMax3D(*hump, indices, v1, v2);

				float h = v2(2) - v1(2);
				float w = std::fabs(wall.cloudPlane_(0, 0)) < 0.1 ? (v2(0) - v1(0)) : (v2(1) - v1(1));
				float hwr = h / w;
				if (hwr < 1.f) hwr = 1.f / hwr;
				float szr = CNT / (w * h);
				LOG(INFO) << ll::unsafe_format("%d points in %.3f x %.3f area: (%.3f, %.3f)", CNT, w, h, hwr, szr);

				if (szr < 300.f || hwr>9.f) continue;
				salientClouds.push_back(geo::getSubSet(hump, indices.indices));

				Salient s;
				s.bbp1_ = geo::P_(v1.block<3, 1>(0, 0));
				s.bbp2_ = geo::P_(v2.block<3, 1>(0, 0));
				s.height_ = ll::max_by([&](int i) {
					const auto& p = hump->points[i];
					return dis_to_plane(p);
				}, indices.indices).second;
				s.area_ = indices.indices.size() / den;

				salients.emplace_back(s);
			}
		}
	}

	{
		std::stringstream ss;
		for (const auto& s : salients) ss << to_string(s) << "\n";
		LOG(INFO) << salients.size() << " salients detected:\n" << ss.str();
	}

#if VISUALIZATION_ENABLED
	SimpleViewer viewer;
	viewer.addCloud(sparsed, 100., 100., 100.);
	// viewer.addCloud(hump, 100., 100., 100.);

	for (const auto& s : salients) {
		double r{ geo::random() }, g{ geo::random() }, b{ geo::random() };
		viewer.addBox(s.bbp1_, s.bbp2_, r, g, b);
	}

	for (const auto& sc : salientClouds) viewer.addCloud(sc);

	viewer.show();
#endif

	return salients;
}

float WhitewashDesigner::getSolidArea(const CloudItem& wall) const {
	LL_ASSERT(wall.type_ == CLOUD_WALL_E);

	const auto& border = wall.cadBorder_;

	auto calc_area = [](const std::vector<seg_pair_t>& segs)->float {
		if (segs.size() == 4)
			return (segs[0].second - segs[0].first).norm() * (segs[1].second - segs[1].first).norm();
		LOG(INFO) << "seg size!=4, need to calc polygon.";
		return 0.f;
	};

	float area = calc_area(wall.cadBorder_.front());
	for (std::size_t i = 1; i < wall.cadBorder_.size(); ++i)
		area -= calc_area(wall.cadBorder_[i]);

	return area;
}

void WhitewashDesigner::getWallConstraintPair(const std::vector<seg_pair_t>& rootCloudBorder,
	const std::vector<seg_pair_t>& rootCADBorder) {
	LL_ASSERT(walls_.size() == rootCADBorder.size() && "should setup walls first.");

	LOG(INFO) << "getWallConstraintPair";

	const double calcLengthTh = 0.01;
	const Eigen::Vector3d& horizenSeg = rootCADBorder.front().first - rootCADBorder.front().second;
	std::vector<std::size_t> vecVerticalIndex;
	std::vector<std::size_t> vecHorizenIndex;
	groupDirectionIndex(horizenSeg, rootCADBorder, vecVerticalIndex, vecHorizenIndex);

	std::map<std::pair<std::size_t, std::size_t>, double> mapWallPair;
	auto getWallPair = [&](const std::vector<std::size_t>& calcIndex) {

		for (std::size_t i = 0; i < calcIndex.size(); i++) {
			seg_pair_t toSeg = rootCADBorder[calcIndex[i]];
			if ((toSeg.first - toSeg.second).norm() < calcLengthTh)
				continue;


			for (int j = calcIndex.size() - 1; j >= i + 1; j--) {
				seg_pair_t calcSeg = rootCADBorder[calcIndex[j]];
				if ((calcSeg.first - calcSeg.second).norm() < calcLengthTh)
					continue;

				bool hasOverlap;
				Eigen::Vector3d s1Pt, e1Pt, s2Pt, e2Pt;
				std::tie(hasOverlap, s1Pt, e1Pt, s2Pt, e2Pt) = calcOverlap(toSeg, calcSeg);


				if (!hasOverlap || (s1Pt - e1Pt).norm() < calcLengthTh
					|| (s2Pt - e2Pt).norm() < calcLengthTh)
					continue;

				std::size_t optIndex, indexOther;
				int dir;
				std::tie(optIndex, indexOther, dir) = getWallGrowAxisAndDir(s1Pt, e1Pt);

				mapWallPair[std::make_pair(calcIndex[i], calcIndex[j])] = s2Pt[indexOther] - s1Pt[indexOther];
			}
		}
	};

	getWallPair(vecVerticalIndex);
	getWallPair(vecHorizenIndex);

	std::vector<double> vecWallLength(rootCADBorder.size(), 0.0);
	std::vector<WallConstraint> vecConstraint;
	for (auto& pairIdx : mapWallPair) {
		LOG(INFO) << pairIdx.first.first << ", " << pairIdx.first.second << ": " << pairIdx.second;

		WallConstraint tmp;
		tmp.i_ = pairIdx.first.first;
		tmp.j_ = pairIdx.first.second;
		tmp.expectedDistance_ = pairIdx.second;

		// auto& segI = rootCloudBorder[tmp.i_];
		// auto& segJ = rootCloudBorder[tmp.j_];

		// std::size_t optIndex, indexOther;
		// int dir;
		// std::tie(optIndex, indexOther, dir) = getWallGrowAxisAndDir(segI.first, segI.second);

		// tmp.lowBound_ = segI.first[indexOther];
		// tmp.highBound_ = segJ.second[indexOther];
		// vecConstraint.emplace_back(tmp);
		if (tmp.expectedDistance_ < 0.) {
			std::swap(tmp.i_, tmp.j_);
			tmp.expectedDistance_ = -tmp.expectedDistance_;
		}

		addConstraint(tmp);

		// vecWallLength[tmp.i_] = (segI.first - segI.second).norm();
		// vecWallLength[tmp.j_] = (segJ.first - segJ.second).norm();

		// LOG(INFO) << "WallIdx: " << tmp.i_ << " - " << tmp.j_ << " expect_dis:" << tmp.expectedDistance_
		//	<< " posI: " << tmp.lowBound_ << " posJ: " << tmp.highBound_ << " dis:" << (tmp.highBound_ - tmp.lowBound_);

		// LOG(INFO) << vecWallLength[tmp.i_] << " " << vecWallLength[tmp.j_];
	}

}


}

