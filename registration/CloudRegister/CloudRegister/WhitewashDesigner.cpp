#include "WhitewashDesigner.h"
#include "ll.hpp"
#include "LPWrapper.h"
#include "glog/logging.h"

namespace CloudReg {

void WhitewashDesigner::setupWalls(std::vector<Wall>&& walls) {
	walls_.swap(walls);

	LL_ASSERT(std::all_of(walls_.begin(), walls_.end(), [](const Wall& w) { 
		return w.salientAreas_.size() == w.salientHeights_.size();
	}));
}

void WhitewashDesigner::addConstraint(const WallConstraint& wc) {
	LL_ASSERT(wc.i_ >= 0 && wc.i_ < walls_.size() && wc.j_ >= 0 && wc.j_ < walls_.size());
	LL_ASSERT(walls_[wc.i_].pos_ < walls_[wc.j_].pos_ && "wrong order.");

	constraints_.emplace_back(wc);
}

bool WhitewashDesigner::solve() {
	LL_ASSERT(!walls_.empty() && !constraints_.empty());

	// prepare
	for (auto& wall : walls_) {
		wall.maxSalientHeight_ = wall.salientHeights_.empty() ? 0. :
			*std::max_element(wall.salientHeights_.begin(), wall.salientHeights_.end());
	}

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
		d1[ii + 2] = 0.;

		c[ii] = wall.length_;
		c[ii + 1] = wall.length_ * 2.;
		c[ii + 2] = wall.length_ * 4.;
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

#define P(a) #a<< ":\n"<< a <<"\n"
	LOG(INFO) << "problem constructed: \n"
		<< P(A) << P(b) << P(c) << P(d1) << P(d2);
#undef P


	auto result = LPWrapper::Solve(c, A, b, d1, d2);
	LOG(INFO) << result.to_string();

	if (!result.isOptimal()) {
		LOG(WARNING) << "failed find optimal.";
		return false;
	}

	// fill result.

	return true;
}

}

