#include "BaseType.h"
#include "GeometryUtils.h"
#include "ll.hpp"

namespace CloudReg {

Eigen::Vector3f Salient::getBoxMin() const { return geo::V_(bbp1_); }
Eigen::Vector3f Salient::getBoxMax() const { return geo::V_(bbp2_); }

std::vector<float> Wall::getSalientChipingHeight() const {
	std::vector<float> hs;
	if (salients_.empty()) return hs;

	hs.reserve(salients_.size());
	if (wallChipping_ > 0.f) {
		// wall chipped, then all salients need to be chipped.
		for (const auto& s : salients_) hs.emplace_back(s.height_ + wallChipping_);
	} else if (salientChipping_ > 0.f) {
		// check if any salient chipped
		float fh = maxSalientHeight_ - salientChipping_;
		for (const auto& s : salients_) hs.emplace_back(std::max(0.f, s.height_ - fh));
	} else {
		// just paint
		hs.resize(salients_.size(), 0.f);
	}

	return hs;
}

}
