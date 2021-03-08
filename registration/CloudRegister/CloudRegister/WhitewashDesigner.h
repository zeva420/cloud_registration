#pragma once

#include "BaseType.h"
#include "CalcMeasureHelper.h"

namespace CloudReg {

class WhitewashDesigner {
public:
	struct Salient {
		float height_;
		float area_;
		Point bbp1_, bbp2_; // bounding box
	};
	struct Wall {
		double length_;
		double pos_;
		std::vector<Salient> salients_;
		double maxSalientHeight_{0.};

		double minWallPaintThickness_{ 0.05f };
		double minSalientPaintThickness_{ 0.02f };
	};
	struct WallConstraint {
		std::size_t i_, j_; //! the order matters.
		double expectedDistance_; // pos_j - pos_i
		double lowBound_{ -0.1 }, highBound_{ 0.1 }; // actualDistance - expectedDistance_ \in [lowBound_, highBound_]
	};

	struct WallGuide {
		double paintThickness_;
		double salientChipping{ 0. };
		double wallChipping_{ 0. };
	};

	void inputData(const CloudItem& floor, const vecItems_t& walls);

	void getWallConstraintPair(const std::vector<seg_pair_t>& rootBorder,
		const std::vector<seg_pair_t>& rootCADBorder);

	void setupWalls(std::vector<Wall>&& walls);
	void addConstraint(const WallConstraint& wc);

	bool solve();

	const std::vector<WallGuide>& getGuidelines() const { return guidelines_; }

private:
	std::vector<Wall> walls_;
	std::vector<WallConstraint> constraints_;
	std::vector<WallGuide> guidelines_;

	std::vector<Salient> detectSalients(const CloudItem& wall) const;
};

}
