#pragma once

#include "BaseType.h"
#include "CalcMeasureHelper.h"

namespace CloudReg {

class WhitewashDesigner {
public:
	struct Wall {
		double length_;
		double pos_;
		std::vector<double> salientHeights_, salientAreas_;
		double maxSalientHeight_;

		double minWallPaintThickness_;
		double minSalientPaintThickness_;
	};
	struct WallConstraint {
		std::size_t i_, j_; //! the order matters.
		double expectedDistance_; // pos_j - pos_i
		double lowBound_, highBound_;
	};

	struct WallGuide {
		double paintThickness_;
		double salientChipping{ 0. };
		double wallChipping_{ 0. };
	};

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
};

}
