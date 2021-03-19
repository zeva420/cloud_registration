#pragma once

#include "BaseType.h"
#include "CalcMeasureHelper.h"

namespace CloudReg {

class WhitewashDesigner {
public:
	struct WallConstraint {
		std::size_t i_, j_; //! the order matters.
		double expectedDistance_; // pos_j - pos_i
		double lowBound_{ -0.1 }, highBound_{ 0.1 }; // actualDistance - expectedDistance_ \in [lowBound_, highBound_]
	};

	struct ConfigParams {
		double minSalientArea_{ 0. }; // low than this would not be considered
		double maxSalientHeight_{ 1. }; // salient higher than this must be clipped

		double minWallPaintThickness_{ 0.05f };
		double minSalientPaintThickness_{ 0.02f };
	};

	void inputData(const CloudItem& floor, const vecItems_t& walls);

	void config(const ConfigParams& cp) { params_ = cp; }

	bool solve();

	const std::vector<Wall>& getWalls() const { return walls_; }

private:
	ConfigParams params_;
	std::vector<Wall> walls_;
	std::vector<WallConstraint> constraints_;

	std::vector<Salient> detectSalients(const CloudItem& wall) const;


	void getWallConstraintPair(const std::vector<seg_pair_t>& rootBorder,
		const std::vector<seg_pair_t>& rootCADBorder);

	void setupWalls(std::vector<Wall>&& walls);
	void addConstraint(const WallConstraint& wc);

	float getSolidArea(const CloudItem& wall) const;
};

}
