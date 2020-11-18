#pragma once

#include "BaseType.h"

namespace CloudReg {
enum ModelItemType {
	ITEM_HOLE_E = 1,
	ITEM_BEAM_E,
	ITEM_BOTTOM_E,
	ITEM_WALL_E,
	ITEM_TOP_E,
	ITEM_MAX_E
};
struct ModelItem {
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		ModelItem(const ModelItemType& type) :itemtype_(type) {

	}

	ModelItem(const ModelItem& other) :itemtype_(other.itemtype_) {
		parentIndex_ = other.parentIndex_;
		points_ = other.points_;
		segments_ = other.segments_;
		highRange_ = other.highRange_;

	}


	using ItemPair_t = std::pair<Eigen::Vector3d, Eigen::Vector3d>;
	void buildSegment() {
		segments_.clear();

		if (points_.size() < 2) return;

		for (std::size_t i = 1; i < points_.size(); i++) {
			segments_.emplace_back(ItemPair_t(points_[i - 1], points_[i]));

		}
		segments_.emplace_back(ItemPair_t(points_.back(), points_.front()));
	}

	std::string toString() const;
	
	Eigen::vector<Eigen::Vector3d> points_;
	Eigen::vector<ItemPair_t> segments_;
	const ModelItemType itemtype_;
	std::size_t parentIndex_ = 9999;
	std::pair<double, double> highRange_ = std::make_pair(0, 0);

};
class CADModel {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	CADModel();
	~CADModel();

	bool initCAD(const std::string& fileName);

	// getters
	bool containModels(ModelItemType type) const{ return mapModelItem_.find(type) != mapModelItem_.end(); }

	const std::vector<ModelItem>& getTypedModelItems(ModelItemType type) const;

	std::vector<ModelItem> tryGetTypedModelItems(ModelItemType type) const;

	std::vector<ModelItem> filterModelItems(std::function<bool (const ModelItem&)> func) const;

	// test
	PointCloud::Ptr genTestCloud() const;
	std::string toString() const;

private:
	using vecItems_t = std::vector<ModelItem>;

	bool savePCD(const std::string& name, std::vector<ModelItem>& vec_item);

	std::map<ModelItemType, vecItems_t> mapModelItem_;
	Eigen::Vector3d centerPt_ = Eigen::Vector3d(0, 0, 0);
};
}
