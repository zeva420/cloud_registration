#pragma once

#include "BaseType.h"

namespace CloudReg {
enum ModelItemType {
	ITEM_HOLE_E = 0,
	ITEM_BEAM_E,
	ITEM_BOTTOM_E,
	ITEM_WALL_E,
	ITEM_TOP_E,
	ITEM_MAX_E
};
// debug func
inline std::string toModelItemName(ModelItemType type) {
	const std::array<std::string, ITEM_MAX_E> TYPENAME{
		"hole", "beam", "bottom", "wall", "top"
	};
	return TYPENAME[type];
}

struct ModelItem {
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		ModelItem(const ModelItemType& type) :itemtype_(type) {

	}

	ModelItem(const ModelItem& other){
		itemtype_ = other.itemtype_;
		parentIndex_ = other.parentIndex_;
		points_ = other.points_;
		segments_ = other.segments_;
		highRange_ = other.highRange_;
		area_ = other.area_;

	}



	void buildSegment() {
		segments_.clear();

		if (points_.size() < 2) return;

		for (std::size_t i = 1; i < points_.size(); i++) {
			segments_.emplace_back(seg_pair_t(points_[i - 1], points_[i]));

		}
		segments_.emplace_back(seg_pair_t(points_.back(), points_.front()));
	}

	std::string toString() const;
	
	Eigen::vector<Eigen::Vector3d> points_;
	std::vector<seg_pair_t> segments_;
	ModelItemType itemtype_;
	std::size_t parentIndex_ = 9999;
	std::pair<double, double> highRange_ = std::make_pair(0, 0);
	double area_ = 0.0; //only hole has this value

};
class CADModel {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	CADModel();
	~CADModel();

	bool initCAD(const std::string& fileName, const bool changeCADOrder);

	// getters
	bool containModels(ModelItemType type) const{ return mapModelItem_.find(type) != mapModelItem_.end(); }

	const std::vector<ModelItem>& getTypedModelItems(ModelItemType type) const;

	std::vector<ModelItem> tryGetTypedModelItems(ModelItemType type) const;

	std::vector<ModelItem> filterModelItems(std::function<bool (const ModelItem&)> func) const;

	// test
	PointCloud::Ptr genTestFrameCloud() const;
	std::map<ModelItemType, std::vector<PointCloud::Ptr>> genFragCloud(double delta=0.01) const;
	std::string toString() const;

	void scaleModel(const double scale);

	using vecItems_t = std::vector<ModelItem>;

	std::tuple<bool, ModelItem> genBottom(const std::vector<std::string>& vecSubStr);
	std::tuple<bool, ModelItem> genWall(const std::vector<std::string>& vecSubStr);
	std::tuple<bool, ModelItem> genHole(const std::vector<std::string>& vecSubStr);
	std::tuple<bool, ModelItem, ModelItem> genBeam(const std::vector<std::string>& vecSubStr);
	std::tuple<bool, ModelItem> genTop();

	void getAxis(const std::pair<Eigen::Vector3d, Eigen::Vector3d>& segment, 
		std::size_t& other_axis_index,double& start_axis, double& other_axis, bool& operate);

	void getAxis_clockwise(const std::pair<Eigen::Vector3d, Eigen::Vector3d>& segment,
		std::size_t& other_axis_index, double& start_axis, double& other_axis, bool& operate);

private:
	bool savePCD(const std::string& name, std::vector<ModelItem>& vec_item);
	void reSortWall();
	void cutWallByBeam();

	bool changeCADOrder_ = false;
	std::map<ModelItemType, vecItems_t> mapModelItem_;
	Eigen::Vector3d centerPt_ = Eigen::Vector3d(0, 0, 0);

	//warn: axis plane, need to in range.
	PointCloud::Ptr InterpolateShape(const Eigen::vector<Eigen::Vector3d>& points, 
		const std::vector<Eigen::vector<Eigen::Vector3d> >& holes, double delta) const;
};
}
