#include "CADModel.h"

#include "funHelper.h"

namespace CloudReg {

// debug func
inline std::string toString(ModelItemType type) {
	const std::array<std::string, ITEM_MAX_E> TYPENAME{
		"hole", "beam", "bottom", "wall","top", "--"
	};
	return TYPENAME[type - ITEM_HOLE_E];
}

std::string ModelItem::toString() const {
	std::stringstream ss;
	ss << "ModelItem{ type: " << CloudReg::toString(itemtype_) << ", "
		<< points_.size() << " points(" << segments_.size() << " segments), parent: " << parentIndex_ << ", "
		<< "height range: " << highRange_.first << " -> " << highRange_.second << "}";
	return ss.str();
}

CADModel::CADModel() {}

CADModel::~CADModel() {
	mapModelItem_.clear();
}



bool CADModel::initCAD(const std::string& fileName) {
	mapModelItem_.clear();


	mapModelItem_[ITEM_BOTTOM_E] = std::move(std::vector<ModelItem>());
	mapModelItem_[ITEM_WALL_E] = std::move(std::vector<ModelItem>());
	mapModelItem_[ITEM_HOLE_E] = std::move(std::vector<ModelItem>());
	mapModelItem_[ITEM_BEAM_E] = std::move(std::vector<ModelItem>());
	mapModelItem_[ITEM_TOP_E] = std::move(std::vector<ModelItem>());

	std::vector<ModelItem> vec_item;
	double maxX = 0, maxY = 0, minX = 0, minY = 0;
	std::ifstream in(fileName);
	if (in) {
		std::string line;
		while (getline(in, line)) {
			std::vector<std::string> vecSubStr = splitByCharacter(line, ',');
			if (vecSubStr.empty()) {
				LOG(ERROR) << "splite string error";
				return false;
			}

			if (vecSubStr[0] == "B") {
				ModelItem item(ITEM_BOTTOM_E);
				std::size_t number = atol(vecSubStr[1].c_str());
				if (number < 1) {
					LOG(ERROR) << "invaild number in ITEM_BOTTOM_E";
					return false;
				}

				
				for (std::size_t i = 2; i < number * 2 + 1; ) {
					Eigen::Vector3d point;
					
					point[0] = atol(vecSubStr[i].c_str());
					point[1] = atol(vecSubStr[i + 1].c_str());
					point[2] = 0;													

					item.points_.emplace_back(point);

					if (point[0] > maxX) maxX = point[0];
					if (point[1] > maxY) maxY = point[1];
					if (point[0] < minX) minX = point[0];
					if (point[1] < minY) minY = point[1];

					i += 2;

				}

				centerPt_[0] = (maxX + minX) / 2;
				centerPt_[1] = (maxY + minY) / 2;
				
				
				item.buildSegment();

				vec_item.emplace_back(item);
				mapModelItem_[ITEM_BOTTOM_E].emplace_back(item);
				const ModelItem& bottom = mapModelItem_[ITEM_BOTTOM_E].front();


			} 
			else if (vecSubStr[0] == "W") {
				if (mapModelItem_[ITEM_BOTTOM_E].empty()) {
					LOG(ERROR) << "not found ITEM_BOTTOM_E";
					return false;
				}
				const ModelItem& bottom = mapModelItem_[ITEM_BOTTOM_E].front();

				ModelItem item(ITEM_WALL_E);
				std::size_t parent = atol(vecSubStr[1].c_str());
				std::size_t high = atol(vecSubStr[2].c_str());



				Eigen::Vector3d ptA = bottom.segments_[parent].first;
				Eigen::Vector3d ptAA = ptA + Eigen::Vector3d(0, 0, high);
				Eigen::Vector3d ptB = bottom.segments_[parent].second;
				Eigen::Vector3d ptBB = ptB + Eigen::Vector3d(0, 0, high);


				item.points_.emplace_back(ptAA);
				item.points_.emplace_back(ptA);
				item.points_.emplace_back(ptB);
				item.points_.emplace_back(ptBB);

				item.highRange_ = std::make_pair(0, high);
				item.buildSegment();
				vec_item.emplace_back(item);
				mapModelItem_[ITEM_WALL_E].emplace_back(item);

			} 
			else if (vecSubStr[0] == "H") {
				if (mapModelItem_[ITEM_BOTTOM_E].empty() || mapModelItem_[ITEM_WALL_E].empty()) {
					LOG(ERROR) << "not found ITEM_BOTTOM_E";
					return false;
				}

				const std::size_t parent = mapModelItem_[ITEM_WALL_E].size() - 1;
				const ModelItem& botton = mapModelItem_[ITEM_BOTTOM_E].back();

				ModelItem item(ITEM_HOLE_E);
				std::size_t number = atol(vecSubStr[1].c_str());
				if (number < 1) {
					LOG(ERROR) << "invaild number in ITEM_HOLE_E";
					return false;
				}


				const auto& segment = botton.segments_[parent];

				std::size_t other_axis_index = 1;
				double other_axis = segment.second[1];
				double start_axis = segment.second[0];
				bool operate = segment.second[0] < segment.first[0] ? true : false;

				if (segment.first[0] == segment.second[0]) {
					other_axis_index = 0;
					other_axis = segment.second[0];
					start_axis = segment.second[1];
					operate = segment.second[1] < segment.first[1] ? true : false;

				}

				double minZ = 999999, maxZ = 0;
				Eigen::vector<Eigen::Vector2d> vecOriPts;
				for (std::size_t i = 2; i < number * 2 + 1; ) {
					Eigen::Vector3d point;

					if (operate)
						point[1 - other_axis_index] = start_axis + atol(vecSubStr[i].c_str());
					else
						point[1 - other_axis_index] = start_axis - atol(vecSubStr[i].c_str());


					point[other_axis_index] = other_axis;
					point[2] = atol(vecSubStr[i + 1].c_str());
					
					vecOriPts.emplace_back(Eigen::Vector2d(atol(vecSubStr[i].c_str()), atol(vecSubStr[i+1].c_str())));
					item.points_.emplace_back(point);

					item.points_.emplace_back(point);

					if (point[2] > maxZ) maxZ = point[2];
					if (point[2] < minZ) minZ = point[2];

					i += 2;


				}

				item.parentIndex_ = parent;
				item.highRange_ = std::make_pair(minZ, maxZ);
				item.area_ = calcArea(vecOriPts);
				

				item.buildSegment();
				vec_item.emplace_back(item);
				mapModelItem_[ITEM_HOLE_E].emplace_back(item);

			} 
			else if (vecSubStr[0] == "L") {
				ModelItem itemFirst(ITEM_BEAM_E);
				ModelItem itemSecond(ITEM_BEAM_E);

				if (mapModelItem_[ITEM_BOTTOM_E].empty() || mapModelItem_[ITEM_WALL_E].empty()) {
					LOG(ERROR) << "not found ITEM_BOTTOM_E";
					return false;
				}

				const std::size_t parent = atol(vecSubStr[1].c_str());
				const double thick = atol(vecSubStr[2].c_str());
				const ModelItem& botton = mapModelItem_[ITEM_BOTTOM_E].back();

				std::size_t number = atol(vecSubStr[3].c_str());
				if (number < 1) {
					LOG(ERROR) << "invaild number in ITEM_BEAM_E";
					return false;
				}

				const auto& segment = botton.segments_[parent];

				std::size_t other_axis_index = 1;
				double other_axis = segment.second[1];
				double start_axis = segment.second[0];
				bool operate = segment.second[0] < segment.first[0] ? true : false;

				if (segment.first[0] == segment.second[0]) {
					other_axis_index = 0;
					other_axis = segment.second[0];
					start_axis = segment.second[1];
					operate = segment.second[1] < segment.first[1] ? true : false;

				}

				const double other_axis_raw = other_axis;
				if (centerPt_[other_axis_index] > other_axis)
					other_axis += thick;
				else
					other_axis -= thick;



				double minZ = 999999, maxZ = 0;
				for (std::size_t i = 4; i < number * 2 + 3; ) {
					Eigen::Vector3d point;
					if (operate)
						point[1 - other_axis_index] = start_axis + atol(vecSubStr[i].c_str());
					else
						point[1 - other_axis_index] = start_axis - atol(vecSubStr[i].c_str());
					point[other_axis_index] = other_axis;
					point[2] = atol(vecSubStr[i + 1].c_str());
					itemFirst.points_.emplace_back(point);

					if (point[2] > maxZ) maxZ = point[2];
					if (point[2] < minZ) minZ = point[2];

					if (itemFirst.points_.size() == 1 
						|| itemFirst.points_.size() == 4)
					{
						itemFirst.points_.emplace_back(point);
						point[other_axis_index] = other_axis_raw;
						itemFirst.points_.emplace_back(point);

					}

					i += 2;

				}

				itemFirst.parentIndex_ = parent;
				itemFirst.highRange_ = std::make_pair(minZ, maxZ);
				itemFirst.buildSegment();
				vec_item.emplace_back(itemFirst);
				mapModelItem_[ITEM_BEAM_E].emplace_back(itemFirst);

				itemSecond.parentIndex_ = parent;
				itemSecond.buildSegment();
				vec_item.emplace_back(itemSecond);
				mapModelItem_[ITEM_BEAM_E].emplace_back(itemSecond);

				
			}

		}
	} else {
		LOG(ERROR) << "no such file" << fileName;
		return false;
	}

	reSortWall();

	//build top
	const std::vector<ModelItem>& vecWall = mapModelItem_[ITEM_WALL_E];
	const std::vector<ModelItem>& vecBeam = mapModelItem_[ITEM_BEAM_E];
	if (vecWall.empty()) return false;
	Eigen::vector<Eigen::Vector3d> vecPoints;
	for (auto& value : vecWall)
	{
		vecPoints.emplace_back(value.points_[0]);
		vecPoints.emplace_back(value.points_[3]);
	}

	
	for(std::size_t index = 0; index < vecBeam.size(); index+=2)
	{
		auto& beam = vecBeam[index];
		std::size_t parent = beam.parentIndex_;
		if (parent > vecWall.size() - 1)
		{
			LOG(ERROR) << "parent out of range";
			return false;
		}

		vecPoints[parent] = beam.points_[1];
		vecPoints[parent+1] = beam.points_[2];

		if (parent == 0)
		{
			vecPoints[vecPoints.size()-1] = beam.points_[1];
			vecPoints[parent + 2] = beam.points_[2];
		}
		else if (parent == vecWall.size()-1)
		{
			vecPoints[parent - 1] = beam.points_[1];
			vecPoints[0] = beam.points_[2];
		}
		else
		{
			vecPoints[parent-1] = beam.points_[1];
			vecPoints[parent+2] = beam.points_[2];

		}
		
	}

	ModelItem item(ITEM_TOP_E);
	for (std::size_t i = 0; i < vecPoints.size(); i += 2)
	{
		item.points_.emplace_back(vecPoints[i]);
	}
	item.buildSegment();
	vec_item.emplace_back(item);
	mapModelItem_[ITEM_TOP_E].emplace_back(item);


	// scale to meters
	//todo: in place
	scaleModel(0.001);

	savePCD("cad_model.pcd", vec_item);
	return true;

}

const std::vector<ModelItem>& CADModel::getTypedModelItems(ModelItemType type) const {
	return mapModelItem_.at(type);
}

std::vector<ModelItem> CADModel::tryGetTypedModelItems(ModelItemType type) const {
	auto iter = mapModelItem_.find(type);
	return iter == mapModelItem_.end() ? std::vector<ModelItem>() : iter->second;
}

std::vector<ModelItem> CADModel::filterModelItems(std::function<bool (const ModelItem&)> func) const {
	std::vector<ModelItem> models;
	for (const auto& pr : mapModelItem_)
		for (const auto& model : pr.second)
			if (func(model)) models.emplace_back(model);

	return models;
}

bool CADModel::savePCD(const std::string& name, std::vector<ModelItem>& vec_item) {
	Eigen::vector<Eigen::Vector3d> vecPoints;
	for (auto& item : vec_item) {
		for (auto& pt_pair : item.segments_) {
			auto vec_tmp = ininterpolateSeg(pt_pair.first, pt_pair.second, .1);
			vecPoints.insert(vecPoints.end(), vec_tmp.begin(), vec_tmp.end());
		}
	}

	writePCDFile(name, vecPoints);
	return true;
}

PointCloud::Ptr CADModel::genTestCloud() const {
	PointCloud::Ptr cloud(new PointCloud());

	for (auto& pr : mapModelItem_) {
		for (auto& item : pr.second) {
			for (auto& pt_pair : item.segments_) {
				auto vec_tmp = ininterpolateSeg(pt_pair.first, pt_pair.second, 0.01);
				for (const auto& v : vec_tmp)
					cloud->points.emplace_back(v(0, 0), v(1, 0), v(2, 0));
			}
		}
	}
	cloud->height = cloud->points.size();
	cloud->width = 1;
	cloud->is_dense = false;

	return cloud;
}

std::string CADModel::toString() const {
	std::stringstream ss;
	ss << "CADModel{ ";
	for (const auto& pr : mapModelItem_)
		ss << CloudReg::toString(pr.first) << " : " << pr.second.size() << ", ";
	ss << "}";

	return ss.str();
}

void CADModel::reSortWall()
{
	std::size_t maxIndex = 0;
	double maxArea = 0;
	const auto& vecHole = mapModelItem_[ITEM_HOLE_E];
	
	std::vector<double> vecArea(mapModelItem_[ITEM_WALL_E].size(), 0);

	for (auto& hole : vecHole)
	{
		if (hole.parentIndex_ < vecArea.size())
		{
			vecArea[hole.parentIndex_] += hole.area_;
			if (vecArea[hole.parentIndex_] > maxArea)
			{
				maxArea = vecArea[hole.parentIndex_];
				maxIndex = hole.parentIndex_;
			}
		}
	}

	if (maxIndex > 0)
	{
		auto& wall = mapModelItem_[ITEM_WALL_E];
		vecItems_t newWall;
		newWall.insert(newWall.end(),wall.begin()+maxIndex, wall.end());
		newWall.insert(newWall.end(),wall.begin(), wall.begin()+maxIndex);
		wall.swap(newWall);
		

		auto& botton = mapModelItem_[ITEM_BOTTOM_E].front().points_;
		Eigen::vector<Eigen::Vector3d> points;
		points.insert(points.end(), botton.begin() + maxIndex, botton.end());
		points.insert(points.end(), botton.begin(), botton.begin() + maxIndex);
		botton.swap(points);
		mapModelItem_[ITEM_BOTTOM_E].front().buildSegment();

		auto& vecHole = mapModelItem_[ITEM_HOLE_E];
		for (auto& hole : vecHole)
		{
			if (hole.parentIndex_ < vecArea.size())
			{
				if (hole.parentIndex_ >= maxIndex)
					hole.parentIndex_ -= maxIndex;
				else
					hole.parentIndex_ += (wall.size() - maxIndex);

			}
		}

		auto& vecBeam = mapModelItem_[ITEM_BEAM_E];
		for (auto& beam : vecBeam)
		{
			if (beam.parentIndex_ < vecArea.size())
			{
				if (beam.parentIndex_ >= maxIndex)
					beam.parentIndex_ -= maxIndex;
				else
					beam.parentIndex_ += (wall.size() - maxIndex);

			}
		}

	}
	
}

void CADModel::scaleModel(const double scale)
{
	for (auto& pr : mapModelItem_) {
		for (auto& model : pr.second) {
			for (auto& v : model.points_) v = v * scale;
			for (auto& v : model.segments_) v = ModelItem::ItemPair_t(v.first * scale, v.second * scale);
			model.highRange_ = std::make_pair(model.highRange_.first * scale, model.highRange_.second * scale);
		}
	}
}

}