#include "CADModel.h"


#include "funHelper.h"
#include "GeometryUtils.h"

namespace CloudReg {

std::string ModelItem::toString() const {
	std::stringstream ss;
	ss << "ModelItem{ type: " << toModelItemName(itemtype_) << ", "
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
	std::ifstream in(fileName);
	if (!in)
	{
		LOG(ERROR) << "no such file" << fileName;
		return false;
	}
	

	std::string line;
	while (getline(in, line)) {
		std::vector<std::string> vecSubStr = splitByCharacter(line, ',');
		if (vecSubStr.empty()) {
			LOG(ERROR) << "splite string error";
			return false;
		}

		if (vecSubStr[0] == "B") {	
			auto ret = genBottom(vecSubStr);
			if (!std::get<0>(ret)) return false;			
		} 
		else if (vecSubStr[0] == "W") {
			auto ret = genWall(vecSubStr);
			if (!std::get<0>(ret)) return false;			
		} 
		else if (vecSubStr[0] == "H") {
			auto ret = genHole(vecSubStr);
			if (!std::get<0>(ret)) return false;			
		} 
		else if (vecSubStr[0] == "L") {
			auto ret = genBeam(vecSubStr);
			if (!std::get<0>(ret)) return false;
		}

	}
	cutWallByBeam();

	reSortWall();

	//build top
	auto ret = genTop();
	if (!std::get<0>(ret)) return false;
	
	// scale model to meters
	scaleModel(0.001);

#ifdef VISUALIZATION_ENABLED
	for (auto& item : mapModelItem_)
	{
		//if (item.first != ITEM_BOTTOM_E && item.first != ITEM_TOP_E) continue;
		//if (item.first != ITEM_BEAM_E && item.first != ITEM_WALL_E && item.first != ITEM_HOLE_E) continue;
		vec_item.insert(vec_item.end(), item.second.begin(), item.second.end());
	}
	savePCD("cad_model.pcd", vec_item);
#endif
	return true;

}
void CADModel::cutWallByBeam()
{
	//counter clockwise
	auto& vecBeam = mapModelItem_[ITEM_BEAM_E];
	const auto& vecWall = mapModelItem_[ITEM_WALL_E];
	std::vector<int> BeamIdx(vecBeam.size()/2,-1);
	for (std::size_t i = 0; i < vecBeam.size(); i+=2)
	{
		
		if (i + 2 < vecBeam.size() && 
			vecBeam[i+2].parentIndex_ - vecBeam[i].parentIndex_ == 1)
			BeamIdx[i/2] = i;
				

		if (i == vecBeam.size() - 2 && 
			vecBeam[i].parentIndex_ == vecWall.size()-1 && 
			vecBeam[0].parentIndex_ == 0)
		{
			BeamIdx[i/2] = 0;
		}
		
		
	}
	
	auto cutFunA = [](const ModelItem& target, const ModelItem& objA,const ModelItem& objB)
		->Eigen::vector<Eigen::Vector3d> {
		
		Eigen::vector<Eigen::Vector3d> newPoints;
		newPoints.emplace_back(target.points_.front());

		newPoints.emplace_back(objB.points_[objB.points_.size() - 1]);
		newPoints.emplace_back(objB.points_[objB.points_.size() - 2]);
		newPoints.emplace_back(objA.points_[objA.points_.size() - 2]);


		newPoints.insert(newPoints.end(), target.points_.begin()+2, target.points_.end());
		return newPoints;
		
	};


	auto cutFunB = [](const ModelItem& target, const ModelItem& objA, const ModelItem& objB)
		->Eigen::vector<Eigen::Vector3d> {

		Eigen::vector<Eigen::Vector3d> newPoints;
		newPoints.insert(newPoints.end(),target.points_.begin(), target.points_.end()-2);
		
		

		newPoints.emplace_back(objA.points_[1]);
		newPoints.emplace_back(objA.points_[0]);
		newPoints.emplace_back(objB.points_[0]);

		newPoints.emplace_back(target.points_.back());
		return newPoints;

		
	};

	auto cutFunC = [](const ModelItem& target, const ModelItem& obj)
		->Eigen::vector<Eigen::Vector3d> {

		Eigen::vector<Eigen::Vector3d> newPoints;

		newPoints.emplace_back(obj.points_[3]);
		newPoints.emplace_back(obj.points_[2]);
		newPoints.insert(newPoints.end(), target.points_.begin()+2, target.points_.end());		
		return newPoints;
	};

	for (std::size_t i = 1; i < vecBeam.size(); i += 2)
	{
		
		if (BeamIdx[i / 2] >= 0)
		{
			//std::cout << "enter:" << (i - 1) << "---" << (i + 1) / 2 << std::endl;
			auto& beamA1 = vecBeam[i - 1];
			auto& beamA2 = vecBeam[i];

			auto& beamB1 = vecBeam[i + 1];
			auto& beamB2 = vecBeam[i + 2];
			if (beamA1.highRange_.first < beamB1.highRange_.first)
			{
				auto newBeamPt = cutFunA(beamA1, beamB1, beamB2);
				beamA1.points_.swap(newBeamPt);
				beamA1.buildSegment();
			}
			else if (beamA1.highRange_.first > beamB1.highRange_.first)
			{
				auto newBeamPt = cutFunB(beamB1, beamA1, beamA2);
				beamB1.points_.swap(newBeamPt);
				beamB1.buildSegment();
			}
			else
			{
				auto newBeamPtA1 = cutFunC(beamA1, beamB1);
				auto newBeamPtA2 = cutFunC(beamA2, beamB2);

				beamA1.points_.swap(newBeamPtA1);
				beamA1.buildSegment();

				beamA2.points_.swap(newBeamPtA2);
				beamA2.buildSegment();
			}


		}
	}
	
	
	for (std::size_t i = 1; i < vecBeam.size(); i += 2)
	{
			auto& beamA1 = vecBeam[i - 1];
			auto& beamA2 = vecBeam[i];

			auto parent = vecBeam[i].parentIndex_;

			std::size_t preWallIndex = parent > 0 ?
				parent - 1 : mapModelItem_[ITEM_WALL_E].size() - 1;
			std::size_t nextWallIndex = parent != mapModelItem_[ITEM_WALL_E].size() - 1 ?
				parent + 1 : 0;


			//std::cout <<  i << " "<<parent << " nextWallIndex:" << nextWallIndex << " preWallIndex:" << preWallIndex << std::endl;

			{
				auto& wall = mapModelItem_[ITEM_WALL_E][nextWallIndex];
				if (wall.highRange_.second > beamA1.highRange_.first)
				{
					auto newWallPt = cutFunB(wall, beamA1, beamA2);
					wall.points_.swap(newWallPt);
					wall.buildSegment();
				}
				
			}
			
			{
				auto& wall = mapModelItem_[ITEM_WALL_E][preWallIndex];
				if (wall.highRange_.second > beamA1.highRange_.first)
				{
					auto newWallPt = cutFunA(wall, beamA1, beamA2);
					wall.points_.swap(newWallPt);
					wall.buildSegment();
				}
				
				
			}			

	}
		
	
}

std::tuple<bool, ModelItem> CADModel::genTop()
{
	const std::vector<ModelItem>& vecWall = mapModelItem_[ITEM_WALL_E];
	const std::vector<ModelItem>& vecBeam = mapModelItem_[ITEM_BEAM_E];
	if (vecWall.empty()) return std::tuple<bool, ModelItem>(false, ModelItem(ITEM_TOP_E));

	Eigen::vector<Eigen::Vector3d> vecPoints;
	for (auto& value : vecWall)
	{
		vecPoints.emplace_back(value.points_[1]);
		vecPoints.emplace_back(value.points_[2]);
	}

	//wall is clockwise
	for (std::size_t index = 0; index < vecBeam.size(); index += 2)
	{
		auto& beam = vecBeam[index];
		std::size_t parent = beam.parentIndex_;
		if (parent > vecWall.size() - 1)
		{
			LOG(ERROR) << "parent out of range";
			return std::tuple<bool, ModelItem>(false, ModelItem(ITEM_MAX_E));
		}

		auto getTopPt = [](const ModelItem& beam)->Eigen::vector<Eigen::Vector3d> {
			Eigen::vector<Eigen::Vector3d> vecPt;
			double topHigh = beam.highRange_.second;
			for (auto& pt : beam.points_)
			{
				if (std::fabs(pt[2] - topHigh) < 0.00001)
					vecPt.emplace_back(pt);
			}
			return vecPt;
		};

		auto vecPt = getTopPt(beam);
		//self modify
		vecPoints[parent*2] = vecPt[0];
		vecPoints[parent*2 + 1] = vecPt[1];

		//neighbor modify
		if (parent == 0)
		{
			vecPoints[vecPoints.size() - 1] = vecPt[0];
			vecPoints[parent*2 + 2] = vecPt[1];;
		}
		else if (parent == vecWall.size() - 1)
		{
			vecPoints[parent*2 - 1] = vecPt[0];;
			vecPoints[0] = vecPt[1];;
		}
		else
		{
			vecPoints[parent*2 - 1] = vecPt[0];;
			vecPoints[parent*2 + 2] = vecPt[1];;

		}

	}

	ModelItem item(ITEM_TOP_E);
	for (std::size_t i = 0; i < vecPoints.size(); i += 2)
		item.points_.emplace_back(vecPoints[i]);
	

	item.buildSegment();
	mapModelItem_[ITEM_TOP_E].emplace_back(item);

	return std::tuple<bool, ModelItem>(true,item);
}

std::tuple<bool, ModelItem, ModelItem> CADModel::genBeam(const std::vector<std::string>& vecSubStr)
{
	

	if (mapModelItem_[ITEM_BOTTOM_E].empty() || mapModelItem_[ITEM_WALL_E].empty()) {
		LOG(ERROR) << "not found ITEM_BOTTOM_E";
		return std::tuple<bool, ModelItem, ModelItem>(false, ModelItem(ITEM_MAX_E), ModelItem(ITEM_MAX_E));
	}

	ModelItem itemFirst(ITEM_BEAM_E);
	ModelItem itemSecond(ITEM_BEAM_E);


	const std::size_t parent = atol(vecSubStr[1].c_str());
	const double thick = atol(vecSubStr[2].c_str());
	const ModelItem& botton = mapModelItem_[ITEM_BOTTOM_E].back();

	std::size_t number = atol(vecSubStr[3].c_str());
	if (number < 1) {
		LOG(ERROR) << "invaild number in ITEM_BEAM_E";
		return std::tuple<bool, ModelItem, ModelItem>(false, ModelItem(ITEM_MAX_E), ModelItem(ITEM_MAX_E));
	}

	const auto& segment = botton.segments_[parent];

	std::size_t other_axis_index = 999;
	double other_axis, start_axis = 0.0;
	bool operate = false;

	getAxis(segment, other_axis_index, start_axis, other_axis, operate);

	const double other_axis_raw = other_axis;
	if (centerPt_[other_axis_index] > other_axis)
		other_axis += thick;
	else
		other_axis -= thick;


	/* first 1---2  second 1 --- 2
	*		 |   |		   |     |
	*		 0---3   (wall)0 --- 3  */
	
	double minZ = 999999, maxZ = 0;
	for (std::size_t i = 4; i < number * 2 + 3; i += 2) {
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

		if (itemFirst.points_.size() == 4)
		{
			itemSecond.points_.emplace_back(point);
			point[other_axis_index] = other_axis_raw;
			itemSecond.points_.emplace_back(point);

		}else if (itemFirst.points_.size() == 1)
		{
			auto tmpPt = point;
			point[other_axis_index] = other_axis_raw;
			itemSecond.points_.emplace_back(point);
			itemSecond.points_.emplace_back(tmpPt);

		}


	}

	//wall	
	auto& wall = mapModelItem_[ITEM_WALL_E][parent];
	wall.points_[1][2] = itemFirst.points_[0][2];
	wall.points_[2][2] = itemFirst.points_[0][2];
	wall.highRange_.second = minZ;
	wall.buildSegment();



	itemFirst.parentIndex_ = parent;
	itemFirst.highRange_ = std::make_pair(minZ, maxZ);
	itemFirst.buildSegment();
	mapModelItem_[ITEM_BEAM_E].emplace_back(itemFirst);

	itemSecond.parentIndex_ = parent;
	itemSecond.buildSegment();
	mapModelItem_[ITEM_BEAM_E].emplace_back(itemSecond);

	return std::tuple<bool, ModelItem,ModelItem>(true, itemFirst, itemSecond);
}

std::tuple<bool, ModelItem> CADModel::genHole(const std::vector<std::string>& vecSubStr)
{
	if (mapModelItem_[ITEM_BOTTOM_E].empty() || mapModelItem_[ITEM_WALL_E].empty()) {
		LOG(ERROR) << "not found ITEM_BOTTOM_E";
		return std::tuple<bool, ModelItem>(false, ModelItem(ITEM_MAX_E));
	}

	const std::size_t parent = mapModelItem_[ITEM_WALL_E].size() - 1;
	const ModelItem& botton = mapModelItem_[ITEM_BOTTOM_E].back();

	ModelItem item(ITEM_HOLE_E);
	std::size_t number = atol(vecSubStr[1].c_str());
	if (number < 1) {
		LOG(ERROR) << "invaild number in ITEM_HOLE_E";
		return std::tuple<bool, ModelItem>(false, ModelItem(ITEM_MAX_E));
	}


	const auto& segment = botton.segments_[parent];
	std::size_t other_axis_index = 999;
	double other_axis, start_axis = 0.0;
	bool operate = false;

	getAxis(segment, other_axis_index, start_axis, other_axis, operate);
	

	double minZ = 999999, maxZ = 0;
	Eigen::vector<Eigen::Vector2d> vecOriPts;
	for (std::size_t i = 2; i < number * 2 + 1; i += 2) {
		Eigen::Vector3d point;

		if (operate)
			point[1 - other_axis_index] = start_axis + atol(vecSubStr[i].c_str());
		else
			point[1 - other_axis_index] = start_axis - atol(vecSubStr[i].c_str());


		point[other_axis_index] = other_axis;
		point[2] = atol(vecSubStr[i + 1].c_str());

		vecOriPts.emplace_back(Eigen::Vector2d(atol(vecSubStr[i].c_str()), atol(vecSubStr[i + 1].c_str())));
		item.points_.emplace_back(point);

		if (point[2] > maxZ) maxZ = point[2];
		if (point[2] < minZ) minZ = point[2];

	}

	item.parentIndex_ = parent;
	item.highRange_ = std::make_pair(minZ, maxZ);
	item.area_ = calcArea(vecOriPts);
	


	item.buildSegment();
	mapModelItem_[ITEM_HOLE_E].emplace_back(item);

	return std::tuple<bool, ModelItem>(true, item);
}

std::tuple<bool, ModelItem> CADModel::genWall(const std::vector<std::string>& vecSubStr)
{	
	if (mapModelItem_[ITEM_BOTTOM_E].empty()) {
		LOG(ERROR) << "not found ITEM_BOTTOM_E";
		return std::tuple<bool, ModelItem>(false, ModelItem(ITEM_MAX_E));
	}
	
	const ModelItem& bottom = mapModelItem_[ITEM_BOTTOM_E].front();
	ModelItem item(ITEM_WALL_E);
	std::size_t parent = atol(vecSubStr[1].c_str());
	std::size_t high = atol(vecSubStr[2].c_str());



	Eigen::Vector3d ptA = bottom.segments_[parent].first;
	Eigen::Vector3d ptAA = ptA + Eigen::Vector3d(0, 0, high);
	Eigen::Vector3d ptB = bottom.segments_[parent].second;
	Eigen::Vector3d ptBB = ptB + Eigen::Vector3d(0, 0, high);

	/* 1---2
	*  |   |
	*  0---3 */

	item.points_.emplace_back(ptB);
	item.points_.emplace_back(ptBB);
	item.points_.emplace_back(ptAA);
	item.points_.emplace_back(ptA);
	
	item.highRange_ = std::make_pair(0, high);
	item.buildSegment();
	mapModelItem_[ITEM_WALL_E].emplace_back(item);

	return std::tuple<bool, ModelItem>(true, item);
}

std::tuple<bool, ModelItem> CADModel::genBottom(const std::vector<std::string>& vecSubStr)
{
	ModelItem item(ITEM_BOTTOM_E);
	std::size_t number = atol(vecSubStr[1].c_str());
	if (number < 1) {
		LOG(ERROR) << "invaild number in ITEM_BOTTOM_E";
		return std::tuple<bool, ModelItem>(false, ModelItem(ITEM_MAX_E));
	}

	double maxX = 0, maxY = 0, minX = 0, minY = 0;
	for (std::size_t i = 2; i < number * 2 + 1; i += 2) {
		Eigen::Vector3d point;

		point[0] = atol(vecSubStr[i].c_str());
		point[1] = atol(vecSubStr[i + 1].c_str());
		point[2] = 0;

		item.points_.emplace_back(point);

		if (point[0] > maxX) maxX = point[0];
		if (point[1] > maxY) maxY = point[1];
		if (point[0] < minX) minX = point[0];
		if (point[1] < minY) minY = point[1];

	}

	centerPt_[0] = (maxX + minX) / 2;
	centerPt_[1] = (maxY + minY) / 2;


	item.buildSegment();
	mapModelItem_[ITEM_BOTTOM_E].emplace_back(item);

	return std::tuple<bool, ModelItem>(true, item);

}

void  CADModel::getAxis(const std::pair<Eigen::Vector3d, Eigen::Vector3d>& segment,
	std::size_t& other_axis_index, double& start_axis, double& other_axis, bool& operate)
{
	other_axis_index = 1;
	other_axis = segment.second[1];
	start_axis = segment.second[0];
	operate = segment.second[0] < segment.first[0] ? true : false;

	if (segment.first[0] == segment.second[0]) {
		other_axis_index = 0;
		other_axis = segment.second[0];
		start_axis = segment.second[1];
		operate = segment.second[1] < segment.first[1] ? true : false;

	}
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
			auto vec_tmp = ininterpolateSeg(pt_pair.first, pt_pair.second, .01);
			vecPoints.insert(vecPoints.end(), vec_tmp.begin(), vec_tmp.end());
		}
	}

	writePCDFile(name, vecPoints);
	return true;
}

PointCloud::Ptr CADModel::genTestFrameCloud() const {
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

std::map<ModelItemType, std::vector<PointCloud::Ptr>> CADModel::genFragCloud(double delta) const{
	std::map<ModelItemType, std::vector<PointCloud::Ptr>> mapCADPoint;

	// walls
	const auto& walls = getTypedModelItems(ITEM_WALL_E);
	std::vector<std::vector<ModelItem>> wallHoles(walls.size());
	if(containModels(ITEM_HOLE_E)) {
		const auto& holes = getTypedModelItems(ITEM_HOLE_E);
		for(const auto& hole: holes) {
			if(hole.parentIndex_>=0 && hole.parentIndex_<walls.size()) wallHoles[hole.parentIndex_].push_back(hole);
			else LOG(WARNING)<<"invalid hole parent index: "<< hole.parentIndex_;
		}
	}

	for(auto pr: ll::zip(walls.begin(), walls.end(), wallHoles.begin(), wallHoles.end())){
		const auto& wall = *pr.first;
		const auto& holes = *pr.second;

		auto cloud = InterpolateShape(wall.points_, ll::mapf([](const ModelItem& mi) { return mi.points_; }, holes), delta);
		mapCADPoint[ITEM_WALL_E].emplace_back(cloud);
	}

	// floor, roof & beams
	std::vector<Eigen::vector<Eigen::Vector3d>> noholes;
	for (const auto& mi : getTypedModelItems(ITEM_TOP_E))
		mapCADPoint[ITEM_TOP_E].emplace_back(InterpolateShape(mi.points_, noholes, delta));

	for (const auto& mi : getTypedModelItems(ITEM_BOTTOM_E))
		mapCADPoint[ITEM_BOTTOM_E].emplace_back(InterpolateShape(mi.points_, noholes, delta));

	for(const auto& mi: getTypedModelItems(ITEM_BEAM_E))
		mapCADPoint[ITEM_BEAM_E].emplace_back(InterpolateShape(mi.points_, noholes, delta));

	return mapCADPoint;
}

std::string CADModel::toString() const {
	std::stringstream ss;
	ss << "CADModel{ ";
	for (const auto& pr : mapModelItem_)
		ss << CloudReg::toModelItemName(pr.first) << " : " << pr.second.size() << ", ";
	ss << "}";

	return ss.str();
}

void CADModel::reSortWall()
{
	

	std::size_t maxIndex = 0;
	double maxArea = 0;
	const auto& vecHole = mapModelItem_[ITEM_HOLE_E];
	

	for (auto& hole : vecHole)
	{
		if (hole.parentIndex_ < mapModelItem_[ITEM_WALL_E].size())
		{
			if (hole.area_ > maxArea)
			{
				maxArea = hole.area_;
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

		auto function = [&](ModelItem& value) {
			if (value.parentIndex_ < mapModelItem_[ITEM_WALL_E].size()) {
				if (value.parentIndex_ >= maxIndex) value.parentIndex_ -= maxIndex;
				else value.parentIndex_ += (wall.size() - maxIndex);
			}
		};

		auto& vecHole = mapModelItem_[ITEM_HOLE_E];
		std::for_each(vecHole.begin(), vecHole.end(),function);

		auto& vecBeam = mapModelItem_[ITEM_BEAM_E];
		std::for_each(vecBeam.begin(), vecBeam.end(), function);

	}

	//resort by clockwise
	{
		auto& wall = mapModelItem_[ITEM_WALL_E];
		vecItems_t newWall{wall.front()};
		newWall.insert(newWall.end(), wall.rbegin(), wall.rend()-1);
		wall.swap(newWall);

		auto& bottom = mapModelItem_[ITEM_BOTTOM_E].front().points_;

		bottom.clear();
		for(auto& value : wall)
			bottom.emplace_back(value.points_[0]);
		mapModelItem_[ITEM_BOTTOM_E].front().buildSegment();

		auto function = [&](ModelItem& value) {
			if (value.parentIndex_ < mapModelItem_[ITEM_WALL_E].size()
				&& value.parentIndex_ > 0) {
					value.parentIndex_ = mapModelItem_[ITEM_WALL_E].size() - value.parentIndex_;
				
			}
		};

		auto& vecHole = mapModelItem_[ITEM_HOLE_E];
		std::for_each(vecHole.begin(), vecHole.end(), function);

		auto& vecBeam = mapModelItem_[ITEM_BEAM_E];
		std::for_each(vecBeam.begin(), vecBeam.end(), function);

	}
	
}

void CADModel::scaleModel(const double scale)
{
	for (auto& pr : mapModelItem_) {
		for (auto& model : pr.second) {
			for (auto& v : model.points_) v = v * scale;
			for (auto& v : model.segments_) v = seg_pair_t(v.first * scale, v.second * scale);
			model.highRange_ = std::make_pair(model.highRange_.first * scale, model.highRange_.second * scale);
		}
	}
}

PointCloud::Ptr CADModel::InterpolateShape(const Eigen::vector<Eigen::Vector3d>& points,
	const std::vector<Eigen::vector<Eigen::Vector3d> >& holes, double delta) const{
	PointCloud::Ptr cloud(new PointCloud());

	Eigen::Vector3d minp = Eigen::Vector3d::Ones() * 100.;
	Eigen::Vector3d maxp = Eigen::Vector3d::Ones() * -100.;
	for (const Eigen::Vector3d& p : points) {
		for (int i = 0; i < 3; ++i) {
			if (p(i) < minp(i)) minp(i) = p(i);
			if (p(i) > maxp(i)) maxp(i) = p(i);
		}
	}

	// this is not very readable, may refactor someday. or, just write more duplicate codes...
	int dim1{0}, dim2{1}, dimfix{2};
	Eigen::Vector3d dp = maxp- minp;
	if(dp(0)<1e-6) std::swap(dimfix, dim1);
	else if(dp(1)<1e-6) std::swap(dimfix, dim2);

	double fixvalue = points.front()[dimfix];

	auto intersect_contour = [](const Eigen::Vector3d& s, const Eigen::Vector3d& e, 
		const Eigen::vector<Eigen::Vector3d>& points){
		Eigen::Vector3d se = e - s;

		Eigen::vector<Eigen::Vector3d> intersections;
		for (std::size_t i = 0; i < points.size(); ++i) {
			Eigen::Vector3d s2 = points[i];
			Eigen::Vector3d e2 = points[(i + 1) % points.size()];
			Eigen::Vector3d s2e2 = e2 - s2;

			Eigen::Vector3d c1 = (s2 - s).cross(se);
			Eigen::Vector3d c2 = s2e2.cross(se);
			double t = -c2.transpose().dot(c1) / c2.squaredNorm();

			if (t > 0. && t < 1.) intersections.emplace_back(s2 + s2e2 * t);
		}

		return intersections;
	};

	for(double dim1s = minp[dim1]; dim1s< maxp[dim1]; dim1s+=delta){
		Eigen::Vector3d s, e;
		s[dim1] = dim1s;
		s[dim2] = minp[dim2];
		s[dimfix] = fixvalue;
		e[dim1] = dim1s;
		e[dim2] = maxp[dim2];
		e[dimfix] = fixvalue;

		Eigen::vector<Eigen::Vector3d> intersections = intersect_contour(s, e, points);
		for(const auto& outline: holes){
			Eigen::vector<Eigen::Vector3d> ins =  intersect_contour(s, e, outline);
			intersections.insert(intersections.end(), ins.begin(), ins.end());
		}
		
		if(!intersections.empty() && intersections.size()%2!=0)
			LOG(INFO)<<"intersections issue, failed to interpolate a segment, size: "<< intersections.size();

		Eigen::Vector3d se = e - s;
		std::sort(intersections.begin(), intersections.end(), [&s, &se](const Eigen::Vector3d& v1, const Eigen::Vector3d& v2){
			return (v1-s).dot(se) < (v2-s).dot(se);
		});

		for(std::size_t i=0; i<intersections.size(); i+=2){
			for(const Eigen::Vector3d& v: ininterpolateSeg(intersections[i], intersections[i+1], delta))
				cloud->points.emplace_back(v(0), v(1), v(2));
		}
	}

	cloud->width = 1;
	cloud->height = cloud->points.size();
	cloud->is_dense = false;
	return cloud;
}

}
