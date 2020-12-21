#include "CADModel.h"

#include "funHelper.h"
#include "GeometryUtils.h"

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
		newPoints.emplace_back(target.points_[2]);
		newPoints.emplace_back(target.points_[3]);

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

std::map<ModelItemType, std::vector<PointCloud::Ptr>> CADModel::genTestFragCloud(double delta) const{
	

	auto sample_segment = [&](const Eigen::Vector3d& a, const Eigen::Vector3d& b, PointCloud::Ptr pCloud){
		for(const auto& v: ininterpolateSeg(a, b, delta)) pCloud->points.emplace_back(v(0), v(1), v(2));
	};

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
	
	for(std::size_t i=0; i<walls.size(); ++i){

		PointCloud::Ptr pCloud(new PointCloud());
		

		const auto& wall = walls[i];
		const auto& holes = wallHoles[i];

		LOG_IF(WARNING, wall.points_.size()!=4)<<"??? the wall is not a rectangle.";
		Eigen::Vector3d a = wall.points_[0];
		Eigen::Vector3d b = wall.points_[3];
		Eigen::Vector3d c = wall.points_[2];
		// 0	-- 3(c)
		// |	   |
		// 1(a) -- 2(b)

		for(float z= b(2); z< c(2); z+= delta){
			// test each hole
			Eigen::vector<Eigen::Vector3d> segends;
			for(const auto& hole: holes){
				if(z<hole.highRange_.first || z> hole.highRange_.second) continue;

				Eigen::vector<Eigen::Vector3d> ips;
				for (const auto& seg : hole.segments_) {
					auto sii = geo::zIntersectSegment(seg.first, seg.second, z);
					if (sii.valid()) { //todo: use margin check
						ips.emplace_back(sii.point_);
					}
				}

				if(ips.size()%2!=0)
					LOG(WARNING) << "intersection points is odd! may intersected with ends."; //todo: may just ignore this line.
				else segends.insert(segends.end(), ips.begin(), ips.end());
			}

			Eigen::Vector3d s(a(0), a(1), z), e(b(0), b(1), z);

			if(segends.empty()) sample_segment(s, e, pCloud);
			else {
				// need a simple sort.
				Eigen::Vector3d n = (e-s).normalized();
				segends.emplace_back(s);
				segends.emplace_back(e);
				
				// sort by proj len
				std::sort(segends.begin(), segends.end(), [&](const Eigen::Vector3d& p1, const Eigen::Vector3d& p2){
					return (p1-s).dot(n)< (p2-s).dot(n); 
				});

				for(std::size_t i=0; i< segends.size(); i+=2) sample_segment(segends[i], segends[i+1], pCloud);
			}
		}

		pCloud->height = pCloud->points.size();
		pCloud->width = 1;
		pCloud->is_dense = false;
		mapCADPoint[ITEM_WALL_E].emplace_back(pCloud);
	}

	// roof & floor
	auto y_intersect_segment = [](const Eigen::Vector3d& a, const Eigen::Vector3d& b, double y){
		geo::SegmentIntersectInfo sii;

		double dy = b(1)- a(1);
		if(std::fabs(dy)<1e-6) return sii;

		sii.lambda_ = (y-a(1))/ dy;
		if(sii.valid()) sii.point_ = a * (1 - sii.lambda_) + b * sii.lambda_;

		return sii;
	};

	auto sample_hor_shape = [&](const ModelItem& mi, PointCloud::Ptr pCloud){
		// aabb
		float z = mi.points_.front()(2);
		float x1 = ll::min_by([](const Eigen::Vector3d& v){ return v(0); }, mi.points_).second;
		float x2 = ll::max_by([](const Eigen::Vector3d& v) { return v(0); }, mi.points_).second;
		float y1 = ll::min_by([](const Eigen::Vector3d& v) { return v(1); }, mi.points_).second;
		float y2 = ll::max_by([](const Eigen::Vector3d& v) { return v(1); }, mi.points_).second;
		for(float y=y1; y<=y2; y+=delta){
			Eigen::vector<Eigen::Vector3d> ips;
			for (const auto& seg : mi.segments_) {
				auto sii = y_intersect_segment(seg.first, seg.second, y);
				if (sii.valid()) { //todo: use margin check
					ips.emplace_back(sii.point_);
				}
			}

			if (ips.size() % 2 != 0)
				LOG(WARNING) << "intersection points is odd! may intersected with ends.";
			else{
				// sort by x
				std::sort(ips.begin(), ips.end(), [&](const Eigen::Vector3d& p1, const Eigen::Vector3d& p2) { return p1(0)<p2(0); });

				for (std::size_t i = 0; i < ips.size(); i += 2) sample_segment(ips[i], ips[i + 1], pCloud);
			}
		}
	};

	for (const auto& mi : getTypedModelItems(ITEM_TOP_E))
	{
		PointCloud::Ptr pCloud(new PointCloud());
		sample_hor_shape(mi, pCloud);
		pCloud->height = pCloud->points.size();
		pCloud->width = 1;
		pCloud->is_dense = false;
		mapCADPoint[ITEM_TOP_E].emplace_back(pCloud);
	}
	for (const auto& mi : getTypedModelItems(ITEM_BOTTOM_E))
	{
		PointCloud::Ptr pCloud(new PointCloud());
		sample_hor_shape(mi, pCloud);
		pCloud->height = pCloud->points.size();
		pCloud->width = 1;
		pCloud->is_dense = false;
		mapCADPoint[ITEM_BOTTOM_E].emplace_back(pCloud);

	}

	

	return mapCADPoint;
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

		auto& botton = mapModelItem_[ITEM_BOTTOM_E].front().points_;
		Eigen::vector<Eigen::Vector3d> points{ botton.front()};
		points.insert(points.end(), botton.rbegin(), botton.rend()-1);
		botton.swap(points);
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
			for (auto& v : model.segments_) v = ModelItem::ItemPair_t(v.first * scale, v.second * scale);
			model.highRange_ = std::make_pair(model.highRange_.first * scale, model.highRange_.second * scale);
		}
	}
}

}