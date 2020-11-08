#include "CADModel.h"

#include "funHelper.h"

namespace CloudReg
{
	CADModel::CADModel()
	{
	}

	CADModel::~CADModel()
	{
		mapModelItem_.clear();
	}



	bool CADModel::initCAD(const std::string& fileName)
	{
		mapModelItem_.clear();


		mapModelItem_[ITEM_BOTTOM_E] = std::move(std::vector<ModelItem>());
		mapModelItem_[ITEM_WALL_E] = std::move(std::vector<ModelItem>());
		mapModelItem_[ITEM_HOLE_E] = std::move(std::vector<ModelItem>());
		mapModelItem_[ITEM_BEAM_E] = std::move(std::vector<ModelItem>());

		std::vector<ModelItem> vec_item;
		std::ifstream in(fileName);
		if (in)
		{
			std::string line;
			while (getline(in, line))
			{
				std::vector<std::string> vecSubStr = splitByCharacter(line, ',');
				if (vecSubStr.empty())
				{
					LOG(ERROR) << "splite string error";
					return false;
				}

				if (vecSubStr[0] == "B")
				{
					ModelItem item(ITEM_BOTTOM_E);
					std::size_t number = atol(vecSubStr[1].c_str());
					if (number < 1) 
					{
						LOG(ERROR) << "invaild number in ITEM_BOTTOM_E";
						return false;
					}
						

					for (std::size_t i = 2; i < number * 2 + 1; )
					{
						Eigen::Vector3d point;
						point[0] = atol(vecSubStr[i].c_str());
						point[1] = atol(vecSubStr[i + 1].c_str());
						point[2] = 0;
						item.points_.emplace_back(point);
						i += 2;

					}
					item.buildSegment();

					vec_item.emplace_back(item);
					mapModelItem_[ITEM_BOTTOM_E].emplace_back(item);
					const ModelItem& bottom = mapModelItem_[ITEM_BOTTOM_E].front();
					

				}
				else if (vecSubStr[0] == "W")
				{
					if (mapModelItem_[ITEM_BOTTOM_E].empty())
					{
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
					

					item.buildSegment();
					vec_item.emplace_back(item);
					mapModelItem_[ITEM_WALL_E].emplace_back(item);

				}
				else if (vecSubStr[0] == "H")
				{
					if (mapModelItem_[ITEM_BOTTOM_E].empty() || mapModelItem_[ITEM_WALL_E].empty())
					{
						LOG(ERROR) << "not found ITEM_BOTTOM_E";
						return false;
					}

					const std::size_t parent = mapModelItem_[ITEM_WALL_E].size() - 1;
					const ModelItem& botton = mapModelItem_[ITEM_BOTTOM_E].back();
					
					ModelItem item(ITEM_HOLE_E);
					std::size_t number = atol(vecSubStr[1].c_str());
					if (number < 1)
					{
						LOG(ERROR) << "invaild number in ITEM_HOLE_E";
						return false;
					}

					double other_axis_index = 1;
					const auto& segment = botton.segments_[parent];

					double other_axis = segment.first[1];
					double start_axis = segment.first[0];
					bool operate = segment.first[0] < segment.second[0] ? true : false;

					if (segment.first[0] == segment.second[0])
					{
						other_axis_index = 0;
						other_axis = segment.first[0];
						start_axis = segment.first[1];
						operate = segment.first[1] < segment.second[1] ? true : false;

					}
					
					

					for (std::size_t i = 2; i < number * 2 + 1; )
					{
						Eigen::Vector3d point;

						if (operate)
							point[1 - other_axis_index] = start_axis + atol(vecSubStr[i].c_str());
						else
							point[1 - other_axis_index] = start_axis - atol(vecSubStr[i].c_str());

				
						point[other_axis_index] = other_axis;
						point[2] = atol(vecSubStr[i + 1].c_str());
						item.points_.emplace_back(point);
						i += 2;


					}

					item.parentIndex_ = parent;
					item.buildSegment();
					vec_item.emplace_back(item);
					mapModelItem_[ITEM_HOLE_E].emplace_back(item);
					
				}
				else if (vecSubStr[0] == "L")
				{
					ModelItem item(ITEM_BEAM_E);

					if (mapModelItem_[ITEM_BOTTOM_E].empty() || mapModelItem_[ITEM_WALL_E].empty())
					{
						LOG(ERROR) << "not found ITEM_BOTTOM_E";
						return false;
					}

					const std::size_t parent = atol(vecSubStr[1].c_str());
					const double thick = atol(vecSubStr[2].c_str());
					const ModelItem& botton = mapModelItem_[ITEM_BOTTOM_E].back();

					std::size_t number = atol(vecSubStr[3].c_str());
					if (number < 1)
					{
						LOG(ERROR) << "invaild number in ITEM_BEAM_E";
						return false;
					}

					double other_axis_index = 1;
					const auto& segment = botton.segments_[parent];

					double other_axis = segment.first[1];
					double start_axis = segment.first[0];
					bool operate = segment.first[0] < segment.second[0] ? true : false;

					if (segment.first[0] == segment.second[0])
					{
						other_axis_index = 0;
						other_axis = segment.first[0];
						start_axis = segment.first[1];
						operate = segment.first[1] < segment.second[1] ? true : false;

					}


					const auto& segment_other = ((parent + 1) < botton.segments_.size())
						? botton.segments_[parent + 1] : botton.segments_[0];

					Eigen::Vector3d AB = segment.first - segment.second;
					Eigen::Vector3d CB = segment_other.second - segment_other.first;
					
					

					for (std::size_t i = 4; i < number * 2 + 3; )
					{
						Eigen::Vector3d point;
						if (operate)
							point[1 - other_axis_index] = start_axis + atol(vecSubStr[i].c_str());
						else
							point[1 - other_axis_index] = start_axis - atol(vecSubStr[i].c_str());
						point[other_axis_index] = other_axis;
						point[2] = atol(vecSubStr[i + 1].c_str());
						item.points_.emplace_back(point);
						i += 2;

					}
					
					item.parentIndex_ = parent;
					item.buildSegment();
					vec_item.emplace_back(item);
					mapModelItem_[ITEM_BEAM_E].emplace_back(item);
				}





			}
		}
		else
		{
			LOG(INFO) << "no such file" << fileName;
			return false;
		}
		savePCD("cad_model.pcd", vec_item);
		return true;

	}

	

	bool CADModel::savePCD(const std::string& name, std::vector<ModelItem>& vec_item)
	{
		Eigen::vector<Eigen::Vector3d> vecPoints;
		for (auto& item : vec_item)
		{
			for (auto& pt_pair : item.segments_)
			{
				auto vec_tmp = ininterpolateSeg(pt_pair.first, pt_pair.second, 1.0);
				vecPoints.insert(vecPoints.end(), vec_tmp.begin(), vec_tmp.end());

			}
		}

		writePCDFile(name, vecPoints);
		return true;
	}
}