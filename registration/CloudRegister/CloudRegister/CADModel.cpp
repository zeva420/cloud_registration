#include "CADModel.h"

#include "funHelper.h"

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

	std::ifstream in(fileName);
	if (in)
	{
		std::string line;
		while (getline(in, line))
		{
			std::vector<std::string> vecSubStr = splitByCharacter(line,',');
			if (vecSubStr.empty()) continue;

			if (vecSubStr[0] == "B")
			{
				ModelItem item(ITEM_BOTTOM_E);
				std::size_t number = atol(vecSubStr[1].c_str());
				if (number < 1) continue;

				for (std::size_t i = 2; i < number * 2 + 1; )
				{
					Eigen::Vector3d point;
					point[0] = atol(vecSubStr[i].c_str());
					point[1] = atol(vecSubStr[i + 1].c_str());
					point[2] = 0;
					//std::cout << point << std::endl;
					item.points_.emplace_back(point);
					i += 2;

				}
				item.buildSegment();		
				
				mapModelItem_[ITEM_BOTTOM_E].emplace_back(item);
				const ModelItem& bottom = mapModelItem_[ITEM_BOTTOM_E].front();
				savePCD("test.pcd", item);

			}
			else if(vecSubStr[0] == "W")
			{
				ModelItem item(ITEM_WALL_E);
				std::size_t parent = atol(vecSubStr[1].c_str());
				std::size_t high = atol(vecSubStr[2].c_str());

				if (mapModelItem_[ITEM_BOTTOM_E].empty())
				{
					return false;
				}

				const ModelItem& bottom = mapModelItem_[ITEM_BOTTOM_E].front();
				
				
				Eigen::Vector3d ptA = bottom.segments_[parent].first;
				Eigen::Vector3d ptAA = ptA + Eigen::Vector3d(0,0,high);
				Eigen::Vector3d ptB = bottom.segments_[parent].second;
				Eigen::Vector3d ptBB = ptA + Eigen::Vector3d(0, 0, high);

				
				item.points_.emplace_back(ptAA);
				item.points_.emplace_back(ptA);
				item.points_.emplace_back(ptB);
				item.points_.emplace_back(ptBB);

				item.buildSegment();
				mapModelItem_[ITEM_WALL_E].emplace_back(item);
				
			}
			else if (vecSubStr[0] == "H")
			{
				ModelItem item(ITEM_HOLE_E);
				item.buildSegment();
				mapModelItem_[ITEM_HOLE_E].emplace_back(item);
			}
			else if (vecSubStr[0] == "L")
			{

			}

			
			
			
			
		}
	}
	else
	{
		LOG(INFO) << "no such file" << fileName;
		return false;
	}
	return true;

}

bool CADModel::savePCD(const std::string& name, ModelItem& item)
{
	Eigen::vector<Eigen::Vector3d> vecPoints;
	for (auto& pt_pair : item.segments_)
	{
		std::size_t index = 0;
		const Eigen::Vector3d& pt1 = pt_pair.first;
		const Eigen::Vector3d& pt2 = pt_pair.second;
		std::size_t diff = fabs(pt1[0] - pt2[0]);
		
		if (pt_pair.first[0] == pt_pair.second[0])
		{
			index = 1;
			diff = fabs(pt1[1] - pt2[1]);

		}
		bool opt = pt2[index] - pt1[index] > 0 ? true : false;

		vecPoints.emplace_back(pt1);
		Eigen::Vector3d tmp = pt1;
		for (std::size_t i = 0; i < diff; i++)
		{
			if(opt)
				tmp[index] += 1;
			else
				tmp[index] -= 1;

			vecPoints.emplace_back(tmp);
		}
		vecPoints.emplace_back(pt2);
	}
	writePCDFile(name, vecPoints);
	return true;
}
