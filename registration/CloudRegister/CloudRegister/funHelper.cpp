#include "funHelper.h"

namespace CloudReg
{

	std::vector<std::string> splitByCharacter(const std::string& strtem, const char a)
	{
		std::vector<std::string> strvec;

		std::string::size_type pos1, pos2;
		pos2 = strtem.find(a);
		pos1 = 0;
		while (std::string::npos != pos2)
		{
			strvec.push_back(strtem.substr(pos1, pos2 - pos1));

			pos1 = pos2 + 1;
			pos2 = strtem.find(a, pos1);
		}
		strvec.push_back(strtem.substr(pos1));
		return strvec;
	}

	bool writePCDFile(const std::string& name, Eigen::vector<Eigen::Vector3d>& vecCloud)
	{
		if (vecCloud.empty()) return false;

		pcl::PointCloud<pcl::PointXYZ> cloud;
		cloud.width = vecCloud.size();
		cloud.height = 1;
		cloud.is_dense = false;
		cloud.points.resize(cloud.width * cloud.height);

		for (size_t i = 0; i < cloud.points.size(); ++i)
		{
			cloud.points[i].x = vecCloud[i][0];
			cloud.points[i].y = vecCloud[i][1];
			cloud.points[i].z = vecCloud[i][2];
		}

		pcl::io::savePCDFile(name, cloud);

		return true;
	}
}