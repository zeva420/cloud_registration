#include "funHelper.h"

namespace CloudReg
{
	Eigen::vector<Eigen::Vector3d> ininterpolateSeg(const Eigen::Vector3d& sPoint, const Eigen::Vector3d& ePoint, const double step)
	{
		Eigen::vector<Eigen::Vector3d> value;


		const std::size_t number = (sPoint - ePoint).norm() / step;

		const double r_N = 1.0 / number;
		const Eigen::Vector3d step_AB = (ePoint - sPoint) * r_N;

		value.emplace_back(sPoint);
		for (std::size_t i = 1; i < number; i++)
		{
			value.emplace_back(Eigen::Vector3d(step_AB * i + sPoint));
		}
		value.emplace_back(ePoint);
		// LOG(INFO) << "ininterpolateSeg number:" << value.size();
		return value;
	}

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

	double calcArea(const Eigen::vector<Eigen::Vector2d>& vecPts)
	{
		if (vecPts.size() < 3) return 0.0;

		double area = 0.0;
		for (std::size_t i = 2; i < vecPts.size();i++)
		{
			Eigen::Vector2d OB = vecPts[0] - vecPts[i - 1];
			Eigen::Vector2d OC = vecPts[0] - vecPts[i];
			area += OB[0] * OC[1] - OC[0] * OB[1];
			
		}
		return fabs(0.5* area);
	}
}