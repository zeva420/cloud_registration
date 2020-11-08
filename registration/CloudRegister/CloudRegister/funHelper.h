#pragma once

#include "BaseType.h"


namespace CloudReg 
{
	Eigen::vector<Eigen::Vector3d> ininterpolateSeg(const Eigen::Vector3d& sPoint, 
		const Eigen::Vector3d& ePoint, const double step);
	std::vector<std::string> splitByCharacter(const std::string& strtem, const char a);
	bool writePCDFile(const std::string& name, Eigen::vector<Eigen::Vector3d>& vecCloud);
}


