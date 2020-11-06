#pragma once

#include "BaseType.h"


namespace CloudReg 
{

	std::vector<std::string> splitByCharacter(const std::string& strtem, const char a);
	bool writePCDFile(const std::string& name, Eigen::vector<Eigen::Vector3d>& vecCloud);
}


