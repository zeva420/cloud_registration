#pragma once

#include "BaseType.h"

namespace CloudReg
{
	class CloudRefine
	{
	public:
		CloudRefine() {};
		~CloudRefine() {};

		bool run(std::map<CloudItemType, vecItems_t>& ret);
	};
}

