#pragma once
#pragma once

#include "BaseType.h"
#include <pcl/visualization/pcl_visualizer.h>

namespace CloudReg {

class SimpleViewer {
public:
	// 255
	void addCloud(PointCloud::Ptr cloud);
	void addCloud(PointCloud::Ptr cloud, double r, double g, double b);

	void addSegment(const Point& s, const Point& e);
	void addSegment(const Point& s, const Point& e, double r, double g, double b);

	void addPoint(const Point& p);
	void addPoint(const Point& p, double r, double g, double b);

	void addLine(const Point& s, const Point& e, double r, double g, double b);

	// would block the thread
	void show();

private:
	pcl::visualization::PCLVisualizer viewer_;
	std::size_t id_{ 0 }; // avoid name setting...

	std::string genName() { return std::to_string(id_++); }
};

}