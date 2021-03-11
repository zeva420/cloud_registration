#pragma once
#pragma once

#include "BaseType.h"

#ifdef VISUALIZATION_ENABLED
#include <pcl/visualization/pcl_visualizer.h>
#endif

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

	void addBox(const Point& min, const Point& max, double r, double g, double b);

	void addZPlane(float z, float halfsize = 2.f);

#ifdef VISUALIZATION_ENABLED
	pcl::visualization::PCLVisualizer& viewer() { return viewer_; }
#endif
	// would block the thread
	void show();

private:
#ifdef VISUALIZATION_ENABLED
	pcl::visualization::PCLVisualizer viewer_;
#endif
	std::size_t id_{ 0 }; // avoid name setting...

	std::string genName() { return std::to_string(id_++); }
};

}

