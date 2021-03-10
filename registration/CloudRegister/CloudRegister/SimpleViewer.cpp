#include "SimpleViewer.h"
#include "GeometryUtils.h"


namespace CloudReg {

void SimpleViewer::addCloud(PointCloud::Ptr cloud, double r, double g, double b) {
#ifdef VISUALIZATION_ENABLED
	pcl::visualization::PointCloudColorHandlerCustom<Point> color(cloud, r, g, b);
	viewer_.addPointCloud(cloud, color, genName());
#endif
}
void SimpleViewer::addCloud(PointCloud::Ptr cloud) {
	double r{ geo::random() * 255. }, g{ geo::random() * 255. }, b{ geo::random() * 255. };
	addCloud(cloud, r, g, b);
}

void SimpleViewer::addSegment(const Point& s, const Point& e, double r, double g, double b) {
#ifdef VISUALIZATION_ENABLED
	viewer_.addSphere(s, 0.02f, r, g, b, genName());
	viewer_.addSphere(e, 0.02f, r, g, b, genName());
	viewer_.addLine(s, e, r, g, b, genName());
#endif
}
void SimpleViewer::addSegment(const Point& s, const Point& e) {
	double r{ geo::random() * 255. }, g{ geo::random() * 255. }, b{ geo::random() * 255. };
	addSegment(s, e, r, g, b);
}

void SimpleViewer::addPoint(const Point& p) {
	double r{ geo::random() * 255. }, g{ geo::random() * 255. }, b{ geo::random() * 255. };
	addPoint(p, r, g, b);
}
void SimpleViewer::addPoint(const Point& p, double r, double g, double b) {
#ifdef VISUALIZATION_ENABLED
	viewer_.addSphere(p, 0.05f, r, g, b, genName());
#endif
}

void SimpleViewer::addLine(const Point& s, const Point& e, double r, double g, double b) {
#ifdef VISUALIZATION_ENABLED
	viewer_.addLine(s, e, r, g, b, genName());
#endif
}

void SimpleViewer::addBox(const Point& min, const Point& max, double r, double g, double b) {
#ifdef VISUALIZATION_ENABLED
	double x1{ min.x }, y1{ min.y }, z1{ min.z }, x2{ max.x }, y2{ max.y }, z2{ max.z };
#define ADD_LINE(i1, i2, i3, j1, j2, j3) addLine(Point(x##i1, y##i2, z##i3), Point(x##j1, y##j2, z##j3), r, g, b)
	ADD_LINE(1, 1, 1, 2, 1, 1);
	ADD_LINE(2, 1, 1, 2, 2, 1);
	ADD_LINE(2, 2, 1, 1, 2, 1);
	ADD_LINE(1, 2, 1, 1, 1, 1);

	ADD_LINE(1, 1, 1, 1, 1, 2);
	ADD_LINE(2, 1, 1, 2, 1, 2);
	ADD_LINE(2, 2, 1, 2, 2, 2);
	ADD_LINE(1, 2, 1, 1, 2, 2);

	ADD_LINE(1, 1, 2, 2, 1, 2);
	ADD_LINE(2, 1, 2, 2, 2, 2);
	ADD_LINE(2, 2, 2, 1, 2, 2);
	ADD_LINE(1, 2, 2, 1, 1, 2);
#undef ADD_LINE
#endif
}

void SimpleViewer::show() {
#ifdef VISUALIZATION_ENABLED
	while (!viewer_.wasStopped()) viewer_.spinOnce(33);
	viewer_.close();
#endif
}

}

