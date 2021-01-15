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

void SimpleViewer::show() {
#ifdef VISUALIZATION_ENABLED
	while (!viewer_.wasStopped()) viewer_.spinOnce(33);
	viewer_.close();
#endif
}

}

