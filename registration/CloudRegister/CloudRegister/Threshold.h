#pragma once
#include "SingletonTemplate.h"
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

namespace CloudReg
{

class threshold
{
public:
    threshold();
    ~threshold();

    bool init(const std::string &config_path);
    void uninit();
    bool initAndEcho(std::string config_path =
                         boost::filesystem::initial_path<boost::filesystem::path>().string());

	double get_corner_dropLength() const { return corner_dropLength_;}
	double get_corner_rulerWidth() const { return corner_calcWidth_; }
	double get_corner_rulerLength() const { return corner_calcLength_; }
	
	double get_height_calcHalfPara() const { return height_calcHalfPara_; }
	
	double get_bay_depth_calcParaZ() const { return bay_depth_calcParaZ_; }
	double get_bay_depth_calcHalfPara() const { return bay_depth_calcHalfPara_; }

	double get_segment_cloud_size_ratio() const { return segment_cloud_size_ratio_; }
	double get_segment_cloud_growth_angle() const { return segment_cloud_growth_angle_; }
	std::size_t get_segment_cloud_plane_number() const { return segment_cloud_plane_number_; }

private:
	void echo_config_options() const;

	double corner_dropLength_ = 0.01;
	double corner_calcWidth_ = 0.005;
	double corner_calcLength_ = 0.13;

	double height_calcHalfPara_ = 0.005;

	double bay_depth_calcParaZ_ = 0.5;
	double bay_depth_calcHalfPara_ = 0.005;

	double segment_cloud_size_ratio_ = 0.25;
	double segment_cloud_growth_angle_ = 10.0;
	std::size_t segment_cloud_plane_number_ = 500;
};

typedef CloudReg::rdSingleton<threshold> TheThreshold;

} 

