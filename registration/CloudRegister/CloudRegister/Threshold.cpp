
#include "Threshold.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <string>
#include <exception>
#include <iostream>

#include "glog/logging.h"

namespace CloudReg
{
using namespace boost::property_tree;

threshold::threshold()
{
}

threshold::~threshold()
{
}

/**
 * @brief initialize threshold config. if not config, use default value
 *
 * @param config_path:
 * @return bool:
 */
bool threshold::init(const std::string &config_path)
{
    std::string ini_file_name = config_path;
    ini_file_name.append("/config.ini");
    boost::filesystem::path targetPath = boost::filesystem::path(ini_file_name);

    if (!targetPath.empty() )
    {
        if ( ! boost::filesystem::exists(targetPath) )
        {
            LOG(WARNING) << "threshold file is not exist: " << ini_file_name;
            return false;
        }
    }

    LOG(INFO) << "config file path:" << config_path;

    try
    {
        ptree pt;
        ini_parser::read_ini(ini_file_name.c_str(), pt);

        corner_dropLength_ = pt.get<double>("threshold.corner_dropLength", corner_dropLength_);
		corner_calcWidth_ = pt.get<double>("threshold.corner_calcWidth", corner_calcWidth_);
		corner_calcLength_ = pt.get<double>("threshold.corner_calcLength", corner_calcLength_);

		height_calcHalfPara_ = pt.get<double>("threshold.height_calcHalfPara", height_calcHalfPara_);
		bay_depth_calcParaZ_ = pt.get<double>("threshold.bay_depth_calcParaZ", bay_depth_calcParaZ_);
		bay_depth_calcHalfPara_ = pt.get<double>("threshold.bay_depth_calcHalfPara", bay_depth_calcHalfPara_);

		segment_cloud_size_ratio_ = pt.get<double>("threshold.segment_cloud_size_ratio", segment_cloud_size_ratio_);
    }
    catch (const ini_parser_error &e)
    {
        auto what_data = e.what();
        LOG(WARNING) << "threshold,ini parser err data:" << what_data ;
        return false;
    }
    catch (const ptree_error &e)
    {
        auto what_data = e.what();
        LOG(WARNING) << "threshold,err data:" << what_data ;
        return false;
    }

    return true;
}

void threshold::uninit()
{
    //do nothing
}

/**
 * @brief display the config value
 *
 */
void threshold::echo_config_options() const
{
    LOG(INFO) <<"****************************************************************" ;
    LOG(INFO) <<"*******************echo_config_options begin********************" ;
	LOG(INFO) << "corner_dropLength: " << corner_dropLength_ << " corner_calcWidth: " << corner_calcWidth_
	<< " corner_calcLength: " << corner_calcLength_;
    
	LOG(INFO) << "*******************echo_config_options end********************";
    LOG(INFO) << "**************************************************************";
}

/**
 * @brief initialize and echo
 *
 * @param config_path:
 * @return bool:
 */
bool threshold::initAndEcho(std::string config_path)
{
    if (! init(config_path))
    {
        LOG(WARNING) << "init server config failed! Use default para";
        echo_config_options();
        return false;
    }
    else
    {
        echo_config_options();
        return true;
    }
}

}
