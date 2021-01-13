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

private:
	void echo_config_options() const;
};

typedef CloudReg::rdSingleton<threshold> TheThreshold;

} 

