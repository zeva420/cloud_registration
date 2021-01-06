#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <boost/filesystem.hpp>


class FileHelper
{
public:
    /*********************************common function**********************************/
    static int getFilenamesFromdir(const std::string &dir,
                                   std::vector<std::string> &filenames,
                                   const std::string &suffix = ".pcd")
    {
        namespace fs = boost::filesystem;

        fs::path path(dir);

        if (!fs::exists(dir))
        {
            return -1;
        }

        fs::directory_iterator end_iter;

        for (fs::directory_iterator iter(path); iter != end_iter; ++iter)
        {
            if (fs::is_regular_file(iter->status()))
            {
                std::string fileExt = fs::extension(iter->path());

                if (fileExt == suffix)
                {
                    filenames.push_back(iter->path().string());
                }
            }

        }

        return (int)filenames.size();
    };

    static bool assurePathExists(const std::string &strPath)
    {
        bool bRet = false;

        boost::filesystem::path targetPath = boost::filesystem::path(strPath);

        if (!targetPath.empty() )
        {
            if ( boost::filesystem::exists(targetPath) )
            {
                bRet = true;
            }
            else
            {
                bRet = boost::filesystem::create_directories(targetPath);
            }
        }
        else
        {
            // empty string represents current directory
            bRet = true;
        }

        return bRet;
    }

};

