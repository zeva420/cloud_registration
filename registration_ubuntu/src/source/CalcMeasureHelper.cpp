#include "CalcHoleMeasure.h"

#include "funHelper.h"
#include <pcl/common/common.h>
#include <pcl/filters/crop_box.h>
#include <pcl/filters/passthrough.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>

namespace CloudReg
{

	void writePCDFile(const std::string& name, const PointCloud::Ptr pCloud, std::vector<seg_pair_t>& border)
	{

		pcl::PointCloud<pcl::PointXYZRGB>::Ptr pCloudRGB(new pcl::PointCloud<pcl::PointXYZRGB>());

		for(auto& pt : pCloud->points)
		{
			pcl::PointXYZRGB p2;
			p2.x = pt.x;	
			p2.y = pt.y;	
			p2.z = pt.z;	
			p2.r = 0;
			p2.g = 255;
			p2.b = 0;

			pCloudRGB->push_back(p2);
		}

		for(auto& seg : border)
		{
			auto vecPts = ininterpolateSeg(seg.first,seg.second,0.01);
			for(auto& pt : vecPts)
			{
				pcl::PointXYZRGB p2;
				p2.x = pt[0];	
				p2.y = pt[1];	
				p2.z = pt[2];	
				p2.r = 255;
				p2.g = 0;
				p2.b = 0;

				pCloudRGB->push_back(p2);
			}
		}
		pcl::io::savePCDFile(name, *pCloudRGB);

	}

	Eigen::Vector3d calcPerpendicular(const Eigen::Vector3d &pt,
			const Eigen::Vector3d &start, const Eigen::Vector3d &end)
	{
		Eigen::Vector3d retVal(0,0,0);
		Eigen::Vector3d segAB = start - end; 
		Eigen::Vector3d segPA = pt - start; 
		
		if(segAB.norm() < EPS_FLOAT_DOUBLE)
		{
			retVal = start;
			return retVal;
		}

		double u = segAB.dot(segPA)/segAB.squaredNorm();
		retVal = start + u * segAB;
		return retVal;
	}



	void groupDirection(const Eigen::Vector3d& horizenSeg, const std::vector<seg_pair_t>& border, 
			std::vector<seg_pair_t>& vecVertical, std::vector<seg_pair_t>& vecHorizen)
	{
		for(auto& seg : border)
		{
			auto curSeg = seg.first - seg.second;
			double cos = horizenSeg.dot(curSeg)/(horizenSeg.norm() * curSeg.norm());
			if (fabs(cos) < 0.0001) 
				vecVertical.emplace_back(seg);
			else 
				vecHorizen.emplace_back(seg);
		}

		std::sort(vecHorizen.begin(), vecHorizen.end(), [&](const seg_pair_t& left, const seg_pair_t& right){
				return left.first[2] < right.first[2];});
	}

	PointCloud::Ptr filerCloudByRange(pcl::PointCloud<pcl::PointXYZ>::Ptr pCloud,
			const pcl::PointXYZ& min, const pcl::PointXYZ& max)
	{
		pcl::CropBox<pcl::PointXYZ> clipper;
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_final(new pcl::PointCloud<pcl::PointXYZ>());
		
		clipper.setMin(Eigen::Vector4f(min.x,min.y,min.z,1.0));
		clipper.setMax(Eigen::Vector4f(max.x,max.y,max.z,1.0));
        clipper.setNegative(false);
        clipper.setInputCloud(pCloud);
        clipper.filter(*cloud_final);
		LOG(INFO) << "inPut: " << pCloud->points.size() << " outPut: " << cloud_final->points.size();
		return cloud_final;
	}


	Eigen::vector<Eigen::Vector3d> getNearestPt(const Eigen::vector<Eigen::Vector3d>& vecCalcPt, 
			const PointCloud::Ptr pCloud, const double maxDist)
	{
		Eigen::vector<Eigen::Vector3d> vecCalcRaw(vecCalcPt.size(),Eigen::Vector3d(0,0,0));
		std::vector<double> vecDist(vecCalcPt.size(), 9999.0);
		for(auto& pt : pCloud->points)
		{
			Eigen::Vector3d ePt(pt.x, pt.y, pt.z);
			for(std::size_t i = 0; i < vecCalcPt.size(); i++)
			{
				double tmp = (ePt - vecCalcPt[i]).squaredNorm();
				//LOG(INFO) << tmp;
				if (tmp < maxDist && tmp < vecDist[i])
				//if (tmp < vecDist[i])
				{
					vecDist[i] = tmp;
					vecCalcRaw[i] = ePt;
				}
			}
		}

		for(std::size_t i = 0; i < vecCalcPt.size(); i++)
		{
			if (vecCalcRaw[i].norm() < EPS_FLOAT_DOUBLE)
			{
				vecCalcRaw[i] = vecCalcPt[i];
				LOG(WARNING) << "can't getNearestPt, use border instead";
			}
		}

		return vecCalcRaw;
	}

}
