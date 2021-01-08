#include "CalcMeasureHelper.h"

#include "funHelper.h"
#include <pcl/common/common.h>
#include <pcl/filters/crop_box.h>
#include <pcl/filters/passthrough.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>

namespace CloudReg
{

	void writePCDFile(const std::string& name, const std::vector<seg_pair_t>& segA, const std::vector<seg_pair_t>& segB)
	{
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr pCloudRGB(new pcl::PointCloud<pcl::PointXYZRGB>());
		for(auto& seg : segA)
		{
			auto vecPts = ininterpolateSeg(seg.first,seg.second,0.01);
			for(auto& pt : vecPts)
			{
				pcl::PointXYZRGB p2;
				p2.x = pt[0];	
				p2.y = pt[1];	
				p2.z = pt[2];	
				p2.r = 0;
				p2.g = 255;
				p2.b = 0;

				pCloudRGB->push_back(p2);
			}
		}
		
		for(auto& seg : segB)
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
			// std::cout << "cor " << seg.first[0] << " " << seg.first[1] << " " << seg.first[2] << std::endl;
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

	void groupDirectionIndex(const Eigen::Vector3d& horizenSeg, const std::vector<seg_pair_t>& border, 
			std::vector<std::size_t>& vecVertical, std::vector<std::size_t>& vecHorizen)
	{
		for(std::size_t i = 0 ; i < border.size(); i++)
		{
			auto& seg = border[i];
			auto curSeg = seg.first - seg.second;
			double cos = horizenSeg.dot(curSeg)/(horizenSeg.norm() * curSeg.norm());
			if (fabs(cos) < 0.0001) 
				vecVertical.emplace_back(i);
			else 
				vecHorizen.emplace_back(i);
		}

	}

	bool isRootInSeg(const seg_pair_t& seg, const Eigen::Vector3d& p)
	{

		auto segAB = seg.first - seg.second;
		auto segAP = seg.first - p;
		double dotA = segAB.dot(segAP);

		auto segBA = seg.second - seg.first;
		auto segBP = seg.second - p;
		double dotB = segBA.dot(segBP);
		//LOG(INFO) << dotA << " " << dotB;
		if (dotA < 0.0 || dotB < 0.0)
			return false;
		
		return true;
	}
	
	std::vector<Eigen::Vector3d> createRulerBox(seg_pair_t ruler, int thicknessDir, double thickness, double width)
	{
		if ((ruler.first - ruler.second).norm() < 1e-8)
		{
			LOG(ERROR) << "ruler length is 0";
			std::vector<Eigen::Vector3d> rulerPoints;
			return rulerPoints;
		}

		//thickness Direction
		Eigen::Vector3d thickn(0,0,0);
		thickn[thicknessDir] = 1;
		//ruler direction
		Eigen::Vector3d rulerab = ruler.second - ruler.first;
		Eigen::Vector3d rulern = rulerab.normalized();
		//width direction
		Eigen::Vector3d widthn = thickn.cross(rulern);

		//first 4 points
		Eigen::Vector3d pta = ruler.first + widthn * (width / 2);
		Eigen::Vector3d ptb = ruler.first - widthn * (width / 2);

		Eigen::Vector3d pt1 = pta + thickn * (thickness / 2);
		Eigen::Vector3d pt2 = pta - thickn * (thickness / 2);
		Eigen::Vector3d pt3 = ptb + thickn * (thickness / 2);
		Eigen::Vector3d pt4 = ptb - thickn * (thickness / 2);

		//last 4 points
		Eigen::Vector3d pt5 = pt1 + rulerab;
		Eigen::Vector3d pt6 = pt2 + rulerab;
		Eigen::Vector3d pt7 = pt3 + rulerab;
		Eigen::Vector3d pt8 = pt4 + rulerab;

		std::vector<Eigen::Vector3d> rulerPoints = {pt1, pt2, pt3, pt4, pt5, pt6, pt7, pt8};
		return rulerPoints;
	}

	bool calIntersection(seg_pair_t line1, seg_pair_t line2, Eigen::Vector3d& intersec)
	{
		Eigen::Vector3d s1e1 = line1.second - line1.first;
		Eigen::Vector3d s2e2 = line2.second - line2.first;
		Eigen::Vector3d c12 = s1e1.cross(s2e2);
		if (c12.norm() <= 0 + 1e-6)
		{
			LOG(INFO) << "Parallel line";
			return false;
		}

		Eigen::Vector3d s1s2 = line2.first - line1.first;
		Eigen::Vector3d cross1 = s1s2.cross(s1e1);
		Eigen::Vector3d cross2 = -s2e2.cross(s1e1);
		double scale = cross2.dot(cross1)/cross2.squaredNorm();
		if (scale < 0 || scale > 1)
		{
			LOG(INFO) << "There is no intersec";
			return false;
		}

		intersec = line2.first + scale*s2e2;
		return true;
	}

	bool calRuler3d(const std::vector<seg_pair_t>& wallBorder, const std::vector<seg_pair_t>& holeBorder, 
					 const seg_pair_t& rotateLine, const Eigen::Vector3d& P0,
					 const float& theta, seg_pair_t& ruler)
	{
		if (wallBorder.empty())
		{
			LOG(ERROR) << "empty input for ruler";
			return false;
		}
		
		if (P0 != rotateLine.first && P0 != rotateLine.second)
		{
			LOG(ERROR) << "P0 is not on the rotateLine for ruler";
			return false;
		}
		
		if (theta == 0)
		{
			ruler = rotateLine;
			return true;
		}

		//step0: Calculate normal vector
		seg_pair_t firstSeg = wallBorder[0];
		seg_pair_t bottomSeg = wallBorder.back();
		Eigen::Vector3d firstV = firstSeg.first - firstSeg.second;
		Eigen::Vector3d bottomV = bottomSeg.first - bottomSeg.second;
		Eigen::Vector3d nV = bottomV.cross(firstV);
		nV = nV.normalized();

		//step1: create R
		float radian = theta * geo::PI/180;
		Eigen::AngleAxisd t_V(radian, nV);
		Eigen::Matrix3d R = t_V.matrix();

		//step3: rotate
		Eigen::Vector3d rotateLineV = rotateLine.second - rotateLine.first;
		Eigen::Vector3d Q = R * rotateLineV + P0;
		seg_pair_t rotateQ = std::make_pair(P0, Q);

		//step4: Calculate the intersection point
		Eigen::Vector3d new3d;
		double rulerLength = 1000000.0;
		std::vector<seg_pair_t> allBorder = wallBorder;
		allBorder.insert(allBorder.end(), holeBorder.begin(), holeBorder.end());
		for (std::size_t i = 0; i < allBorder.size(); i++)
		{
			seg_pair_t line = allBorder[i];
			if (line == rotateLine)
				continue;

			Eigen::Vector3d s1e1 = line.second - line.first;
			Eigen::Vector3d s2e2 = rotateQ.second - rotateQ.first;
			Eigen::Vector3d c12 = s1e1.cross(s2e2);
			if (c12.norm() <= 0 + 1e-6 && ((line.first - P0).norm() < 1e-3 || (line.second - P0).norm() < 1e-3))
			{
				LOG(INFO) << "find the Vertical line"; // only for theta is 90
				ruler = line;
				LOG(INFO) <<"find the ruler end point is " << line.first[0] << " " <<line.first[1] 
					<< " " << line.first[2];
				LOG(INFO) <<"find the ruler end point is " << line.second[0] << " " <<line.second[1] 
					<< " " << line.second[2];
				return true;
			}
			
			Eigen::Vector3d intersec;
			if (calIntersection(rotateQ,line , intersec))
			{
				double length = (intersec - P0).norm();
				if (length <= 0 + 1e-6)
					continue;
				if (length < rulerLength)
				{
					rulerLength = length;
					new3d = intersec;
				}
			}
		}

		if (std::fabs(1000000.0 - rulerLength)< 1e-8)
		{
			LOG(ERROR) << "can not find ruler, please check theta";
			return false;
		}
		LOG(INFO) <<"find the ruler end point is " << new3d[0] << " " <<new3d[1] << " " << new3d[2];
		ruler = std::make_pair(P0, new3d);
		return true;
	}

}
