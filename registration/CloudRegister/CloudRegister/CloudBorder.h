#pragma once

#include "BaseType.h"

namespace CloudReg 
{
class CloudBorder
{
public:
    CloudBorder()
	{
	}

    ~CloudBorder()
    {
    }

	std::map<size_t, std::map<size_t, std::vector<Eigen::Vector3d>>> calcWallNodes(
			const std::string &name, 
			const PointCloud::Ptr cloud,
			const Eigen::Vector4d &cloudPlane,
			const std::vector<std::vector<seg_pair_t>> &cadBorder,
			const Eigen::Vector4d &cadPlane,
			const std::vector<seg_pair_t> &cloudOuterSegs);

private:		
	void prepareData(
			const std::string &name, 
			const PointCloud::Ptr cloud,
			const Eigen::Vector4d &cloudPlane,
			const std::vector<std::vector<seg_pair_t>> &cadBorder,
			const Eigen::Vector4d &cadPlane,
			const std::vector<seg_pair_t> &cloudOuterSegs);

	void getOverlapOfOuterSegAndHoleSeg(
			std::vector<seg_pair_t> &candidateSegs);

	bool detectLinesInCloudWall(
				std::vector<Eigen::VectorXf> &detectLineCoeffs,
				std::vector<PointCloud::Ptr> &detectLinePoints);

	int judgeCloudInWhichSideOfLine(const std::string &fileName,
			const Eigen::VectorXf &line, const PointCloud::Ptr inliers);

	std::pair<bool, int> inSameSideOfLine(
			const Eigen::Vector4d &cloudPlane, const Eigen::VectorXf &line, 
			const PointCloud::Ptr inliers, const float rate_Th);

	bool tooCloseToLine(const Eigen::VectorXf &line, 
			const PointCloud::Ptr inliers, const float rate_Th);

	bool removeOuterSideLines(
		const std::vector<Eigen::VectorXf> &detectLineCoeffs,
		const std::vector<PointCloud::Ptr> &detectLinePoints,
		std::vector<Eigen::VectorXf> &leftLineCoeffs,
		std::vector<PointCloud::Ptr> &leftLinePoints);

	void addCandidateSegs(const std::vector<seg_pair_t> &candidateSegs,
				std::vector<Eigen::VectorXf> &lineCoeffs,
				std::vector<PointCloud::Ptr> &linePoints,
				std::set<int> &candidateSegIdxs);

	bool matchLineToCadHoleSeg(
		const std::vector<Eigen::VectorXf> &lineCoeffs,
		const std::vector<PointCloud::Ptr> &linePoints,
		const std::set<int> &candidateSegIdxs,
		std::map<size_t, std::map<size_t, std::vector<size_t>>> &mapHole2Lines);

	void getIntersectionNodeOfLines(
		const std::map<size_t, std::map<size_t, std::vector<size_t>>> &mapHole2Lines,
		const std::vector<Eigen::VectorXf> &lineCoeffs,
		const std::vector<PointCloud::Ptr> &linePoints,
		std::map<size_t, std::map<size_t, std::vector<Eigen::Vector3d>>> &mapHole2Nodes);

	double dist_to_seg(const Eigen::Vector3d &point, const seg_pair_t &seg);

	double aveDistToSeg(const PointCloud::Ptr inputCloud, const seg_pair_t &seg);

	int findNextSeg(const std::vector<seg_pair_t> &vecSegs,
						const std::vector<bool> &vecValid,
						const seg_pair_t &currtSeg,  seg_pair_t &nextSeg);

	bool anticlockwiseSortSegVec(const Eigen::Vector3d &planeNormal, 
								const std::vector<seg_pair_t> &vecSegs,
								std::vector<seg_pair_t> &sortSegs, 
								std::vector<int> &matchIdxs);

	PointCloud::Ptr getArea( 
		const seg_pair_t &seg, PointCloud::Ptr inputCloud, 
		const Eigen::Vector4d &cloudPlane, 
		const double calcLength, bool left);

	int judge_flag(const PointCloud::Ptr pCloudLeft, 
							const PointCloud::Ptr pCloudRight);

	PointCloud::Ptr pointsOverlapWithSeg(const PointCloud::Ptr inputCloud, 
							const seg_pair_t &seg);

	void savePcdFileOfLines(const std::string &file_name,
				const std::vector<PointCloud::Ptr> &linePoints);

	void savePcdFileOfBoundPoints(const std::string &file_name,
				const PointCloud::Ptr boundPoints);

	void savePcdFileOfBelongLines(const std::string &file_name,
		const std::vector<PointCloud::Ptr> &linePoints,
		const std::vector<size_t> &lineIdxs, const seg_pair_t &seg);

	void savePcdFileOfBothSideAreaOfLines(const std::string &file_name,
			const PointCloud::Ptr &pCloudLeft, const PointCloud::Ptr &pCloudRight,
			const seg_pair_t &seg);

private:
	std::string name_; 
	PointCloud::Ptr cloud_;
	Eigen::Vector4d cloudPlane_;
	std::vector<std::vector<seg_pair_t>> cadBorder_;
	Eigen::Vector4d cadPlane_;
	std::vector<seg_pair_t> cloudOuterSegs_;

	std::vector<seg_pair_t> cadOuterSegs_;
	std::vector<seg_pair_t> cadHoleSegs_;
	std::map<size_t, std::vector<std::pair<size_t, size_t>>> cadSegNeighbours_;

	PointCloud::Ptr cloudSampling_;
};
}


