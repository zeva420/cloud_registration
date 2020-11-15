#ifndef EDGE_CONSTRAINT_H
#define EDGE_CONSTRAINT_H


#include "g2o/core/base_unary_edge.h"
#include "g2o/core/base_binary_edge.h"
#include "g2o/core/base_multi_edge.h"
#include "g2o/types/sba/types_sba.h"
#include "g2o/types/sba/types_six_dof_expmap.h"
#include "g2o/types/sim3/types_seven_dof_expmap.h"

namespace g2o
{
  
class EdgePtToPlaneDist : public g2o::BaseUnaryEdge<1, double, g2o::VertexSE3Expmap>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	EdgePtToPlaneDist(Eigen::Vector3d pt, Eigen::Vector4d plane)
			: g2o::BaseUnaryEdge<1, double, g2o::VertexSE3Expmap>(),
				pt_(pt), 
                plane_(plane)
	{

	}

    double calcPt2PlaneDist(const Eigen::Vector3d &point, const Eigen::Vector4d &plane)
    {
        //plane: Ax+By+Cz+D=0,  point to plane dist = ||Ax+By+Cz+D||/(sqrt(A*A+B*B+C*C))
        Eigen::Vector3d n  = plane.block<3,1>(0,0);
        double d = plane(3);
        double dis = std::abs(n.dot(point) + d) / n.norm();
        return dis;
    }

    void computeError()
    {
		const g2o::VertexSE3Expmap* v = static_cast<const g2o::VertexSE3Expmap*>(_vertices[0]);
		g2o::SE3Quat T = v->estimate();

        Eigen::Vector3d newPt = T.rotation().toRotationMatrix() * pt_ + T.translation();
        double d = calcPt2PlaneDist(newPt, plane_);
		_error[0] = d;
    }

    virtual bool read(std::istream& is);
    virtual bool write(std::ostream& os) const;

private:
	Eigen::Vector3d pt_;
    Eigen::Vector4d plane_;
};

} //namespace g2o

#endif // EDGE_CONSTRAINT_H
