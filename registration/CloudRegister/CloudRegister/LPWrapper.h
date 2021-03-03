#pragma once

#include "BaseType.h"

namespace CloudReg {

// min <c, x>
// s.t. Ax <= b 
// d1 <= x <= d2

class LPWrapper {
public:
	struct OptResult {
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

		Eigen::VectorXd x_;
		double obj_;
		int ret_;

		bool isOptimal() const { return ret_ == 0; }

		std::string to_string() const;
	};

	static OptResult Solve(const Eigen::VectorXd& c, const Eigen::MatrixXd& A, const Eigen::VectorXd& b,
		const Eigen::VectorXd& d1, const Eigen::VectorXd& d2);
};

}
