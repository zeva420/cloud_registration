#include "LPWrapper.h"
#include "ll.hpp"
#include "glog/logging.h"
#include "lp_lib.h"

namespace CloudReg {

std::string LPWrapper::OptResult::to_string() const {
	if (!isOptimal_) return "{\nfailed\n}";
	std::stringstream ss;
	ss << "{\noptimal with objective: " << obj_ << "\nx: " << x_.transpose() << "\n}";
	return ss.str();
}

LPWrapper::OptResult LPWrapper::Solve(const Eigen::VectorXd& c, const Eigen::MatrixXd& A, const Eigen::VectorXd& b,
	const Eigen::VectorXd& d1, const Eigen::VectorXd& d2) {
	const int DIM = c.rows();

	LL_ASSERT(A.cols() == DIM && A.rows() == b.rows() && d1.rows() == DIM && d2.rows() == DIM &&
		"wrong dimension.");

	lprec* lp = make_lp(0, DIM);

	int* colno = new int[DIM];
	std::iota(colno, colno + DIM, 1);
	double* row = new double[DIM];

	// le
	set_add_rowmode(lp, TRUE);
	for (int r = 0; r < A.rows(); ++r) {
		for (int ci = 0; ci < A.cols(); ++ci)
			row[ci] = A(r, ci);
		add_constraintex(lp, DIM, row, colno, LE, b(r, 0));
	}
	set_add_rowmode(lp, FALSE);

	// bd
	for (int ci = 0; ci < DIM; ++ci) set_bounds(lp, ci + 1, d1(ci, 0), d2(ci, 0));

	// obj
	for (int ci = 0; ci < DIM; ++ci) row[ci] = c(ci, 0);
	set_obj_fnex(lp, DIM, row, colno);
	set_minim(lp);

	OptResult result;
	result.isOptimal_ = solve(lp) == OPTIMAL;
	if (result.isOptimal_) {
		result.x_.resize(DIM);
		get_variables(lp, result.x_.data());
		result.obj_ = get_objective(lp);
	}

	delete[] colno;
	delete[] row;
	delete_lp(lp);

	return result;
}

}
