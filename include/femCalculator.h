#pragma once

#include <vector>
#include <Eigen/dense>
#include <Eigen/sparse>

#include "tpmsTetGenerator.h"
namespace tpmsfem {

	class femCalculator {
	public:
	femCalculator(){}
	~femCalculator() {
		unConstrainedDisplacements.resize(0);
		unConstrainedForces.resize(0);
		dofIds.resize(0);
		dofRules.resize(0);
		K11.resize(0, 0);
		K12.resize(0, 0);
		K21.resize(0, 0);
		K22.resize(0, 0);
		hatf.resize(0);
		barf.resize(0);
	}

	void init(tpmsTetGenerator *tets, const double E, const double mu, const std::vector<int> &pointrules, const Eigen::MatrixXd displacements, const Eigen::MatrixXd forces);
	void compute(Eigen::MatrixXd& points);


	public:
		tpmsTetGenerator *tetMesh;

		int unConstrainedNum = 0;

		Eigen::VectorXd unConstrainedDisplacements;
		Eigen::VectorXd unConstrainedForces;
		Eigen::VectorXi dofRules;
		Eigen::VectorXi dofIds;
		Eigen::VectorXd hatf, barf, baru,hatu;
		Eigen::SparseMatrix<double> K11, K12, K21, K22;
	};

}