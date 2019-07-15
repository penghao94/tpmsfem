#pragma once

#include <vector>
#include <string>
#include <Eigen/dense>
#include <TPMSGenerator.h>
namespace tpmsfem {
	class tpmsTetGenerator {
	public:
		tpmsTetGenerator(tpmsgen::TPMSGenerator<double, double> *tg) :tg_(tg){}
		~tpmsTetGenerator()	{
			points.resize(0, 0);
			facets.resize(0, 0);
			elements.resize(0, 0);
			volumes.resize(0);
		}

		bool tpmsToTetrahedron(const std::string& switches,double levelset, Eigen::MatrixXd &V, Eigen::MatrixXi& F);

		static void preComputation(tpmsTetGenerator * tets, Eigen::MatrixXd& D, std::vector<Eigen::MatrixXd>& Kelements);

		int pointSize() { return points.rows(); }
		int elementSize() { return elements.rows(); }
		int facetSize() { return facets.rows(); }

		Eigen::MatrixXd points; // #V by 3 vertex position list
		Eigen::MatrixXi elements; //#T by 4 list of tet face indices
		Eigen::MatrixXi facets;
		Eigen::VectorXd volumes; //#T by 1 list of volume for each tetrahedron
		double scale = 0.004;

	private:
		tpmsgen::TPMSGenerator<double,double> *tg_;
		
	};
}