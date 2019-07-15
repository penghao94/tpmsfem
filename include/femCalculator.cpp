#include "femCalculator.h"
#include<iostream>
void tpmsfem::femCalculator::init(tpmsTetGenerator * tets, const double E, const double mu, const std::vector<int>& pointrules, const Eigen::MatrixXd displacements, const Eigen::MatrixXd forces)
{
	tetMesh = tets;

	int nv = tetMesh->pointSize();
	int ne = tetMesh->elementSize();

	dofIds.setZero(nv * 3);

	//point rules set -1 if points are fixed points; 0 if the force is 0 or gravity and 1 if the force is loaded
	std::vector<int> fixIndex, freeIndex, handleIndex;
	for (int i = 0; i < pointrules.size(); i++) {
		switch (pointrules[i])
		{
		default:
			break;
		case -1:
			fixIndex.push_back(i);
			break;
		case 0:
			freeIndex.push_back(i);
			break;
		case 1:
			handleIndex.push_back(i);
			break;
		}
	}

	// Initialize the user-defined force and displacement at non-fixed vertex
	// For fixed vertex, we can simply removed it to make the optimization problem unconstrainted
	//Set global index while ignoring the fix vertexes
	int index = 0;
	for (int i = 0; i < pointrules.size(); i++)
		dofIds(i) = pointrules[i] >= 0 ? index++ : -1;

	unConstrainedNum = index;
	dofRules.resize(unConstrainedNum);
	unConstrainedDisplacements.resize(unConstrainedNum);
	unConstrainedForces.resize(unConstrainedNum);

	//for free points
	int freeDofNum = freeIndex.size();
	hatf.resize(freeDofNum);
	for (int i = 0; i < freeIndex.size(); i++) {
		int uid = dofIds(freeIndex[i]);
		dofRules(uid) = i;
		unConstrainedForces(uid) = forces(freeIndex[i] / 3, freeIndex[i] % 3);
		unConstrainedDisplacements(uid) = displacements(freeIndex[i] / 3, freeIndex[i] % 3);
		hatf(i) = forces(freeIndex[i] / 3, freeIndex[i] % 3);
	}

	//for handle points
	int handleDofNum = handleIndex.size();
	barf.resize(handleDofNum);

	for (int i = 0; i < handleIndex.size(); i++) {
		int uid = dofIds(handleIndex[i]);
		dofRules(uid) = unConstrainedNum + i;
		unConstrainedForces(uid) = forces(handleIndex[i] / 3, handleIndex[i] % 3);
		unConstrainedDisplacements(uid) = displacements(handleIndex[i] / 3, handleIndex[i] % 3);
		barf(i) = forces(handleIndex[i] / 3, handleIndex[i] % 3);
	}


	//Ke values
	const double Lambda = mu / (1 + mu) / (1 - 2 * mu);// lame coefficient
	const double Mu = 1.0 / (2.0*(1.0 + mu));//lame cofficient

	// Linear tensor of elasticity D
	const auto computeLinearTensor = [](double lambda, double mu)->Eigen::MatrixXd {
		Eigen::MatrixXd tmp(6, 6);
		tmp << lambda + 2 * mu, lambda, lambda, 0, 0, 0,
			lambda, lambda + 2 * mu, lambda, 0, 0, 0,
			lambda, lambda, lambda + 2 * mu, 0, 0, 0,
			0, 0, 0, mu, 0, 0,
			0, 0, 0, 0, mu, 0,
			0, 0, 0, 0, 0, mu;
		return tmp;
	};
	Eigen::MatrixXd De = computeLinearTensor(Lambda,Mu);

	std::vector<Eigen::MatrixXd> Kelements;
	tpmsTetGenerator::preComputation(tetMesh,De,Kelements);

	//   freeDofNum handleDofNum 
	//   [ K11(E)   K12(E) ] freeDofNum
	//K= [				   ] 
	//   [ K21(E)   K22(E) ] handleDofNum

	K11.resize(freeDofNum, freeDofNum);
	K12.resize(freeDofNum, handleDofNum);
	K21.resize(handleDofNum, freeDofNum);
	K22.resize(handleDofNum, handleDofNum);


	const auto allocate=[&](const int num, int& k, int & index)->void{
		k = num / unConstrainedNum;
		index = num % unConstrainedNum;
	};

	int r = 0, c = 0, specA = 0, specB = 0;
	std::vector<Eigen::Triplet<double>> triplets11, triplets12, triplets21, triplets22;
	for (int i = 0; i < ne; i++) {

		for (int ii = 0; ii < 4; ii++) {
			for (int jj = 0; jj < 3; jj++) {
				int vidA = tetMesh->elements(i, ii);
				int gidA = dofIds(3 * vidA + jj);
				if (gidA < 0)continue;

				allocate(dofRules(gidA), r, specA);

				for (int m = 0; m < 4; m++) {
					for (int n = 0; n < 3; n++) {
						int vidB = tetMesh->elements(i, m);
						int gidB = dofIds(3 * vidB + n);
						if (gidB < 0) continue;

						allocate(dofRules(gidB), c, specB);
						switch (2 * r + c)
						{
						case 0:
							triplets11.emplace_back(specA, specB, Kelements[i](3 * ii + jj, 3 * m + n)*E);
							break;
						case 1:
							triplets12.emplace_back(specA, specB, Kelements[i](3 * ii + jj, 3 * m + n)*E);
							break;
						case 2:
							triplets21.emplace_back(specA, specB, Kelements[i](3 * ii + jj, 3 * m + n)*E);
							break;
						case 3:
							triplets22.emplace_back(specA, specB, Kelements[i](3 * ii + jj, 3 * m + n)*E);
							break;
						default:
							break;
						}

					}
				}
			}
		}
	}
	K11.setFromTriplets(triplets11.begin(), triplets11.end());
	K12.setFromTriplets(triplets12.begin(), triplets12.end());
	K21.setFromTriplets(triplets21.begin(), triplets21.end());
	K22.setFromTriplets(triplets22.begin(), triplets22.end());
	
}

void tpmsfem::femCalculator::compute(Eigen::MatrixXd & points)
{
	int nv = tetMesh->pointSize();
	int ne = tetMesh->elementSize();

	K11.makeCompressed();

	Eigen::SparseLU<Eigen::SparseMatrix<double>> solverK11;

	solverK11.compute(K11);

	//wave(K)
	Eigen::SparseMatrix<double> waveK;
	Eigen::SparseMatrix<double> inverseK11K12;
	inverseK11K12 = solverK11.solve(K12);
	waveK = K22 - K21 * inverseK11K12;

	//bar(u)

	Eigen::SparseLU<Eigen::SparseMatrix<double>> solverWaveK;
	solverWaveK.compute(waveK);

	baru = solverWaveK.solve(barf);

	hatu = solverK11.solve(hatf - K12 * baru);

	Eigen::VectorXd u(unConstrainedNum);

	for (int i = 0; i < unConstrainedNum; i++) {
		u(i) = dofRules(i) < unConstrainedNum ? hatu(dofRules(i)) : baru(dofRules(i) % unConstrainedNum);

		int nv = tetMesh->pointSize();

		points = tetMesh->points;

		for (int i = 0; i < nv; i++) {
			for (int j = 0; j < 3; j++) {
				int vid = dofIds(3 * i + j);

				if (vid < 0)
					continue;
				else
					points(i, j) += u(vid);
			}
		}
	}

}
