#include<iostream>
#include<fstream>
#include<vector>
#include "TPMSGenerator.h"
#include "tpmsTetGenerator.h"
#include "femCalculator.h"
#include<Eigen/dense>


bool writeOBJ(const std::string & str, const Eigen::MatrixXd & V, const Eigen::MatrixXi & F,Eigen::MatrixXd &C)
{
	using namespace std;
	using namespace Eigen;
	assert(V.cols() == 3 && "V should have 3 columns");
	Eigen::MatrixXi FF(F.rows() * 4, 3);

	for (int i = 0; i < F.rows(); i++) {
		FF.row(4 * i) << F(i, 0), F(i, 2), F(i, 1);
		FF.row(4 * i + 1) << F(i, 0), F(i, 3), F(i, 2);
		FF.row(4 * i + 2) << F(i, 0), F(i, 1), F(i, 3);
		FF.row(4 * i + 3) << F(i, 1), F(i, 2), F(i, 3);
	}
	ofstream s(str);
	if (!s.is_open())
	{
		fprintf(stderr, "IOError: writeOBJ() could not open %s\n", str.c_str());
		return false;
	}
	s <<
		V.format(IOFormat(FullPrecision, DontAlignCols, " ", "\n", "v ", "", "", "\n")) <<
		(FF.array() + 1).format(IOFormat(FullPrecision, DontAlignCols, " ", "\n", "f ", "", "", "\n"));
	s.close();
	string strc = str.substr(0, str.length() - 4);
	ofstream s1(strc + ".col");
	s1 << C.rows() << " " << C.cols() << "\n";
	s1 << C.format(IOFormat(FullPrecision, DontAlignCols, " ", "\n", "", "", "", "\n"));
	s1.close();
	return true;
}

bool writeLOG(const std::string & str, const Eigen::VectorXd &param, int index = 0) {
	using namespace std;
	using namespace Eigen;

	ofstream csv;

	if (index == 0)
		csv.open(str);
	else
		csv.open(str, ios::app);

	if (!csv.is_open())
	{
		fprintf(stderr, "IOError: writeCSV() could not open %s\n", str.c_str());
		return false;
	}

	csv << param.format(IOFormat(FullPrecision, DontAlignCols, "", ",", "", "", "", "\n"));
	csv.close();
}

bool writePOLY(const std::string & str, const Eigen::MatrixXd & V, const Eigen::MatrixXi & E)
{
	using namespace std;
	using namespace Eigen;
	assert(V.cols() == 3 && "V should have 3 columns");
	ofstream sv(str + "_node.csv");
	if (!sv.is_open())
	{
		fprintf(stderr, "IOError: writeOBJ() could not open %s\n", str.c_str());
		return false;
	}
	for (int i = 0; i < V.rows(); i++)
		sv<< V(i, 0) << "," << V(i, 1) << "," << V(i, 2) << "\n";
	sv.close();
	ofstream se(str + "_ele.csv");
	if (!se.is_open())
	{
		fprintf(stderr, "IOError: writePOLY() could not open %s\n", str.c_str());
		return false;
	}
	auto EPLUS = (E.array() + 1).matrix();
	for (int i = 0; i < EPLUS.rows(); i++)
		se << EPLUS(i, 0) << "," << EPLUS(i, 1) << "," << EPLUS(i, 2) << "," << EPLUS(i, 3) << "\n";
	se.close();
	return true;
}


int main(){

	const auto &f = [](double x, double y, double z)->double {
		return 10 * (cos(x)*sin(y) + cos(y)*sin(z) + cos(z)*sin(x))
			- 0.5*(cos(2 * x)*cos(2 * y) + cos(2 * y)*cos(2 * z) + cos(2 * z)*cos(2 * x));
		//return cos(x) + cos(y) + cos(z);
	};

	//std::vector<double>Vmin = { 0,0,0 };
	//std::vector<double>Vmax = { 2 * PI,2 * PI,2 * PI };

	double x = 12.5, y = 2.5, z = 17.5;
	std::vector<double>Vmin = { (x - 10)*0.1*PI,(y - 10)*0.1*PI,(z - 10)*0.1*PI };
	std::vector<double>Vmax = { (x + 10)*0.1*PI,(y + 10)*0.1*PI,(z + 10)*0.1*PI };

	tpmsgen::TPMSGenerator<double, double> *tg=new tpmsgen::TPMSGenerator<double, double>(f, Vmin, Vmax, 15);

	tpmsfem::tpmsTetGenerator* tets=new tpmsfem::tpmsTetGenerator(tg);

	tpmsfem::femCalculator *fem=new tpmsfem::femCalculator();

	double E = 1.59e7;
	double mu = 0.45;
	for (int r = 0; r < 100; r++) {
		double levelset = 14 - r * 0.09;
		
		Eigen::MatrixXd V;
		Eigen::MatrixXi F;
		Eigen::MatrixXd C;
		tets->scale = 0.004;
		tets->tpmsToTetrahedron("pq1.2/10Y", levelset, V, F);
		C.setOnes(tets->points.rows(), 3);

		Eigen::VectorXd genE(1);
		Eigen::VectorXd genMu(1);

		for (int c = 0; c < 1; c++) {
			

			double maxZ = tets->points.col(2).maxCoeff();
			double minZ = tets->points.col(2).minCoeff();

			int nv = tets->pointSize();

			std::vector<int> pointrules(nv * 3, 0);

			std::vector<int> fixPoints;
			std::vector<int> specPoints;


			for (int i = 0; i < nv; i++) {
				if (abs(tets->points(i, 2) - minZ) <= 0.005*abs(maxZ - minZ)) {
					fixPoints.push_back(i);
					pointrules[3 * i] = -1;
					pointrules[3 * i + 1] = -1;
					pointrules[3 * i + 2] = -1;
					C.row(i) << 0, 0, 1;
				}

				if (abs(tets->points(i, 2) - maxZ) <= 0.005*abs(maxZ - minZ)) {
					specPoints.push_back(i);
					pointrules[3 * i + 2] = 1;
					C.row(i) << 1, 0, 0;
				}
			}

			Eigen::MatrixXd spec(specPoints.size(), 3);
			for (int i = 0; i < specPoints.size(); i++) spec.row(i) = tets->points.row(specPoints[i]);

			double minY = spec.col(1).minCoeff();
			double maxY = spec.col(1).maxCoeff();

			Eigen::MatrixXd disps(nv, 3);
			Eigen::MatrixXd forces(nv, 3);

			disps.setConstant(0.0);
			forces.setConstant(0.0);

			double load= tets->scale*(c + 1) / 20*sqrt(c+1);

			for (int i : specPoints) {
				disps.row(i) << 0, 0, 0;
				forces.row(i) << 0, 0, load;
			}

			Eigen::MatrixXd V_new;
			
			fem->init(tets, E, mu, pointrules, disps, forces);
			fem->compute(V_new);
			double avgu = fem->baru.sum() / fem->baru.rows();
			genE(c) = load*fem->baru.rows() / avgu;

			double rate = abs(avgu)/ (maxZ - minZ);

			genMu(c) = (1 - sqrt(1 / (1 + rate))) / rate;
			writeOBJ("../data/result.obj", V_new, tets->elements,C);
			std::cout << r << "\t" << c <<"\t"<<genE(c)<<"\t"<<genMu(c)<<std::endl;			
		}

		writeLOG("../data/E_size_"+std::to_string(tets->scale)+".csv", genE, r);
		writeLOG("../data/Mu_size_"+ std::to_string(tets->scale) + ".csv", genMu, r);
	}

    return 0;
}