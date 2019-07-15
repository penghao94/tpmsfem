#include "tpmsTetGenerator.h"

#include <cassert>
#include <iostream>

#ifndef TETLIBRARY
#define TETLIBRARY 
#endif
#include "tetgen.h"

bool tpmsfem::tpmsTetGenerator::tpmsToTetrahedron(const std::string & switches,double levelset,Eigen::MatrixXd &V, Eigen::MatrixXi& F)
{
	//generate a tpms triangle mesh by tpms generator
	
	

	tg_->makeLevelSet(0x3F, levelset, V, F);

	// divide tpms mesh into tetrahedron mesh

	tetgenio in;

	//assign tetgen input

	in.firstnumber = 0;
	in.numberofpoints = V.rows();
	in.pointlist = new double[in.numberofpoints * 3];

	for (int i = 0; i < V.rows(); i++) {
		in.pointlist[i * 3 + 0] = V(i, 0);
		in.pointlist[i * 3 + 1] = V(i, 1);
		in.pointlist[i * 3 + 2] = V(i, 2);
	}

	in.numberoffacets = F.rows();
	in.facetlist = new tetgenio::facet[in.numberoffacets];
	in.facetmarkerlist = new int[in.numberoffacets];

	//loop over face
	for (int i = 0; i < F.rows(); i++) {
		in.facetmarkerlist[i] = i;
		tetgenio::facet * f = &in.facetlist[i];
		f->numberofpolygons = 1;
		f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
		f->numberofholes = 0;
		f->holelist = NULL;
		tetgenio::polygon * p = &f->polygonlist[0];
		p->numberofvertices = F.cols();
		p->vertexlist = new int[p->numberofvertices];
		// loop around face
		for (int j = 0; j < F.cols(); j++)
		{
			p->vertexlist[j] = F(i, j);
		}
	}

	//parser tetgen out

	tetgenio out;
	try
	{
		char * cswitches = new char[switches.size() + 1];
		strcpy(cswitches, switches.c_str());
		::tetrahedralize(cswitches, &in, &out);
		delete[] cswitches;
	}
	catch (int e)
	{
		std::cerr << "^" << __FUNCTION__ << ": TETGEN CRASHED...KABOOM!!" << std::endl;
		return false;
	}
	if (out.numberoftetrahedra == 0)
	{
		std::cerr << "^" << __FUNCTION__ << ": Tetgen failed to create tets" << std::endl;
		return false;
	}

	if (out.pointlist == NULL)
	{
		std::cerr << "^tetgenio_to_tetmesh Error: point list is NULL\n" << std::endl;
		return false;
	}


	points.resize(0, 0);
	points = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(out.pointlist, out.numberofpoints, 3);
	points = points/(2*PI)*scale;
	elements.resize(0, 0);
	elements= Eigen::Map<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(out.tetrahedronlist, out.numberoftetrahedra, 4);

	facets.resize(0, 0);
	facets = Eigen::Map<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(out.trifacelist, out.numberoftrifaces, 3);

	volumes.resize(out.numberoftetrahedra);
	for (int i = 0; i < out.numberoftetrahedra; i++) {
		Eigen::Vector3d a = points.row(elements(i, 0));
		Eigen::Vector3d b = points.row(elements(i, 1));
		Eigen::Vector3d c = points.row(elements(i, 2));
		Eigen::Vector3d d = points.row(elements(i, 3));
		volumes(i) = 1.0 / 6.0 * std::abs((d - a).dot((b - a).cross(c - a)));
	}
	return true;
}

void tpmsfem::tpmsTetGenerator::preComputation(tpmsTetGenerator * tets, Eigen::MatrixXd& De, std::vector<Eigen::MatrixXd>& Kelements)
{
	int nv = tets->pointSize();
	int ne = tets->elementSize();

	//percompute lame cofficient for each tet
	Kelements.resize(ne);
	for (int i = 0; i < ne; i++) {
		double MInv[16];
		Eigen::MatrixXd M;
		M.setOnes(4, 4);
		for (int v = 0; v < 4; v++) {
			int vid = tets->elements(i, v);
			M.block(0, v, 3, 1) = tets->points.row(vid).transpose();
		}

		Eigen::MatrixXd temp;
		temp.noalias() = M.inverse();

		for (int ii = 0; ii < 4; ii++)
			for (int jj = 0; jj < 4; jj++)
				MInv[4 * ii + jj] = temp(ii, jj);


		Eigen::MatrixXd Be(6, 12);

		Be << MInv[0], 0, 0, MInv[4], 0, 0, MInv[8], 0, 0, MInv[12], 0, 0,
			0, MInv[1], 0, 0, MInv[5], 0, 0, MInv[9], 0, 0, MInv[13], 0,
			0, 0, MInv[2], 0, 0, MInv[6], 0, 0, MInv[10], 0, 0, MInv[14],
			MInv[1], MInv[0], 0, MInv[5], MInv[4], 0, MInv[9], MInv[8], 0, MInv[13], MInv[12], 0,
			0, MInv[2], MInv[1], 0, MInv[6], MInv[5], 0, MInv[10], MInv[9], 0, MInv[14], MInv[13],
			MInv[2], 0, MInv[0], MInv[6], 0, MInv[4], MInv[10], 0, MInv[8], MInv[14], 0, MInv[12];



		//Ke = Be' * D * Be * Ve
		Kelements[i].noalias() = Be.transpose()*De*Be*tets->volumes(i);

	}
}
