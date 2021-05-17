// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <cstdlib>
#include <stdlib.h>
#include <time.h>
#include <exception>
#include <chrono>
#include <stdio.h>
#include <filesystem>

#include "../dune/heat-DG/IncludesDUNE.hh"
#include "../dune/heat-DG/heat/include_problem.hh"
//TODO: Change problem name to: BM-phasechange-0d

int main(int argc, char **argv)
{
	try
	{
		// Maybe initialize MPI
		 Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

		if (argc != 2)
		{
			if (helper.rank() == 0)
			{
				std::cout << "usage: ./ex2 <input_file.ini> " << std::endl;
			}
			return 1;
		}
		std::string PATH = "/home/peiravim/dune/heat-DG/dune/heat-DG/heat/";
		char input[80];
	    sscanf(argv[1],"%39s", input);
	    std::string input_file = "/home/peiravim/dune/heat-DG/dune/heat-DG/heat/inputs/";
	    input_file += input;
	    std::cout<< "input file: " << input_file << std::endl ;

	    Dune::ParameterTree ptree;
	    Dune::ParameterTreeParser ptreeparser;
	    ptreeparser.readINITree(input_file,ptree);
	    ptreeparser.readOptions(argc,argv,ptree);

		std::string fileName = ptree.get("output.file_name",(std::string)"test");
		std::string pathName = ptree.get("output.path_name",(std::string)"test");
		pathName += "outputs/";
		pathName += fileName ;
		if (helper.rank() == 0)
		{
				std::filesystem::create_directory(pathName);
		}

		
		/**************************************************************************************************/
		// MESH
	    MeshParameters<Dune::ParameterTree> mesh(ptree);
	    const int dim = mesh.dimension;

		/*____________________________________________*/

		
		Dune::FieldVector<double, dim> L(0.0); // L represents the right top node of the rectangular/cuboidal domain
		L[0] = mesh.X_length;
		if (dim == 1)
		{
			L[0] = mesh.Z_length;
		}
		if (dim == 2)
		{
			L[1] = mesh.Z_length;
		}
		else if (dim == 3)
		{
			L[1] = mesh.Y_length;
			L[2] = mesh.Z_length;
		}
		std::array<int, dim> N(Dune::filledArray<dim, int>(1));
		N[0] = mesh.X_cells;
		if (dim == 1)
			N[0] = mesh.Z_cells;
		if (dim == 2)
			N[1] = mesh.Z_cells;
		else if (dim == 3)
		{
			N[1] = mesh.Y_cells;
			N[2] = mesh.Z_cells;
		}
#ifdef YASP
		typedef Dune::YaspGrid<dim> Grid;
		std::bitset<dim> periodic(false);

		int overlap = 1;
		std::shared_ptr<Grid> grid = std::shared_ptr<Grid>(new Grid(L, N, periodic, overlap, Dune::MPIHelper::getCollectiveCommunication()));
		grid->refineOptions(false); // keep overlap in cells

#elif defined(ALUGRID)
		typedef Dune::ALUGrid<dim, dim, Dune::cube, Dune::nonconforming> Grid;
		auto ll = Dune::FieldVector<Grid::ctype, dim>{{0, L[1]}};
		auto ur = Dune::FieldVector<Grid::ctype, dim>{{L[0], 0}};
		std::array<unsigned int, dim> elements;
		elements[0] = N[0];
		elements[1] = N[1];
		// std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createSimplexGrid(ll, ur, elements);
		std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createCubeGrid(ll, ur, elements); // load balance the grid
		// std::string filename = ptree.get("grid.alugrid.name",
        //                                  (std::string)"grid.msh");
		// auto grid_file = PATH;
		// grid_file += "grids/";
		// grid_file += filename;
        // Dune::GridFactory<Grid> factory;
        // Dune::GmshReader<Grid>::read(factory,grid_file,false,false);
        // std::shared_ptr<Grid> grid(factory.createGrid());

#elif defined(UG) 
		typedef  Dune::UGGrid<dim> Grid;
		auto ll = Dune::FieldVector<Grid::ctype, dim>{{0, L[1]}};
		auto ur = Dune::FieldVector<Grid::ctype, dim>{{L[0], 0}};
		std::array<unsigned int, dim> elements;
		elements[0] = N[0];
		elements[1] = N[1];
		// // std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createSimplexGrid(ll, ur, elements);
		std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createCubeGrid(ll, ur, elements); // load balance the grid
		// std::string filename = ptree.get("grid.alugrid.name",
        //                                  (std::string)"grid.msh");
		// auto grid_file = PATH;
		// grid_file += "grids/";
		// grid_file += filename;
        // Dune::GridFactory<Grid> factory;
        // Dune::GmshReader<Grid>::read(factory,grid_file,false,false);
        // std::shared_ptr<Grid> grid(factory.createGrid());

#endif

		typedef Grid::LeafGridView GV;
		GV gv = grid->leafGridView();
		// grid->globalRefine(1);
  		// Transfer partitioning from ParMETIS to our grid
  		grid->loadBalance();

		//grid->globalRefine(1);
		// Dune::VTKWriter<GV> vtkWriter(gv);
  		// vtkWriter.write(std::string("gridviews"));
		//exit(0);
		driver(gv, ptree, helper);
		// driver_Sh(gv, ptree, helper);

	}
	catch (Dune::Exception &e)
	{
		std::cerr << "Dune reported error: " << e << std::endl;
	}
	catch (...)
	{
		std::cerr << "Unknown exception thrown!" << std::endl;
	}
}
