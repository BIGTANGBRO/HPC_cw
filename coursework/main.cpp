#include <iostream>
#include <cmath>
#include <fstream>
#include <cblas.h>
#include <cstdlib>
#include "omp.h"
#include "ReactionDiffusion.hpp"
#include <time.h>
#include <boost/program_options.hpp>

namespace po = boost::program_options;
//main program, author:Jiaxuan Tang


void printMatrix(double *H, int &Nx, int &Ny){
	for (int i = 0;i<Ny;i++){
		for (int j = 0;j<Nx;j++){
			std::cout << H[j*Ny + i] << " ";
		}
		std::cout << std::endl;
	}
}

void printVector(double *H, int &Nx){
	for (int i = 0;i<Nx;i++){
		std::cout << H[i] << " ";
	}
}

int main(int argc, char *argv[])	
{
	//read from the command line
	po::options_description opts(
		"Available options.");
	opts.add_options()
		("dt", po::value<double>()->default_value(0.001),
		 "dt value.")
		("T", po::value<int>()->default_value(100),
		 "T value.")
		("Nx", po::value<int>()->default_value(101),
		 "Nx value.")
		("Ny", po::value<int>()->default_value(101),
		 "Ny value.")
		("a", po::value<double>()->default_value(0.75),
		 "a value.")
		("b", po::value<double>()->default_value(0.06),
		 "b value.")
		("mu1", po::value<double>()->default_value(5.0),
		 "mu1 value.")
		("mu2", po::value<double>()->default_value(0.0),
		 "mu2 value.")
		("eps", po::value<double>()->default_value(50.0),
		 "eps value.");
	
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, opts), vm);
	po::notify(vm);
	if (vm.count("Help")){
		std::cout << "opts" << std::endl;
	}
	
	//variables
	double dt = vm["dt"].as<double>();
	int T = vm["T"].as<int>();
	int Nx = vm["Nx"].as<int>();
	int Ny = vm["Ny"].as<int>();
	double a = vm["a"].as<double>();
	double b = vm["b"].as<double>();
	double mu1 = vm["mu1"].as<double>();
	double mu2 = vm["mu2"].as<double>();
	double eps = vm["eps"].as<double>();
	
	clock_t start, end;
	start = clock();
	
	//calculation process
	ReactionDiffusion reaction;
	reaction.SetParameters(Nx, Ny, T, dt, a, b, mu1, mu2, eps);
	reaction.SetInitialConditions();
	reaction.TimeIntegrations();
	reaction.writeInTxt();
	
	//time calculation for serial
	end = clock();
	std::cout << "Time = " << double(end-start)/CLOCKS_PER_SEC << "s" << std::endl;
	std::cout << "The process finished" << std::endl;
	return 0;
}