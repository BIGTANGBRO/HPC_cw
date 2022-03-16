#include <iostream>
#include <cmath>
#include <fstream>
#include <cblas.h>
#include <cstdlib>
#include <ReactionDiffusion.hpp>
using namespace std;

void printMatrix(double *H, int &Nx, int &Ny){
	for (int i = 0;i<Ny;i++){
		for (int j = 0;j<Nx;j++){
			cout << H[j*Ny + i] << " ";
		}
		cout << endl;
	}
}

void printVector(double *H, int &Nx){
	for (int i = 0;i<Nx;i++){
		cout << H[i] << " ";
	}
}

int main()	
{
	ReactionDiffusion reaction;
	reaction.SetInitialConditions();
	reaction.TimeIntegrations();
	reaction.writeInTxt();
	cout << "The process finished" << endl;
	return 0;
}