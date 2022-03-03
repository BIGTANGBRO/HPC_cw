#include <iostream>
#include <cmath>
#include <cblas.h>
#include <cstdlib>
using namespace std;

void fillMatrixUInitial(int &Nx, int &Ny, double *matrix){
	//Nx is row
	//Ny is column
	for (int i = 0;i < Nx;i++){
		for (int j = 0;j < Ny;j++){
			if (j > Ny){
				matrix[j*Nx + i] = 1;
			}else{
				matrxi[j*Nx + i] = 0;
			}
		}
			
	}
}

void printMatrix(){
	
}

int main()
{
	//data initialization
	double arg;
	int T;
	int Nx;
	int Nx;
	double a;
	double b;
	double mu1;
	double mu2;
	double eps;
	
	T = 100;
	dt = 0.001;
	Nx = 101;
	Ny = 101;
	a = 0.75;
	b = 0.06;
	mu1 = 5.0;
	mu2 = 5.0;
	
	double *matrix = new double[Nx*Ny];
	
	
	return 0;
}
