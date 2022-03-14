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
			if (j > Ny/2){
				matrix[j*Nx + i] = 1;
			}else{
				matrix[j*Nx + i] = 0;
			}
		}
			
	}
}

void fillMatrixVInitial(int &Nx, int &Ny, double *matrix, double &a){
	//Nx is row
	//Ny is column
	for (int i = 0;i < Nx;i++){
		for (int j = 0;j < Ny;j++){
			if (i < Nx/2){
				matrix[j*Nx + i] = a/2.0;
			}else{
				matrix[j*Nx + i] = 0;
			}
		}
	}			
}

void fillMatrix1(int &Nx, int &Ny, double *A, double &mu, double &h, double &dt){
	//Nx is row
	//Ny is column
	int j;
	for (int i = 1;i < Ny-1;i++){
		j = i-1;
		A[(0+j)*Nx + i] = 1.0 * h*h * dt / mu;
		A[(1+j)*Nx + i] = 4.0 * h*h * dt / mu;
		A[(2+j)*Nx + i] = 1.0 * h*h * dt / mu;
	}
	
	A[0*Nx + 0] = -1.0 * h*h * dt / mu;
	A[1*Nx + 0] = 1.0 * h*h * dt / mu;
	
	A[(Nx - 2)*Nx + (Nx-1)] = -1.0 * h*h * dt / mu;
	A[(Nx - 1)*Nx + (Nx-1)] = 1.0 * h*h * dt / mu;
}

void fillMatrix2(int &Nx, int &Ny, double *B, double &mu, double &h, double &dt){
	int j;
	for (int i = 1;i < Ny-1;i++){
		j = i-1;
		B[(0+j)*Nx + i] = 1.0 * h*h * dt / mu;
		B[(2+j)*Nx + i] = 1.0 * h*h * dt / mu;
	}
}

void printMatrix(double *H, int &Nx, int &Ny){
	for (int i = 0;i<Nx;i++){
		for (int j = 0;j<Ny;j++){
			cout << H[j*Nx + i] << " ";
		}
		cout << endl;
	}
}

int main()
{
	//data initialization
	double arg;
	double dt;
	int T;
	int Nx;
	int Ny;
	double a;
	double b;
	double mu1;
	double mu2;
	double eps;
	double h = 0.5;
	
	T = 100;
	dt = 0.001; 
	Nx = 101;
	Ny = 101;
	a = 0.75;
	b = 0.06;
	mu1 = 5.0;
	mu2 = 0.0;

	//initialize
	double *uInit = new double[Nx*Ny];
	double *vInit = new double[Nx*Ny];
	double *m1 = new double[Nx*Ny];
	double *m2 = new double[Nx*Ny];
	
	//fill the matrix
	fillMatrixUInitial(Nx, Ny, uInit);
	fillMatrixVInitial(Nx, Ny, vInit, a);
	fillMatrix1(Nx, Ny, m1, mu1, h, dt);
	fillMatrix2(Nx, Ny, m2, mu1, h, dt);
		
	//then do the iteration
	//solve for 1 iter
	
	double *eqn1 = new double[Nx*Ny];
	
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, Nx, Ny, Nx, 1.0, m1, Nx, uInit, Nx, 0.0, eqn1, Nx);
	printMatrix(eqn1, Nx, Ny);
	
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, Nx, Ny, Nx, 1.0, m2, Nx, uInit, Nx, 0.0, eqn1, Nx);

	
	
	return 0;
}