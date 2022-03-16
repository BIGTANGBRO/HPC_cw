#include <iostream>
#include <cmath>
#include <cblas.h>
#include <cstdlib>
using namespace std;

void fillMatrixUInitial(int &Nx, int &Ny, double *matrix){
	//Nx is column
	//Ny is row
	for (int i = 0;i < Ny;i++){
		for (int j = 0;j < Nx;j++){
			if (j > Nx/2){
				matrix[j*Ny + i] = 1;
			}else{
				matrix[j*Ny + i] = 0;
			}
		}
			
	}
}

void fillMatrixVInitial(int &Nx, int &Ny, double *matrix, double &a){
	//Nx is column 
	//Ny is row
	for (int i = 0;i < Ny;i++){
		for (int j = 0;j < Nx;j++){
			if (i < Ny/2){
				matrix[j*Ny + i] = a/2.0;
			}else{
				matrix[j*Ny + i] = 0;
			}
		}
	}			
}

void fillMatrix1(int &Nx, int &Ny, double *A, double &mu, double &h, double &dt){
	//Nx is row
	//Ny is column
	int j;
	for (int i = 1;i < Nx-1;i++){
		j = i-1;
		A[(0+j)*Ny + i] = 1.0 * mu / (h*h);
		A[(1+j)*Ny + i] = -2.0 * mu / (h*h);
		A[(2+j)*Ny + i] = 1.0 * mu / (h*h);
	}
	
	A[0*Ny + 0] = -1.0 * mu / (h*h);
	A[1*Ny + 0] = 1.0 * mu / (h*h);
	
	A[(Ny - 2)*Ny + (Ny-1)] = 1.0 * mu / (h*h);
	A[(Ny - 1)*Ny + (Ny-1)] = -1.0 * mu / (h*h);
}

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

void fillCol(int &col, int &Nx, int &Ny, double *Matrix, double *colH){
	for (int i = 0;i< Ny;i++){
		Matrix[col*Ny + i] = colH[i];
	}
}

void fillRow(int &row, int &Nx, int &Ny, double *Matrix, double *rowH){
	for (int i = 0;i< Nx;i++){
		Matrix[i*Ny + row] = rowH[i];
	}
}


void getCol(int &col, int &Nx, int &Ny, double *Matrix, double *colH){
	for (int i = 0;i< Ny;i++){
		colH[i] = Matrix[col*Ny + i];
	}
}

void getRow(int &row, int &Nx, int &Ny, double *Matrix, double *rowH){
	for (int i = 0;i< Nx;i++){
		rowH[i] = Matrix[i*Ny + row];
	}
}

double *getf1(int &Nx, int &Ny, double *u, double *v, double &eps, double &a, double &b){
	double *f1 = new double[Ny*Nx];
	
	for (int i = 0;i<Ny;i++){
		for (int j = 0;j<Nx;j++){
			f1[j*Ny+i] = eps*u[j*Ny + i]*(1.0-u[j*Ny+i])*(u[j*Ny+i]-(v[j*Ny+i]+b)/a);
		}
	}
	return f1;
}

double *getf2(int &Nx, int &Ny, double *u, double *v){
	double *f2 = new double[Ny*Nx];
	for (int i = 0;i<Ny;i++){
		for (int j = 0;j<Nx;j++){
			f2[j*Ny+i] = u[j*Ny + i] * u[j*Ny + i] * u[j*Ny + i] - v[j*Ny+i];
		}
	}
	return f2;
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
	double h = 1.0;
	
	T = 10;
	dt = 0.001; 
	Nx = 51;
	Ny = 51;
	a = 0.75;
	b = 0.06;
	mu1 = 5.0;
	mu2 = 0.0;
	eps = 50.0;
	
	//initialize
	double *u = new double[Ny*Nx];
	double *v = new double[Ny*Nx];
	double *A1 = new double[Ny*Nx];
	double *A2 = new double[Ny*Nx];
	
	//fill the matrix
	fillMatrixUInitial(Nx, Ny, u);
	fillMatrixVInitial(Nx, Ny, v, a);
	//for U
	fillMatrix1(Nx, Ny, A1, mu1, h, dt);
	//for V
	fillMatrix1(Ny, Nx, A2, mu2, h, dt);
	
	//then do the iteration
	for (int i = 0;i<T/dt;i++){
		double *tempU1 = new double[Ny*Nx];
		double *tempV1 = new double[Ny*Nx];
		
		for (int j = 0;j<Nx;j++){
			double *tempUCol1 = new double[Ny];			
			double *tempVCol1 = new double[Ny];

			double *uCol = new double[Ny];
			double *vCol = new double[Ny];
		
			getCol(j, Nx, Ny, u, uCol);
			getCol(j, Nx, Ny, v, vCol);
		
			cblas_dgemv(CblasColMajor, CblasNoTrans, Ny, Nx, 1.0, A1, Ny, uCol, 1, 0.0, tempUCol1 ,1);
			fillCol(j, Nx, Ny, tempU1, tempUCol1);		
		
			cblas_dgemv(CblasColMajor, CblasNoTrans, Ny, Nx, 1.0, A2, Ny, vCol, 1, 0.0, tempVCol1 ,1);
			fillCol(j, Nx, Ny, tempV1, tempVCol1);			
		}

		
		double *tempU2 = new double[Ny*Nx];
		double *tempV2 = new double[Ny*Nx];
		for (int j = 0;j<Ny;j++){
			double *tempURow2 = new double[Nx];			
			double *tempVRow2 = new double[Nx];

			double *uRow = new double[Nx];
			double *vRow = new double[Nx];
		
			getRow(j, Nx, Ny, u, uRow);
			getRow(j, Nx, Ny, v, vRow);

			cblas_dgemv(CblasColMajor, CblasTrans, Ny, Nx, 1.0, A1, Ny, uRow, 1, 0.0, tempURow2, 1);
			fillRow(j, Nx, Ny, tempU2, tempURow2);		
		
			cblas_dgemv(CblasColMajor, CblasTrans, Ny, Nx, 1.0, A2, Ny, vRow, 1, 0.0, tempVRow2, 1);
			fillRow(j, Nx, Ny, tempV2, tempVRow2);	
		}
	
		double *f1 = getf1(Nx, Ny, u, v, eps, a, b);
		double *f2 = getf2(Nx, Ny, u, v);
		
		printMatrix(f1,Nx,Ny);
		cblas_daxpy(Nx*Ny,dt,tempU1,1,u,1);
		cblas_daxpy(Nx*Ny,dt,tempU2,1,u,1);

		cblas_daxpy(Nx*Ny,dt,tempV1,1,v,1);
		cblas_daxpy(Nx*Ny,dt,tempV2,1,v,1);
	
		cblas_daxpy(Nx*Ny,dt,f1,1,u,1);
		cblas_daxpy(Nx*Ny,dt,f2,1,v,1);		
		printMatrix(u,Nx,Ny);
		cout << endl;
	
	}
	cout << u[1] << endl;

	return 0;
}