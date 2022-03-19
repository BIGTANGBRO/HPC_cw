#include "ReactionDiffusion.hpp"
#include <cblas.h>
#include <cstdlib>
#include <fstream>
#include "omp.h"
using namespace std;
//author: Jiaxuan Tang

ReactionDiffusion::ReactionDiffusion()
{
	this ->Nx = 101;
	this ->Ny = 101;
	this ->u = new double[this->Nx*this->Ny];
	this ->v = new double[this->Nx*this->Ny];
	this ->T = 100;
	this ->dt = 0.001;
	this ->a = 0.75;
	this ->b = 0.06;
	this ->mu1 = 5.0;
	this ->mu2 = 0.0;
	this ->eps = 50.0;
	this ->h = 1.0;
}

void ReactionDiffusion::SetParameters(int &Nx, int &Ny, int &T, double &dt, double &a, double &b, double &mu1, double &mu2, double &eps){
	this -> Nx = Nx;
	this -> Ny = Nx;
	this -> u = new double[this->Nx*this->Ny];
	this -> v = new double[this->Nx*this->Ny];
	this -> T = T;
	this -> dt = dt;
	this -> a = a;
	this -> b = b;
	this -> mu1 = mu1;
	this -> mu2 = mu2;
	this -> eps = eps;
}

void ReactionDiffusion::SetInitialConditions(){
	//Nx is column
	//Ny is row
	for (int i = 0;i < Ny;i++){
		for (int j = 0;j < Nx;j++){
			if (j > Nx/2){
				this->u[j*Ny + i] = 1;
			}else{
				this->u[j*Ny + i] = 0;
			}
		}
	}
	
	for (int i = 0;i < Ny;i++){
		for (int j = 0;j < Nx;j++){
			if (i < Ny/2){
				this->v[j*Ny + i] = this->a/2.0;
			}else{
				this->v[j*Ny + i] = 0;
			}
		}
	}	
}

//fill the big matrix1
double* ReactionDiffusion::fillMatrixNy(double &mu){
	double *A = new double[this->Ny * this->Ny];
	int i;
	
	//#pragma omp parallel for
	for (i = 1;i < this->Ny-1;i++){
		//A[(i-1)*Ny + i] = 1.0 * mu / (h*h);
		A[(i)*Ny + i] = -2.0 * mu / (h*h);
		A[(1+i)*Ny + i] = 1.0 * mu / (h*h);
	}
	
	A[0*Ny + 0] = -1.0 * mu / (h*h);
	A[1*Ny + 0] = 1.0 * mu / (h*h);
	
	//A[(Ny - 2)*Ny + (Ny-1)] = 1.0 * mu / (h*h);
	A[(Ny - 1)*Ny + (Ny-1)] = -1.0 * mu / (h*h);
	
	return A;
}

//fill the big matrix2
double* ReactionDiffusion::fillMatrixNx(double &mu){
	double *B = new double[this->Nx * this->Nx];
	int i;
	
	//#pragma omp parallel for
	for (i = 1;i < this->Nx-1;i++){
		//B[(i-1)*Nx + i] = 1.0 * mu / (h*h);
		B[(i)*Nx + i] = -2.0 * mu / (h*h);
		B[(i+1)*Nx + i] = 1.0 * mu / (h*h);
	}
	
	B[0*Nx + 0] = -1.0 * mu / (h*h);
	B[1*Nx + 0] = 1.0 * mu / (h*h);
	
	//B[(Nx - 2)*Nx + (Nx-1)] = 1.0 * mu / (h*h);
	B[(Nx - 1)*Nx + (Nx-1)] = -1.0 * mu / (h*h);
	
	return B;
}

//f1 calculation
double* ReactionDiffusion::getf1(){
	double *f1 = new double[Ny*Nx];
	int i = 0;
	int j = 0;
	
	//for each element
	//#pragma omp parallel for
	for (i = 0;i<Ny;i++){
		for (j = 0;j<Nx;j++){
			f1[j*Ny+i] = eps*u[j*Ny + i]*(1.0-u[j*Ny+i])*(u[j*Ny+i]-(v[j*Ny+i]+b)/a);
		}
	}
	
	return f1;
} 

//f2 calculation
double* ReactionDiffusion::getf2(){
	double *f2 = new double[Ny*Nx];
	int i = 0;
	int j = 0;
	
	//#pragma omp parallel for
	for (i = 0;i<Ny;i++){
		for (j = 0;j<Nx;j++){
			f2[j*Ny+i] = u[j*Ny + i] * u[j*Ny + i] * u[j*Ny + i] - v[j*Ny+i];
		}
	}
	
	return f2;
}

void ReactionDiffusion::TimeIntegrations(){
	//initilization
	double *A1 = fillMatrixNy(this->mu1);
	double *A2 = fillMatrixNy(this->mu2);
	double *B1 = fillMatrixNx(this->mu1);
	double *B2 = fillMatrixNx(this->mu2);
	
	//Do the iteration
	for (int i = 0;i<T/dt;i++){
		double *tempU1 = new double[Ny*Nx];
		double *tempV1 = new double[Ny*Nx];
		
		double *tempU2 = new double[Ny*Nx];
		double *tempV2 = new double[Ny*Nx];
		
		//using symmetrical matrix, calculate u in x and y directions
		#pragma omp parallel
		{
			#pragma omp single
			{
				#pragma omp task
				cblas_dsymm(CblasColMajor, CblasLeft, CblasUpper, Ny, Nx, 1.0, A1, Ny, u, Ny, 0.0, tempU1, Ny);
				//cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, Ny, Nx, Ny, 1.0, A1, Ny, u, Ny, 0.0, tempU1, Ny);
								
				#pragma omp task
				cblas_dsymm(CblasColMajor, CblasLeft, CblasUpper, Ny, Nx, 1.0, A2, Ny, v, Ny, 0.0, tempV1, Ny);
				//cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, Ny, Nx, Ny, 1.0, A2, Ny, v, Ny, 0.0, tempV1, Ny);

				#pragma omp task
				cblas_dsymm(CblasColMajor, CblasRight, CblasUpper, Ny, Nx, 1.0, B1, Nx, u, Ny, 0.0, tempU2, Ny);
				//cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, Ny, Nx, Nx, 1.0, u, Ny, B1, Nx, 0.0, tempU2, Ny);
				
				#pragma omp task
				cblas_dsymm(CblasColMajor, CblasRight, CblasUpper, Ny, Nx, 1.0, B2, Nx, v, Ny, 0.0, tempV2, Ny);
				//cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, Ny, Nx, Nx, 1.0, v, Ny, B2, Nx, 0.0, tempV2, Ny);
			}
		}
		
		double *f1 = this->getf1();
		double *f2 = this->getf2();
		
		//add together to get u and v for each iteration
		#pragma omp parallel
		{
			#pragma omp single
			{
				#pragma omp task
				cblas_daxpy(Nx*Ny,dt,tempU1,1,u,1);
				#pragma omp task
				cblas_daxpy(Nx*Ny,dt,tempU2,1,u,1);
				
				#pragma omp task
				cblas_daxpy(Nx*Ny,dt,tempV1,1,v,1);
				#pragma omp task
				cblas_daxpy(Nx*Ny,dt,tempV2,1,v,1);
				
				#pragma omp task
				cblas_daxpy(Nx*Ny,dt,f1,1,u,1);
				#pragma omp task
				cblas_daxpy(Nx*Ny,dt,f2,1,v,1);	
				
			}
		}
	}	
}

void ReactionDiffusion::writeInTxt(){
	ofstream outfile;
	//write in local directory
	outfile.open("output.txt");
	for (int i = 0;i < Ny;i++){
		for (int j = 0;j < Nx;j++){
			outfile << "x" << j << " " << "y" << i << " " << this->u[j*Ny + i] << " " << this->v[j*Ny + i] << endl;
		}
		outfile << endl;
	}
	
	outfile.close();
}

ReactionDiffusion::~ReactionDiffusion()
{
}

