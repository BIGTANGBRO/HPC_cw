#include "ReactionDiffusion.hpp"
#include <cblas.h>
#include <cstdlib>
#include <fstream>
using namespace std;

ReactionDiffusion::ReactionDiffusion()
{
	this ->Nx = 50;
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

void ReactionDiffusion::SetParameters(int &Nx, int &Ny, int &T, double &dt, double &a, double &b, double &mu1, double &mu2, double &eps, double &h){
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
	this -> h = h;
}

void ReactionDiffusion::SetInitialConditions(){
	//Nx is column
	//Ny is row
	for (int i = 0;i < Ny;i++){
		for (int j = 0;j < Nx;j++){
			if (j > Nx/2){
				u[j*Ny + i] = 1;
			}else{
				u[j*Ny + i] = 0;
			}
		}
	}
	
	for (int i = 0;i < Ny;i++){
		for (int j = 0;j < Nx;j++){
			if (i < Ny/2){
				v[j*Ny + i] = this->a/2.0;
			}else{
				v[j*Ny + i] = 0;
			}
		}
	}	
}

double* ReactionDiffusion::fillMatrixNy(double &mu){
	double *A = new double[this->Ny * this->Ny];
	int j;
	for (int i = 1;i < this->Ny-1;i++){
		j = i-1;
		A[(0+j)*Ny + i] = 1.0 * mu / (h*h);
		A[(1+j)*Ny + i] = -2.0 * mu / (h*h);
		A[(2+j)*Ny + i] = 1.0 * mu / (h*h);
	}
	
	A[0*Ny + 0] = -1.0 * mu / (h*h);
	A[1*Ny + 0] = 1.0 * mu / (h*h);
	
	A[(Ny - 2)*Ny + (Ny-1)] = 1.0 * mu / (h*h);
	A[(Ny - 1)*Ny + (Ny-1)] = -1.0 * mu / (h*h);
	
	return A;
}

double* ReactionDiffusion::fillMatrixNx(double &mu){
	double *B = new double[this->Nx * this->Nx];
	int j;
	for (int i = 1;i < this->Nx-1;i++){
		j = i-1;
		B[(0+j)*Nx + i] = 1.0 * mu / (h*h);
		B[(1+j)*Nx + i] = -2.0 * mu / (h*h);
		B[(2+j)*Nx + i] = 1.0 * mu / (h*h);
	}
	
	B[0*Nx + 0] = -1.0 * mu / (h*h);
	B[1*Nx + 0] = 1.0 * mu / (h*h);
	
	B[(Nx - 2)*Nx + (Nx-1)] = 1.0 * mu / (h*h);
	B[(Nx - 1)*Nx + (Nx-1)] = -1.0 * mu / (h*h);
	
	return B;
}

double* ReactionDiffusion::getf1(){
	double *f1 = new double[Ny*Nx];
	for (int i = 0;i<Ny;i++){
		for (int j = 0;j<Nx;j++){
			f1[j*Ny+i] = eps*u[j*Ny + i]*(1.0-u[j*Ny+i])*(u[j*Ny+i]-(v[j*Ny+i]+b)/a);
		}
	}
	return f1;
} 

double* ReactionDiffusion::getf2(){
	double *f2 = new double[Ny*Nx];
	for (int i = 0;i<Ny;i++){
		for (int j = 0;j<Nx;j++){
			f2[j*Ny+i] = u[j*Ny + i] * u[j*Ny + i] * u[j*Ny + i] - v[j*Ny+i];
		}
	}
	return f2;
}

void ReactionDiffusion::TimeIntegrations(){
	double *A1 = fillMatrixNy(this->mu1);
	double *A2 = fillMatrixNy(this->mu2);
	double *B1 = fillMatrixNx(this->mu1);
	double *B2 = fillMatrixNx(this->mu2);
	
	//then do the iteration
	for (int i = 0;i<T/dt;i++){
		double *tempU1 = new double[Ny*Nx];
		double *tempV1 = new double[Ny*Nx];
		
		double *tempU2 = new double[Ny*Nx];
		double *tempV2 = new double[Ny*Nx];
		
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, Ny, Nx, Ny, 1.0, A1, Ny, u, Ny, 0.0, tempU1, Ny);
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, Ny, Nx, Ny, 1.0, A2, Ny, v, Ny, 0.0, tempV1, Ny);
		
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, Ny, Nx, Nx, 1.0, u, Ny, B1, Nx, 0.0, tempU2, Ny);
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, Ny, Nx, Nx, 1.0, v, Ny, B2, Nx, 0.0, tempV2, Ny);
	
		double *f1 = this->getf1();
		double *f2 = this->getf2();
		
		cblas_daxpy(Nx*Ny,dt,tempU1,1,u,1);
		cblas_daxpy(Nx*Ny,dt,tempU2,1,u,1);

		cblas_daxpy(Nx*Ny,dt,tempV1,1,v,1);
		cblas_daxpy(Nx*Ny,dt,tempV2,1,v,1);
	
		cblas_daxpy(Nx*Ny,dt,f1,1,u,1);
		cblas_daxpy(Nx*Ny,dt,f2,1,v,1);	
	}
	
}

void ReactionDiffusion::writeInTxt(){
	ofstream outfile;
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

