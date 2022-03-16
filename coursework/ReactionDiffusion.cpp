#include "ReactionDiffusion.hpp"

ReactionDiffusion::ReactionDiffusion(double *u, double *v)
{
	this->u = u;
	this->v = v;
}

void ReactionDiffusion::SetInitialConditions(int &Nx, int &Ny, double &a){
	for (int i = 0;i < Nx;i++){
		for (int j = 0;j < Ny;j++){
			if (j > Ny/2){
				this->u[j*Nx + i] = 1;
			}else{
				this->u[j*Nx + i] = 0;
			}
		}
	}
	
	for (int i = 0;i < Nx;i++){
		for (int j = 0;j < Ny;j++){
			if (i < Nx/2){
				this->v[j*Nx + i] = a/2.0;
			}else{
				this->v[j*Nx + i] = 0;
			}
		}
	}
}

void ReactionDiffusion::SetParameters(){

}

void ReactionDiffusion::TimeIntegrations(){
	
}

ReactionDiffusion::~ReactionDiffusion()
{
}

