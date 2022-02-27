#include "ReactionDiffusion.hpp"

ReactionDiffusion::ReactionDiffusion(double &u, double &v)
{
	this -> u = u;
	this -> v = v;
}

void ReactionDiffusion::SetInitialConditions(){
	
}

void ReactionDiffusion::SetParameters(double &u, double &v){
	this -> u = u;
	this -> v = v;
}

void ReactionDiffusion::TimeIntegrations(){
	
}

ReactionDiffusion::~ReactionDiffusion()
{
}

