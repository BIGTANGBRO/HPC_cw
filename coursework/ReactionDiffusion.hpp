#ifndef REACTIONDIFFUSION_HPP
#define REACTIONDIFFUSION_HPP

class ReactionDiffusion
{
//declaration
public:
	ReactionDiffusion();
	ReactionDiffusion(double *u, double *v);
	~ReactionDiffusion();
	
	void SetParameters();
	void SetInitialConditions(int &Nx, int &Ny, double &a);
	void TimeIntegrations(); 
	
private:
	double *u;
	double *v;

};

#endif // REACTIONDIFFUSION_HPP
