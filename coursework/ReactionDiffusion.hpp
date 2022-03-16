#ifndef REACTIONDIFFUSION_HPP
#define REACTIONDIFFUSION_HPP

class ReactionDiffusion
{
//declaration
public:
	ReactionDiffusion();
	~ReactionDiffusion();
		
	//methods
	void SetParameters(int &Nx, int &Ny, int &T, double &dt, double &a, double &b, double &mu1, double &mu2, double &eps, double &h);
	void SetInitialConditions();

	void TimeIntegrations(); 
	void writeInTxt();


private:
	int Nx;
	int Ny;
	double *u;
	double *v;
	int T;
	double dt;
	double a;
	double b;
	double mu1;
	double mu2;
	double eps;
	double h;
	
	double *fillMatrixNy(double &mu);
	double *fillMatrixNx(double &mu);
	double *getf1();
	double *getf2();
};

#endif // REACTIONDIFFUSION_HPP
