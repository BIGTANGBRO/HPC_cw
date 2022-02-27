#ifndef REACTIONDIFFUSION_HPP
#define REACTIONDIFFUSION_HPP

class ReactionDiffusion
{
//declaration
public:
	ReactionDiffusion();
	~ReactionDiffusion();
	void SetParameters(double &u, double &v);
	void SetInitialConditions();
	void TimeIntegrations(); 
	
private:
	double u;
	double v;

};

#endif // REACTIONDIFFUSION_HPP
