#ifndef LMOPTIM_CLASS
#define LMOPTIM_CLASS

#include "Sweep.h"

class LMOptim
{
public:
	//variables
<<<<<<< HEAD
	// Numerical Jacobian
	const double delta = pow(10, -3);
	// Levenberg-Marquart
	const double BiSq_threshold = 10; //weighting func: w(x) = (1-(x/BiSq_threshold)^2)^2
	const int maxIterLM = 200; // Max iteration
	const double lambda_ini = 10; // Initial lambda
	const double lambda_scale = 10; // Lambda scale
	const double lambda_max = 10e6;// Max lambda
=======
	double delta = 0;
>>>>>>> 3c624b338322549276388e3a8f766ac207eaccbc

	//functions
	LMOptim();
	~LMOptim();
    double Distance2EdgePlane(LoamPt &pt, Sweep &OldSweep, VectorXd &EstTransform, int EnPflag);
	MatrixXd GetJacobian(VectorXd &DistanceVec, MatrixXd &W, Sweep &OldSweep, Sweep &NewSweep, VectorXd EstTransform);
	VectorXd GetDistanceVec(Sweep &OldSweep, Sweep &NewSweep, VectorXd EstTransform);
	VectorXd TransformEstimate(Sweep &OldSweep, Sweep &NewSweep);

};


#endif // !LMOPTIM_CLASS
