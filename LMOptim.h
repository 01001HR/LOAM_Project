#ifndef LMOPTIM_CLASS
#define LMOPTIM_CLASS



#include "Sweep.h"

class LMOptim
{
public:
	//variables


	//functions
	LMOptim();
	~LMOptim();
	static double Distance2EdgePlane(LoamPt &pt, Sweep &NewSweep, Sweep &OldSweep, VectorXd EstTransform, int EnPflag);
	MatrixXd GetJacobian(VectorXd DistanceVector, Sweep &OldSweep, Sweep &NewSweep, VectorXd EstTransform);
	VectorXd TransformEstimate(Sweep &OldSweep, Sweep &NewSweep);

};


#endif // !LMOPTIM_CLASS
