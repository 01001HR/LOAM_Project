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
	static double Distance2EdgePlane(LoamPt &pt, Sweep &OldSweep, VectorXd EstTransform, int EnPflag);
	MatrixXd GetJacobian(VectorXd DistanceVec, Sweep &OldSweep, Sweep &NewSweep, VectorXd EstTransform);
	VectorXd LMOptim::GetDistanceVec(Sweep &OldSweep, Sweep &NewSweep, VectorXd EstTransform);
	VectorXd TransformEstimate(Sweep &OldSweep, Sweep &NewSweep);

};


#endif // !LMOPTIM_CLASS
