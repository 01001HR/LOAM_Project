#include "LMOptim.h"

LMOptim::LMOptim() {

}

LMOptim::~LMOptim() {

}

double LMOptim::Distance2EdgePlane(LoamPt &pt, Sweep &NewSweep, Sweep &OldSweep, VectorXd EstTransform, int EnPflag) {
	// EnPflag = 1: Edge | EnPflag = 2: Plane
	Vector3d xi = pt.xyz;
	Vector3d xj = OldSweep.ptCloud[pt.nearPt1[0]][pt.nearPt1[1]].xyz;
	Vector3d xl = OldSweep.ptCloud[pt.nearPt2[0]][pt.nearPt2[1]].xyz;
	Vector3d T_trans, T_rot, omega, xi_hat;
	Matrix3d eye3 = Matrix3d::Identity(), omega_hat, R;
	double Interp, theta, Distance;

	Interp = (NewSweep.timeStamps[pt.sliceID] - NewSweep.tStart) / (NewSweep.tCur - NewSweep.tStart);

	double Tmax = abs(EstTransform(1)*Interp);
	for (int Idx = 1; Idx < 6; Idx++) {
		if (abs(EstTransform(Idx)*Interp) > Tmax) {
			Tmax = abs(EstTransform(Idx)*Interp);
		}
	}
	if (Tmax < pow(10, -5)) {
		xi_hat = xi;
	}
	else {
		T_rot << EstTransform(3), EstTransform(4), EstTransform(5);
		T_rot = T_rot*Interp;
        theta = T_rot.norm();
		T_trans << EstTransform(0), EstTransform(1), EstTransform(2);
		omega << EstTransform(3) / theta, EstTransform(4) / theta, EstTransform(5) / theta;
		omega_hat << 0, -omega(2), omega(1),
			omega(2), 0, -omega(0),
			-omega(1), omega(0), 0;
		R = eye3 + omega_hat*sin(theta) + omega_hat*omega_hat*(1 - cos(theta));
		R.transposeInPlace();
		xi_hat = R*(xi - T_trans);
	}
	if (EnPflag == 1) {
		Distance = ((xi_hat - xj).cross(xi_hat - xl)).norm() / (xj - xl).norm();
	}
	else if (EnPflag == 2) {
		Vector3d xm = OldSweep.ptCloud[pt.nearPt3[0]][pt.nearPt3[1]].xyz;
		Distance = abs((xi_hat - xj).dot((xj - xl).cross(xj - xm))) / ((xj - xl).cross(xj - xm)).norm();
	}
	else {
		cout << "Edge or Plane indicator not provided!" << endl;
	}

	return Distance;
}

MatrixXd LMOptim::GetJacobian(VectorXd DistanceVectorEig, Sweep &OldSweep, Sweep &NewSweep,
	VectorXd EstTransform) {

	vector<vector<double>> Jacobian_Full;
	vector<double> JacobianRow(6), DistanceVector;
	VectorXd InterpTransform_Delta, InterpTransform;
	double Delta = pow(10, -6), OldDistance;
	int EnPflag = 1; //edge distance
	for (int i = 0; i < NewSweep.edgePts.size(); i++) {
		InterpTransform = EstTransform*(NewSweep.timeStamps[i] - NewSweep.tStart) / (NewSweep.tEnd - NewSweep.tStart);
		for (auto &entry : NewSweep.edgePts[i]) {
			OldDistance = Distance2EdgePlane(NewSweep.ptCloud[i][entry], NewSweep, OldSweep, InterpTransform, EnPflag);
			for (int col = 0; col < 6; col++) {
				InterpTransform_Delta = InterpTransform;
				InterpTransform_Delta(col) = InterpTransform(col) + Delta;
				JacobianRow[col] = (Distance2EdgePlane(NewSweep.ptCloud[i][entry], NewSweep, OldSweep, InterpTransform_Delta, EnPflag) -
					OldDistance) / Delta;
			}
			Jacobian_Full.push_back(JacobianRow);
			DistanceVector.push_back(OldDistance);
		}
	}
	EnPflag = 2; //plane distance
	for (int i = 0; i < NewSweep.planePts.size(); i++) {
		InterpTransform = EstTransform*(NewSweep.timeStamps[i] - NewSweep.tStart) / (NewSweep.tEnd - NewSweep.tStart);
		for (auto &entry : NewSweep.planePts[i]) {
			Distance2EdgePlane(NewSweep.ptCloud[i][entry], NewSweep, OldSweep, InterpTransform, EnPflag);
			for (int col = 0; col < 6; col++) {
				InterpTransform_Delta = InterpTransform;
				InterpTransform_Delta(col) = InterpTransform(col) + Delta;
				JacobianRow[col] = (Distance2EdgePlane(NewSweep.ptCloud[i][entry], NewSweep, OldSweep, InterpTransform_Delta, EnPflag) -
					OldDistance) / Delta;
			}
			Jacobian_Full.push_back(JacobianRow);
			DistanceVector.push_back(OldDistance);
		}
	}
	MatrixXd JacobianFullEigen(Jacobian_Full.size(), 6);
	for (int m = 0; m < Jacobian_Full.size(); m++) {
		for (int n = 0; n < 6; n++) {
			JacobianFullEigen(m, n) = Jacobian_Full[m][n];
		}
		DistanceVectorEig(m) = DistanceVector[m];
	}
	return JacobianFullEigen;
}

VectorXd LMOptim::TransformEstimate(Sweep &OldSweep, Sweep &NewSweep) {
	double lambda = 1;
	double lambda_scale = 10;
	VectorXd TransformInitial;


	return TransformInitial;


}
