#include "LMOptim.h"

LMOptim::LMOptim() {

}

LMOptim::~LMOptim() {

}

double LMOptim::Distance2EdgePlane(LoamPt &pt, Sweep &OldSweep, VectorXd EstTransform, int EnPflag) {
	// EnPflag = 1: Edge | EnPflag = 2: Plane
	Vector3d xi = pt.xyz;
	Vector3d xj = OldSweep.ptCloud[pt.nearPt1[0]][pt.nearPt1[1]].xyz;
	Vector3d xl = OldSweep.ptCloud[pt.nearPt2[0]][pt.nearPt2[1]].xyz;
	Vector3d T_trans, T_rot, omega, xi_hat;
	Matrix3d eye3 = Matrix3d::Identity(), omega_hat, R;
	double theta, Distance;
	double Tmax = abs(EstTransform(1));
	for (int Idx = 1; Idx < 6; Idx++) {
		if (abs(EstTransform(Idx)) > Tmax) {
Tmax = abs(EstTransform(Idx));
		}
	}
	if (Tmax < pow(10, -5)) {
		xi_hat = xi;
	}
	else {
		T_rot << EstTransform(3), EstTransform(4), EstTransform(5);
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

MatrixXd LMOptim::GetJacobian(VectorXd DistanceVec, Sweep &OldSweep, Sweep &NewSweep, VectorXd EstTransform) {
	MatrixXd Jacobian = MatrixXd::Zero(NewSweep.numEdges + NewSweep.numPlanes, 6);
	VectorXd InterpTransform_Delta, InterpTransform;
	DistanceVec = VectorXd::Zero(NewSweep.numEdges + NewSweep.numPlanes);
	double Delta = pow(10, -6); //numerical Jacobian step-size
	double OldDistance;
	int cnt = 0;
	int SliceID;
	//edge distance
	for (auto &slices : NewSweep.edgePts) {
		SliceID = slices.first;
		InterpTransform = EstTransform*(NewSweep.timeStamps[SliceID] - NewSweep.tStart) / (NewSweep.tCur - NewSweep.tStart);
		for (auto &entry : slices.second) {
			OldDistance = Distance2EdgePlane(NewSweep.ptCloud[SliceID][entry], OldSweep, InterpTransform, (int)1);
			for (int col = 0; col < 6; col++) {
				InterpTransform_Delta = InterpTransform;
				InterpTransform_Delta(col) = InterpTransform(col) + Delta;
				Jacobian(cnt, col) = (Distance2EdgePlane(NewSweep.ptCloud[SliceID][entry], OldSweep, InterpTransform_Delta, (int)1) -
					OldDistance) / Delta;
			}
			DistanceVec(cnt) = OldDistance;
			cnt++;
		}
	}
    //plane distance
	for (auto &slices : NewSweep.planePts) {
		SliceID = slices.first;
		InterpTransform = EstTransform*(NewSweep.timeStamps[SliceID] - NewSweep.tStart) / (NewSweep.tCur - NewSweep.tStart);
		for (auto &entry : slices.second) {
			OldDistance = Distance2EdgePlane(NewSweep.ptCloud[SliceID][entry], OldSweep, InterpTransform, (int)2);
			for (int col = 0; col < 6; col++) {
				InterpTransform_Delta = InterpTransform;
				InterpTransform_Delta(col) = InterpTransform(col) + Delta;
				Jacobian(cnt, col) = (Distance2EdgePlane(NewSweep.ptCloud[SliceID][entry], NewSweep, InterpTransform_Delta, (int)2) -
					OldDistance) / Delta;
			}
			DistanceVec(cnt) = OldDistance;
			cnt++;
		}
	}
	return Jacobian;
}

VectorXd LMOptim::GetDistanceVec( Sweep &OldSweep, Sweep &NewSweep, VectorXd EstTransform) {
	int SliceID;
	VectorXd InterpTransform;
	VectorXd DistanceVec = VectorXd::Zero(NewSweep.numEdges + NewSweep.numPlanes);
	int cnt = 0;
	for (auto &slices : NewSweep.edgePts) {
		SliceID = slices.first;
		InterpTransform = EstTransform*(NewSweep.timeStamps[SliceID] - NewSweep.tStart) / (NewSweep.tCur - NewSweep.tStart);
		for (auto &ptID : slices.second) {
			DistanceVec(cnt) = Distance2EdgePlane(NewSweep.ptCloud[SliceID][ptID], OldSweep, InterpTransform, (int)1);
			cnt++;
		}
	}
	for (auto &slices : NewSweep.planePts) {
		SliceID = slices.first;
		InterpTransform = EstTransform*(NewSweep.timeStamps[SliceID] - NewSweep.tStart) / (NewSweep.tCur - NewSweep.tStart);
		for (auto &ptID : slices.second) {
			DistanceVec(cnt) = Distance2EdgePlane(NewSweep.ptCloud[SliceID][ptID], OldSweep, InterpTransform, (int)2);
			cnt++;
		}
	}
	return DistanceVec;
}

VectorXd LMOptim::TransformEstimate(Sweep &OldSweep, Sweep &NewSweep) {
	bool iterate1 = 1, iterate2;
	int cnt = 0;
	int MaxIter = 10; // Max iteration number
	double lambda = 1, lambda_scale = 10;
	VectorXd OldTransform, NewTransform, OldDistanceVec, NewDistanceVec;
	MatrixXd Jacobian, JTJ, JTJ_Diag;
	OldTransform << 0, 0, 0, 0, 0, 0;
	// Find feature points in NewSweep
	for (int sliceIdx = 0; sliceIdx < NewSweep.ptCloud.size(); sliceIdx++) {
		NewSweep.FindEdges(sliceIdx);
	}
	// Jacobian and distance vector
	while (iterate1) {
		Jacobian = GetJacobian(OldDistanceVec, OldSweep, NewSweep, OldTransform);
		JTJ = Jacobian.transpose()*Jacobian;
		JTJ_Diag = MatrixXd::Zero(Jacobian.rows(), Jacobian.rows());
		for (int diagIdx = 0; diagIdx < Jacobian.rows(); diagIdx++) {
			JTJ_Diag(diagIdx, diagIdx) = JTJ(diagIdx, diagIdx);
		}
		// LM steps
		iterate2 = 1;
		while (iterate2) {
			NewTransform = OldTransform - (JTJ + lambda*JTJ_Diag).inverse()*Jacobian.transpose()*OldDistanceVec;
			// Compare
			NewDistanceVec = GetDistanceVec(OldSweep, NewSweep, NewTransform);
			OldDistanceVec = GetDistanceVec(OldSweep, NewSweep, OldTransform);
			if (NewDistanceVec.transpose()*NewDistanceVec < OldDistanceVec.transpose()*OldDistanceVec) {
				// improved, take new transform, tune lambda
				OldTransform = NewTransform;
				lambda = lambda / lambda_scale;
				cnt++;
				iterate2 = 0;
			}
			else {
				// not improved, keep old transform, tune lambda
				lambda = lambda*lambda_scale;
				cnt++;
			}
			if (0) {
				// convergence check
				cout << "Transformation estimation converged!" << endl;
				iterate1 = 0;
				iterate2 = 0;
			}
			if (cnt == MaxIter) {
				// max iteration check
				cout << "Max iteration reached! Optimality not guaranteed!" << endl;
				iterate1 = 0;
				iterate2 = 0;
			}
		}
	}
	return NewTransform;
}
