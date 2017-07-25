#include "LMOptim.h"

LMOptim::LMOptim() {

}

LMOptim::~LMOptim() {

}

double LMOptim::Distance2EdgePlane(LoamPt &pt, Sweep &OldSweep, VectorXd EstTransform, int EnPflag) {
	// EnPflag = 1: Edge | EnPflag = 2: Plane
	Sweep TEMP;
	Vector3d xi = pt.xyz;
	Vector3d xj = OldSweep.ptCloud[pt.nearPt1[0]][pt.nearPt1[1]].xyz;
	Vector3d xl = OldSweep.ptCloud[pt.nearPt2[0]][pt.nearPt2[1]].xyz;
	Vector3d T_trans(3), T_rot(3), omega(3), xi_hat;
	Matrix3d eye3 = Matrix3d::Identity(), omega_hat, R;
	double theta, Distance;
	double Tmax = abs(EstTransform(0)), Rmax = abs(EstTransform(3));
	int Idx;
	// check translation
	for (Idx = 1; Idx < 3; Idx++) {
		if (abs(EstTransform(Idx)) > Tmax) {
            Tmax = abs(EstTransform(Idx));
		}
	}
	if (Tmax < pow(10, -3)) {
		T_trans << 0, 0, 0; // no translation
	}
	else {
		T_trans << EstTransform(0), EstTransform(1), EstTransform(2);
	}
	// check rotation
	for (Idx = 4; Idx < 6; Idx++) {
		if (abs(EstTransform(Idx)) > Rmax) {
			Rmax = abs(EstTransform(Idx));
		}
	}
	if (Rmax < pow(10, -3)) {
		R = Matrix3d::Identity(); // no rotation
	}
	else {
		T_rot << EstTransform(3), EstTransform(4), EstTransform(5);
		theta = T_rot.norm();
		omega << EstTransform(3) / theta, EstTransform(4) / theta, EstTransform(5) / theta;
		omega_hat << 0, -omega(2), omega(1),
			omega(2), 0, -omega(0),
			-omega(1), omega(0), 0;
		R = eye3 + omega_hat*sin(theta) + omega_hat*omega_hat*(1 - cos(theta));
	}
    // 
	xi_hat = R.inverse()*(xi - T_trans);
	if (EnPflag == 1) {
		//Vector3d a = xi_hat - xj;
		//Vector3d b = xj - xl;
		//Vector3d c = a.cross(b);
		//double d = c.norm();
		//double e = b.norm();
		//Distance = d / e;	
		//Distance = ((xi_hat - xj).cross(xi_hat - xl)).norm() / (xj - xl).norm();
		Distance = TEMP.Dist2Line(xi_hat, xj, xl);
	}
	else if (EnPflag == 2) {
		Vector3d xm = OldSweep.ptCloud[pt.nearPt3[0]][pt.nearPt3[1]].xyz;
		//Distance = abs((xi_hat - xj).dot((xj - xl).cross(xj - xm))) / ((xj - xl).cross(xj - xm)).norm();
		Distance = TEMP.Dist2Plane(xi_hat, xj, xl, xm);
	}
	else {
		cout << "Edge or Plane indicator not provided!" << endl;
	}

	return Distance;
}

MatrixXd LMOptim::GetJacobian(VectorXd DistanceVec,MatrixXd &W, Sweep &OldSweep, Sweep &NewSweep, VectorXd EstTransform) {
	MatrixXd Jacobian = MatrixXd::Zero(NewSweep.numEdges + NewSweep.numPlanes, 6);
	VectorXd InterpTransform_Delta, InterpTransform;
	DistanceVec = VectorXd::Zero(NewSweep.numEdges + NewSweep.numPlanes);
	W = MatrixXd::Zero(NewSweep.numEdges + NewSweep.numPlanes,  NewSweep.numEdges + NewSweep.numPlanes);
	double Delta = pow(10, -6); //numerical Jacobian step-size
	double BiSq_threshold = 10; //weighting func: w(x) = (1-(x/BiSq_threshold)^2)^2
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
			if (OldDistance < BiSq_threshold) {
				W(cnt, cnt) = pow((1 - pow(OldDistance / BiSq_threshold, 2)), 2);
			}
			cout <<"Jacobian row: "<< cnt << endl;
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
				Jacobian(cnt, col) = (Distance2EdgePlane(NewSweep.ptCloud[SliceID][entry], OldSweep, InterpTransform_Delta, (int)2) -
					OldDistance) / Delta;
			}
			DistanceVec(cnt) = OldDistance;
			if (OldDistance < BiSq_threshold) {
				W(cnt, cnt) = pow((1 - pow(OldDistance / BiSq_threshold, 2)), 2);
			}
			cout <<"Jacobian row: "<< cnt << endl;
			cnt++;
		}
	}
	
	return Jacobian;
}

VectorXd LMOptim::GetDistanceVec(Sweep &OldSweep, Sweep &NewSweep, VectorXd EstTransform) {
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
	bool iterate1 = 1, iterate2, converged;
	int cnt = 0;
	int n = 3;
	int MaxIter = 10; // Max iteration number
	double lambda = 1, lambda_scale = 10; // Lambda for LM
	double convergence_threshold_residual = 10, convergence_threshold_diffs = 0.1; // threshold for LM convergence
	double sum_diff = 0;
	VectorXd prev_diffs = convergence_threshold_diffs*VectorXd::Zero(n);
	VectorXd OldTransform(6), NewTransform, OldDistanceVec, NewDistanceVec;
	MatrixXd Jacobian, JTWJ, JTWJ_Diag, W;
	OldTransform << 0, 0, 0, 0, 0, 0;
	// Find feature points in NewSweep
	cout << "Start looking for correspondences" << endl;
	NewSweep.FindAllEdges();
	OldSweep.FindAllEdges();
	for (int sliceIdx = 0; sliceIdx < NewSweep.ptCloud.size(); sliceIdx++) {
		NewSweep.FindCorrespondences(sliceIdx, OldSweep);		
	}
	// Iteration
	cout << "Start LM iteration" << endl;
	while (iterate1) {
		Jacobian = GetJacobian(OldDistanceVec, W, OldSweep, NewSweep, OldTransform);
		JTWJ = Jacobian.transpose()*W*Jacobian;
		JTWJ_Diag = MatrixXd::Zero(6, 6);
		for (int diagIdx = 0; diagIdx < 6; diagIdx++) {
			JTWJ_Diag(diagIdx, diagIdx) = JTWJ(diagIdx, diagIdx);
		}
		// LM steps
		iterate2 = 1;
		while (iterate2) {
			NewTransform = OldTransform - (JTWJ + lambda*JTWJ_Diag).inverse()*Jacobian.transpose()*W*OldDistanceVec;
			// Compare
			NewDistanceVec = GetDistanceVec(OldSweep, NewSweep, NewTransform);
			if (NewDistanceVec.norm() < OldDistanceVec.norm()) {
				// improved: record n previous transform diffs, take new transform, tune lambda
				for (int i = 0; i < n - 1; i++) {
					prev_diffs(i) = prev_diffs(i + 1);
				}
				prev_diffs(n - 1) = (NewDistanceVec - OldDistanceVec).norm();
				lambda = lambda / lambda_scale;
				iterate2 = 0;
				OldTransform = NewTransform;
			}
			else {
				// not improved, keep old transform, tune lambda
				lambda = lambda*lambda_scale;
			}
			cout << "LM iteration " << cnt << endl;
			cnt++;
			// convergence check: 1. Residual (d->0); 2. Variance of T in previous n iteration
			if (NewDistanceVec.norm() < convergence_threshold_residual) {
				converged = 1;
			}
			else if (prev_diffs.sum() < convergence_threshold_diffs) {
				converged = 1;
			}
			if (converged) {
				iterate1 = 0;
				iterate2 = 0;
			}
			// max iteration check
			if (cnt == MaxIter) {
				cout << "Max iteration reached! Optimality not guaranteed!" << endl;
				iterate1 = 0;
				iterate2 = 0;
			}
		}
	}

	return NewTransform;
}
