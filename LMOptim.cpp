#include "LMOptim.h"

LMOptim::LMOptim() {

}

LMOptim::~LMOptim() {

}

double LMOptim::Distance2EdgePlane(LoamPt &pt, Sweep &OldSweep, VectorXd &EstTransform, int EnPflag) {
	// EnPflag = 1: Edge | EnPflag = 2: Plane
	Vector3d xi = pt.xyz;
	Vector3d xj = OldSweep.ptCloud[pt.nearPt1[0]][pt.nearPt1[1]].xyz;
	Vector3d xl = OldSweep.ptCloud[pt.nearPt2[0]][pt.nearPt2[1]].xyz;
	Vector3d xm;
	Vector3d T_trans, T_rot, omega, xi_hat;
	Matrix3d eye3 = Matrix3d::Identity(), omega_hat, R;
	double theta, d;
	int Idx;
	T_trans << EstTransform(0), EstTransform(1), EstTransform(2);
	T_rot << EstTransform(3), EstTransform(4), EstTransform(5);
	theta = T_rot.norm();
	// check rotation
	if (theta < 10e-6) {
		R = Matrix3d::Identity(); // no rotation
	}
	else {
		omega = T_rot/theta;
		omega_hat << 0, -omega(2), omega(1),
			omega(2), 0, -omega(0),
			-omega(1), omega(0), 0;
		R = eye3 + omega_hat*sin(theta) + omega_hat*omega_hat*(1 - cos(theta));
	}
	xi_hat = R.inverse()*(xi - T_trans);
	// Distance to edge
	if (EnPflag == 1) {
		if ((xj - xl).norm() < 10e-6) {
			d = (xi_hat - xj).norm();
		}
		else {
			d = ((xi_hat - xj).cross(xi_hat - xl)).norm() / (xj - xl).norm();
		}
	}
	// Distance to plane (Needs exception check)
	else if (EnPflag == 2) {





		xm = OldSweep.ptCloud[pt.nearPt3[0]][pt.nearPt3[1]].xyz;
		d = abs((xi_hat - xj).dot((xj - xl).cross(xj - xm))) / ((xj - xl).cross(xj - xm)).norm();
	}
	//else {
	//	//cout << "Edge or Plane indicator not provided!" << endl;
	//}
	return d;
}

MatrixXd LMOptim::GetJacobian(VectorXd &DistanceVec,MatrixXd &W, Sweep &OldSweep, Sweep &NewSweep, VectorXd EstTransform) {
	//int numPts = NewSweep.numEdges + NewSweep.numPlanes;
	int numPts = NewSweep.numEdges + NewSweep.numPlanes;
	VectorXd InterpTransform_Delta, InterpTransform, distanceVecDel, EstTransformDel;
	MatrixXd Jacobian = MatrixXd::Zero(numPts, 6);
	DistanceVec = VectorXd::Zero(numPts);
	W = MatrixXd::Zero(numPts, numPts);
	int cnt = 0, dim, i;
	DistanceVec = GetDistanceVec(OldSweep, NewSweep, EstTransform);
	for (dim = 0; dim < 6; dim++) {
		EstTransformDel = EstTransform;
		EstTransformDel(dim) = EstTransformDel(dim) + delta;
		distanceVecDel = GetDistanceVec(OldSweep, NewSweep, EstTransformDel);
		for (i = 0; i < numPts; i++) {
			Jacobian(i, dim) = (distanceVecDel(i) - DistanceVec(i)) / delta;
			if (dim == 0) {
				// Set weights
				if (DistanceVec(i) >= BiSq_threshold) {
					W(i, i) = 0;
				}
				else {
					W(i, i) = pow((1 - pow(DistanceVec(i) / BiSq_threshold, 2)), 2);
				}
			}
		}
	}
	return Jacobian;
}

VectorXd LMOptim::GetDistanceVec(Sweep &OldSweep, Sweep &NewSweep, VectorXd EstTransform) {
	
	int numPts = NewSweep.numEdges + NewSweep.numPlanes;
	//int numPts = NewSweep.numEdges;
	VectorXd InterpTransform;
	VectorXd DistanceVec = VectorXd::Zero(numPts);
	int cnt = 0;
	double tRatio;
	//Edge points
	for (auto &slice : NewSweep.edgePts)
	{
		InterpTransform = (NewSweep.timeStamps[slice.first] - NewSweep.tStart) / (NewSweep.tCur - NewSweep.tStart)*EstTransform;
		for (auto &idx : slice.second)
		{
			DistanceVec(cnt) = Distance2EdgePlane(NewSweep.ptCloud[slice.first][idx], OldSweep, InterpTransform, (int)1);
			cnt++;
		}
	}

	for (auto &slice : NewSweep.planePts)
	{
		tRatio = (NewSweep.timeStamps[slice.first] + 1 - NewSweep.tStart) / (NewSweep.tCur - NewSweep.tStart);
		for (auto &idx : slice.second)
		{
			auto &x = NewSweep.ptCloud[slice.first][idx];
			DistanceVec(cnt) = NewSweep.Dist2Plane(BackTransform(x.xyz, EstTransform, tRatio), OldSweep.ptCloud[x.nearPt1[0]][x.nearPt1[1]].xyz, OldSweep.ptCloud[x.nearPt2[0]][x.nearPt2[1]].xyz, OldSweep.ptCloud[x.nearPt3[0]][x.nearPt3[1]].xyz);
			cnt++;
		}
	}
	return DistanceVec;
}

VectorXd LMOptim::TransformEstimate(Sweep &OldSweep, Sweep &NewSweep) {
	bool iterate1 = 1, iterate2, converged = 0;
	int cnt = 0, diagIdx;
	int numPts = NewSweep.numEdges + NewSweep.numPlanes;
	double lambda = lambda_ini; // Initial lambda
	double convergence_threshold_residual = 0.002*numPts; // Convergence criteria: Levenberg-Marquart
	VectorXd oldTransform(6), newTransform, oldDistanceVec, newDistanceVec;
	MatrixXd Jacobian, JTWJ, JTWJ_Diag, W;
	oldTransform << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
	// Find feature points in Sweep
	cout << "Start looking for correspondences" << endl;
	NewSweep.FindAllEdges();
	OldSweep.FindAllEdges();
	for (int sliceIdx = 0; sliceIdx < NewSweep.ptCloud.size(); sliceIdx++) {
		NewSweep.FindCorrespondences(sliceIdx, OldSweep);		
	}
	// Iteration
	cout << "Correspondences found. Start LM iteration" << endl;
	while (iterate1) {
		Jacobian = GetJacobian(oldDistanceVec, W, OldSweep, NewSweep, oldTransform);
		JTWJ = Jacobian.transpose()*W*Jacobian;
		/*cout << "J: " << endl << Jacobian << endl;
		cout << "JTWJ: " << endl << JTWJ << endl;*/
		JTWJ_Diag = MatrixXd::Zero(6, 6);
		for (diagIdx = 0; diagIdx < 6; diagIdx++) {
			JTWJ_Diag(diagIdx, diagIdx) = JTWJ(diagIdx, diagIdx);
		}
		iterate2 = 1;
		while (iterate2) {
			newTransform = oldTransform - (JTWJ + lambda*JTWJ_Diag).inverse()*Jacobian.transpose()*W*oldDistanceVec;
			newDistanceVec = GetDistanceVec(OldSweep, NewSweep, newTransform);
			if (newDistanceVec.norm() < oldDistanceVec.norm()) {
				// improved, keep new transform, tune down lambda, break iter2
				lambda = lambda / lambda_scale;
				if (lambda >= lambda_max) {
					lambda = lambda_max;
				}
				oldTransform = newTransform;
				iterate2 = 0;
			}
			else {
				// not improved, keep old transform, tune up lambda
				lambda = lambda*lambda_scale;
			}
			cout << "LM iteration " << cnt << endl;
			cout << "Old T " << oldTransform << endl;
			cout << "New T " << newTransform << endl;
			cout << "Distance norm " << newDistanceVec.norm() << endl;
			cout << "Lambda " << lambda << endl;
			cnt++;
			cout << cnt << endl;
			// convergence check (residual)
			if (newDistanceVec.norm() < convergence_threshold_residual) {
				cout << "Residual converged!" << endl;
				iterate1 = 0;
				iterate2 = 0;
			}
			// max iteration check
			if (cnt == maxIterLM) {
				cout << "Max iteration reached! Optimality not guaranteed!" << endl;
				iterate1 = 0;
				iterate2 = 0;
			}
		}
	}
	cout << "Estimated T: " << newTransform << endl;
	cout << "Number of iteration " << cnt << endl;
	cout << "Distance residual " << newDistanceVec.norm() << endl;
	getchar();
	return newTransform;
}
