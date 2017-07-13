#include "Sweep.h"

Sweep::Sweep()
{

}

Sweep::~Sweep()
{

}

Sweep::Sweep(std::vector<std::vector<double>> &inputSlice)
{
	sweepID = NULL;
	AddSlice(sweepID, numSlices + 1, (double)NULL, inputSlice);
}

void Sweep::AddSlice(int sweepNumber, int sliceNumber, double timeStamp, std::vector<std::vector<double>> &inputSlice)
{
	std::vector<LoamPt> slice;
	LoamPt tempPt;
	numSlices++;
	for (auto &xyzPt : inputSlice)
	{
		slice.push_back(LoamPt(xyzPt, sweepNumber, sliceNumber, timeStamp));
	}
	ptCloud.push_back(slice);
}

void Sweep::AddSlice(std::vector<std::vector<double>> &inputSlice)
{
	std::vector<LoamPt> slice;
	for (auto &xyzPt : inputSlice)
	{
		slice.push_back(LoamPt(xyzPt));
	}
	ptCloud.push_back(slice);
	numSlices++;
}

void Sweep::FindEdges(int sliceIdx)
{
	auto &slice = ptCloud[sliceIdx];
	std::vector<std::vector<double>> curvatures; // {curvatureValue, ptIndex}
	int numPts = slice.size();
	curvatures.resize(numPts);
	int firstPt = kernalSize / 2;
	int lastPt = numPts - kernalSize / 2;
	int ptPerRegion = (lastPt - firstPt) / 4;
	Vector3d distVec = { 0,0,0 };

	distVec = kernalSize*slice[firstPt].xyz;

	// find distance vector for initial point
	for (int i = -(kernalSize / 2); i < (kernalSize / 2) + 1; i++)
	{
		distVec -= slice[firstPt + i].xyz;
	}

	curvatures[firstPt] = { distVec.norm() / slice[firstPt].xyz.norm(), (double)firstPt };

	// calculate curvature values for all points after the first point
	for (int i = firstPt + 1; i < lastPt; i++)
	{
		distVec = distVec - kernalSize*(slice[i - 1].xyz - slice[i].xyz) - slice[i + kernalSize / 2].xyz + slice[i - kernalSize / 2 - 1].xyz;
		curvatures[i] = { distVec.norm() / slice[i].xyz.norm(), (double)i };
	}

	// find and record edges/plains via curvature sorting
	for (int reg = 0; reg < regionPerSlice; reg++)
	{
		SortCurvatures(sliceIdx, curvatures, firstPt + ptPerRegion*reg, firstPt + ptPerRegion*(reg + 1));
	}
}

void Sweep::SortCurvatures(int sliceIdx, std::vector<std::vector<double>> &curveVec, int startIdx, int endIdx)
{
	std::vector<std::vector<double>> curvatures;
	curvatures.resize(endIdx - startIdx);

	// get subset of curvature values corresponding to region between startIdx, endIdx
	for (int i = 0; i < endIdx - startIdx; i++)
	{
		curvatures[i] = curveVec[startIdx + i];
	}

	// sort the curvature vector
	MergeSort(curvatures);

	int edges = 0, planes = 0;

	int planeTurn = false;

	while (curvatures.size() > 0 && (edges < edgePerRegion || planes < planePerRegion))
	{
		// try to find edge point
		if (planeTurn == false && edges < edgePerRegion)
		{
			if (FindBestEdgePt(sliceIdx, curvatures) == true)
			{
				edges++;
			}
		}
		// try to find plane point
		else if (planes < planePerRegion)
		{
			if (FindBestPlanePt(sliceIdx, curvatures) == true)
			{
				planes++;
			}
		}
		planeTurn = !planeTurn;
	}
}

bool Sweep::FindBestEdgePt(int sliceIdx, std::vector<std::vector<double>> &curveVec)
{
	// best edges have highest curvature values
	std::vector<double> pt;
	pt = curveVec[curveVec.size() - 1];
	curveVec.erase(curveVec.end() - 1);
	if (EvaluateEdge(sliceIdx, pt) == true) // valid edge point
	{
		// save edgePt's {sliceIdx, ptIdx}
		edgePts[sliceIdx].push_back((int)pt[1]);
		return true;
	}
	return false;
}

bool Sweep::FindBestPlanePt(int sliceIdx, std::vector<std::vector<double>> &curveVec)
{
	// best planes have low curvature values
	std::vector<double> pt;
	pt = curveVec[0];
	curveVec.erase(curveVec.begin());
	if (EvaluatePlane(sliceIdx, pt) == true) // valid plane point
	{
		// save planePt's {sliceIdx, ptIdx}
		planePts[sliceIdx].push_back((int)pt[1]);
		return true;
	}
	return false;
}

bool Sweep::EvaluateEdge(int sliceIdx, std::vector<double> &potentialPt) // Checks to make sure that no nearby points are drastically closer to the sensor
{
	double ptDist = ptCloud[sliceIdx][potentialPt[1]].xyz[0]; // treating ptDist[0] = X, the depth-distance w.r.t. sensor
	for (int i = -(kernalSize / 2); i < 0; i++)
	{
		if ((ptCloud[sliceIdx][potentialPt[1] + i + 1].xyz[0] - ptCloud[sliceIdx][potentialPt[1] + i].xyz[0]) > edgeFindThreshold) // large increase in distance as we approach the potential point == occlusion
		{
			return false;
		}
	}
	for (int i = 1; i < (kernalSize / 2) + 1; i++)
	{
		if ((ptCloud[sliceIdx][potentialPt[1] + i + 1].xyz[0] - ptCloud[sliceIdx][potentialPt[1] + i].xyz[0]) < -edgeFindThreshold) // large decrease in distance as we move away from the potential point == occlusion
		{
			return false;
		}
	}
	// no nearby pts indicate occlusion, pt is valid!
	return true;
}

bool Sweep::EvaluatePlane(int sliceIdx, std::vector<double> &potentialPt) // Checks to make sure that no nearby points are drastically closer to the sensor
{
	Vector3d distVec;
	Vector3d Xi = ptCloud[sliceIdx][potentialPt[1]].xyz;
	distVec = kernalSize*Xi;
	for (int i = -kernalSize / 2; i < kernalSize / 2; i++)
	{
		distVec -= ptCloud[sliceIdx][potentialPt[1] + i].xyz;
	}

	if (distVec.dot(Xi) / (distVec.norm()*Xi.norm()) >= 0.5) // a*b = |a||b|cos(theta) ---> dot(a,b)/(|a||b|) = cos(theta)
	{
		return true;
	}
	else
	{
		return false;
	}
}

double Sweep::Distance2(LoamPt &pt, Sweep &OldSweep, VectorXd EstTransform, int &EnPflag) {
	// EnPflag = 1: Edge | EnPflag = 2: Plane
	Vector3d xi = pt.xyz;
	Vector3d xj = OldSweep.ptCloud[pt.nearPt1[0]][pt.nearPt1[1]].xyz;
	Vector3d xl = OldSweep.ptCloud[pt.nearPt2[0]][pt.nearPt2[1]].xyz;
	Vector3d T_trans, T_rot, omega, xi_hat;
	Matrix3d eye3, omega_hat, R;
	double Distance;

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
		double theta = T_rot.norm();
		T_trans << EstTransform(0), EstTransform(1), EstTransform(2);
		omega << EstTransform(3) / theta, EstTransform(4) / theta, EstTransform(5) / theta;
		eye3 << 1, 0, 0, 0, 1, 0, 0, 0, 1;
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

MatrixXd Sweep::GetJacobian(VectorXd DistanceVectorEig, Sweep &OldSweep, Sweep &NewSweep,
	VectorXd EstTransform) {

	vector<vector<double>> Jacobian_Full;
	vector<double> JacobianRow(6), DistanceVector;
	VectorXd EstTransform_Delta;
	double Delta = pow(10, -6), OldDistance;
	int EnPflag = 1; //edge distance
	for (int i = 0; i < NewSweep.edgePts.size(); i++) {
		for (auto &entry : edgePts[i]) {
			OldDistance = Distance2(NewSweep.ptCloud[i][entry], OldSweep, EstTransform, EnPflag);
			for (int col = 0; col < 6; col++) {
				EstTransform_Delta = EstTransform;
				EstTransform_Delta(col) = EstTransform(col) + Delta;
				JacobianRow[col] = (Distance2(NewSweep.ptCloud[i][entry], OldSweep, EstTransform_Delta, EnPflag) -
					OldDistance) / Delta;
			}
			Jacobian_Full.push_back(JacobianRow);
			DistanceVector.push_back(OldDistance);
		}
	}
	EnPflag = 2; //plane distance
	for (int i = 0; i < NewSweep.planePts.size(); i++) {
		for (auto &entry : planePts[i]) {
			Distance2(NewSweep.ptCloud[i][entry], OldSweep, EstTransform, EnPflag);
			for (int col = 0; col < 6; col++) {
				EstTransform_Delta = EstTransform;
				EstTransform_Delta(col) = EstTransform(col) + Delta;
				JacobianRow[col] = (Distance2(NewSweep.ptCloud[i][entry], OldSweep, EstTransform_Delta, EnPflag) -
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

VectorXd LMoptimization(Sweep &OldSweep, Sweep &NewSweep) {
	double lambda = 1;
	double lambda_scale = 10;
	VectorXd TransformInitial;


	return TransformInitial;


}

void Sweep::FindCorrespondences(int sliceNumber, Sweep &OldSweep)
{
	// Find all edge-point correspondences
	for (auto &elem : edgePts[sliceNumber])
	{
		FindNearestLine(ptCloud[sliceNumber][elem], OldSweep);
	}

	// Find all plane-point correspondences
	for (auto &elem : planePts[sliceNumber])
	{
		FindNearestPlane(ptCloud[sliceNumber][elem], OldSweep);
	}
}

void Sweep::FindNearestLine(LoamPt &curEdgePt, Sweep &OldSweep)
{
	// Back-transform the point to the beginning of the current sweep.
	Vector3d Xi_ = BackTransform(curEdgePt.xyz, transform, (timeStamps[curEdgePt.sliceID] - tStart)/(tCur - tStart));

	std::vector<int> bestPt = { 0,0 };
	double bestDist = 10e10, tempDist;
	std::vector<double> distances;

	// Find the nearest edgepoint located in the +/- n-neighboring slices of the previous sweep
	for (int i = curEdgePt.sliceID - 2 ; curEdgePt.sliceID + i < 3; i++)
	{
		for (auto &oldIdx : OldSweep.edgePts[i%maxNumSweeps])
		{
			tempDist = (OldSweep.ptCloud[i%maxNumSweeps][oldIdx].xyz - Xi_).norm();
			if ((tempDist < bestDist))
			{
				bestDist = tempDist;
				bestPt = { i%maxNumSweeps, oldIdx };
			}
		}
	}
	curEdgePt.nearPt1 = bestPt;
	bestDist = 10e10;

	// Find edgepoint closest to nearPt1 located in either adjacent slice of the previous sweep
	for (int i = curEdgePt.nearPt1[0] - 1; i < curEdgePt.nearPt1[0] + 2; i += 2)
	{
		for (auto &oldIdx : OldSweep.edgePts[i%maxNumSweeps])
		{
			tempDist = (OldSweep.ptCloud[i%maxNumSweeps][oldIdx].xyz - OldSweep.ptCloud[curEdgePt.nearPt1[0]][curEdgePt.nearPt1[1]].xyz).norm();
			if (tempDist < bestDist)
			{
				bestDist = tempDist;
				bestPt = { i%maxNumSweeps, oldIdx };
			}
		}
	}
	curEdgePt.nearPt2 = bestPt;
	int x;
	curEdgePt.dist = Distance2(curEdgePt, OldSweep, transform, x); // this needs to be changed
}

void Sweep::FindNearestPlane(LoamPt &curEdgePt, Sweep &OldSweep)
{
	// Back-transform the point to the beginning of the current sweep.
	Vector3d Xi_ = BackTransform(curEdgePt.xyz, transform, (timeStamps[curEdgePt.sliceID] - tStart) / (tCur - tStart));
}