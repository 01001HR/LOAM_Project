#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <math.h>
#include <Eigen\Core>
#include <random>
#include <fstream>
#include <iterator>
#include <map>
#include "LinearAlgebraFns.h"

using namespace Eigen;
using namespace std;

class LoamPt
{
public:
	//std::vector<double> xyz;
	Vector3d xyz;
	int sweepID;
	int sliceID;
	std::vector<int> nearPt1, nearPt2, nearPt3; // these will be len 2 vectors with pt = {sliceIndex, pointIndex}
	int filled;
	double timeStamp;

	LoamPt();
	~LoamPt();
	LoamPt(const std::vector<double> &xyzInput); // no time given
	LoamPt(const std::vector<double> &xyzInput, const double time);
	LoamPt(const double x, const double y, const double z, const double time);
	LoamPt(const LoamPt &otherPt);
	LoamPt &operator = (const LoamPt &otherPt);
	//std::vector<double> Minus(const LoamPt &otherPt);
	//std::vector<double> Minus(const std::vector<double> &otherPt);

	//static inline int Size(const std::vector<double> &xyz);

	bool SetXYZ(const std::vector<double> &newXYZ);

	inline const double GetX();
	inline const double GetY();
	inline const double GetZ();
	inline const double GetTime();

	bool Distance(double &dist, LoamPt &otherPt);
	Vector3d Transform(Matrix4d &xformMatrix4x4);
	//Vector3d Transform(const std::vector<double[4]> &xformMatrix4x4); // xformMatrix is stored as stacked row-vectors
	void TransformSelf(Matrix4d &xformMatrix4x4); // xformMatrix is stored as stacked row-vectors
};

LoamPt::LoamPt()
{
	filled = 0;
}

LoamPt::~LoamPt()
{
}

LoamPt::LoamPt(const LoamPt &otherPt) // copy constructor
{
	if (otherPt.filled == 1)
	{
		xyz = otherPt.xyz;
		timeStamp = otherPt.timeStamp;
		filled = 1;
	}
	else
	{
		std::cout << "You are either trying to self assign, or the incoming pt has not been initialized/filled" << std::endl;
	}
}

LoamPt &LoamPt::operator = (const LoamPt &otherPt) // copy operator
{
	if (this != &otherPt && otherPt.filled == 1) // if the incoming instance is a different pt and has values
	{
		xyz = otherPt.xyz;
		timeStamp = otherPt.timeStamp;
		filled = 1;
	}
	else
	{
		std::cout << "You are either trying to self assign, or the incoming pt has not been initialized/filled" << std::endl;
	}
	return *this;
}

LoamPt::LoamPt(const std::vector<double> &xyzInput)
{
	if (SetXYZ(xyzInput) == true)
	{
		timeStamp = NULL;
		filled = 1;
		std::cout << "Warning, the incoming datapoints have no timestamps" << std::endl;
	}
}

LoamPt::LoamPt(const std::vector<double> &xyzInput, double time) // vector input constructor
{
	std::cout << "In here" << std::endl;
	if (SetXYZ(xyzInput) == true)
	{
		timeStamp = time;
		filled = 1;
	}
}

LoamPt::LoamPt(double x, double y, double z, double time) // individual value constructor
{
	xyz = { x,y,z };
	timeStamp = time;
	filled = 1;
}

bool LoamPt::Distance(double &dist, LoamPt &otherPt)
{
	if ((filled == 1) && (otherPt.filled == 1) && (this != &otherPt)) // if both points are filled
	{
		// calculate pt2pt distance, return success
		dist = (xyz-otherPt.xyz).norm();
		return true;
	}
	else
	{
		// return failure
		return false;
	}
}

bool LoamPt::SetXYZ(const std::vector<double> &newXYZ)
{
	//int size = Size(newXYZ);
	//std::cout << size << std::endl;
	if (newXYZ.size() == 3)
	{
		// point is set to new value and declared filled
		xyz = { newXYZ[0], newXYZ[1], newXYZ[2] };
		filled = 1;
		return true;
	}
	else
	{
		// point stays the same, filled declaration remains what it was before fill attempt
		std::cout << "The incoming xyz point has " << newXYZ.size() << " values rather than 3" << std::endl;
		return false;
	}
}

inline const double LoamPt::GetX()
{
	if (filled == 1) return xyz[0];
}

inline const double LoamPt::GetY()
{
	if (filled == 1) return xyz[1];
}

inline const double LoamPt::GetZ()
{
	if (filled == 1) return xyz[2];
}

inline const double LoamPt::GetTime()
{
	if (filled == 1) return timeStamp;
}

Vector3d LoamPt::Transform(Matrix4d &xformMatrix4x4)
{
	Vector4d augVec = { xyz[0], xyz[1], xyz[2], 1 }, newVec;
	newVec = xformMatrix4x4*augVec;
	return { newVec[0], newVec[1], newVec[2] };
}


void LoamPt::TransformSelf(Matrix4d &xformMatrix4x4)
{
	// assumes the 4x4 matrix is filled/correct
	Vector4d augVec = { xyz[0], xyz[1], xyz[2], 1 };
	if (filled == 1)
	{
		augVec = xformMatrix4x4*augVec;
		xyz = { augVec[0], augVec[1], augVec[2] };
	}
	else
	{
		std::cout << "Tried transforming an unfilled point" << std::endl;
	}
}


class Sweep
{
public:
	std::vector<std::vector<LoamPt>> ptCloud;			// Cloud of LoamPts where ptCloud[i] = slice of LoamPts, and ptCloud[i][j] = individual LoamPt
	std::map<int, std::vector<int>> edgePts, planePts;	// 2D Vector containing the slice/pt indices of pts declared as edge/plane points
	std::vector<double> timeStamps;						// Vector of time values, where timeStamps[i] is the timeStamp corresponding to slice i of the ptCloud
	double tStart, tEnd;
	VectorXd transform;
	int sweepID, numSlices = -1, kernalSize = 11, regionPerSlice = 4, edgePerRegion = 2, planePerRegion = 4, edgeFindThreshold = 3;
	Sweep();
	~Sweep();
	Sweep(std::vector<std::vector<double>> &inputSlice);
	void AddSlice(int sweepNumber, int sliceNumber, std::vector<std::vector<double>> &inputSlice);
	void AddSlice(std::vector<std::vector<double>> &inputSlice);
	void FindAllEdges();
	void FindEdges(int sliceIndex);
	void SortCurvatures(int sliceIdx, std::vector<std::vector<double>> &curveVec, int startIdx, int endIdx);
	bool FindBestEdgePt(int sliceIdx, std::vector<std::vector<double>> &curveVec);
	bool EvaluateEdge(int sliceIdx, std::vector<double> &potentialPt);
	bool FindBestPlanePt(int sliceIdx, std::vector<std::vector<double>> &curveVec);
	bool EvaluatePlane(int sliceIdx, std::vector<double> &potentialPt);
	double Distance2(LoamPt &pt, Sweep &OldSweep, VectorXd EstTransform, int &EnPflag);
	MatrixXd GetJacobian(Sweep &OldSweep);
};

Sweep::Sweep()
{

}

Sweep::~Sweep()
{

}

Sweep::Sweep(std::vector<std::vector<double>> &inputSlice)
{
	sweepID = NULL;
	AddSlice(sweepID, numSlices+1, inputSlice);
}

void Sweep::AddSlice(int sweepNumber, int sliceNumber, std::vector<std::vector<double>> &inputSlice)
{
	std::vector<LoamPt> slice;
	for (auto &xyzPt : inputSlice)
	{
		slice.push_back(LoamPt(xyzPt));
	}
	ptCloud.push_back(slice);
	numSlices++;
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
		distVec = distVec - kernalSize*(slice[i-1].xyz - slice[i].xyz) - slice[i + kernalSize / 2].xyz + slice[i - kernalSize / 2 - 1].xyz;
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
	for (int i = 0; i < endIdx-startIdx; i++)
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
	pt = curveVec[curveVec.size()-1];
	curveVec.erase(curveVec.end()-1);
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
	for (int i = -(kernalSize/2); i < 0; i++)
	{
		if ((ptCloud[sliceIdx][potentialPt[1] + i + 1].xyz[0] - ptCloud[sliceIdx][potentialPt[1]+i].xyz[0]) > edgeFindThreshold) // large increase in distance as we approach the potential point == occlusion
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

//double Sweep::Distance2(LoamPt &pt, Sweep &OldSweep, VectorXd EstTransform, int &EnPflag) {
//	// EnPflag = 1: Edge | EnPflag = 2: Plane
//	Vector3d xi = pt.xyz;
//	Vector3d xj = OldSweep.ptCloud[pt.nearPt1[0]][pt.nearPt1[1]].xyz;
//	Vector3d xl = OldSweep.ptCloud[pt.nearPt2[0]][pt.nearPt2[1]].xyz;
//	Vector3d T_trans, T_rot, omega, xi_hat;
//	//VectorXd EstTransform = OldSweep.transform;
//	Matrix3d eye3, omega_hat, R;
//	double theta = T_rot.norm(), Distance;
//
//	T_trans << EstTransform(0), EstTransform(1), EstTransform(2);
//	T_rot << EstTransform(3), EstTransform(4), EstTransform(5);
//	omega << EstTransform(3)/theta, EstTransform(4)/theta, EstTransform(5)/theta;
//	eye3 << 1, 0, 0, 0, 1, 0, 0, 0, 1;
//	omega_hat << 0, -omega(2), omega(1),
//		         omega(2), 0, -omega(0),
//	         	-omega(1), omega(0), 0;
//	R = eye3 + omega_hat*sin(theta) + omega_hat*omega_hat*(1 - cos(theta));
//	R.transposeInPlace();
//	xi_hat = R*(xi - T_trans);
//
//	if (EnPflag == 1) {
//		Distance = ((xi_hat - xj).cross(xi_hat - xl)).norm() / (xj - xl).norm();
//	}
//	else if (EnPflag == 2) {
//		Vector3d xm = OldSweep.ptCloud[pt.nearPt3[0]][pt.nearPt3[1]].xyz;
//		Distance = abs((xi_hat - xj).dot((xj - xl).cross(xj - xm))) / ((xj - xl).cross(xj - xm)).norm();
//	}
//	else {
//		cout << "Edge or Plane indicator not provided!" << endl;
//	}
//
//	return Distance;
//}

//MatrixXd Sweep::GetJacobian(Sweep &OldSweep) {
//	MatrixXd JacobianFull, JacobianRow(1, 6);
//	VectorXd EstTransform = OldSweep.transform, EstTransform_Delta;
//	double Delta = pow(10, -6);
//	int cnt = 0;
//	int EnPflag = 1;
//	for (int i = 0; i < OldSweep.edgePts.size(); i++) {
//		for (auto &entry : edgePts[i]) {
//			for (int col = 0; col < 6; col++) {
//				EstTransform_Delta = EstTransform;
//				EstTransform_Delta(col) = EstTransform(col) + Delta;
//				JacobianRow(0, col) = (Distance2(OldSweep.ptCloud[i][entry], OldSweep, EstTransform_Delta, EnPflag) - 
//					Distance2(OldSweep.ptCloud[i][entry], OldSweep, EstTransform, EnPflag)) / Delta;
//			}
//			JacobianFull << JacobianFull, JacobianRow;
//		}
//		
//	}
//
//	return JacobianFull;
//}


int main(void)
{

	LoamPt myNewPt = LoamPt();
	std::vector<double> pt1 = { 1,1,1 }, pt2 = { 0, 1, 0 }, pt3, pt4;

	LoamPt filledPt = LoamPt({ 1,1,1 }, 20);
	myNewPt.SetXYZ(pt2);


	std::vector<std::vector<double>> newSlicePoints;

	for (double i = 0.0; i < 10.0; i++)
	{
		newSlicePoints.push_back({ i,i,i });
	}

	std::cout << "Size of this point = " << pt1.size() << std::endl;




	Sweep NewSweep;

	for (int i = 0; i < 5; i++)
	{
		NewSweep.AddSlice(0, i, newSlicePoints);
	}


	Vector3d xyz2 = NewSweep.ptCloud[0][6].xyz;

	std::cout << "There are " << NewSweep.ptCloud[0].size() << " pts in the first slice of our ptcloud" << std::endl;

	double c = Mult(std::vector<double>{ 1, 1, 1 }, std::vector<double>{ 1 });

	std::vector<double> &d = Mult(std::vector<double>{ 1, 1, 1 }, 2.0);
	std::vector<double> &e = Mult(2.0, std::vector<double>{ 1, 1, 1 });

	double a = 5;

	double f = Mult(d, e);

	std::vector<double> g = Divide(d, 3.0);

	double h = Dist(g);

	double hh = Dist(e, g);

	std::vector<double> aa = Divide(e, 3.0);

	std::vector<std::vector<double>> unsortedVec, sortedVec;

	for (int i = 0; i < 16; i++)
	{
		unsortedVec.push_back({ (double)rand(), (double)i });
	}

	sortedVec = unsortedVec;

	MergeSort(sortedVec);

	std::vector<std::vector<double>> testSlice;

	for (double i = 0; i < 100; i++)
	{
		std::cout << (int)i / 25 << std::endl;
		switch ((int)i / 25)
		{
		case 0:
			testSlice.push_back({ 50.0,0.0,i - 50 });
			continue;
		case 1:
			testSlice.push_back({ 50.0 + i - 25, 0, i - 50 });
			continue;
		case 2:
			testSlice.push_back({ 50.0 + 100.0, 0.0, i - 50 });
			continue;
		case 3:
			testSlice.push_back({ 50.0 + 50.0, 0.0, i - 50 });
			continue;
		}
	}

	Sweep testSweep(testSlice);

	for (auto &elem : testSlice)
	{
		elem[0] += (double)rand() / RAND_MAX * 2; // add a bit of noise
	}

	testSweep.FindEdges(testSweep.numSlices);

	testSweep.AddSlice(testSlice);

	testSweep.FindEdges(testSweep.numSlices);

	return 0;


}