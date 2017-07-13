#ifndef LOAM_PT_CLASS
#define LOAM_PT_CLASS

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <math.h>

#ifdef WIN32
#include <Eigen\Core>
#include <Eigen\Dense>
#else
#include <Eigen/Core>
#include <Eigen/Dense>
#endif

#include <random>
#include <fstream>
#include <iterator>
#include <map>


using namespace Eigen;
using namespace std;

class LoamPt
{
public:
	Vector3d xyz;
	int sweepID = 0;
	int sliceID = 0;
	std::vector<int> nearPt1, nearPt2, nearPt3; // these will be len 2 vectors with pt = {sliceIndex, pointIndex}
	int filled = 0;
	double timeStamp = 0.0;

	LoamPt();
	~LoamPt();
	LoamPt(const std::vector<double> &xyzInput); // no time given
	LoamPt(const std::vector<double> &xyzInput, double time);
	LoamPt(const std::vector<double> &xyzInput, int sweepInput, int sliceInput, double time);
	LoamPt(const double x, const double y, const double z, const double time);
	LoamPt(const LoamPt &otherPt);
	LoamPt &operator = (const LoamPt &otherPt);

	bool SetXYZ(const std::vector<double> &newXYZ);

	inline const double GetX();
	inline const double GetY();
	inline const double GetZ();
	inline const double GetTime();

	bool Distance(double &dist, LoamPt &otherPt);
	Vector3d Transform(Matrix4d &xformMatrix4x4);
	void TransformSelf(Matrix4d &xformMatrix4x4); // xformMatrix is stored as stacked row-vectors
	Vector3d BackTransform(VectorXd &tVec);
};

#endif // !1

