#ifndef SWEEP_CLASS
#define SWEEP_CLASS

#include "LoamPt.h"
#include "LinearAlgebraFns.h"

class Sweep
{
public:
	std::vector<std::vector<LoamPt>> ptCloud;			// Cloud of LoamPts where ptCloud[i] = slice of LoamPts, and ptCloud[i][j] = individual LoamPt
	std::map<int, std::vector<int>> edgePts, planePts;	// 2D Vector containing the slice/pt indices of pts declared as edge/plane points
	std::vector<double> timeStamps;						// Vector of time values, where timeStamps[i] is the timeStamp corresponding to slice i of the ptCloud
	double tStart, tEnd, tCur;
	VectorXd transform;
	int sweepID, numSlices = -1, kernalSize = 11, regionPerSlice = 4, edgePerRegion = 2, planePerRegion = 4, edgeFindThreshold = 3, maxNumSlices = 720;
	int numEdges = 0, numPlanes = 0;
	Sweep();
	~Sweep();
	Sweep(std::vector<std::vector<double>> &inputSlice);
	void AddSlice(int sweepNumber, int sliceNumber, double timeStamp, std::vector<std::vector<double>> &inputSlice);
	void AddSlice(std::vector<std::vector<double>> &inputSlice);
	void FindAllEdges();
	void FindEdges(int sliceIndex);
	void SortCurvatures(int sliceIdx, std::vector<std::vector<double>> &curveVec, int startIdx, int endIdx);
	bool FindBestEdgePt(int sliceIdx, std::vector<std::vector<double>> &curveVec);
	bool EvaluateEdge(int sliceIdx, std::vector<double> &potentialPt);
	bool FindBestPlanePt(int sliceIdx, std::vector<std::vector<double>> &curveVec);
	bool EvaluatePlane(int sliceIdx, std::vector<double> &potentialPt);
	double Dist2Line(Vector3d &x_, Vector3d &x1, Vector3d &x2);
	double Dist2Plane(Vector3d &x_, Vector3d &x1, Vector3d &x2, Vector3d &x3);
	void FindCorrespondences(int sliceNumber, Sweep &OldSweep);
	void FindNearestLine(LoamPt &curEdgePt, Sweep &OldSweep);
	void FindNearestPlane(LoamPt &curPlanePt, Sweep &OldSweep);
	void TransformAll(VectorXd transform);
};

#endif // !SWEEP_CLASS