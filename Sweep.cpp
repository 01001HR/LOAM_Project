#include "Sweep.h"
#include "LMOptim.h"

Sweep::Sweep()
{

}

Sweep::~Sweep()
{

}

Sweep::Sweep(std::vector<std::vector<double>> &inputSlice)
{
	sweepID = NULL;
	transform.resize(6);
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

void Sweep::FindAllEdges(void)
{
	// Loop through all current slices and find best feature points
	for (int i = 0; i < ptCloud.size(); i++)
	{
		std::cout << "Finding Edges in i = " << i << std::endl;
		FindEdges(ptCloud[i][0].sliceID);
	}
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
	MergeSort(curvatures,0);

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
		numEdges++;
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
		numPlanes++;
		return true;
	}
	return false;
}

bool Sweep::EvaluateEdge(int sliceIdx, std::vector<double> &potentialPt) // Checks to make sure that no nearby points are drastically closer to the sensor
{
	double ptDist = ptCloud[sliceIdx][potentialPt[1]].xyz[0]; // treating ptDist[0] = X, the depth-distance w.r.t. sensor
	for (int i = -(kernalSize / 2); i < 0; i++)
	{
		if (potentialPt[1] + i < 0) {
			continue;
		}
		else if ((ptCloud[sliceIdx][potentialPt[1] + i + 1].xyz[0] - ptCloud[sliceIdx][potentialPt[1] + i].xyz[0]) > edgeFindThreshold) // large increase in distance as we approach the potential point == occlusion
		{
			return false;
		}
	}
	for (int i = 1; i < (kernalSize / 2) + 1; i++)
	{
		if (potentialPt[1] + i + 1 >= ptCloud[sliceIdx].size())
		{
			break;
		}
		else if ((ptCloud[sliceIdx][potentialPt[1] + i + 1].xyz[0] - ptCloud[sliceIdx][potentialPt[1] + i].xyz[0]) < -edgeFindThreshold) // large decrease in distance as we move away from the potential point == occlusion
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


double Sweep::Dist2Line(Vector3d &x_, Vector3d &x1, Vector3d &x2)
{
	// Takes a back-transformed point x_ and returns the orthogonal distance between the point and line formed by nearby pts x1,x2
	if ((x2 - x1).norm() < 10e-6)
	{
		return (x_ - x1).norm();
	}
	return ((x_ - x1).cross(x_ - x2)).norm() / (x2 - x1).norm();
}

double Sweep::Dist2Plane(Vector3d &x_, Vector3d &x1, Vector3d &x2, Vector3d &x3)
{
	// Takes a back-transformed point x_ and returns the orthogonal distance between the point and line formed by nearby pts x1,x2,x3
	if ((x2-x1).norm() < 10e-6)
	{
		return Dist2Line(x_, x1, x3);
	}
	else if ((x3 - x1).norm() < 10e-6)
	{
		return Dist2Line(x_, x1, x2);
	}
	else
	{
		return ((x_ - x1).cross((x1 - x2).cross(x1 - x3))).norm() / ((x1 - x2).cross(x1 - x3)).norm();
	}
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
	Vector3d x_ = BackTransform(curEdgePt.xyz, transform, (timeStamps[curEdgePt.sliceID] - tStart)/(tCur - tStart));

	std::vector<int> bestPt = { 0,0 };
	double bestDist = 10e10, tempDist;
	std::vector<double> distances;

	// Find the nearest edgepoint located in the +/- n-neighboring slices of the previous sweep
	int j;
	for (int i = curEdgePt.sliceID - 2; i < curEdgePt.sliceID + 3; i++)
	{
		j = i % maxNumSlices;
		for (auto &oldIdx : OldSweep.edgePts[j])
		{
			tempDist = (OldSweep.ptCloud[j][oldIdx].xyz - x_).norm();
			if (tempDist < bestDist)
			{
				bestDist = tempDist;
				bestPt = { j, oldIdx };
			}
		}
	}
	curEdgePt.nearPt1 = bestPt;
	bestDist = 10e10;

	// Find edgepoint2 closest to nearPt1 located in either adjacent slice of the previous sweep
	for (auto &i : { (curEdgePt.nearPt1[0] - 1) % maxNumSlices, (curEdgePt.nearPt1[0] + 1) % maxNumSlices })
	{
		for (auto &oldIdx : OldSweep.edgePts[i])
		{
			tempDist = (OldSweep.ptCloud[i][oldIdx].xyz - OldSweep.ptCloud[curEdgePt.nearPt1[0]][curEdgePt.nearPt1[1]].xyz).norm();
			if (tempDist < bestDist)
			{
				bestDist = tempDist;
				bestPt = { i, oldIdx };
			}
		}
	}
	curEdgePt.nearPt2 = bestPt;
	curEdgePt.dist = Dist2Line(x_, ptCloud[curEdgePt.nearPt1[0]][curEdgePt.nearPt1[1]].xyz, ptCloud[curEdgePt.nearPt2[0]][curEdgePt.nearPt2[1]].xyz); // this needs to be changed
}

void Sweep::FindNearestPlane(LoamPt &curEdgePt, Sweep &OldSweep)
{
	// Back-transform the point to the beginning of the current sweep.
	Vector3d x_ = BackTransform(curEdgePt.xyz, transform, (timeStamps[curEdgePt.sliceID] - tStart) / (tCur - tStart));

	std::vector<int> bestPt = { 0,0 };
	double bestDist = 10e10, tempDist;
	std::vector<double> distances;

	// Find the nearest planePoint located in the +/- n-neighboring slices of the previous sweep
	int j;
	for (int i = curEdgePt.sliceID - 2; curEdgePt.sliceID + i < 3; i++)
	{
		j = i % maxNumSlices;
		for (auto &oldIdx : OldSweep.edgePts[j])
		{
			tempDist = (OldSweep.ptCloud[j][oldIdx].xyz - x_).norm();
			if (tempDist < bestDist)
			{
				bestDist = tempDist;
				bestPt = { j, oldIdx };
			}
		}
	}
	curEdgePt.nearPt1 = bestPt;

	bestDist = 10e10;
	// Find the closest planePoint to nearPt1 in the same slice of the previous sweep
	for (auto &oldIdx : OldSweep.edgePts[curEdgePt.nearPt1[0]])
	{
		if (oldIdx != curEdgePt.nearPt1[1])
		{
			tempDist = (OldSweep.ptCloud[curEdgePt.nearPt1[0]][curEdgePt.nearPt1[1]].xyz - OldSweep.ptCloud[curEdgePt.nearPt1[0]][oldIdx].xyz).norm();
			if (tempDist < bestDist)
			{
				bestDist = tempDist;
				bestPt = { curEdgePt.nearPt1[0], oldIdx };
			}
		}
	}
	curEdgePt.nearPt3 = bestPt;

	bestDist = 10e10;
	// Find the closest to planePoint to nearPt1 located in either adjacent slice of the previous sweep
	for (auto &i : { (curEdgePt.nearPt1[0] - 1) % maxNumSlices, (curEdgePt.nearPt1[0] + 1) % maxNumSlices })
	{
		for (auto &oldIdx : OldSweep.edgePts[i])
		{
			tempDist = (OldSweep.ptCloud[i][oldIdx].xyz - OldSweep.ptCloud[curEdgePt.nearPt1[0]][curEdgePt.nearPt1[1]].xyz).norm();
			if (tempDist < bestDist)
			{
				bestDist = tempDist;
				bestPt = { i, oldIdx };
			}
		}
	}
	curEdgePt.nearPt2 = bestPt;
	curEdgePt.dist = Dist2Plane(x_, ptCloud[curEdgePt.nearPt1[0]][curEdgePt.nearPt1[1]].xyz, ptCloud[curEdgePt.nearPt2[0]][curEdgePt.nearPt2[1]].xyz, ptCloud[curEdgePt.nearPt3[0]][curEdgePt.nearPt3[1]].xyz); // this needs to be changed
}

void Sweep::TransformAll(VectorXd transform)
{
	for (int i = 0; i < ptCloud.size(); i++)
	{
		for (auto &pt : ptCloud[i])
		{
			pt.xyz = ForwardTransform(pt.xyz, transform, (timeStamps[i] - tStart)/(tCur-tStart));
		}
	}
}
