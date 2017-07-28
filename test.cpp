//#include "Sweep.h"
#include "LMOptim.h"
//#include <xmmintrin.h>
//#include <pmmintrin.h>
//#define _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
//#define _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);

int main(void)
{
	std::map<int, std::vector<double>> newMap;
	newMap[1] = { 1,1,1 };

	std::vector<double> tt = newMap[1], b = newMap[2];

	std::vector<std::vector<double>> output = ParseBinary<float, double>((float)1, (double)1, "velo_test.bin");

	std::vector<std::vector<std::vector<double>>> orderedPoints;

	orderedPoints = OrganizePoints(output);

	Sweep oldTestSweep, newTestSweep;

	newTestSweep.transform.resize(6);
	newTestSweep.transform << 1.0, -0.5, 1.0, -0.1, 0.0, 0.2;

	int i = 0;
	for (auto &slice : orderedPoints)
	{
		oldTestSweep.AddSlice(0, i, i, slice);
		oldTestSweep.timeStamps.push_back(i);
		oldTestSweep.tCur = i;
		i++;
	}
	oldTestSweep.tStart = 0;
	oldTestSweep.tEnd = i;
	oldTestSweep.tCur = i;

	std::cout << "OldSweepFilled" << std::endl;

	newTestSweep.tStart = i;
	i = 0;
	for (auto &slice : orderedPoints)
	{
		newTestSweep.AddSlice(1, i, 719 + i, slice);
		newTestSweep.timeStamps.push_back(719+i);
		newTestSweep.tCur = 719+i;
		i++;
		//if (i == 19) {
		//	break;
		//}
	}
	newTestSweep.tEnd = 719+i;
	newTestSweep.tCur = 719+i;

	std::cout << "NewSweepFilled" << std::endl;

	std::cout << "Transforming New sweep" << std::endl;
	newTestSweep.TransformAll(newTestSweep.transform);

	std::cout << "Finding All Edges" << std::endl;
	oldTestSweep.FindAllEdges();
	newTestSweep.FindAllEdges();

	std::cout << "Finding correspondences" << std::endl;
	for (int i = 0; i < newTestSweep.ptCloud.size(); i++)
	{
		std::cout << i << std::endl;
		newTestSweep.FindCorrespondences(i, oldTestSweep);
	}
	std::cout << "Done with correspondences" << std::endl;

	std::cout << "EdgePts:" << std::endl;
	for (auto &edgeKeyVal : newTestSweep.edgePts)
	{
		for (auto &idx : edgeKeyVal.second)
		{
			std::cout << " | "  << newTestSweep.ptCloud[edgeKeyVal.first][idx].dist << " | ";
		}
	}

	std::cout << "PlanePts:" << std::endl;
	for (auto &planeKeyVal : newTestSweep.planePts)
	{
		for (auto &idx : planeKeyVal.second)
		{
			std::cout << " | " << newTestSweep.ptCloud[planeKeyVal.first][idx].dist << " | ";
		}
	}

	Vector3d p = { 1, 2, 3 }, p_, p2_, pf, pf2;

	VectorXd twist(6);
	twist <<  10, 10, 10, 1, 1, 1;

	p_ = BackTransform(p, twist, 1.0);

	p2_ = BackTransform(p, twist, 0.5);

	pf = ForwardTransform(p, twist, 1.0);
	pf2 = ForwardTransform(p, twist, 0.5);

	LoamPt myNewPt = LoamPt();
	std::vector<double> pt1 = { 1,1,1 }, pt2 = { 0, 1, 0 }, pt3, pt4;

	LoamPt filledPt = LoamPt({ 1,1,1 }, 20);
	myNewPt.SetXYZ(pt2);


	std::vector<std::vector<double>> newSlicePoints;

	for (double i = 0.0; i < 10.0; i++)
	{
		newSlicePoints.push_back({ i,i,i });
	}

	//std::cout << "Size of this point = " << pt1.size() << std::endl;

	Sweep NewSweep;

	for (int i = 0; i < 5; i++)
	{
		NewSweep.AddSlice(1, 0, i, newSlicePoints);
	}

	Vector3d xyz2 = NewSweep.ptCloud[0][6].xyz;

	std::vector<std::vector<double>> unsortedVec, sortedVec;

	for (int i = 0; i < 16; i++)
	{
		unsortedVec.push_back({ (double)rand(), (double)i });
	}

	sortedVec = unsortedVec;

	MergeSort(sortedVec,0);

	std::vector<std::vector<double>> testSlice;

	for (double i = 0; i < 100; i++)
	{
		//std::cout << (int)i / 25 << std::endl;
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

	Sweep testSweep;

	testSweep.AddSlice(1, 0, 0.0, testSlice);

	for (auto &elem : testSlice)
	{
		elem[0] += (double)rand() / RAND_MAX * 2; // add a bit of noise
	}

	testSweep.FindEdges(testSweep.numSlices);

	testSweep.AddSlice(1, 1, 0.0, testSlice);

	testSweep.FindEdges(testSweep.numSlices);

	Matrix3d I = Matrix3d::Identity();

	LMOptim optimizer;

	VectorXd estT = optimizer.TransformEstimate(oldTestSweep, newTestSweep);


	return 0;
}