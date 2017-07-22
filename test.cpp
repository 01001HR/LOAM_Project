#include "Sweep.h"

//#include <xmmintrin.h>
//#include <pmmintrin.h>
//#define _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
//#define _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);

int main(void)
{
	std::map<int, std::vector<double>> newMap;
	newMap[1] = { 1,1,1 };

	std::vector<double> tt = newMap[1], b = newMap[2];

	int key;
	for (auto &keyValPair : newMap)
	{
		key = keyValPair.first;
		std::cout << key << std::endl;
	}

	for (auto &elem : tt)
	{
		std::cout << elem << std::endl;
	}

	for (auto &elem : b)
	{
		std::cout << elem << std::endl;
	}

	std::vector<std::vector<double>> output = ParseBinary<float, double>((float)1, (double)1, "velo_test.bin");

	std::vector<std::vector<std::vector<double>>> orderedPoints;

	orderedPoints = OrganizePoints(output);

	Sweep oldTestSweep, newTestSweep;

	newTestSweep.transform.resize(6);
	newTestSweep.transform << 1.0, 1.0, 1.0, 0.0, 0.0, M_PI/2.0;

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
	for (auto &slice : orderedPoints)
	{
		newTestSweep.AddSlice(0, i, i, slice);
		newTestSweep.timeStamps.push_back(i);
		newTestSweep.tCur = i;
		i++;
	}
	newTestSweep.tEnd = i;
	newTestSweep.tCur = i;

	std::cout << "NewSweepFilled" << std::endl;

	newTestSweep.TransformAll(newTestSweep.transform);

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

	std::cout << "Size of this point = " << pt1.size() << std::endl;

	Sweep NewSweep;

	for (int i = 0; i < 5; i++)
	{
		NewSweep.AddSlice(1, 0, i, newSlicePoints);
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

	MergeSort(sortedVec,0);

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

	return 0;
}