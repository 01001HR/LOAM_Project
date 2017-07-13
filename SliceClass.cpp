#include "Sweep.h"

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

	Sweep testSweep;

	testSweep.AddSlice(1, 0, 0.0, testSlice);

	for (auto &elem : testSlice)
	{
		elem[0] += (double)rand() / RAND_MAX * 2; // add a bit of noise
	}

	testSweep.FindEdges(testSweep.numSlices);

	testSweep.AddSlice(1, 1, 0.0, testSlice);

	testSweep.FindEdges(testSweep.numSlices);

	getchar();

	return 0;
}