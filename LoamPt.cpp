#include "LoamPt.h"

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
		sweepID = otherPt.sweepID;
		sliceID = otherPt.sliceID;
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

LoamPt::LoamPt(const std::vector<double> &xyzInput, int sweepInput, int sliceInput, double time) // vector input constructor
{
	std::cout << "In here" << std::endl;
	if (SetXYZ(xyzInput) == true)
	{
		timeStamp = time;
		sweepID = sweepInput;
		sliceID = sliceInput;
		filled = 1;
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
		dist = (xyz - otherPt.xyz).norm();
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
    else return 0.0;
}

inline const double LoamPt::GetY()
{
	if (filled == 1) return xyz[1];
    else return 0.0;
}

inline const double LoamPt::GetZ()
{
	if (filled == 1) return xyz[2];
    else return 0.0;
}

inline const double LoamPt::GetTime()
{
	if (filled == 1) return timeStamp;
    else return 0.0;
}

Vector3d LoamPt::Transform(Matrix4d &xformMatrix4x4)
{
	Vector4d augVec = { xyz[0], xyz[1], xyz[2], 1 }, newVec;
	newVec = xformMatrix4x4*augVec;
	return{ newVec[0], newVec[1], newVec[2] };
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