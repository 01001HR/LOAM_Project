#ifndef LIN_ALG_FNS
#define LIN_ALG_FNS

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <math.h>
#include <Eigen\Core>
#include <random>
#include <fstream>
#include <iterator>

using namespace Eigen;

template <class numType>
inline numType Mult(const std::vector<numType> &v1, const std::vector<numType> &v2)
{
	numType product = (numType)0;
	if (v1.size() == v2.size() && v1.size() > 1)
	{
		for (int i = 0; i < fmin(v1.size(), v2.size()); i++)
		{
			product += v1[i] * v2[i];
		}
		return product;
	}

	else
	{
		std::cout << "Input vector sizes do not match" << std::endl;
	}
}

template <class numType>
inline std::vector<numType> Mult(std::vector<numType> v1, const numType v2)
{
	if (v1.size() >= 1)
	{
		for (auto &elem : v1)
		{
			elem *= v2;
		}
		return v1;
	}
	else
	{
		std::cout << "Input vector is not size >= 1" << std::endl;
	}
}

template <class numType>
inline std::vector<numType> Mult(const numType v1, std::vector<numType> v2)
{
	if (v2.size() >= 1)
	{
		for (auto &elem : v2)
		{
			elem *= v1;
		}
		return v2;
	}
	else
	{
		std::cout << "Input vector is not size >= 1" << std::endl;
	}
}

template <class numType>
inline std::vector<numType> Add(std::vector<numType> v1, const std::vector<numType> &v2)
{
	if (v1.size() == v2.size() && v1.size() > 1)
	{
		for (int i = 0; i < v1.size(); i++)
		{
			v1[i] += v2[i];
		}
		return v1;
	}
	else
	{
		std::cout << "Input vector sizes do not match" << std::endl;
	}
}

template <class numType>
inline std::vector<numType> Add(std::vector<numType> v1, const numType v2)
{
	if (v1.size() >= 1)
	{
		for (auto &elem : v1)
		{
			elem += v2;
		}
		return v1;
	}
	else
	{
		std::cout << "Input vector is not size >= 1" << std::endl;
	}
}

template <class numType>
inline std::vector<numType> Add(const numType v1, std::vector<numType> v2)
{
	if (v2.size() >= 1)
	{
		for (auto &elem : v2)
		{
			elem += v1;
		}
		return v2;
	}
	else
	{
		std::cout << "Input vector is not size >= 1" << std::endl;
	}
}

template <class numType>
inline std::vector<numType> Minus(std::vector<numType> v1, const std::vector<numType> &v2)
{
	std::cout << v1.size() << std::endl;
	std::cout << "vec subtract" << v1[0] << v1[1] << v1[2] << v2[0] << v2[1] << v2[2] << std::endl;
	if (v1.size() == v2.size() && v1.size() > 1)
	{
		for (int i = 0; i < v1.size(); i++)
		{
			v1[i] -= v2[i];
			std::cout << "vec subtract" << v1[0] << v1[1] << v1[2] << std::endl;
		}
		return v1;
	}
	else
	{
		std::cout << "Input vector sizes do not match" << std::endl;
	}
}

template <class numType>
inline std::vector<numType> Minus(std::vector<numType> v1, const numType v2)
{
	if (v1.size() >= 1)
	{
		for (auto &elem : v1)
		{
			elem -= v2;
		}
		return v1;
	}
	else
	{
		std::cout << "Input vector is not size >= 1" << std::endl;
	}
}

template <class numType>
inline std::vector<numType> Minus(const numType v1, std::vector<numType> v2)
{
	if (v2.size() >= 1)
	{
		for (auto &elem : v2)
		{
			elem -= v1;
		}
		return v2;
	}
	else
	{
		std::cout << "Input vector is not size >= 1" << std::endl;
	}
}

template <class numType>
inline std::vector<numType> Divide(std::vector<numType> v1, const numType v2)
{
	return Mult(v1, (numType)1 / v2);
}

template <class numType>
inline std::vector<numType> Divide(const numType v1, std::vector<numType> v2)
{
	return Mult((numType)1 / v1, v2);
}

template <class numType>
inline numType Dist(const std::vector<numType> &v1)
{
	numType dist = (numType)0;
	for (auto &elem : v1)
	{
		dist += pow(elem, 2);
	}
	return sqrt(dist);
}

template <class numType>
inline numType Dist(std::vector<numType> v1, const std::vector<numType> &v2)
{
	return Dist(Minus(v1, v2));
}

inline void Merge(std::vector<std::vector<double>> &incomingArray, int idx1, int end1, int idx2, int end2)
{
	std::vector<std::vector<double>> sortedArray; // initial allocation, numbers will be changed during sorting
	int firstIdx = idx1, lastIdx = end2;
	sortedArray.resize(lastIdx - firstIdx, { 0,0 });
	int i = 0;
	for (;;)
	{
		if (idx1 < end1 && idx2 < end2)
		{
			if (incomingArray[idx1][0] <= incomingArray[idx2][0])
			{
				sortedArray[i] = std::vector<double>(incomingArray[idx1]);
				idx1++;
			}
			else
			{
				sortedArray[i] = std::vector<double>(incomingArray[idx2]);
				idx2++;
			}
		}
		else if (idx1 < end1)
		{
			sortedArray[i] = std::vector<double>(incomingArray[idx1]);
			idx1++;
		}
		else if (idx2 < end2)
		{
			sortedArray[i] = std::vector<double>(incomingArray[idx2]);
			idx2++;
		}
		else
		{
			for (int j = 0; j < lastIdx - firstIdx; j++)
			{
				incomingArray[firstIdx + j] = std::vector<double>(sortedArray[j]);
			}
			return;
		}
		i++;
	}
}

inline void MergeSort(std::vector<std::vector<double>> &vec)
{
	int size = 1;
	int len = vec.size();
	int start1, end1, start2, end2;
	for (auto &elem : vec)
	{
		std::cout << " " << elem[0] << " " << elem[1];
	}
	std::cout << std::endl;
	while (size < vec.size())
	{
		for (int i = 0; i < len - size; i += size * 2)
		{
			start1 = i;
			end1 = start1 + size;
			start2 = end1;
			end2 = (int)fmin(start2 + size, len);
			Merge(vec, start1, end1, start2, end2);
			for (auto &elem : vec)
			{
				std::cout << " " << elem[0] << " " << elem[1];
			}
			std::cout << std::endl;
		}
		size *= 2;
	}
}

#endif