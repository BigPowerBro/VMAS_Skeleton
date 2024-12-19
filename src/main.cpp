#pragma once
#include <iostream>
#include "VMAS_Skeleton.h"

int main()
{
	std::string fpath = "../data/";
	std::string filename = fpath+ std::string("cuboid.obj");

	Mesh mesh;
	if (!OpenMesh::IO::read_mesh(mesh, filename))
	{
		std::cout << "无法读取三角网格" << std::endl;
	}

	VMAS_Skeleton ske(mesh, 0.2);
}