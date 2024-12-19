#pragma once
#include <iostream>
#include "VMAS_Skeleton.h"

int main()
{
	std::string fpath = "C:/MyProgram/SDK/project/c++/dip/VMAS_Skeleton/data/";
	std::string filename = fpath+ std::string("cube.obj");

	Mesh mesh;
	if (!OpenMesh::IO::read_mesh(mesh, filename))
	{
		std::cout << "无法读取三角网格" << std::endl;
	}

	VMAS_Skeleton ske(mesh, 0.2);
}