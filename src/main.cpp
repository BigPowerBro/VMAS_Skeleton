#pragma once
#include <iostream>
#include "VMAS_Skeleton.h"

int main()
{
	std::string filename = "D:\\study\\DIP\\pre\\obj\\bunny.obj";

	Mesh mesh;
	if (!OpenMesh::IO::read_mesh(mesh, filename))
	{
		std::cout << "�޷���ȡ��������" << std::endl;
	}

	VMAS_Skeleton ske(mesh);
}