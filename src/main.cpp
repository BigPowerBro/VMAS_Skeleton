#pragma once
#include <iostream>
//#include "VMAS_Skeleton.h"
#include "VSfast.h"

int main()
{
	std::string fpath = "../data/";
	std::string filename = fpath+ std::string("bunny.obj");

	Mesh mesh;
	if (!OpenMesh::IO::read_mesh(mesh, filename))
	{
		std::cout << "无法读取三角网格" << std::endl;
	}

	//VMAS_Skeleton ske(mesh, 0.2);
	VSfast vsf;
	vsf.init(mesh, 0.2, 10);
	vsf.run();
	vsf.write_color_obj(fpath+std::string("test.obj"));
}