#pragma once
#include <iostream>
//#include "VMAS_Skeleton.h"
#include "VSfast.h"

void normalizeMesh(Mesh& mesh) {
    // 1. 计算每个坐标轴的最小值和最大值
    Mesh::Point min_pt(std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max());
    Mesh::Point max_pt(std::numeric_limits<float>::lowest(),
        std::numeric_limits<float>::lowest(),
        std::numeric_limits<float>::lowest());

    for (auto v : mesh.vertices()) {
        const Mesh::Point& p = mesh.point(v);
        min_pt[0] = std::min(min_pt[0], p[0]);
        min_pt[1] = std::min(min_pt[1], p[1]);
        min_pt[2] = std::min(min_pt[2], p[2]);
        max_pt[0] = std::max(max_pt[0], p[0]);
        max_pt[1] = std::max(max_pt[1], p[1]);
        max_pt[2] = std::max(max_pt[2], p[2]);
    }

    // 2. 计算包围盒的对角线长度
    Mesh::Point diagonal = max_pt - min_pt;
    float diagonal_length = diagonal.length();

    // 如果对角线长度为零（意味着网格没有有效的尺寸），则不进行缩放
    if (diagonal_length > 0) {
        float scale_factor = 1.0f / diagonal_length;  // 计算缩放因子

        // 3. 对每个顶点进行缩放，使网格整体归一化
        for (auto v : mesh.vertices()) {
            Mesh::Point& p = mesh.point(v);
            // 先平移到原点，再按比例缩放，最后恢复到归一化后的范围
            p = (p - min_pt) * scale_factor;
        }
    }
    else {
        std::cerr << "警告：网格的对角线长度为零，无法进行有效的缩放。" << std::endl;
    }
}

int main()
{
	std::string fpath = "../data/";
	std::string filename = fpath+ std::string("lamp.obj");

	Mesh mesh;
	if (!OpenMesh::IO::read_mesh(mesh, filename))
	{
		std::cout << "无法读取三角网格" << std::endl;
	}

    normalizeMesh(mesh);

	//VMAS_Skeleton ske(mesh, 0.2);
	VSfast vsf;
	vsf.init(mesh, 0.2, 30);
	vsf.run();
	vsf.write_color_obj(fpath+std::string("test"));
}

