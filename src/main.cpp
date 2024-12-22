#pragma once
#include <iostream>
//#include "VMAS_Skeleton.h"
#include "VSfast.h"

void normalizeMesh(Mesh& mesh) {
    // 1. ����ÿ�����������Сֵ�����ֵ
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

    // 2. �����Χ�еĶԽ��߳���
    Mesh::Point diagonal = max_pt - min_pt;
    float diagonal_length = diagonal.length();

    // ����Խ��߳���Ϊ�㣨��ζ������û����Ч�ĳߴ磩���򲻽�������
    if (diagonal_length > 0) {
        float scale_factor = 1.0f / diagonal_length;  // ������������

        // 3. ��ÿ������������ţ�ʹ���������һ��
        for (auto v : mesh.vertices()) {
            Mesh::Point& p = mesh.point(v);
            // ��ƽ�Ƶ�ԭ�㣬�ٰ��������ţ����ָ�����һ����ķ�Χ
            p = (p - min_pt) * scale_factor;
        }
    }
    else {
        std::cerr << "���棺����ĶԽ��߳���Ϊ�㣬�޷�������Ч�����š�" << std::endl;
    }
}

int main()
{
	std::string fpath = "../data/";
	std::string filename = fpath+ std::string("lamp.obj");

	Mesh mesh;
	if (!OpenMesh::IO::read_mesh(mesh, filename))
	{
		std::cout << "�޷���ȡ��������" << std::endl;
	}

    normalizeMesh(mesh);

	//VMAS_Skeleton ske(mesh, 0.2);
	VSfast vsf;
	vsf.init(mesh, 0.2, 30);
	vsf.run();
	vsf.write_color_obj(fpath+std::string("test"));
}

