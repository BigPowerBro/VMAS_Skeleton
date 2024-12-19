#pragma once
#include <iostream>
// Eigen
#include <Eigen\Dense>
#include <Eigen\Core>

#include "MeshDef.h"

typedef Eigen::Vector3f Point3D;
typedef Eigen::Vector3f Vector3D;

struct Sphere
{
	Point3D q;
	float r;
};

struct Plane
{
	Point3D p;
	Vector3D n;
};

class VMAS_Skeleton
{
public:
	VMAS_Skeleton(Mesh& mesh) { Init(mesh); };
	~VMAS_Skeleton() {};

public:
	void Init(Mesh& mesh);
private:
	// 距离函数
	float dps(Point3D p, Sphere s);
	float dpns(Plane a, Sphere s);
	float Qpns(Plane a, Sphere s);
	float Dvs(const int v, const Sphere s);  // Dvs
	float Qvs(const int v, const Sphere s);  // Qvs
	float sign(float a);

	// 能量计算
	float E_SQEM(const std::vector<int> cluster, const Sphere sphere);
	float E_euclidean(const std::vector<int> cluster, const Sphere sphere);

	/// 更新球体的位置的半径
	void update_sphere();
	// 给定一族点 计算对应球体
	Sphere update_single_sphere(const std::vector<int>& cluster, const Sphere init_sphere, float lambda, const float eps);

private:
	Mesh input_mesh;  //输入网格
	std::vector<Sphere> sphere_list; //存放所有 sphere
	std::vector<std::vector<int>> clusters; //所有类别，每个类别的所有顶点
};

