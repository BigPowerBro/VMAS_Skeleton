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

struct SphereClass
{
	int id;
	Sphere s;
	std::vector<int> cluster;
	float E;
};

class VMAS_Skeleton
{
public:
	VMAS_Skeleton(Mesh& mesh, float lambda) { Init(mesh, lambda); };
	~VMAS_Skeleton() {};

public:
	void Init(Mesh& mesh,float lambda);
private:
	// ���뺯��
	float dps(Point3D p, Sphere s);
	float dpns(Plane a, Sphere s);
	float Qpns(Plane a, Sphere s);
	float Dvs(const int v, const Sphere s);  // Dvs
	float Qvs(const int v, const Sphere s);  // Qvs
	float sign(float a);

	// ��������
	float E_c(const std::vector<int> cluster, const Sphere sphere);
	float E_SQEM(const std::vector<int> cluster, const Sphere sphere);
	float E_euclidean(const std::vector<int> cluster, const Sphere sphere);

	/// ���������λ�õİ뾶
	void update_sphere();
	// ����һ��� �����Ӧ����
	Sphere update_single_sphere(const std::vector<int>& cluster, const Sphere init_sphere, float lambda, const float eps);

	// �������㷨 
	Sphere ShrinkingBall(const int v);
	float ComputeRadius(const Point3D p, const Point3D q, const Vector3D n);
	Point3D CalClosestPoint(const Point3D p);



private:
	Mesh input_mesh;  //��������
	float lambda;
	
	std::vector<SphereClass> sphere_classes;  
	int *vertices_type;
};

