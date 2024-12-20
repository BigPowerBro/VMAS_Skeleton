#ifndef VSFAST_H
#define VSFAST_H
#include<Eigen/dense>
#include<vector>
#include <MeshDef.h>
#include<iostream>
#include<map>
#include<set>
struct Sphere
{
	Eigen::Vector4d s;//����+�뾶
	std::vector<int> cluster;//��ĵ����
	double E;//����

	Sphere(Eigen::Vector4d s0, std::vector<int> cluster0, double E0) : s(s0), cluster(cluster0), E(E0) {};
};

class VSfast
{
public:
	void init(Mesh mesh,double lambda,double threshold);//��ʼ��
	void updata_spheres_s_E();//�������ĺͰ뾶
	void updata_spheres_cluster();//������ĵ����
	void split_spheres();//������

	void correction_spheres();//����ȥģ�͵�����С��ģ���ڲ�

	void run();

public:
	void cal_spheres_adjacency();//��������ڽӹ�ϵ

private:
	// ���������ĺ���
	void cal_QEM_matrix(const std::vector<int>& cluster, Eigen::Matrix4d& A, Eigen::Vector4d& b, double& c);
	void cal_ps_as(const std::vector<int>& cluster, Eigen::MatrixXd& ps, Eigen::VectorXd& as);
	double E_SQEM(const Eigen::Matrix4d A, const Eigen::Vector4d b, const double c, const Eigen::Vector4d s);
	double E_euclidean(const Eigen::MatrixXd ps, const Eigen::VectorXd as, const Eigen::Vector4d s);
	//����һ�������λ�úͰ뾶 �Լ� ����
	void updata_single_sphere(std::shared_ptr<Sphere> sphere, const double eps);
	double line_search(double a, double b, std::function<double(double)> func);  
	
	/// �������㷨
	Eigen::Vector4d shringking_ball(const int v);
	Eigen::Vector4d cal_closest_point(const Eigen::Vector4d c);
	double compute_radius(const Eigen::Vector4d p, const Eigen::Vector4d q, const Eigen::Vector4d n);

private:
	Mesh mesh;
	Eigen::MatrixXd point_pos; // n*4������ w=0
	Eigen::MatrixXd point_n;// n*4�ķ�����,w=1
	Eigen::VectorXd point_area;//n*1ÿ������Χ�����ε�����֮һ������
	std::vector<std::shared_ptr<Sphere>> spheres;//���ָ��
	Eigen::MatrixXi spheres_adjacency;//����ڽӹ�ϵ
	double lambda;//ϵ��
	double threshold;//������ֵ
	
private:
};






#endif // !VSFAST_H