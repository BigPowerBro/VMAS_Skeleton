#ifndef VSFAST_H
#define VSFAST_H
#include<Eigen/dense>
#include<vector>
#include <MeshDef.h>
#include<iostream>
#include<map>
#include<set>
#include <Eigen/Sparse>
#include<string>
#include<fstream>
#include "../ANN/ANN.h"

class ANNClass
{
public:
	//ANNClass() {};
	//ANNClass(std::vector<Eigen::Vector3d>& vectorP);
	ANNClass(Eigen::MatrixXd LPR);
	Eigen::Vector3d AnnSearch(Eigen::Vector3d vIn);
	//Eigen::Vector3d AnnSearch(Eigen::Vector3d vIn, int& id);

public:
	//ann����
	//3ά����
	int dim = 3;
	//���ݸ���
	int	maxPts;
	//�ٽ���
	int	k = 1;
	//����������ĵ�
	ANNidxArray nnIdx = new ANNidx[k];
	//distance
	ANNdistArray dists = new ANNdist[k]; 
	//kd��
	ANNkd_tree* kdTree;
	//������������
	ANNpoint queryPt = annAllocPt(dim);
	//��������
	double eps = 0.01;
	//data
	ANNpointArray dataPts;

};


struct Sphere
{
	int id;
	Eigen::Vector4d s;//����+�뾶
	std::vector<int> cluster;//��ĵ����
	double E;//����

	Sphere(int id0, Eigen::Vector4d s0, std::vector<int> cluster0, double E0) : id(id0), s(s0), cluster(cluster0), E(E0) {};
};

class VSfast
{
public:
	void init(Mesh mesh,double lambda, int max_sphere_ball);//��ʼ��
	void update_spheres_s_E();//�������ĺͰ뾶
	void update_spheres_cluster();//������ĵ����
	void update_spheres();
	void split_spheres();//������
	

	void correction_spheres();//����ȥģ�͵�����С��ģ���ڲ�

	void run();

public:
	void cal_spheres_adjacency();//��������ڽӹ�ϵ
	void write_color_obj(std::string fname);
private:
	// ���������ĺ���
	void cal_QEM_matrix(const std::vector<int>& cluster, Eigen::Matrix4d& A, Eigen::Vector4d& b, double& c);
	void cal_ps_as(const std::vector<int>& cluster, Eigen::MatrixXd& ps, Eigen::VectorXd& as);
	double Dvs(const int v, const Eigen::Vector4d s);
	double Qvs(const int v, const Eigen::Vector4d s);
	double E_SQEM(const Eigen::Matrix4d A, const Eigen::Vector4d b, const double c, const Eigen::Vector4d s);
	double E_euclidean(const Eigen::MatrixXd ps, const Eigen::VectorXd as, const Eigen::Vector4d s);
	//����һ�������λ�úͰ뾶 �Լ� ����
	void update_single_sphere(std::shared_ptr<Sphere> sphere, const double eps);
	double line_search(double a, double b, std::function<double(double)> func);  
	
	/// �������㷨
	Eigen::Vector4d shringking_ball(Eigen::Vector4d p, const Eigen::Vector4d n);
	double compute_radius(const Eigen::Vector3d p, const Eigen::Vector3d q, const Eigen::Vector3d n);
	Eigen::Vector3d cal_closest_point(const Eigen::Vector3d c, const Eigen::Vector3d p);
	Eigen::Vector4d shringking_ball_ann(const int v);
	

	

private:
	Mesh mesh;
	Mesh ske_mesh;
	Eigen::MatrixXd point_pos; // n*4������ w=0
	Eigen::MatrixXd point_n;// n*4�ķ�����,w=1
	Eigen::VectorXd point_area;//n*1ÿ������Χ�����ε�����֮һ������
	std::vector<std::shared_ptr<Sphere>> spheres;//���ָ��
	Eigen::MatrixXi spheres_adjacency;//����ڽӹ�ϵ,��������
	std::map<int,std::set<int>> spheres_adjacency_map;//����ڽӹ�ϵ,��map��ʾ
	double lambda;//ϵ��
	double threshold; //������ֵ1 ����ֹͣ�����Ż�
	int max_sphere_num;  // ����Ҫ���ٸ���
	
private:
	Eigen::SparseMatrix<int> point_adjacency;
	std::shared_ptr<ANNClass> annclass;
};






#endif 