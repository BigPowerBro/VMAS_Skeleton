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
	Eigen::Vector4d s;//球心+半径
	std::vector<int> cluster;//球的点的类
	double E;//能量

	Sphere(Eigen::Vector4d s0, std::vector<int> cluster0, double E0) : s(s0), cluster(cluster0), E(E0) {};
};

class VSfast
{
public:
	void init(Mesh mesh,double lambda,double threshold);//初始化
	void updata_spheres_s_E();//更新球心和半径
	void updata_spheres_cluster();//更新球的点的类
	void split_spheres();//分裂球

	void correction_spheres();//将出去模型的球缩小到模型内部

	void run();

public:
	void cal_spheres_adjacency();//计算球的邻接关系

private:
	// 计算能量的函数
	void cal_QEM_matrix(const std::vector<int>& cluster, Eigen::Matrix4d& A, Eigen::Vector4d& b, double& c);
	void cal_ps_as(const std::vector<int>& cluster, Eigen::MatrixXd& ps, Eigen::VectorXd& as);
	double E_SQEM(const Eigen::Matrix4d A, const Eigen::Vector4d b, const double c, const Eigen::Vector4d s);
	double E_euclidean(const Eigen::MatrixXd ps, const Eigen::VectorXd as, const Eigen::Vector4d s);
	//更新一个球体的位置和半径 以及 能量
	void updata_single_sphere(std::shared_ptr<Sphere> sphere, const double eps);
	double line_search(double a, double b, std::function<double(double)> func);  
	
	/// 收缩球算法
	Eigen::Vector4d shringking_ball(const int v);
	Eigen::Vector4d cal_closest_point(const Eigen::Vector4d c);
	double compute_radius(const Eigen::Vector4d p, const Eigen::Vector4d q, const Eigen::Vector4d n);

private:
	Mesh mesh;
	Eigen::MatrixXd point_pos; // n*4的向量 w=0
	Eigen::MatrixXd point_n;// n*4的法向量,w=1
	Eigen::VectorXd point_area;//n*1每个点周围三角形的三分之一面积求和
	std::vector<std::shared_ptr<Sphere>> spheres;//球的指针
	Eigen::MatrixXi spheres_adjacency;//球的邻接关系
	double lambda;//系数
	double threshold;//能量阈值
	
private:
};






#endif // !VSFAST_H