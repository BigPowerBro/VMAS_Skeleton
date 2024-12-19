#ifndef VSFAST_H
#define VSFAST_H
#include<Eigen/dense>
#include<vector>
#include <MeshDef.h>
struct Sphere
{
	Eigen::Vector4d s;//球心+半径
	std::vector<int> cluster;//球的点的类
	double E;//能量

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
	Mesh mesh;
	Eigen::MatrixXd point_pos; // n*3的向量
	Eigen::MatrixXd point_n;//n*4的法向量,w=1
	Eigen::MatrixXd point_area;//n*1每个点周围三角形的三分之一面积求和
	std::vector<std::shared_ptr<Sphere>> spheres;//球的指针
	Eigen::MatrixXi spheres_adjacency;//球的邻接关系
	double lambda;//系数
	double threshold;//能量阈值
};






#endif // !VSFAST_H