#ifndef VSFAST_H
#define VSFAST_H
#include<Eigen/dense>
#include<vector>
#include <MeshDef.h>
struct Sphere
{
	Eigen::Vector4d s;//����+�뾶
	std::vector<int> cluster;//��ĵ����
	double E;//����

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
	Mesh mesh;
	Eigen::MatrixXd point_pos; // n*3������
	Eigen::MatrixXd point_n;//n*4�ķ�����,w=1
	Eigen::MatrixXd point_area;//n*1ÿ������Χ�����ε�����֮һ������
	std::vector<std::shared_ptr<Sphere>> spheres;//���ָ��
	Eigen::MatrixXi spheres_adjacency;//����ڽӹ�ϵ
	double lambda;//ϵ��
	double threshold;//������ֵ
};






#endif // !VSFAST_H