#include<VSfast.h>

void VSfast::init(Mesh mesh, double lambda, double threshold)
{
	this->mesh = mesh;
	this->lambda = lambda; 
	this->threshold = threshold;
	int count_v = mesh.n_vertices();
	point_pos = Eigen::MatrixXd(count_v, 3);
	point_n = Eigen::MatrixXd(count_v, 4);
	point_area = Eigen::MatrixXd(count_v, 1);

	

	for (auto v:mesh.vertices())
	{
		auto v_p = mesh.point(v);
		point_pos.row(v.idx()) = Eigen::Vector3d(v_p[0], v_p[1], v_p[2]);
		auto v_n = mesh.calc_vertex_normal(v);
		point_n.row(v.idx()) = Eigen::Vector4d(v_n[0], v_n[1], v_n[2],1);
		//std::cout << point_pos.row(v.idx()) << std::endl;
		
		double face_num = 0; double face_area = 0;
		for (auto f : v.faces())
		{
			face_num++;
			face_area+=mesh.calc_face_area(f);
		}
		point_area(v.idx(),0) = face_area / face_num ;
	}
}
#include<VSfast.h>

double VSfast::E_SQEM(const std::vector<int>& cluster, const Eigen::Vector4d s)
{
	return 0.5 * s.transpose() * A * s - b.dot(s);
}

double VSfast::E_euclidean(const std::vector<int>& cluster, const Eigen::Vector4d s)
{

	return 0.0;
}
