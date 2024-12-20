#include<VSfast.h>

void VSfast::init(Mesh mesh, double lambda)
{
	this->mesh = mesh;
	this->lambda = lambda; 
	this->threshold1 = 1e-3;
	this->threshold2 = 1.0;
	int count_v = mesh.n_vertices();
	std::vector<int> init_cluster;
	point_pos = Eigen::MatrixXd(count_v, 4);
	point_n = Eigen::MatrixXd(count_v, 4);
	point_area = Eigen::VectorXd(count_v, 1);
	for (auto v:mesh.vertices())
	{
		auto v_p = mesh.point(v);
		point_pos.row(v.idx()) = Eigen::Vector4d(v_p[0], v_p[1], v_p[2], 0);
		auto v_n = mesh.calc_vertex_normal(v);
		point_n.row(v.idx()) = Eigen::Vector4d(v_n[0], v_n[1], v_n[2],1);
		//std::cout << point_pos.row(v.idx()) << std::endl;
		
		double face_area = 0;
		for (auto f : v.faces())
		{
			face_area+=mesh.calc_face_area(f);
		}
		point_area(v.idx()) = face_area / 3;

		init_cluster.push_back(v.idx());
	}
	
	//添加邻接矩阵

	point_adjacency = Eigen::SparseMatrix<int>(count_v, count_v);
	for (auto e : mesh.edges())
	{
		point_adjacency.insert(e.v0().idx(), e.v1().idx()) = 1;
		point_adjacency.insert(e.v1().idx(), e.v0().idx()) = 1;
	}

	Eigen::Vector4d init_s = {1,1,1,1};
	std::shared_ptr<Sphere> s_ptr = std::make_shared<Sphere>(init_s, init_cluster, 0);
	updata_single_sphere(s_ptr, 1e-5);
}

void VSfast::updata_spheres_s_E()
{
	for (auto sphere : this->spheres)
	{
		updata_single_sphere(sphere, 1e-5);
	}

	std::sort(spheres.begin(), spheres.end(), [](const std::shared_ptr<Sphere>& s1, const std::shared_ptr<Sphere>& s2) {
		return s1->E > s2->E;  // 按能量降序排序
		});
}

void VSfast::updata_spheres_cluster()
{
	int n = spheres.size();
	int m = point_pos.rows();

	for (auto sphere_ptr : spheres) sphere_ptr->cluster.clear();
	
	for (int i = 0; i < m; i++)
	{
		double Evs = LONG_MAX;
		int s_id = 0;
		for (int j = 0; j < n; j++)
		{
			double Evs_tmp = Qvs(i, spheres[j]->s) + lambda * Dvs(i, spheres[j]->s);
			if (Evs_tmp < Evs)
			{
				Evs = Evs_tmp;
				s_id = j;
			}
		}
		spheres[s_id]->cluster.push_back(i);
	}
}

void VSfast::split_spheres()
{
	int n = spheres.size();
	cal_spheres_adjacency();
	std::set<int> deleted_sphere_i;
	for (int i = 0; i < n; i++)
	{
		if (spheres[i]->E > threshold2)
		{
			if (deleted_sphere_i.find(i) != deleted_sphere_i.end())
			{
				int max_id = spheres[i]->cluster[0];
				double max_e = 0;
				for (auto v : spheres[i]->cluster)
				{
					double energy=Dvs(v, spheres[i]->s) + Qvs(v, spheres[i]->s);
					if (max_e < energy)
					{
						max_id = v;
						max_e = energy;
					}
				}
				//分裂,并且把其相邻的sphere加入deleted_sphere_i
				
				//分裂
				spheres.push_back(std::make_shared<Sphere>( shringking_ball(max_id),std::vector<int>(),0));

				
				//加入
				Eigen::ArrayXi row = spheres_adjacency.row(i);
				for (int j = 0; j < row(j); j++)
				{
					if (row(j) > 0)
					{
						deleted_sphere_i.insert(j);
					}
				}
				
			}
		}
	}
}

void VSfast::correction_spheres()
{

}

void VSfast::run()
{
}

void VSfast::cal_spheres_adjacency()
{
	std::vector<int> vid_to_color(point_n.rows());
	int i = 0; 
	for (auto sph : spheres)
	{
		sph->id = i;
		for (auto v : sph->cluster)
		{
			vid_to_color[v] = i;
		}
		i++;
	}
	spheres_adjacency = Eigen::MatrixXi(i, i);
	spheres_adjacency.setZero();
	for (auto e : mesh.edges())
	{
		int v0id = e.v0().idx();
		int v1id = e.v1().idx();
		if (vid_to_color[v0id] != vid_to_color[v1id])
		{
			spheres_adjacency(vid_to_color[v0id], vid_to_color[v1id]) = 1;

			spheres_adjacency(vid_to_color[v1id], vid_to_color[v0id]) = 1;
		}
	}
}

void VSfast::cal_QEM_matrix(const std::vector<int>& cluster, Eigen::Matrix4d& A, Eigen::Vector4d& b, double& c)
{
	A = Eigen::Matrix4d::Zero();
	b = Eigen::Vector4d::Zero();
	c = 0;
	for (int v : cluster)
	{
		auto vh = mesh.vertex_handle(v);
		auto p = point_pos.row(v).transpose();
		for (auto f = mesh.vf_begin(vh); f != mesh.vf_end(vh); f++)
		{
			auto fn = mesh.calc_face_normal(*f);
			double area = mesh.calc_face_area(*f);
			Eigen::Vector4d n = { fn[0],fn[1],fn[2],1 };

			A += area / 3.0 * 2 * n * n.transpose();
			b += area / 3.0 * 2 * n.dot(p) * n;
			c += area / 3.0 * n.dot(p) * n.dot(p);
		}
	}
	std::cout << "c: " << c << std::endl;
}

void VSfast::cal_ps_as(const std::vector<int>& cluster, Eigen::MatrixXd& ps, Eigen::VectorXd& as)
{
	int m = cluster.size();
	ps.resize(m, 4);
	as.resize(m);
	int i = 0;
	for (auto v:cluster)
	{
		ps.row(i) = point_pos.row(v);
		as(i) = point_area(v);
		i++;
	}
}

double VSfast::Dvs(const int v, const Eigen::Vector4d s)
{
	Eigen::Vector4d p = point_pos.row(v).transpose();
	Eigen::Vector4d n = point_n.row(v).transpose();
	Eigen::Vector4d q = { s.x(),s.y(),s.z(),0 };
	double r = s.w();

	double dvs = (p - q).norm() - r;

	return point_area(v) * dvs * dvs;
}

double VSfast::Qvs(const int v, const Eigen::Vector4d s)
{
	Eigen::Matrix4d A;
	Eigen::Vector4d b;
	double c;

	cal_QEM_matrix({ v }, A, b, c);

	return 0.5 * s.transpose() * A * s - b.dot(s) + c;
}

double VSfast::E_SQEM(const Eigen::Matrix4d A, const Eigen::Vector4d b, const double c, const Eigen::Vector4d s)
{
	return 0.5 * s.transpose() * A * s - b.dot(s) + c;
}

double VSfast::E_euclidean(const Eigen::MatrixXd ps,const Eigen::VectorXd as, const Eigen::Vector4d s)
{
	Eigen::Vector4d q = { s.x(),s.y(),s.z(), 0};
	double r = s.w();

	Eigen::MatrixXd vpq = ps.rowwise() - q.transpose();
	Eigen::VectorXd d = vpq.rowwise().norm();
	d.array() -= r;
	d = d.array().square();
	d = d.array() * as.array();
	
	return d.sum();
}

void VSfast::updata_single_sphere(std::shared_ptr<Sphere> sphere, const double eps)
{
	
	// 1.计算 A b 以及 ps as
	Eigen::Matrix4d A;
	Eigen::Vector4d b;
	double c;
	cal_QEM_matrix(sphere->cluster, A, b, c);
	Eigen::MatrixXd ps;
	Eigen::VectorXd as;
	cal_ps_as(sphere->cluster, ps, as);
	
	// 2.梯度下降
	Eigen::Vector4d s = sphere->s;
	Eigen::Vector4d q = { s.x(),s.y(),s.z(),0 };
	double r = s.w();
	double zero_eps = 1e-8;
	
	int LOOP = 1000;
	for (int i = 0; i < LOOP; i++)
	{
		// 2.1 计算梯度 g
		
		Eigen::MatrixXd vpq = -(ps.rowwise() - q.transpose());
		Eigen::VectorXd vpq_norm = vpq.rowwise().norm();
		Eigen::VectorXd d = vpq_norm.array() - r;
		
		Eigen::MatrixXd gis = vpq;
		gis.col(3) = -vpq_norm;
		for (int i = 0; i < 4; i++) gis.col(i) = gis.col(i).array() / (vpq_norm.array() + zero_eps);
		for (int i = 0; i < 4; i++) gis.col(i) = gis.col(i).array() * d.array() * as.array();
		gis *= 2 * lambda;


		
		Eigen::Vector4d g = -1 * (A * s - b + gis.colwise().sum().transpose()).normalized();
		// 2.2 线搜索算法
		auto func_E = [this, A, b, c, ps, as, s, g] (float h) {
			Eigen::Vector4d s_new = s + h * g;
			return E_SQEM(A, b, c, s_new) + lambda * E_euclidean(ps, as, s_new);
			};

		double h = line_search(0, 1, func_E);

		std::cout << "s: " << s.x() << " " << s.y() << " " << s.z() << " " << s.w() << std::endl;
		std::cout << "g: " << g.x() << " " << g.y() << " " << g.z() << " " << g.w() << std::endl;
		std::cout << "E1: " << E_SQEM(A, b, c, s) << "E2: " << E_euclidean(ps, as, s) << std::endl;
		std::cout << "h: " << h << std::endl;

		// 2.3 计算新的s
		if (h * g.norm() < eps) break;
		s += h * g;
		
	}

	sphere->s = s;
	sphere->E = (E_SQEM(A, b, c, s) + lambda * E_euclidean(ps, as, s)) / as.sum();
	
}

double VSfast::line_search(double a, double b, std::function<double(double)> func)
{
	double PHI = (1 + std::sqrt(5)) / 2.0;
	double c = b - (b - a) / PHI;
	double d = a + (b - a) / PHI;

	double tol = 1e-5;

	while (std::abs(b - a) > tol)
	{
		if (func(c) < func(d))
		{
			b = d;
		}
		else
		{
			a = c;
		}
		c = b - (b - a) / PHI;
		d = a + (b - a) / PHI;
	}
	return (a + b) / 2.0;
}

Eigen::Vector4d VSfast::shringking_ball(const int v)
{
	
	//1. 初始化
	Eigen::Vector4d p = point_pos.row(v).transpose();
	Eigen::Vector4d n = point_n.row(v).transpose();
	Eigen::Vector4d q = point_pos.row((v + 1) % point_pos.rows()).transpose();
	double r = compute_radius(p, q, n);

	// 2.迭代
	double r_new = 0.0;
	Eigen::Vector4d c; //圆心
	double tol = 1e-8;
	while (abs(r - r_new) > tol)
	{
		c = p - r * n;
		q = cal_closest_point(c);
		r_new = compute_radius(p, q, n);
	}

	return Eigen::Vector4d(c.x(), c.y(), c.z(), r);
	
}

Eigen::Vector4d VSfast::cal_closest_point(const Eigen::Vector4d c)
{
	Eigen::MatrixXd vpc = this->point_pos.rowwise() - c.transpose();
	Eigen::VectorXd ds = vpc.rowwise().norm();
	int min_index = 0;
	ds.minCoeff(&min_index);
	return point_pos.row(min_index).transpose();
}

double VSfast::compute_radius(const Eigen::Vector4d p, const Eigen::Vector4d q, const Eigen::Vector4d n)
{
	float d = (p - q).norm();
	float theta = std::acos(n.dot(p - q) / d);
	return d / std::cos(theta) / 2.0;
}
