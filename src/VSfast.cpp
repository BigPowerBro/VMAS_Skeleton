#include<VSfast.h>
#include <cmath>
// 创建球体网格
void generateSphere(const Eigen::Vector4d& sphere, std::vector<Eigen::Vector3d>& vertices, std::vector<std::vector<int>>& faces, int latitudeCount = 20, int longitudeCount = 20) {
	double R = sphere[3];  // 半径
	Eigen::Vector3d center(sphere[0], sphere[1], sphere[2]);  // 球心
	int vertexCount = 0;

	for (int i = 0; i <= latitudeCount; ++i) {
		double theta = M_PI * i / latitudeCount;  // 纬度角
		for (int j = 0; j <= longitudeCount; ++j) {
			double phi = 2 * M_PI * j / longitudeCount;  // 经度角

			// 球面上的点
			double x = R * sin(theta) * cos(phi) + center[0];
			double y = R * sin(theta) * sin(phi) + center[1];
			double z = R * cos(theta) + center[2];

			vertices.push_back(Eigen::Vector3d(x, y, z));
			++vertexCount;

			// 创建面（三角形面片）
			if (i < latitudeCount && j < longitudeCount) {
				int current = i * (longitudeCount + 1) + j;
				int next = (i + 1) * (longitudeCount + 1) + j;
				int nextLong = (i + 1) * (longitudeCount + 1) + (j + 1) % (longitudeCount + 1);
				int currentLong = i * (longitudeCount + 1) + (j + 1) % (longitudeCount + 1);

				// 第一个三角形
				faces.push_back({ current, next, nextLong });
				// 第二个三角形
				faces.push_back({ current, nextLong, currentLong });
			}
		}
	}
}

// 保存到 .obj 文件
void saveToObj(const std::vector<Eigen::Vector4d>& spheres, const std::string& filename) {
	std::ofstream file(filename);

	if (!file.is_open()) {
		std::cerr << "无法打开文件 " << filename << std::endl;
		return;
	}

	int vertexOffset = 0;
	int faceOffset = 0;

	// 为每个球体生成网格
	for (size_t i = 0; i < spheres.size(); ++i) {
		// 写入新的物体
		file << "o Sphere" << i + 1 << std::endl;

		// 为每个球体生成顶点和面
		std::vector<Eigen::Vector3d> vertices;
		std::vector<std::vector<int>> faces;
		generateSphere(spheres[i], vertices, faces);

		// 写入顶点
		for (const auto& vertex : vertices) {
			file << "v " << vertex[0] << " " << vertex[1] << " " << vertex[2] << std::endl;
		}

		// 写入面
		for (const auto& face : faces) {
			file << "f " << face[0] + 1 + vertexOffset << " " << face[1] + 1 + vertexOffset << " " << face[2] + 1 + vertexOffset << std::endl;
		}

		// 更新偏移量
		vertexOffset += vertices.size();
		faceOffset += faces.size();
	}

	file.close();
	std::cout << "保存到 " << filename << " 完成。" << std::endl;
}


void VSfast::init(Mesh mesh, double lambda, int max_sphere_num)
{
	this->mesh = mesh;
	this->lambda = lambda; 
	this->threshold = 1e-3;
	this->max_sphere_num = max_sphere_num;

	int count_v = mesh.n_vertices();
	std::vector<int> init_cluster;
	point_pos = Eigen::MatrixXd(count_v, 4);
	point_n = Eigen::MatrixXd(count_v, 4);
	point_area = Eigen::VectorXd(count_v, 1);
	mesh.request_vertex_normals();
	mesh.request_vertex_status();
	for (auto v:mesh.vertices())
	{
		auto v_p = mesh.point(v);
		point_pos.row(v.idx()) = Eigen::Vector4d(v_p[0], v_p[1], v_p[2], 0);
		auto v_n = mesh.calc_normal(v);
		//std::cout << v_n << std::endl;
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
	std::shared_ptr<Sphere> s_ptr = std::make_shared<Sphere>(0, init_s, init_cluster, 0);
	spheres.push_back(s_ptr);
	update_single_sphere(s_ptr, 1e-5);
	annclass = std::make_shared<ANNClass>(point_pos.leftCols(3));
}

void VSfast::update_spheres_s_E()
{
	for (auto sphere : this->spheres)
	{
		update_single_sphere(sphere, 1e-5);
		std::cout << sphere->E << " ";
	}
	std::cout << std::endl;

	std::sort(spheres.begin(), spheres.end(), [](const std::shared_ptr<Sphere>& s1, const std::shared_ptr<Sphere>& s2) {
		return s1->E > s2->E;  // 按能量降序排序
		});
}

void VSfast::update_spheres_cluster()
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

void VSfast::update_spheres()
{
	clock_t start, end;
	double E = 0;
	double E_new = -10000;
	do {
		start = clock();
		update_spheres_cluster();
		end = clock();
		std::cout << "聚类用时：" << (double)(end - start) / CLOCKS_PER_SEC << std::endl;

		start = clock();
		update_spheres_s_E();
		end = clock();
		std::cout << "优化能量用时：" << (double)(end - start) / CLOCKS_PER_SEC << std::endl;

		while (true)
		{
			if (isnan(spheres[0]->E)) spheres.erase(spheres.begin());
			else break;
		}

		E_new = E;
		E = 0;
		for (auto sphere : spheres) E += sphere->E / sphere->cluster.size();

	} while (abs(E - E_new) > this->threshold);

	
}

void VSfast::split_spheres()
{
	int n = spheres.size();
	cal_spheres_adjacency();

	std::set<int> deleted_sphere_i;
	double threshold = spheres[n / 2]->E;
	for (int i = 0; i < n; i++)
	{
		if (spheres[i]->E >= threshold)
		{
			if (deleted_sphere_i.find(i) == deleted_sphere_i.end())
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
				spheres.push_back(std::make_shared<Sphere>(spheres.size(), shringking_ball_ann(max_id), std::vector<int>(), 0));
				
				//加入
				Eigen::ArrayXi row = spheres_adjacency.row(i);
				for (int j = 0; j < row.rows(); j++)
				{
					if (row(j) > 0)
					{
						deleted_sphere_i.insert(j);
					}
				}
				
			}
		}
	}
	std::cout << "spheres.size(): " << spheres.size() << std::endl;
}

void VSfast::correction_spheres()
{
	for (auto sphere : spheres)
	{
		auto s = sphere->s;
		auto p3 = annclass->AnnSearch(s.head<3>());
		Eigen::Vector4d p = { p3.x(),p3.y(),p3.z(), 0 };
		
		Eigen::VectorXd d = (point_pos.rowwise() - p.transpose()).rowwise().norm();
		int index; d.minCoeff(&index);
		Eigen::Vector4d n = point_n.row(index).transpose();
		
		//Eigen::Vector4d n = p - s;
		sphere->s = shringking_ball(p, n);
	}
}

void VSfast::run()
{
	int i = 0;
	while (spheres.size() < max_sphere_num) // 当球的个数大于max_sphere_num时 结束
	{
		split_spheres();
		clock_t start, end;
		double E = 0;
		double E_new = -10000;
		do {
			start = clock();
			update_spheres_cluster();
			end = clock();
			std::cout << "聚类用时：" << (double)(end - start) / CLOCKS_PER_SEC << std::endl;

			start = clock();
			update_spheres_s_E();
			end = clock();
			std::cout << "优化能量用时：" << (double)(end - start) / CLOCKS_PER_SEC << std::endl;

			while (true)
			{
				if (isnan(spheres[0]->E)) spheres.erase(spheres.begin());
				else break;
			}

			E_new = E;
			E = 0;
			for (auto sphere : spheres) E = std::max(E,sphere->E);

			std::vector<Eigen::Vector4d> S_vec;
			for (auto sphere : spheres) S_vec.push_back(sphere->s);
			saveToObj(S_vec, "sphere" + std::to_string(i) + ".obj");
			i++;

		} while (abs(E - E_new) > this->threshold);
		correction_spheres();
		
		
	}

	
	/*
	Eigen::Vector4d p = point_pos.row(1200).transpose();
	Eigen::Vector4d n = point_n.row(1200).transpose();
	auto S = shringking_ball(p, n);
	spheres.clear();
	spheres.push_back(std::make_shared<Sphere>(0, S, std::vector<int>(), 0));
	*/

}

void VSfast::cal_spheres_adjacency()
{
	spheres_adjacency_map.clear();
	std::vector<int> vid_to_color(point_n.rows());
	int i = 0; 
	for (auto sph : spheres)
	{
		sph->id = i;
		for (auto v : sph->cluster)
		{
			vid_to_color[v] = i;
		}
		spheres_adjacency_map[i] = std::set<int>();
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

			spheres_adjacency_map[vid_to_color[v0id]].insert(vid_to_color[v1id]);
			spheres_adjacency_map[vid_to_color[v1id]].insert(vid_to_color[v0id]);
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

void VSfast::update_single_sphere(std::shared_ptr<Sphere> sphere, const double eps)
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


		
		Eigen::Vector4d g = -1 * (A * s - b + gis.colwise().sum().transpose());
		// 2.2 线搜索算法
		auto func_E = [this, A, b, c, ps, as, s, g] (float h) {
			Eigen::Vector4d s_new = s + h * g;
			return E_SQEM(A, b, c, s_new) + lambda * E_euclidean(ps, as, s_new);
			};

		double h = line_search(0, 1, func_E);
		//double h = 0.1;
		//std::cout << "s: " << s.x() << " " << s.y() << " " << s.z() << " " << s.w() << std::endl;
		//std::cout << "g: " << g.x() << " " << g.y() << " " << g.z() << " " << g.w() << std::endl;
		//std::cout << "E1: " << E_SQEM(A, b, c, s) << "E2: " << E_euclidean(ps, as, s) << std::endl;
		//std::cout << "h: " << h << std::endl;

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

	double tol = 1e-4;

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

Eigen::Vector4d VSfast::shringking_ball(Eigen::Vector4d p, const Eigen::Vector4d n)
{
	//1. 初始化
	Eigen::Vector3d p3 = p.head<3>();
	Eigen::Vector3d n3 = n.head<3>();

	double r = 0;
	Eigen::Vector3d q;
	for (int i = 0; i < point_pos.size(); i++)
	{
		q = point_pos.row(i).transpose().head<3>();
		if (q == p3) continue;
		r = compute_radius(p3, q, n3);
		if (r > 0) break;
	}
	
	// 2.迭代
	double r_new = 0.0;
	Eigen::Vector3d c; //圆心
	double tol = 1e-5;
	while (true)
	{
		c = p3 - r * n3;
		q = cal_closest_point(c, p3);
		r_new = compute_radius(p3, q, n3);
		if (abs(r - r_new) <= tol) break;
		r = r_new;
	}
	return Eigen::Vector4d(c.x(), c.y(), c.z(), r);
}

double VSfast::compute_radius(const Eigen::Vector3d p, const Eigen::Vector3d q, const Eigen::Vector3d n)
{
	double d = (p - q).norm();
	return d * d / n.dot(p - q) / 2.0;
}

Eigen::Vector3d VSfast::cal_closest_point(const Eigen::Vector3d c, const Eigen::Vector3d p)
{
	Eigen::	MatrixXd ps = point_pos.block(0, 0, point_pos.rows(), 3);
	Eigen::VectorXd dp = (ps.rowwise() - p.transpose()).rowwise().norm();
	int pi; dp.minCoeff(&pi);
	Eigen::VectorXd dc = (ps.rowwise() - c.transpose()).rowwise().norm();
	dc(pi) = LONG_MAX;
	int ci; dc.minCoeff(&ci);
	return point_pos.row(ci).transpose().head<3>();
}

Eigen::Vector4d VSfast::shringking_ball_ann(const int v)
{
	Eigen::Vector3d v_p = point_pos.row(v).head<3>();
	Eigen::Vector3d v_n = point_n.row(v).head<3>();
	//std::cout << point_n << std::endl;
	double theta = 0.9;
	double r = 1.0;
	Eigen::Vector3d c;
	while (true)
	{
		c = v_p - r * v_n;
		Eigen::Vector3d miniP= annclass->AnnSearch(c);

		if ((miniP - c).norm() <= r)
		{
			r = r * theta;
		}
		else
		{
			break;
		}
	}

	return Eigen::Vector4d(c[0], c[1], c[2], r);
}

void VSfast::write_color_obj(std::string fname)
{
	std::ofstream out(fname+std::string(".obj"));
	std::ofstream outply(fname + std::string(".ply"));

	std::vector<int> vid_to_color(point_n.rows());
	std::vector<Eigen::Vector3d> color_bar;
	int i = 0;
	for (auto sph : spheres)
	{
		sph->id = i;
		for (auto v : sph->cluster)
		{
			vid_to_color[v] = i;
		}
		spheres_adjacency_map[i] = std::set<int>();
		i++;
	}

	for (int j = 0; j < i; j++)
	{
		color_bar.push_back(Eigen::Vector3d(j / double(i), j / double(i), j / double(i)));
	}
	for (auto v : mesh.vertices())
	{
		auto p = mesh.point(v);
		out << "v " << p[0] << " " << p[1] << " " << p[2] << " " << color_bar[vid_to_color[v.idx()]].transpose() << std::endl;
	}
	for (auto f : mesh.faces())
	{
		out << "f ";
		for (auto v : f.vertices())
		{
			out << v.idx()+1 << " ";
		}
		out << std::endl;
	}
	cal_spheres_adjacency();
	int v_size = spheres.size();
	std::vector<std::pair<int, int>> edges;
	for (int i = 0; i < v_size; i++)
	{
		for (auto j: spheres_adjacency_map[i])
		{
			if (j > i)
			{
				edges.push_back(std::pair<int, int>(i, j) );
			}
		}
	}
	outply << "ply " << std::endl;
	outply << "format ascii 1.0 " << std::endl;
	outply << "element vertex " <<v_size<<" " << std::endl;
	outply << "property float x " << std::endl;
	outply << "property float y " << std::endl;
	outply << "property float z " << std::endl;
	outply << "element edge "<< edges.size()<<" " << std::endl;
	//outply << "property list uchar int vertex_indices " << std::endl;
	
		
	outply << "property int vertex1 " << std::endl;
	outply << "property int vertex2 " << std::endl;
	outply << "end_header " << std::endl;
	std::vector<Eigen::Vector4d> spheres_s;
	for (int i = 0; i < v_size; i++)
	{
		outply << spheres[i]->s[0] << " " << spheres[i]->s[1] << " " << spheres[i]->s[2] << std::endl;
		spheres_s.push_back(spheres[i]->s);
	}
	for (int i = 0; i < edges.size(); i++)
	{
		outply <<edges[i].first << " " << edges[i].second << std::endl;
	}

	saveToObj(spheres_s, fname + "shpere.obj");
}

ANNClass::ANNClass(Eigen::MatrixXd LPR)
{
	maxPts = LPR.rows();

	std::cout << "ANN_KDTREE" << " SIZE:" << maxPts << std::endl;
	dataPts = annAllocPts(maxPts, dim);
	for (int i = 0; i < maxPts; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			*(dataPts[i] + j) = LPR(i, j);
			//std::cout << *(dataPts[i] + j) << " ";
			//std::cout << Target_points[i]->pos[j] << " ";
		}
		//std::cout << std::endl;
	}
	kdTree = new ANNkd_tree(					// build search structure
		dataPts,					// the data points
		maxPts,						// number of points
		dim);						// dimension of space
}

Eigen::Vector3d ANNClass::AnnSearch(Eigen::Vector3d vIn)
{
	Eigen::Vector3d& lim_p = vIn;
	//limit_point(point, lim_p);
	queryPt[0] = lim_p[0];
	queryPt[1] = lim_p[1];
	queryPt[2] = lim_p[2];
	kdTree->annkSearch(						// search
		queryPt,						// query point
		k,								// number of near neighbors
		nnIdx,							// nearest neighbors (returned)
		dists,							// distance (returned)
		eps);							// error bound
	Eigen::Vector3d sum;
	sum.setZero();

	for (int i = 0; i < k; i++)
	{
		sum[0] += *(dataPts[nnIdx[i]] + 0);
		sum[1] += *(dataPts[nnIdx[i]] + 1);
		sum[2] += *(dataPts[nnIdx[i]] + 2);
	}
	return  sum / double(k);
}
