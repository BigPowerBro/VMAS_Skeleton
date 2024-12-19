#pragma once
#include "VMAS_Skeleton.h"

void VMAS_Skeleton::Init(Mesh& mesh)
{
    input_mesh = mesh;

    //所有顶点聚为一类
    std::vector<int> init_cluster;
    for (int i = 0; i < input_mesh.n_vertices(); i++) init_cluster.push_back(i);

    // 初始球体
    Sphere init_sphere = { {0,0,0},0 };
    init_sphere = update_single_sphere(init_cluster, init_sphere, 0.2, 1e-5);
    
}

float VMAS_Skeleton::dps(Point3D p, Sphere s)
{
    return (p - s.q).norm() - s.r;
}

float VMAS_Skeleton::dpns(Plane a, Sphere s)
{
    return a.n.dot(a.p - s.q) - s.r;
}

float VMAS_Skeleton::Qpns(Plane a, Sphere s)
{
    return a.n.dot(s.q) * a.n.dot(s.q) + 2 * s.r * a.n.dot(s.q) + s.r * s.r -
        2 * a.n.dot(a.p) * (a.n.dot(s.q) + s.r) + a.n.dot(a.p) * a.n.dot(a.p);
}

float VMAS_Skeleton::Dvs(const int v, const Sphere s)
{
    assert(v >= 0 && v < input_mesh.n_vertices());

    auto v_h = input_mesh.vertex_handle(v);

    auto v_p = input_mesh.point(v_h);
    Point3D p = { (float)v_p[0],(float)v_p[1],(float)v_p[2] };
    float d = dps(p, s);

    float k = 0; //d前系数
    for (auto f = input_mesh.vf_begin(v_h); f != input_mesh.vf_end(v_h); f++)
    {
        k += input_mesh.calc_face_area(*f);
    }
    
    return 1 / 3.0 * k * d * d;
}

float VMAS_Skeleton::Qvs(const int v,const Sphere s)
{
    assert(v >= 0 && v < input_mesh.n_vertices());

    auto v_h = input_mesh.vertex_handle(v);

    auto v_p = input_mesh.point(v_h);
    Point3D p = { (float)v_p[0],(float)v_p[1],(float)v_p[2] };

    float qem = 0;
    for (auto f = input_mesh.vf_begin(v_h); f != input_mesh.vf_end(v_h); f++)
    {
        auto f_n = input_mesh.calc_face_normal(*f).normalize();
        Vector3D n = { (float)f_n[0],(float)f_n[1],(float)f_n[2] };

        float A = input_mesh.calc_face_area(*f);

        Plane plane = { p,n };

        qem += 1 / 3.0 * A * Qpns(plane, s);
    }

    return qem;
}

float VMAS_Skeleton::sign(float a)
{
    return a >= 0 ? 1 : -1;
}

float VMAS_Skeleton::E_SQEM(const std::vector<int> cluster, const Sphere sphere)
{
    float E = 0;
    for (auto v : cluster) E += Qvs(v, sphere);
    return E;
}

float VMAS_Skeleton::E_euclidean(const std::vector<int> cluster, const Sphere sphere)
{
    float E = 0;
    for (auto v : cluster) E += Dvs(v, sphere);
    return E;
}

void VMAS_Skeleton::update_sphere()
{
    assert(sphere_list.size() == clusters.size()); //保证球个数和类别个数一致
}

Sphere VMAS_Skeleton::update_single_sphere(const std::vector<int>& cluster, const Sphere init_sphere, float lambda, const float eps)
{
    int n = cluster.size();
    Eigen::Vector4f s(init_sphere.q.x(), init_sphere.q.y(), init_sphere.q.z(), init_sphere.r);
    Sphere sphere = init_sphere;
    
    
    //1. 计算 A b k
    Eigen::Matrix4f A = Eigen::Matrix4f::Zero(), A_ij;
    Eigen::Vector4f b = Eigen::Vector4f::Zero(), b_ij;
    Eigen::VectorXf k = Eigen::VectorXf::Zero(n);
    for (int ii=0;ii<cluster.size();ii++)
    {
        int v = cluster[ii];
        auto v_h = input_mesh.vertex_handle(v);
        auto v_p = input_mesh.point(v_h);
        Point3D p = { (float)v_p[0],(float)v_p[1],(float)v_p[2] };
        for (auto f_it = input_mesh.vf_begin(v_h); f_it != input_mesh.vf_end(v_h); f_it++)
        {
            auto f_n = input_mesh.calc_face_normal(*f_it).normalized();
            float area = input_mesh.calc_face_area(*f_it);

            Eigen::Vector3f n = { (float)f_n[0],(float)f_n[1],(float)f_n[2] };
            Eigen::Matrix3f ntn = n * n.transpose();
            float ntp = n.dot(p);

            A_ij(3, 3) = 1;
            b_ij(3) = 1;
            for (int i = 0; i < 3; i++)
            {
                A_ij(i, 3) = n(i);
                A_ij(3, i) = n(i);
                b_ij(i) = n(i);
                for (int j = 0; j < 3; j++)
                {
                    A_ij(i, j) = ntn(i, j);
                }
            }
            A_ij *= 2 * area / 3.0;
            b_ij *= 2 * ntp * area / 3.0;

            A += A_ij;
            b += b_ij;

            k(ii) += area / 3.0;
        }
    }
    //std::cout << A << std::endl;
    // 
    // 2.梯度下降法 迭代求解

    int LOOP = 1000;
    float E = LONG_MAX;
    for (int ii = 0; ii < LOOP; ii++)
    {
        //2.1 计算梯度g和下降方向dk=-g
        Eigen::Vector4f g = A * s - b;
        
        for (int i = 0; i < n; i++)
        {
            int v = cluster[i];
            auto v_h = input_mesh.vertex_handle(v);
            auto v_p = input_mesh.point(v_h);
            Point3D p = { (float)v_p[0],(float)v_p[1],(float)v_p[2] };

            float d = dps(p, sphere);
            auto vpq = sphere.q - p;
            float dpq = vpq.norm();
            if (dpq == 0) dpq = 1;
            Eigen::Vector4f gi = { vpq.x(), vpq.y(), vpq.z(), -dpq };
            g += 2 * lambda * sign(d) / dpq * gi;
        }
        
        //g.normalized();
        Eigen::Vector4f dk = -g;
        

        //3.4 线搜索步长 黄金分割法
        float a = 0, b = 1;
        float PHI = (1 + std::sqrt(5)) / 2.0;
        float tol = 1e-5;
        
        float c = b - (b - a) / PHI;
        float d = a + (b - a) / PHI;

        auto func_E = [this, cluster, s, dk, lambda](float h) {
            Eigen::Vector4f s_new = s + h * dk;
            Sphere sphere = { {s_new.x(),s_new.y(),s_new.z()}, s_new.w() };
            return E_SQEM(cluster, sphere) + lambda * E_euclidean(cluster, sphere);
            };

        while (std::abs(b - a) > tol)
        {
            if (func_E(c) < func_E(d)) 
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

        float h = (a + b) / 2.0;

        if (h*dk.norm() < eps) break;
        if (func_E(h) > E) break;


        
        // 3.4 更新
        s += h* dk;
        sphere.q = { s.x(),s.y(),s.z() };
        sphere.r = s.w();
        
        E = E_SQEM(cluster, sphere) + lambda * E_euclidean(cluster, sphere);
        std::cout << "下降方向dk:" << dk.x() << " " << dk.y() << " " << dk.z() << std::endl;
        std::cout << "能量E：" << E << std::endl;
        std::cout << "球体S:" << sphere.q.x() << " " << sphere.q.y() << " " << sphere.q.z() << " " << sphere.r << std::endl;
    }
    return sphere;
}
