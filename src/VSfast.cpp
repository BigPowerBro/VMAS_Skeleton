#include<VSfast.h>

double VSfast::E_SQEM(const std::vector<int>& cluster, const Eigen::Vector4d s)
{
	return 0.5 * s.transpose() * A * s - b.dot(s);
}

double VSfast::E_euclidean(const std::vector<int>& cluster, const Eigen::Vector4d s)
{

	return 0.0;
}
