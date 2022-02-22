#include <functional>
#include "dglsolver.h"

namespace dglsolver
{

void rk4(unsigned nIter,
		double h,
		const std::vector<double>& x0,
		std::vector<std::vector<double> >& res,
		std::function< void(const std::vector<double>&, std::vector<double>&) > rhs)
{
	assert(res.size() >= nIter + 1);

	// Dimension des Problems
	unsigned n = x0.size();

	// Schreibe Startwerte als 0.te Spalte in res
	res[0] = x0;

	// Speicher allokieren
	std::vector<double> w0(n), w1(n), w2(n), w3(n);

	for (unsigned i = 0; i < nIter; i++)
	{
		rhs(res[i], w0);
		for (unsigned j = 0; j < n; j++) w1[j] = res[i][j] + h / 2. * w0[j];
		rhs(w1, w2);
		for (unsigned j = 0; j < n; j++) w1[j] = res[i][j] + h / 2. * w2[j];
		rhs(w1, w3);
		for (unsigned j = 0; j < n; j++) w1[j] = res[i][j] + h * w3[j];
		for (unsigned j = 0; j < n; j++) res[i + 1][j] += res[i][j] + h / 6. * (w0[j] + 2 * w2[j] + 2 * w3[j]);
		rhs(w1, w0);
		for (unsigned j = 0; j < n; j++) res[i + 1][j] += h / 6. * w0[j];
	}
}

void rk4step(double h,
		const std::vector<double>& y0,
		std::vector<double>& y1,
		std::vector<double>& w0,
		std::vector<double>& w1,
		std::vector<double>& w2,
		std::vector<double>& w3,
		std::function< void(const std::vector<double>&, std::vector<double>&) > rhs)
{
	unsigned n = y0.size();
	assert(w0.size() == n);
	assert(w1.size() == n);
	assert(w2.size() == n);
	assert(w3.size() == n);
	assert(y1.size() == n);

	rhs(y0, w0);
	for (unsigned j = 0; j < n; j++) w1[j] = y0[j] + h / 2. * w0[j];
	rhs(w1, w2);
	for (unsigned j = 0; j < n; j++) w1[j] = y0[j] + h / 2. * w2[j];
	rhs(w1, w3);
	for (unsigned j = 0; j < n; j++) w1[j] = y0[j] + h * w3[j];
	for (unsigned j = 0; j < n; j++) y1[j] = y0[j] + h / 6. * (w0[j] + 2 * w2[j] + 2 * w3[j]);
	rhs(w1, w0);
	for (unsigned j = 0; j < n; j++) y1[j] += h / 6. * w0[j];
}

}
