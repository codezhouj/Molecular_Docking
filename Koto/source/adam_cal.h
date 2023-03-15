#ifndef VINA_ADAM_CAL_H
#define VINA_ADAM_CAL_H

#include"common.h"


inline void torsions_sqrt(flv& v) {
	for (sz i = 0; i < v.size(); i++)
		v[i] = sqrt(v[i]);
}

inline flv torsions_sqr(const flv v) {
	flv tmp;
	for (sz i = 0; i < v.size(); i++)
		tmp.push_back(v[i] * v[i]);
	return tmp;
}

void torsions_add(flv& res, flv aux) {
	for (sz i = 0; i < res.size(); i++)
		res[i] += aux[i];
}

inline flv torsions_div_s(const flv v, const fl s) {
	flv tmp;
	for (sz i = 0; i < v.size(); i++)
		tmp.push_back(v[i] / s);
	return tmp;
}

inline void torsions_add_s(flv& v, const fl s) {
	for (sz i = 0; i < v.size(); i++)
		v[i] = v[i] + s;
}

inline flv torsions_add_torsions(const flv v1, const flv v2) {
	flv tmp;
	for (sz i = 0; i < v1.size(); i++)
		tmp.push_back(v1[i] + v2[i]);
	return tmp;
}

inline flv torsions_div_torsions(const flv v1, const flv v2) {
	flv tmp;
	for (sz i = 0; i < v1.size(); i++)
		tmp.push_back(v1[i] / v2[i]);
	return tmp;
}

void torsions_minus(flv& v) {
	for (int i = 0; i < v.size(); i++)
		v[i] = -v[i];
}

inline flv s_mul_torsions(const fl s, const flv v) {
	flv tmp;
	for (sz i = 0; i < v.size(); i++)
		tmp.push_back(s * v[i]);
	return tmp;
}

inline void vec_sqrt(vec& v)
{
	for (sz i = 0; i < v.size(); i++)
		v[i] = sqrt(v[i]);
}

inline vec vec_division(const vec& a, const vec& b) {
	vec tmp;
	for (sz i = 0; i < a.size(); i++)
		tmp[i] = a[i] / b[i];
	return tmp;
}

void vec_minus(vec& a) {
	for (int i = 0; i < a.size(); i++)
		a[i] = -a[i];
}

fl dot(energy_cal& energy, const change& a, const change& b) {
	//sz n = energy.m->ligand_degrees_of_freedom(0);
	sz n = 6;
	fl tmp = 0;
	for (sz i = 0; i < n; i++) {
		tmp += a(i) * b(i);
	}
	return tmp;
}

#endif // !VINA_ADAM_CAL_H

