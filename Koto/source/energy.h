#ifndef VINA_ENERGY_H
#define VINA_ENERGY_H

#include "model.h"


struct energy_cal {
	model* m;//不使用指针
	const precalculate* p;
	const igrid* ig;
	const vec v;
	energy_cal(model* m_, const precalculate* p_, const igrid* ig_, const vec& v_) : m(m_), p(p_), ig(ig_), v(v_) {}
	fl operator()(const conf& c, change& g) {
		const fl tmp = m->eval_deriv(*p, *ig, v, c, g);
		return tmp;
	}
};

#endif
