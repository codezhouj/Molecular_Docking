#ifndef VINA_INDIVIDUAL_H
#define VINA_INDIVIDUAL_H

#include "model.h"
#include "conf.h"
#include "incrementable.h"
#include <algorithm>
#include "random.h"
#include "energy.h"


struct H_element {
	double M_CR;
	double M_F;
};

struct History_memory {
	std::vector<H_element> H;
	int k;
};//history memory

bool MyCompare(const output_type& x, const output_type& y);

struct lshade {
	sz num_saved_mins;//最后out的容量
	int num_steps;//迭代次数
	fl min_rmsd;//评估最后的结果
	lshade() : num_steps(800),num_saved_mins(20), min_rmsd(1.0){ }
	void operator()(model& m, output_container& out, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator) const;
};

struct lshade_aux {
	int po_size;//种群规模
	int po_size_min;//种群最小规模
	int A_size;//外部存档的容量
	int H_size;//成功历史记录集的容量，参数默认值是5
	double p;//决定pbest的参数,参数默认值是0.1
	int dim;//问题的维度，7+degrees_of_freedom
	int NFE;//当前适应度函数评估的次数
	int MAX_NFE;//预估的最大适应度函数评估的次数
	int best_index; //记录当前种群变异最优的个体
	History_memory HM;//成功历史记录集
	std::vector<output_type> A;//外部存档
	std::vector<output_type> nowpopulation;//初始种群
	std::vector<output_type> mutatepopulation;//变异交叉后的种群
	lshade_aux() :  po_size(100), A_size(10), H_size(5), p(0.1), dim(7), po_size_min(10), NFE(0), MAX_NFE(1000) { }
	lshade_aux(const model& m) : H_size(5), p(0.1), NFE(0) {
		dim = m.ligand_degrees_of_freedom(0) + 6;//只考虑只有一个配体的情况
		po_size = dim * 8;//r^N^init的默认值
//		po_size = 50;
		po_size_min = 4;
		nowpopulation.resize(po_size);
		mutatepopulation.resize(po_size);
		MAX_NFE = po_size * 600;
		A_size = po_size;//r_arc最小值
		A.resize(0);
		HM.H.resize(H_size);
		HM.k = 0;
	}
	void mutate_cross1(model& m, energy_cal& energy, rng& generator, const vec& corner1, const vec& corner2, const int& num_steps);
	void mutate_cross2(model& m, energy_cal& energy, rng& generator, const vec& corner1, const vec& corner2);
	void select(model& m, energy_cal& energy, rng& generator);

};








#endif
