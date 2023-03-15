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
	sz num_saved_mins;//���out������
	int num_steps;//��������
	fl min_rmsd;//�������Ľ��
	lshade() : num_steps(800),num_saved_mins(20), min_rmsd(1.0){ }
	void operator()(model& m, output_container& out, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator) const;
};

struct lshade_aux {
	int po_size;//��Ⱥ��ģ
	int po_size_min;//��Ⱥ��С��ģ
	int A_size;//�ⲿ�浵������
	int H_size;//�ɹ���ʷ��¼��������������Ĭ��ֵ��5
	double p;//����pbest�Ĳ���,����Ĭ��ֵ��0.1
	int dim;//�����ά�ȣ�7+degrees_of_freedom
	int NFE;//��ǰ��Ӧ�Ⱥ��������Ĵ���
	int MAX_NFE;//Ԥ���������Ӧ�Ⱥ��������Ĵ���
	int best_index; //��¼��ǰ��Ⱥ�������ŵĸ���
	History_memory HM;//�ɹ���ʷ��¼��
	std::vector<output_type> A;//�ⲿ�浵
	std::vector<output_type> nowpopulation;//��ʼ��Ⱥ
	std::vector<output_type> mutatepopulation;//���콻������Ⱥ
	lshade_aux() :  po_size(100), A_size(10), H_size(5), p(0.1), dim(7), po_size_min(10), NFE(0), MAX_NFE(1000) { }
	lshade_aux(const model& m) : H_size(5), p(0.1), NFE(0) {
		dim = m.ligand_degrees_of_freedom(0) + 6;//ֻ����ֻ��һ����������
		po_size = dim * 8;//r^N^init��Ĭ��ֵ
//		po_size = 50;
		po_size_min = 4;
		nowpopulation.resize(po_size);
		mutatepopulation.resize(po_size);
		MAX_NFE = po_size * 600;
		A_size = po_size;//r_arc��Сֵ
		A.resize(0);
		HM.H.resize(H_size);
		HM.k = 0;
	}
	void mutate_cross1(model& m, energy_cal& energy, rng& generator, const vec& corner1, const vec& corner2, const int& num_steps);
	void mutate_cross2(model& m, energy_cal& energy, rng& generator, const vec& corner1, const vec& corner2);
	void select(model& m, energy_cal& energy, rng& generator);

};








#endif
