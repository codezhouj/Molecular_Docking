#include <cmath>
#include <fstream>
#include "lshade.h"
#include "coords.h"
#include "mutate.h"
#include"tee.h"
#include"adam_cal.h"
#include <sstream>

const fl beta1 = 0.5;

const fl beta2 = 0.999;
const fl sigma = 0.00000001;
const fl lr = 0.01;
static int iters = 0;

bool MyCompare(const output_type& x, const output_type& y) {
	return x.e < y.e;
}


void MyNormalize(vec& angles) {
	while (angles[0] > 2 * pi) {
		angles[0] -= 2 * pi;
	}
	while (angles[0] < 0) {
		angles[0] += 2 * pi;
	}
	while (angles[1] > 2 * pi) {
		angles[1] -= 2 * pi;
	}
	while (angles[1] < -2 * pi) {
		angles[1] += 2 * pi;
	}
	if (angles[1] > pi && angles[1] < 2 * pi) {
		angles[1] = 2 * pi - angles[1];
	}
	if (angles[1] > -2 * pi && angles[1] < -pi) {
		angles[1] = 2 * pi + angles[1];
	}
	if (angles[1] > -pi && angles[1] < 0) {
		angles[1] = -angles[1];
	}
	while (angles[2] > 2 * pi) {
		angles[2] -= 2 * pi;
	}
	while (angles[2] < 0) {
		angles[2] += 2 * pi;
	}
	assert(angles[0] > 0 && angles[0] < 2 * pi);
	assert(angles[1] > 0 && angles[1] < pi);
	assert(angles[2] > 0 && angles[2] < 2 * pi);

}

void quaternion_to_3angles(const output_type& x, vec& three_angle) {
	double theta = 2 * acos(x.c.ligands[0].rigid.orientation.R_component_1());
	if (theta < 0.00000001) {
		three_angle[0] = 0;
		three_angle[1] = 0;
		three_angle[2] = 0;
	}
	else {
		vec quaternion_v;//v_x,v_y,v_z
		quaternion_v[0] = x.c.ligands[0].rigid.orientation.R_component_2() / sin(theta / 2); //v_x
		quaternion_v[1] = x.c.ligands[0].rigid.orientation.R_component_3() / sin(theta / 2); //v_y
		quaternion_v[2] = x.c.ligands[0].rigid.orientation.R_component_4() / sin(theta / 2); //v_z
		fl alpha = acos(quaternion_v[2]);
		fl beta = 0;
		fl sin_beta = quaternion_v[1] / sin(alpha);
		fl cos_bata = quaternion_v[0] / sin(alpha);
		if (sin_beta > 0)
			beta = acos(cos_bata);
		else if (sin_beta < 0)
			beta = 2 * pi - acos(cos_bata);
		three_angle[0] = theta;
		three_angle[1] = alpha;
		three_angle[2] = beta;
	}
	assert(three_angle[0] >= 0 && three_angle[0] <= 2 * pi);
	assert(three_angle[1] >= 0 && three_angle[1] <= pi);
	assert(three_angle[2] >= 0 && three_angle[2] <= 2 * pi);
}
void threeangles_to_quaternion(const vec& three_angle, std::vector<double>& quater) {
	if (fabs(three_angle[0] - 0.0) < 1e-6) {
		quater[0] = 0;
		quater[1] = 0;
		quater[2] = 0;
		quater[3] = 0;
	}
	else {
		quater[0] = three_angle[0] / 2;//二分之theta角
		quater[1] = sin(three_angle[1]) * cos(three_angle[2]);
		quater[2] = sin(three_angle[1]) * sin(three_angle[2]);
		quater[3] = cos(three_angle[1]);
	}
}

void increment_single_position(conf& c, change& g, fl beta1, fl beta2, fl& V_d_position, fl& S_d_position, sz p_index, fl lr, fl sigma, unsigned t) {
	fl d_position = g.ligands[0].rigid.position[p_index];
	V_d_position = (beta1 * V_d_position + (1 - beta1) * d_position) / (1 - pow(beta1, t));
	S_d_position = (beta2 * S_d_position + (1 - beta2) * d_position * d_position) / (1 - pow(beta2, t));
	c.ligands[0].rigid.position[p_index] = c.ligands[0].rigid.position[p_index] - lr * V_d_position / (sqrt(S_d_position) + sigma);
}

void increment_position(conf& c, change& g, fl beta1, fl beta2, vec& V_d_position, vec& S_d_position, fl lr, fl sigma, unsigned t) {
	vec d_position = g.ligands[0].rigid.position;

	V_d_position = beta1 * V_d_position + (1 - beta1) * d_position;
	S_d_position = beta2 * S_d_position + (1 - beta2) * elementwise_product(d_position, d_position);

	vec V_d_position_corr = V_d_position / (1 - pow(beta1, t));
	vec S_d_position_corr = S_d_position / (1 - pow(beta2, t));
	vec_sqrt(S_d_position_corr);
	S_d_position_corr += sigma;

	c.ligands[0].rigid.position = c.ligands[0].rigid.position - lr * vec_division(V_d_position_corr, S_d_position_corr);
}


void increment_orientation(conf& c, change& g, fl beta1, fl beta2, vec& V_d_orientation, vec& S_d_orientation, fl lr, fl sigma, unsigned t) {
	vec d_orientation = g.ligands[0].rigid.orientation;

	V_d_orientation = beta1 * V_d_orientation + (1 - beta1) * d_orientation;
	S_d_orientation = beta2 * S_d_orientation + (1 - beta2) * elementwise_product(d_orientation, d_orientation);

	vec V_d_orientation_corr = V_d_orientation / (1 - pow(beta1, t));
	vec S_d_orientation_corr = S_d_orientation / (1 - pow(beta2, t));
	vec_sqrt(S_d_orientation_corr);
	S_d_orientation_corr += sigma;

	vec rotation = -lr * vec_division(V_d_orientation_corr, S_d_orientation_corr);
	quaternion_increment(c.ligands[0].rigid.orientation, rotation);
}

void increment_single_torsions(conf& c, change& g, fl beta1, fl beta2, sz t_index, flv& V_d_torsions, flv& S_d_torsions, fl lr, fl sigma, unsigned t) {
	fl d_torsions = 0;
	fl torsion_;
	for (int i = t_index; i < V_d_torsions.size(); i++) {
		V_d_torsions[i] = (beta1 * V_d_torsions[i] + (1 - beta1) * g.ligands[0].torsions[i]) / (1 - pow(beta1, t));
		S_d_torsions[i] = (beta2 * S_d_torsions[i] + (1 - beta2) * g.ligands[0].torsions[i] * g.ligands[0].torsions[i]) / (1 - pow(beta2, t));
		torsion_ = lr * V_d_torsions[i] / (sqrt(S_d_torsions[i]) + sigma);
		c.ligands[0].torsions[i] -= normalized_angle(torsion_);
		normalize_angle(c.ligands[0].torsions[i]);
	}
}

void increment_torsions(conf& c, change& g, fl beta1, fl beta2, flv& V_d_torsions, flv& S_d_torsions, fl lr, fl sigma, unsigned t) {
	flv d_torsions = g.ligands[0].torsions;

	V_d_torsions = torsions_add_torsions(s_mul_torsions(beta1, V_d_torsions), s_mul_torsions(1 - beta1, d_torsions));
	S_d_torsions = torsions_add_torsions(s_mul_torsions(beta2, S_d_torsions), s_mul_torsions(1 - beta2, torsions_sqr(d_torsions)));

	flv V_d_torsions_corr = torsions_div_s(V_d_torsions, 1 - pow(beta1, t));
	flv S_d_torsions_corr = torsions_div_s(S_d_torsions, 1 - pow(beta2, t));
	torsions_sqrt(S_d_torsions_corr);
	torsions_add_s(S_d_torsions_corr, sigma);

	flv torsions_ = s_mul_torsions(lr, torsions_div_torsions(V_d_torsions_corr, S_d_torsions_corr));

	for (sz i = 0; i < torsions_.size(); i++) {
		c.ligands[0].torsions[i] -= normalized_angle(torsions_[i]);
		normalize_angle(c.ligands[0].torsions[i]);
	}
}

void increment(conf& x_new, const change& p, const fl alpha) {
	x_new.ligands[0].rigid.position += alpha * p.ligands[0].rigid.position;
	vec rotation;
	rotation = alpha * p.ligands[0].rigid.orientation;
	quaternion_increment(x_new.ligands[0].rigid.orientation, rotation);

	for (sz i = 0; i < x_new.ligands[0].torsions.size(); i++) {
		x_new.ligands[0].torsions[i] += normalized_angle(alpha * p.ligands[0].torsions[i]);
		normalize_angle(x_new.ligands[0].torsions[i]);
	}
}

void eval_d(change& p, change& g, const fl beta1, const fl beta2, vec& V_d_position, vec& S_d_position, vec& V_d_orientation, vec& S_d_orientation, unsigned t) {
	vec d_position = g.ligands[0].rigid.position;
	vec d_orientation = g.ligands[0].rigid.orientation;
	flv d_torsions = g.ligands[0].torsions;  //保存导数
	fl sigma = 0.00000001;

	V_d_position = beta1 * V_d_position + (1 - beta1) * d_position;
	S_d_position = beta2 * S_d_position + (1 - beta2) * elementwise_product(d_position, d_position);
	vec V_d_position_corr = V_d_position / (1 - pow(beta1, t));
	vec S_d_position_corr = S_d_position / (1 - pow(beta2, t));
	vec_sqrt(S_d_position_corr);
	S_d_position_corr += sigma;
	p.ligands[0].rigid.position = vec_division(V_d_position_corr, S_d_position_corr);
	vec_minus(p.ligands[0].rigid.position);

	V_d_orientation = beta1 * V_d_orientation + (1 - beta1) * d_orientation;
	S_d_orientation = beta2 * S_d_orientation + (1 - beta2) * elementwise_product(d_orientation, d_orientation);
	vec V_d_orientation_corr = V_d_orientation / (1 - pow(beta1, t));
	vec S_d_orientation_corr = S_d_orientation / (1 - pow(beta2, t));
	vec_sqrt(S_d_orientation_corr);
	S_d_orientation_corr += sigma;
	p.ligands[0].rigid.orientation = vec_division(V_d_orientation_corr, S_d_orientation_corr);
	vec_minus(p.ligands[0].rigid.orientation);
}


void lshade_aux::mutate_cross1(model& m, energy_cal& energy, rng& generator, const vec& corner1, const vec& corner2, const int& num_steps) {
	conf_size s = m.get_size();
	change g(s);
	fl best_e = max_fl;
	best_index = -1;
	int po_size_cur = nowpopulation.size();
	for (int i = 0; i < po_size_cur; i++) {
		NFE++;
		int r_i = random_int(0, H_size - 1, generator);
		nowpopulation[i].CR = random_normal(HM.H[r_i].M_CR, 0.1, generator);
		if (nowpopulation[i].CR > 1) {
			nowpopulation[i].CR = 1;
		}
		else if (nowpopulation[i].CR < 0) {
			nowpopulation[i].CR = 0;
		}
		nowpopulation[i].F = random_cauchy(HM.H[r_i].M_F, 0.05, generator);//为每个个体随机生成变异因子
		while (nowpopulation[i].F <= 0) {
			nowpopulation[i].F = random_cauchy(HM.H[r_i].M_F, 0.05, generator);
		}
		if (nowpopulation[i].F > 1) {
			nowpopulation[i].F = 1;
		}
		int pbest = random_int(0, p * po_size_cur - 1, generator);//生成最优种群的随机个体pbest
		int index1, index2;//生成不同的x_r1和x_r2个体下标
		index1 = random_int(0, po_size_cur - 1, generator);
		while (i == index1) {//当前个体不能和该个体索引相同
			index1 = random_int(0, po_size_cur - 1, generator);
		}
		index2 = random_int(0, po_size_cur + A.size() - 1, generator);
		bool isA = false;
		while (index1 == index2 || index2 == i) {//r2个体从种群和外部存档A的并集中随机选取
			index2 = random_int(0, po_size_cur + A.size() - 1, generator);
		}
		if (index2 >= po_size_cur) {//超出种群规模的部分就是外部存档A中的元素
			index2 -= po_size_cur;
			isA = true;
		}
		std::vector<double> quater;
		quater.resize(4);
		vec angle_de;
		bool is_mutate_qt = false;
		int dim_j = random_int(0, dim - 1, generator);

		for (int j = 0; j < dim; j++) {
			double p_cr = random_fl(0, 1.0, generator);
			if (p_cr < nowpopulation[i].CR || j == dim_j) { //概率正好需要交叉,且至少一个维度是原个体的数据
				if (isA) {//从外部存档中取个体
					if (j < 3) {
						mutatepopulation[i].c.ligands[0].rigid.position[j] = nowpopulation[i](j) + nowpopulation[i].F * (nowpopulation[pbest](j) - nowpopulation[i](j))
							+ nowpopulation[i].F * (nowpopulation[index1](j) - A[index2](j));
						if (mutatepopulation[i].c.ligands[0].rigid.position[j] < corner1[j]) {
							mutatepopulation[i].c.ligands[0].rigid.position[j] = corner1[j];
						}
						else if (mutatepopulation[i].c.ligands[0].rigid.position[j] > corner2[j]) {
							mutatepopulation[i].c.ligands[0].rigid.position[j] = corner2[j];
						}
					}
					else if (j >= 3 && j < 6) {
						//is_mutate_qt = true;

						mutatepopulation[i].rotor_angle[j - 3] = nowpopulation[i].rotor_angle[j - 3] + nowpopulation[i].F * (nowpopulation[pbest].rotor_angle[j - 3] - nowpopulation[i].rotor_angle[j - 3])
								    + nowpopulation[i].F * (nowpopulation[index1].rotor_angle[j - 3] - A[index2].rotor_angle[j - 3]);
						
						//MyNormalize(angle_de);
						//threeangles_to_quaternion(angle_de, quater);
					}
					else {
						mutatepopulation[i].c.ligands[0].torsions[j - 6] = nowpopulation[i](j) + nowpopulation[i].F * (nowpopulation[pbest](j) - nowpopulation[i](j))
							+ nowpopulation[i].F * (nowpopulation[index1](j) - A[index2](j));
						normalize_angle(mutatepopulation[i].c.ligands[0].torsions[j - 6]);
					}
				}
				else {//从原始种群中取个体
					if (j < 3) {
						vec V_d_position(zero_vec);
						vec S_d_position(zero_vec);
						fl t = 1;
						mutatepopulation[i].c.ligands[0].rigid.position[j] = nowpopulation[i](j) + nowpopulation[i].F * (nowpopulation[pbest](j) - nowpopulation[i](j))
							+ nowpopulation[i].F * (nowpopulation[index1](j) - nowpopulation[index2](j));
						if (mutatepopulation[i].c.ligands[0].rigid.position[j] < corner1[j]) {
							mutatepopulation[i].c.ligands[0].rigid.position[j] = corner1[j];
						}
						else if (mutatepopulation[i].c.ligands[0].rigid.position[j] > corner2[j]) {
							mutatepopulation[i].c.ligands[0].rigid.position[j] = corner2[j];
						}
					}
					else if (j >= 3 && j < 6) {
						mutatepopulation[i].rotor_angle[j - 3] = nowpopulation[i].rotor_angle[j - 3] + nowpopulation[i].F * (nowpopulation[pbest].rotor_angle[j - 3] - nowpopulation[i].rotor_angle[j - 3])
								+ nowpopulation[i].F * (nowpopulation[index1].rotor_angle[j - 3] - nowpopulation[index2].rotor_angle[j - 3]);
						
						//MyNormalize(angle_de);
						//threeangles_to_quaternion(angle_de, quater);
					}
					else {
						mutatepopulation[i].c.ligands[0].torsions[j - 6] = nowpopulation[i](j) + nowpopulation[i].F * (nowpopulation[pbest](j) - nowpopulation[i](j))
							+ nowpopulation[i].F * (nowpopulation[index1](j) - nowpopulation[index2](j));
						normalize_angle(mutatepopulation[i].c.ligands[0].torsions[j - 6]);
					}
				}
			}
			else {//不满足交叉概率，直接继承原个体的基因
				if (j < 3) {
					mutatepopulation[i].c.ligands[0].rigid.position[j] = nowpopulation[i](j);
				}
				else if (j >= 3 && j < 6) {
					//is_mutate_qt = false;
					mutatepopulation[i].rotor_angle[j - 3] = nowpopulation[i].rotor_angle[j - 3];
				}
				else {
					mutatepopulation[i].c.ligands[0].torsions[j - 6] = nowpopulation[i](j);
					normalize_angle(mutatepopulation[i].c.ligands[0].torsions[j - 6]);
				}
			}
		}
		MyNormalize(mutatepopulation[i].rotor_angle);
		threeangles_to_quaternion(mutatepopulation[i].rotor_angle, quater);
		fl cos_theta = cos(quater[0]);
		fl sin_theta = sin(quater[0]);
		qt q(cos_theta, quater[1] * sin_theta, quater[2] * sin_theta, quater[3] * sin_theta);
		mutatepopulation[i].c.ligands[0].rigid.orientation = q;
		mutatepopulation[i].e = energy(mutatepopulation[i].c, g);
		if (mutatepopulation[i].e < best_e) {
			best_e = mutatepopulation[i].e;
			best_index = i;
		}
	}
}

void lshade_aux::select(model& m, energy_cal& energy, rng& generator) {
	std::vector<double> SCR;
	std::vector<double> SF;
	std::vector<double> fitness_improvement;
	int count_zero = 0;
	int po_size_cur = nowpopulation.size();
	for (int i = 0; i < po_size_cur; i++) {
		if (nowpopulation[i].e <= mutatepopulation[i].e) {
			continue;
		}
		else {//变异个体更优
			A.push_back(nowpopulation[i]);//外部归档A保存失败个体
			SCR.push_back(nowpopulation[i].CR);
			SF.push_back(nowpopulation[i].F);
			if (SCR.back() == 0) {
				count_zero++;
			}
			fitness_improvement.push_back(fabs(nowpopulation[i].e - mutatepopulation[i].e));
			nowpopulation[i] = mutatepopulation[i];
		}
	}
	sz total_torsions = m.ligand_degrees_of_freedom(0);
	conf_size s = m.get_size();
	change g(s);

	//std::cout << "iter " << ++iters << std::endl;
	vec V_d_position(zero_vec);
	vec S_d_position(zero_vec);
	vec V_d_orientation(zero_vec);
	vec S_d_orientation(zero_vec);
	flv V_d_torsions(total_torsions, 0);
	flv S_d_torsions(total_torsions, 0);
	int t = 1;
	//std::cout << "energy:" << std::setw(15) << energy(mutatepopulation[best_index].c, g) << std::setw(15);
	for (int k = 0; k < 40; k++) {
		mutatepopulation[best_index].e = energy(mutatepopulation[best_index].c, g);
		increment_position(mutatepopulation[best_index].c, g, beta1, beta2, V_d_position, S_d_position, lr, sigma, t);
		increment_orientation(mutatepopulation[best_index].c, g, beta1, beta2, V_d_orientation, S_d_orientation, lr, sigma, t);
		increment_torsions(mutatepopulation[best_index].c, g, beta1, beta2, V_d_torsions, S_d_torsions, lr, sigma, t);
		t++;
	}
	mutatepopulation[best_index].e = energy(mutatepopulation[best_index].c, g);
	quaternion_to_3angles(mutatepopulation[best_index], mutatepopulation[best_index].rotor_angle);
	//std::cout << "energy:" << std::setw(15) << energy(mutatepopulation[best_index].c, g) << std::endl;
	if (nowpopulation[best_index].e > mutatepopulation[best_index].e) {
		nowpopulation[best_index] = mutatepopulation[best_index];
	}
	std::sort(nowpopulation.begin(), nowpopulation.end(), MyCompare);//排序
	//generation += 1;
	if (!SCR.empty() && !SF.empty()) {
		double sum = 0;
		double mean_WL_CR = 0;
		double mean_WL_F = 0;
		for (int i = 0; i < fitness_improvement.size(); i++) {
			sum += fitness_improvement[i];
		}
		double numerator_F = 0;//F分子
		double denominator_F = 0;//F分母
		double numerator_CR = 0;//CR分子
		double denominator_CR = 0;//CR分母
		for (int i = 0; i < SCR.size(); i++) {
			double Wi = fitness_improvement[i] / sum;
			numerator_F += Wi * SF[i] * SF[i];//分子
			numerator_CR += Wi * SCR[i] * SCR[i];
			denominator_F += Wi * SF[i];//分母
			denominator_CR += Wi * SCR[i];
		}
		mean_WL_F = numerator_F / denominator_F;
		if (SCR.size() == count_zero) {//terminal value
			mean_WL_CR = -1;
		}
		else {
			mean_WL_CR = numerator_CR / denominator_CR;
		}
		HM.H[HM.k].M_F = mean_WL_F;
		HM.H[HM.k].M_CR = mean_WL_CR;
		HM.k++;
		if (HM.k >= H_size) HM.k = 0;
	}
	po_size_cur = round(((po_size_min - po_size) * NFE / MAX_NFE) + po_size);
	nowpopulation.resize(po_size_cur);
	mutatepopulation.resize(po_size_cur);
	int A_size_cur = round(((po_size_min - A_size) * NFE / MAX_NFE) + A_size);//假定A的最小容量为po_size
	while (A.size() > A_size_cur) {//保证A的size不超过A_size
		A.erase(A.begin() + random_int(0, A.size() - 1, generator));
	}
}

path my_make_path(const std::string& str) {
	return path(str);
}

void my_write_all_output(model& m, const output_container& out, sz how_many, const std::string& output_name, const std::vector<std::string>& remarks) {
	if (out.size() < how_many)
		how_many = out.size();
	VINA_CHECK(how_many <= remarks.size());
	ofile f(my_make_path(output_name));
	VINA_FOR(i, how_many) {
		m.set(out[i].c);
		m.write_model(f, i + 1, remarks[i]); // so that model numbers start with 1
	}
}

void lshade::operator()(model& m, output_container& out, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator) const
{
	//output_container out_iter;
	//std::string out_iter_name = "1sqa_iter.pdbqt";
	//std::vector<std::string> remarks;
	//int many = 0;

	vec authentic_v(1000, 1000, 1000);
	vec v(10, 10, 10);
	fl best_e = max_fl;
	conf_size s = m.get_size();
	energy_cal energy(&m, &p, &ig, v);//计算能量
	lshade_aux aux(m);
	change g(s);//导数
	output_type tmp(s, 0);
	for (int i = 0; i < aux.po_size; i++) {//生成Po_Size初始种群
		tmp.c.randomize(corner1, corner2, generator);

		for (int j = 0; j < 3; j++) {
			if (j == 1) {
				tmp.rotor_angle[j] = random_fl(0, pi, generator);
			}
			else {
				tmp.rotor_angle[j] = random_fl(0, 2 * pi, generator);
			}
		}
		std::vector<double> quater;
		quater.resize(4);
		threeangles_to_quaternion(tmp.rotor_angle, quater);
		fl cos_theta = cos(quater[0]);
		fl sin_theta = sin(quater[0]);
		qt q(cos_theta, quater[1] * sin_theta, quater[2] * sin_theta, quater[3] * sin_theta);
		quaternion_normalize(q);//归一化
		tmp.c.ligands[0].rigid.orientation = q;
		tmp.e = energy(tmp.c, g);
		aux.nowpopulation[i] = tmp;
		aux.mutatepopulation[i] = tmp;
	}
	std::sort(aux.nowpopulation.begin(), aux.nowpopulation.end(), MyCompare);
	for (int i = 0; i < aux.H_size; i++) {//初始历史数据集是元素都是0.5
		aux.HM.H[i].M_CR = 0.5;
		aux.HM.H[i].M_F = 0.5;
	}
	for (int i = 0; i < num_steps; i++) {
		if (increment_me)
			++(*increment_me);
		aux.mutate_cross1(m, energy, generator, corner1, corner2, i);

		aux.select(m, energy, generator);//select函数中对nowpopulation进行了排序操作

		/*
		* 每隔20代记录当前种群中最好的构象
		if (i % 20 == 0) {
			many++;
			m.set(aux.nowpopulation[0].c);
			aux.nowpopulation[0].coords = m.get_heavy_atom_movable_coords();
			out_iter.push_back(new output_type(aux.nowpopulation[0]));
			std::string how_many_str = std::to_string(many);
			remarks.push_back("");
		}*/

		//my_write_all_output(m, out_iter, many, out_iter_name, remarks);

		if (aux.nowpopulation[0].e < best_e || out.size() < num_saved_mins) {//排序后能量最低的形变个体排在第一个
			m.set(aux.nowpopulation[0].c);//要重新设定
			aux.nowpopulation[0].coords = m.get_heavy_atom_movable_coords();
			add_to_output_container(out, aux.nowpopulation[0], min_rmsd, num_saved_mins); // 20 - max size
			if (aux.nowpopulation[0].e < best_e)
				best_e = aux.nowpopulation[0].e;
		}
	}
	//std::cout << "best:" << aux.nowpopulation[0].e << std::endl;
	VINA_CHECK(!out.empty());
	VINA_CHECK(out.front().e <= out.back().e); // make sure the sorting worked in the correct order

}
