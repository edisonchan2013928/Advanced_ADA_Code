#include "KAF.h"

double s_kernel_value(model& our_model, double s_dist)
{
	if (s_dist >= our_model.bandwidth_s)
		return 0;

	if (our_model.k_type == 2) //Epanechnikov kernel
		return (1 - (s_dist * s_dist) / (our_model.bandwidth_s* our_model.bandwidth_s));
	if (our_model.k_type == 3) //Quartic kernel
		return (1 - (s_dist * s_dist) / (our_model.bandwidth_s * our_model.bandwidth_s)) * (1 - (s_dist * s_dist) / (our_model.bandwidth_s * our_model.bandwidth_s));
	if (our_model.k_type == 4) //Uniform kernel
		return (1 / our_model.bandwidth_s);
	if (our_model.k_type == 0) //Gaussian kernel
		return exp(-s_dist * s_dist / (our_model.bandwidth_s * our_model.bandwidth_s));

	return -1;
}

double t_kernel_value(model& our_model, double t_dist)
{
	if (t_dist >= our_model.bandwidth_t)
		return 0;

	if (our_model.k_type == 2) //Epanechnikov kernel
		return (1 - (t_dist * t_dist) / (our_model.bandwidth_t * our_model.bandwidth_t));
	if (our_model.k_type == 3) //Quartic kernel
		return (1 - (t_dist * t_dist) / (our_model.bandwidth_t * our_model.bandwidth_t)) * (1 - (t_dist * t_dist) / (our_model.bandwidth_t * our_model.bandwidth_t));
	if (our_model.k_type == 4) //Uniform kernel
		return (1 / our_model.bandwidth_t);
	if (our_model.k_type == 0) //Gaussian kernel
		return exp(-t_dist * t_dist / (our_model.bandwidth_t * our_model.bandwidth_t));

	return -1;
}

double edge_KAF(model& our_model, Lixel& cur_l, int edge_index)
{
	int node_index_pi_c = our_model.edge_set[edge_index].n1;
	int node_index_pi_d = our_model.edge_set[edge_index].n2;
	int node_index_pi_a = our_model.node_index_pi_a;
	int node_index_pi_b = our_model.node_index_pi_b;
	double dist_s_q_p, dist_t_q_p;

	double edge_KA_value = 0;
	Edge& e = our_model.edge_set[edge_index];
	vector<Point>& PS = e.PS;
	if (PS.size() == 0)
		return 0;

	//basic algorithm for edge_KAF
	if (our_model.method == 1 || our_model.method == 2)
	{
		if (cur_l.edge_index == edge_index)
		{
			for (int i = 0; i < (int)our_model.edge_set[edge_index].PS.size(); i++)
			{
				dist_s_q_p = fabs(PS[i].dist_n1 - cur_l.dist_n1);
				dist_t_q_p = fabs(PS[i].time - our_model.time_vec[our_model.cur_t_index]);
				edge_KA_value += s_kernel_value(our_model, dist_s_q_p) * t_kernel_value(our_model, dist_t_q_p);
			}

			return edge_KA_value;
		}

		for (int i = 0; i < (int)our_model.edge_set[edge_index].PS.size(); i++)
		{
			dist_t_q_p = fabs(PS[i].time - our_model.time_vec[our_model.cur_t_index]);

			if (our_model.method == 1)
				dist_s_q_p = min(our_model.sp_node_vec[node_index_pi_c].cur_sp_value + PS[i].dist_n1,
					our_model.sp_node_vec[node_index_pi_d].cur_sp_value + PS[i].dist_n2);

			if (our_model.method == 2)
				dist_s_q_p = min(cur_l.dist_n1 + min(our_model.sp_node_vec_node_a[node_index_pi_c].cur_sp_value + PS[i].dist_n1, our_model.sp_node_vec_node_a[node_index_pi_d].cur_sp_value + PS[i].dist_n2),
					cur_l.dist_n2 + min(our_model.sp_node_vec_node_b[node_index_pi_c].cur_sp_value + PS[i].dist_n1, our_model.sp_node_vec_node_b[node_index_pi_d].cur_sp_value + PS[i].dist_n2));

			edge_KA_value += s_kernel_value(our_model, dist_s_q_p) * t_kernel_value(our_model, dist_t_q_p);
		}
	}

	return edge_KA_value;
}

void NTKDV_basic(model& our_model, Lixel& cur_l)
{
	double edge_KA_value;
	int e_index;

	for (int e = 0; e < our_model.access_edge_list.size(); e++)
	{
		e_index = our_model.access_edge_list[e];

		for (int t = 0; t < our_model.T; t++)
		{
			our_model.cur_t_index = t;
			edge_KA_value = edge_KAF(our_model, cur_l, e_index);
			cur_l.KDE_value_vec[t] += edge_KA_value;
		}
	}
}

