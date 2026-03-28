#include "TIC_phi.h"

bool check_validity(model& our_model)
{
	for (int t = 0; t < our_model.T - 1; t++)
	{
		if (our_model.time_vec[t] + our_model.bandwidth_t >= our_model.time_vec[t + 1] - our_model.bandwidth_t)
			return false;
	}
	return true;
}

void init_TIC_phi(model& our_model)
{
	vector<Point> PS_vec;
	vector<double> dist_n1_vec;
	vector<Bucket> bucket_vec;

	for (int t = 0; t < our_model.T; t++)
		our_model.max_inteval_vec.push_back(our_model.time_vec[t] + our_model.bandwidth_t);

	for (int e = 0; e < our_model.m; e++)
	{
		for (int t = 0; t < our_model.T; t++)
		{
			our_model.edge_set[e].PS_new_vec.push_back(PS_vec);
			our_model.edge_set[e].multi_dist_n1_vec.push_back(dist_n1_vec);
			our_model.edge_set[e].multi_bucket_vec.push_back(bucket_vec);
		}
	}
}

void PS_new_distribution_phi(model& our_model)
{
	Point p;
	int t_index;

	for (int e_index = 0; e_index < our_model.m; e_index++)
	{
		Edge& e = our_model.edge_set[e_index];

		for (int p_index = 0; p_index < e.PS.size(); p_index++)
		{
			vector<double>::iterator position_iterator = lower_bound(our_model.max_inteval_vec.begin(), our_model.max_inteval_vec.end(), e.PS[p_index].time);
			t_index = (int)(position_iterator - our_model.max_inteval_vec.begin());
			
			if (t_index > our_model.T - 1)
				continue;

			if (e.PS[p_index].time < our_model.time_vec[t_index] - our_model.bandwidth_t)
				continue;
			
			p.dist_n1 = e.PS[p_index].dist_n1;
			p.dist_n2 = e.PS[p_index].dist_n2;
			p.edge_index = e.PS[p_index].edge_index;
			p.time = e.PS[p_index].time;

			p.aug_sum_dist_c_p_g_0_deg_0 = 0; p.aug_sum_dist_c_p_g_0_deg_1 = 0; p.aug_sum_dist_c_p_g_0_deg_2 = 0;
			p.aug_sum_dist_c_p_g_1_deg_0 = 0; p.aug_sum_dist_c_p_g_1_deg_1 = 0; p.aug_sum_dist_c_p_g_1_deg_2 = 0;
			p.aug_sum_dist_c_p_g_2_deg_0 = 0; p.aug_sum_dist_c_p_g_2_deg_1 = 0; p.aug_sum_dist_c_p_g_2_deg_2 = 0;

			p.aug_sum_dist_d_p_g_0_deg_0 = 0; p.aug_sum_dist_d_p_g_0_deg_1 = 0; p.aug_sum_dist_d_p_g_0_deg_2 = 0;
			p.aug_sum_dist_d_p_g_1_deg_0 = 0; p.aug_sum_dist_d_p_g_1_deg_1 = 0; p.aug_sum_dist_d_p_g_1_deg_2 = 0;
			p.aug_sum_dist_d_p_g_2_deg_0 = 0; p.aug_sum_dist_d_p_g_2_deg_1 = 0; p.aug_sum_dist_d_p_g_2_deg_2 = 0;

			e.PS_new_vec[t_index].push_back(p);
			e.multi_dist_n1_vec[t_index].push_back(p.dist_n1);
		}
	}
}

void PS_new_phi_clear(model& our_model)
{
	for (int e_index = 0; e_index < our_model.m; e_index++)
	{
		Edge& e = our_model.edge_set[e_index];
		for (int t_index = 0; t_index < our_model.T; t_index++)
		{
			e.PS_new_vec[t_index].clear();
			e.multi_dist_n1_vec[t_index].clear();
		}

		e.PS_new_vec.clear();
		e.multi_dist_n1_vec.clear();
	}
}

double edge_KAF_TIC_phi(model& our_model, Lixel& cur_l, int edge_index)
{
	int node_index_pi_c = our_model.edge_set[edge_index].n1;
	int node_index_pi_d = our_model.edge_set[edge_index].n2;
	int node_index_pi_a = our_model.node_index_pi_a;
	int node_index_pi_b = our_model.node_index_pi_b;
	int left_position, right_position, middle_position;
	double left_value, right_value;
	double dist_q_c, dist_q_d;
	double density_g_0, density_g_1, density_g_2;

	double edge_KA_value = 0;
	Edge& e = our_model.edge_set[edge_index];
	//vector<Point>& PS_new = e.PS_new;
	vector<Point>& PS_new = e.PS_new_vec[our_model.cur_t_index];
	if (PS_new.size() == 0)
		return 0;

	//q is in the same edge
	if (cur_l.edge_index == edge_index)
	{
		int q_position;
		double q_value;
		double a_P_left_g_0_deg_0, a_P_left_g_0_deg_1, a_P_left_g_0_deg_2;
		double a_P_left_g_1_deg_0, a_P_left_g_1_deg_1, a_P_left_g_1_deg_2;
		double a_P_left_g_2_deg_0, a_P_left_g_2_deg_1, a_P_left_g_2_deg_2;
		double a_P_right_g_0_deg_0, a_P_right_g_0_deg_1, a_P_right_g_0_deg_2;
		double a_P_right_g_1_deg_0, a_P_right_g_1_deg_1, a_P_right_g_1_deg_2;
		double a_P_right_g_2_deg_0, a_P_right_g_2_deg_1, a_P_right_g_2_deg_2;

		left_value = cur_l.dist_n1 - our_model.bandwidth_s;
		right_value = cur_l.dist_n1 + our_model.bandwidth_s;
		q_value = cur_l.dist_n1;

		left_position = bucket_find_position_phi(left_value, our_model, e);
		right_position = bucket_find_position_phi(right_value, our_model, e);
		q_position = bucket_find_position_phi(q_value, our_model, e);

		if (right_position == -1) //left, q and right positions are all -1
			return 0;

		if (q_position == -1) //left and q positions are all -1
		{
			a_P_right_g_0_deg_0 = PS_new[right_position].aug_sum_dist_c_p_g_0_deg_0; a_P_right_g_0_deg_1 = PS_new[right_position].aug_sum_dist_c_p_g_0_deg_1; a_P_right_g_0_deg_2 = PS_new[right_position].aug_sum_dist_c_p_g_0_deg_2;
			a_P_right_g_1_deg_0 = PS_new[right_position].aug_sum_dist_c_p_g_1_deg_0; a_P_right_g_1_deg_1 = PS_new[right_position].aug_sum_dist_c_p_g_1_deg_1; a_P_right_g_1_deg_2 = PS_new[right_position].aug_sum_dist_c_p_g_1_deg_2;
			a_P_right_g_2_deg_0 = PS_new[right_position].aug_sum_dist_c_p_g_2_deg_0; a_P_right_g_2_deg_1 = PS_new[right_position].aug_sum_dist_c_p_g_2_deg_1; a_P_right_g_2_deg_2 = PS_new[right_position].aug_sum_dist_c_p_g_2_deg_2;

			density_g_0 = (1 - our_model.gamma_s * cur_l.dist_n1 * cur_l.dist_n1) * a_P_right_g_0_deg_0
				+ 2 * our_model.gamma_s * cur_l.dist_n1 * a_P_right_g_0_deg_1 - our_model.gamma_s * a_P_right_g_0_deg_2;
			density_g_1 = (1 - our_model.gamma_s * cur_l.dist_n1 * cur_l.dist_n1) * a_P_right_g_1_deg_0
				+ 2 * our_model.gamma_s * cur_l.dist_n1 * a_P_right_g_1_deg_1 - our_model.gamma_s * a_P_right_g_1_deg_2;
			density_g_2 = (1 - our_model.gamma_s * cur_l.dist_n1 * cur_l.dist_n1) * a_P_right_g_2_deg_0
				+ 2 * our_model.gamma_s * cur_l.dist_n1 * a_P_right_g_2_deg_1 - our_model.gamma_s * a_P_right_g_2_deg_2;

			edge_KA_value = (1 - our_model.gamma_t * our_model.time_vec[our_model.cur_t_index] * our_model.time_vec[our_model.cur_t_index]) * density_g_0
				+ 2 * our_model.gamma_t * our_model.time_vec[our_model.cur_t_index] * density_g_1 - our_model.gamma_t * density_g_2;

			return edge_KA_value;
		}

		if (left_position == -1) //left position is -1
		{
			a_P_left_g_0_deg_0 = PS_new[q_position].aug_sum_dist_c_p_g_0_deg_0; a_P_left_g_0_deg_1 = PS_new[q_position].aug_sum_dist_c_p_g_0_deg_1; a_P_left_g_0_deg_2 = PS_new[q_position].aug_sum_dist_c_p_g_0_deg_2;
			a_P_left_g_1_deg_0 = PS_new[q_position].aug_sum_dist_c_p_g_1_deg_0; a_P_left_g_1_deg_1 = PS_new[q_position].aug_sum_dist_c_p_g_1_deg_1; a_P_left_g_1_deg_2 = PS_new[q_position].aug_sum_dist_c_p_g_1_deg_2;
			a_P_left_g_2_deg_0 = PS_new[q_position].aug_sum_dist_c_p_g_2_deg_0; a_P_left_g_2_deg_1 = PS_new[q_position].aug_sum_dist_c_p_g_2_deg_1; a_P_left_g_2_deg_2 = PS_new[q_position].aug_sum_dist_c_p_g_2_deg_2;

			a_P_right_g_0_deg_0 = PS_new[right_position].aug_sum_dist_c_p_g_0_deg_0 - PS_new[q_position].aug_sum_dist_c_p_g_0_deg_0; a_P_right_g_0_deg_1 = PS_new[right_position].aug_sum_dist_c_p_g_0_deg_1 - PS_new[q_position].aug_sum_dist_c_p_g_0_deg_1; a_P_right_g_0_deg_2 = PS_new[right_position].aug_sum_dist_c_p_g_0_deg_2 - PS_new[q_position].aug_sum_dist_c_p_g_0_deg_2;
			a_P_right_g_1_deg_0 = PS_new[right_position].aug_sum_dist_c_p_g_1_deg_0 - PS_new[q_position].aug_sum_dist_c_p_g_1_deg_0; a_P_right_g_1_deg_1 = PS_new[right_position].aug_sum_dist_c_p_g_1_deg_1 - PS_new[q_position].aug_sum_dist_c_p_g_1_deg_1; a_P_right_g_1_deg_2 = PS_new[right_position].aug_sum_dist_c_p_g_1_deg_2 - PS_new[q_position].aug_sum_dist_c_p_g_1_deg_2;
			a_P_right_g_2_deg_0 = PS_new[right_position].aug_sum_dist_c_p_g_2_deg_0 - PS_new[q_position].aug_sum_dist_c_p_g_2_deg_0; a_P_right_g_2_deg_1 = PS_new[right_position].aug_sum_dist_c_p_g_2_deg_1 - PS_new[q_position].aug_sum_dist_c_p_g_2_deg_1; a_P_right_g_2_deg_2 = PS_new[right_position].aug_sum_dist_c_p_g_2_deg_2 - PS_new[q_position].aug_sum_dist_c_p_g_2_deg_2;

			density_g_0 = (1 - our_model.gamma_s * cur_l.dist_n1 * cur_l.dist_n1) * a_P_left_g_0_deg_0
				+ 2 * our_model.gamma_s * cur_l.dist_n1 * a_P_left_g_0_deg_1 - our_model.gamma_s * a_P_left_g_0_deg_2
				+ (1 - our_model.gamma_s * cur_l.dist_n1 * cur_l.dist_n1) * a_P_right_g_0_deg_0
				+ 2 * our_model.gamma_s * cur_l.dist_n1 * a_P_right_g_0_deg_1 - our_model.gamma_s * a_P_right_g_0_deg_2;

			density_g_1 = (1 - our_model.gamma_s * cur_l.dist_n1 * cur_l.dist_n1) * a_P_left_g_1_deg_0
				+ 2 * our_model.gamma_s * cur_l.dist_n1 * a_P_left_g_1_deg_1 - our_model.gamma_s * a_P_left_g_1_deg_2
				+ (1 - our_model.gamma_s * cur_l.dist_n1 * cur_l.dist_n1) * a_P_right_g_1_deg_0
				+ 2 * our_model.gamma_s * cur_l.dist_n1 * a_P_right_g_1_deg_1 - our_model.gamma_s * a_P_right_g_1_deg_2;

			density_g_2 = (1 - our_model.gamma_s * cur_l.dist_n1 * cur_l.dist_n1) * a_P_left_g_2_deg_0
				+ 2 * our_model.gamma_s * cur_l.dist_n1 * a_P_left_g_2_deg_1 - our_model.gamma_s * a_P_left_g_2_deg_2
				+ (1 - our_model.gamma_s * cur_l.dist_n1 * cur_l.dist_n1) * a_P_right_g_2_deg_0
				+ 2 * our_model.gamma_s * cur_l.dist_n1 * a_P_right_g_2_deg_1 - our_model.gamma_s * a_P_right_g_2_deg_2;

			edge_KA_value = (1 - our_model.gamma_t * our_model.time_vec[our_model.cur_t_index] * our_model.time_vec[our_model.cur_t_index]) * density_g_0
				+ 2 * our_model.gamma_t * our_model.time_vec[our_model.cur_t_index] * density_g_1 - our_model.gamma_t * density_g_2;

			return edge_KA_value;
		}

		//left position >= 0
		a_P_left_g_0_deg_0 = PS_new[q_position].aug_sum_dist_c_p_g_0_deg_0 - PS_new[left_position].aug_sum_dist_c_p_g_0_deg_0; a_P_left_g_0_deg_1 = PS_new[q_position].aug_sum_dist_c_p_g_0_deg_1 - PS_new[left_position].aug_sum_dist_c_p_g_0_deg_1; a_P_left_g_0_deg_2 = PS_new[q_position].aug_sum_dist_c_p_g_0_deg_2 - PS_new[left_position].aug_sum_dist_c_p_g_0_deg_2;
		a_P_left_g_1_deg_0 = PS_new[q_position].aug_sum_dist_c_p_g_1_deg_0 - PS_new[left_position].aug_sum_dist_c_p_g_1_deg_0; a_P_left_g_1_deg_1 = PS_new[q_position].aug_sum_dist_c_p_g_1_deg_1 - PS_new[left_position].aug_sum_dist_c_p_g_1_deg_1; a_P_left_g_1_deg_2 = PS_new[q_position].aug_sum_dist_c_p_g_1_deg_2 - PS_new[left_position].aug_sum_dist_c_p_g_1_deg_2;
		a_P_left_g_2_deg_0 = PS_new[q_position].aug_sum_dist_c_p_g_2_deg_0 - PS_new[left_position].aug_sum_dist_c_p_g_2_deg_0; a_P_left_g_2_deg_1 = PS_new[q_position].aug_sum_dist_c_p_g_2_deg_1 - PS_new[left_position].aug_sum_dist_c_p_g_2_deg_1; a_P_left_g_2_deg_2 = PS_new[q_position].aug_sum_dist_c_p_g_2_deg_2 - PS_new[left_position].aug_sum_dist_c_p_g_2_deg_2;
		a_P_right_g_0_deg_0 = PS_new[right_position].aug_sum_dist_c_p_g_0_deg_0 - PS_new[q_position].aug_sum_dist_c_p_g_0_deg_0; a_P_right_g_0_deg_1 = PS_new[right_position].aug_sum_dist_c_p_g_0_deg_1 - PS_new[q_position].aug_sum_dist_c_p_g_0_deg_1; a_P_right_g_0_deg_2 = PS_new[right_position].aug_sum_dist_c_p_g_0_deg_2 - PS_new[q_position].aug_sum_dist_c_p_g_0_deg_2;
		a_P_right_g_1_deg_0 = PS_new[right_position].aug_sum_dist_c_p_g_1_deg_0 - PS_new[q_position].aug_sum_dist_c_p_g_1_deg_0; a_P_right_g_1_deg_1 = PS_new[right_position].aug_sum_dist_c_p_g_1_deg_1 - PS_new[q_position].aug_sum_dist_c_p_g_1_deg_1; a_P_right_g_1_deg_2 = PS_new[right_position].aug_sum_dist_c_p_g_1_deg_2 - PS_new[q_position].aug_sum_dist_c_p_g_1_deg_2;
		a_P_right_g_2_deg_0 = PS_new[right_position].aug_sum_dist_c_p_g_2_deg_0 - PS_new[q_position].aug_sum_dist_c_p_g_2_deg_0; a_P_right_g_2_deg_1 = PS_new[right_position].aug_sum_dist_c_p_g_2_deg_1 - PS_new[q_position].aug_sum_dist_c_p_g_2_deg_1; a_P_right_g_2_deg_2 = PS_new[right_position].aug_sum_dist_c_p_g_2_deg_2 - PS_new[q_position].aug_sum_dist_c_p_g_2_deg_2;

		density_g_0 = (1 - our_model.gamma_s * cur_l.dist_n1 * cur_l.dist_n1) * a_P_left_g_0_deg_0
			+ 2 * our_model.gamma_s * cur_l.dist_n1 * a_P_left_g_0_deg_1 - our_model.gamma_s * a_P_left_g_0_deg_2
			+ (1 - our_model.gamma_s * cur_l.dist_n1 * cur_l.dist_n1) * a_P_right_g_0_deg_0
			+ 2 * our_model.gamma_s * cur_l.dist_n1 * a_P_right_g_0_deg_1 - our_model.gamma_s * a_P_right_g_0_deg_2;

		density_g_1 = (1 - our_model.gamma_s * cur_l.dist_n1 * cur_l.dist_n1) * a_P_left_g_1_deg_0
			+ 2 * our_model.gamma_s * cur_l.dist_n1 * a_P_left_g_1_deg_1 - our_model.gamma_s * a_P_left_g_1_deg_2
			+ (1 - our_model.gamma_s * cur_l.dist_n1 * cur_l.dist_n1) * a_P_right_g_1_deg_0
			+ 2 * our_model.gamma_s * cur_l.dist_n1 * a_P_right_g_1_deg_1 - our_model.gamma_s * a_P_right_g_1_deg_2;

		density_g_2 = (1 - our_model.gamma_s * cur_l.dist_n1 * cur_l.dist_n1) * a_P_left_g_2_deg_0
			+ 2 * our_model.gamma_s * cur_l.dist_n1 * a_P_left_g_2_deg_1 - our_model.gamma_s * a_P_left_g_2_deg_2
			+ (1 - our_model.gamma_s * cur_l.dist_n1 * cur_l.dist_n1) * a_P_right_g_2_deg_0
			+ 2 * our_model.gamma_s * cur_l.dist_n1 * a_P_right_g_2_deg_1 - our_model.gamma_s * a_P_right_g_2_deg_2;

		edge_KA_value = (1 - our_model.gamma_t * our_model.time_vec[our_model.cur_t_index] * our_model.time_vec[our_model.cur_t_index]) * density_g_0
			+ 2 * our_model.gamma_t * our_model.time_vec[our_model.cur_t_index] * density_g_1 - our_model.gamma_t * density_g_2;

		return edge_KA_value;
	}

	//q is in different edge
	//compute the distances between q to c and q to d 
	dist_q_c = min(cur_l.dist_n1 + our_model.sp_node_vec_node_a[node_index_pi_c].cur_sp_value,
		cur_l.dist_n2 + our_model.sp_node_vec_node_b[node_index_pi_c].cur_sp_value);
	dist_q_d = min(cur_l.dist_n1 + our_model.sp_node_vec_node_a[node_index_pi_d].cur_sp_value,
		cur_l.dist_n2 + our_model.sp_node_vec_node_b[node_index_pi_d].cur_sp_value);

	//Case 1:
	if (dist_q_c > our_model.bandwidth_s && dist_q_d > our_model.bandwidth_s)
		return 0;

	//Case 2:
	if (dist_q_c <= our_model.bandwidth_s && dist_q_d > our_model.bandwidth_s)
	{
		left_value = our_model.bandwidth_s - dist_q_c;
		left_position = bucket_find_position_phi(left_value, our_model, e);

		if (left_position == -1) //left position is -1
			return 0;
		else
		{
			density_g_0 = (1 - our_model.gamma_s * dist_q_c * dist_q_c) * PS_new[left_position].aug_sum_dist_c_p_g_0_deg_0
				- 2 * our_model.gamma_s * dist_q_c * PS_new[left_position].aug_sum_dist_c_p_g_0_deg_1 - our_model.gamma_s * PS_new[left_position].aug_sum_dist_c_p_g_0_deg_2;
			density_g_1 = (1 - our_model.gamma_s * dist_q_c * dist_q_c) * PS_new[left_position].aug_sum_dist_c_p_g_1_deg_0
				- 2 * our_model.gamma_s * dist_q_c * PS_new[left_position].aug_sum_dist_c_p_g_1_deg_1 - our_model.gamma_s * PS_new[left_position].aug_sum_dist_c_p_g_1_deg_2;
			density_g_2 = (1 - our_model.gamma_s * dist_q_c * dist_q_c) * PS_new[left_position].aug_sum_dist_c_p_g_2_deg_0
				- 2 * our_model.gamma_s * dist_q_c * PS_new[left_position].aug_sum_dist_c_p_g_2_deg_1 - our_model.gamma_s * PS_new[left_position].aug_sum_dist_c_p_g_2_deg_2;

			edge_KA_value = (1 - our_model.gamma_t * our_model.time_vec[our_model.cur_t_index] * our_model.time_vec[our_model.cur_t_index]) * density_g_0
				+ 2 * our_model.gamma_t * our_model.time_vec[our_model.cur_t_index] * density_g_1 - our_model.gamma_t * density_g_2;

			return edge_KA_value;
		}
	}

	//Case 3:
	if (dist_q_c > our_model.bandwidth_s && dist_q_d <= our_model.bandwidth_s)
	{
		right_value = our_model.bandwidth_s - dist_q_d;
		right_position = bucket_find_position_phi(e.length - right_value, our_model, e) + 1;

		//if (right_position == PS.size() - 1)
		if (right_position == PS_new.size())
			return 0;
		else
		{
			density_g_0 = (1 - our_model.gamma_s * dist_q_d * dist_q_d) * PS_new[right_position].aug_sum_dist_d_p_g_0_deg_0
				- 2 * our_model.gamma_s * dist_q_d * PS_new[right_position].aug_sum_dist_d_p_g_0_deg_1 - our_model.gamma_s * PS_new[right_position].aug_sum_dist_d_p_g_0_deg_2;
			density_g_1 = (1 - our_model.gamma_s * dist_q_d * dist_q_d) * PS_new[right_position].aug_sum_dist_d_p_g_1_deg_0
				- 2 * our_model.gamma_s * dist_q_d * PS_new[right_position].aug_sum_dist_d_p_g_1_deg_1 - our_model.gamma_s * PS_new[right_position].aug_sum_dist_d_p_g_1_deg_2;
			density_g_2 = (1 - our_model.gamma_s * dist_q_d * dist_q_d) * PS_new[right_position].aug_sum_dist_d_p_g_2_deg_0
				- 2 * our_model.gamma_s * dist_q_d * PS_new[right_position].aug_sum_dist_d_p_g_2_deg_1 - our_model.gamma_s * PS_new[right_position].aug_sum_dist_d_p_g_2_deg_2;

			//edge_KA_value = (1 - our_model.gamma * dist_q_d * dist_q_d) * (PS.size() - right_position)
			//	- 2 * our_model.gamma * dist_q_d * PS[right_position].aug_sum_dist_d_p_deg_1 - our_model.gamma * PS[right_position].aug_sum_dist_d_p_deg_2;

			edge_KA_value = (1 - our_model.gamma_t * our_model.time_vec[our_model.cur_t_index] * our_model.time_vec[our_model.cur_t_index]) * density_g_0
				+ 2 * our_model.gamma_t * our_model.time_vec[our_model.cur_t_index] * density_g_1 - our_model.gamma_t * density_g_2;

			return edge_KA_value;
		}
	}

	//Case 4:
	if (dist_q_c <= our_model.bandwidth_s && dist_q_d <= our_model.bandwidth_s)
	{
		left_value = our_model.bandwidth_s - dist_q_c;
		left_position = bucket_find_position_phi(left_value, our_model, e);

		right_value = our_model.bandwidth_s - dist_q_d;
		right_position = bucket_find_position_phi(e.length - right_value, our_model, e) + 1;

		if (left_position == -1 || right_position == PS_new.size())
		{
			if (left_position == -1 && right_position == PS_new.size())
				return 0;
			if (right_position == PS_new.size())
			{
				density_g_0 = (1 - our_model.gamma_s * dist_q_c * dist_q_c) * PS_new[left_position].aug_sum_dist_c_p_g_0_deg_0
					- 2 * our_model.gamma_s * dist_q_c * PS_new[left_position].aug_sum_dist_c_p_g_0_deg_1 - our_model.gamma_s * PS_new[left_position].aug_sum_dist_c_p_g_0_deg_2;
				density_g_1 = (1 - our_model.gamma_s * dist_q_c * dist_q_c) * PS_new[left_position].aug_sum_dist_c_p_g_1_deg_0
					- 2 * our_model.gamma_s * dist_q_c * PS_new[left_position].aug_sum_dist_c_p_g_1_deg_1 - our_model.gamma_s * PS_new[left_position].aug_sum_dist_c_p_g_1_deg_2;
				density_g_2 = (1 - our_model.gamma_s * dist_q_c * dist_q_c) * PS_new[left_position].aug_sum_dist_c_p_g_2_deg_0
					- 2 * our_model.gamma_s * dist_q_c * PS_new[left_position].aug_sum_dist_c_p_g_2_deg_1 - our_model.gamma_s * PS_new[left_position].aug_sum_dist_c_p_g_2_deg_2;

				edge_KA_value = (1 - our_model.gamma_t * our_model.time_vec[our_model.cur_t_index] * our_model.time_vec[our_model.cur_t_index]) * density_g_0
					+ 2 * our_model.gamma_t * our_model.time_vec[our_model.cur_t_index] * density_g_1 - our_model.gamma_t * density_g_2;

				return edge_KA_value;
			}

			if (left_position == -1)
			{
				density_g_0 = (1 - our_model.gamma_s * dist_q_d * dist_q_d) * PS_new[right_position].aug_sum_dist_d_p_g_0_deg_0
					- 2 * our_model.gamma_s * dist_q_d * PS_new[right_position].aug_sum_dist_d_p_g_0_deg_1 - our_model.gamma_s * PS_new[right_position].aug_sum_dist_d_p_g_0_deg_2;
				density_g_1 = (1 - our_model.gamma_s * dist_q_d * dist_q_d) * PS_new[right_position].aug_sum_dist_d_p_g_1_deg_0
					- 2 * our_model.gamma_s * dist_q_d * PS_new[right_position].aug_sum_dist_d_p_g_1_deg_1 - our_model.gamma_s * PS_new[right_position].aug_sum_dist_d_p_g_1_deg_2;
				density_g_2 = (1 - our_model.gamma_s * dist_q_d * dist_q_d) * PS_new[right_position].aug_sum_dist_d_p_g_2_deg_0
					- 2 * our_model.gamma_s * dist_q_d * PS_new[right_position].aug_sum_dist_d_p_g_2_deg_1 - our_model.gamma_s * PS_new[right_position].aug_sum_dist_d_p_g_2_deg_2;

				edge_KA_value = (1 - our_model.gamma_t * our_model.time_vec[our_model.cur_t_index] * our_model.time_vec[our_model.cur_t_index]) * density_g_0
					+ 2 * our_model.gamma_t * our_model.time_vec[our_model.cur_t_index] * density_g_1 - our_model.gamma_t * density_g_2;

				return edge_KA_value;
			}
		}
		else
		{
			if (left_position < right_position)
			{
				density_g_0 = (1 - our_model.gamma_s * dist_q_c * dist_q_c) * PS_new[left_position].aug_sum_dist_c_p_g_0_deg_0
					- 2 * our_model.gamma_s * dist_q_c * PS_new[left_position].aug_sum_dist_c_p_g_0_deg_1 - our_model.gamma_s * PS_new[left_position].aug_sum_dist_c_p_g_0_deg_2
					+ (1 - our_model.gamma_s * dist_q_d * dist_q_d) * PS_new[right_position].aug_sum_dist_d_p_g_0_deg_0
					- 2 * our_model.gamma_s * dist_q_d * PS_new[right_position].aug_sum_dist_d_p_g_0_deg_1 - our_model.gamma_s * PS_new[right_position].aug_sum_dist_d_p_g_0_deg_2;

				density_g_1 = (1 - our_model.gamma_s * dist_q_c * dist_q_c) * PS_new[left_position].aug_sum_dist_c_p_g_1_deg_0
					- 2 * our_model.gamma_s * dist_q_c * PS_new[left_position].aug_sum_dist_c_p_g_1_deg_1 - our_model.gamma_s * PS_new[left_position].aug_sum_dist_c_p_g_1_deg_2
					+ (1 - our_model.gamma_s * dist_q_d * dist_q_d) * PS_new[right_position].aug_sum_dist_d_p_g_1_deg_0
					- 2 * our_model.gamma_s * dist_q_d * PS_new[right_position].aug_sum_dist_d_p_g_1_deg_1 - our_model.gamma_s * PS_new[right_position].aug_sum_dist_d_p_g_1_deg_2;

				density_g_2 = (1 - our_model.gamma_s * dist_q_c * dist_q_c) * PS_new[left_position].aug_sum_dist_c_p_g_2_deg_0
					- 2 * our_model.gamma_s * dist_q_c * PS_new[left_position].aug_sum_dist_c_p_g_2_deg_1 - our_model.gamma_s * PS_new[left_position].aug_sum_dist_c_p_g_2_deg_2
					+ (1 - our_model.gamma_s * dist_q_d * dist_q_d) * PS_new[right_position].aug_sum_dist_d_p_g_2_deg_0
					- 2 * our_model.gamma_s * dist_q_d * PS_new[right_position].aug_sum_dist_d_p_g_2_deg_1 - our_model.gamma_s * PS_new[right_position].aug_sum_dist_d_p_g_2_deg_2;

				edge_KA_value = (1 - our_model.gamma_t * our_model.time_vec[our_model.cur_t_index] * our_model.time_vec[our_model.cur_t_index]) * density_g_0
					+ 2 * our_model.gamma_t * our_model.time_vec[our_model.cur_t_index] * density_g_1 - our_model.gamma_t * density_g_2;

				return edge_KA_value;
			}
			else
			{
				middle_position = bucket_find_position_phi((dist_q_d - dist_q_c + e.length) / 2.0, our_model, e);

				if (middle_position > -1 && middle_position < (int)PS_new.size() - 1)
				{
					density_g_0 = (1 - our_model.gamma_s * dist_q_c * dist_q_c) * PS_new[middle_position].aug_sum_dist_c_p_g_0_deg_0
						- 2 * our_model.gamma_s * dist_q_c * PS_new[middle_position].aug_sum_dist_c_p_g_0_deg_1 - our_model.gamma_s * PS_new[middle_position].aug_sum_dist_c_p_g_0_deg_2
						+ (1 - our_model.gamma_s * dist_q_d * dist_q_d) * PS_new[middle_position + 1].aug_sum_dist_d_p_g_0_deg_0
						- 2 * our_model.gamma_s * dist_q_d * PS_new[middle_position + 1].aug_sum_dist_d_p_g_0_deg_1 - our_model.gamma_s * PS_new[middle_position + 1].aug_sum_dist_d_p_g_0_deg_2;

					density_g_1 = (1 - our_model.gamma_s * dist_q_c * dist_q_c) * PS_new[middle_position].aug_sum_dist_c_p_g_1_deg_0
						- 2 * our_model.gamma_s * dist_q_c * PS_new[middle_position].aug_sum_dist_c_p_g_1_deg_1 - our_model.gamma_s * PS_new[middle_position].aug_sum_dist_c_p_g_1_deg_2
						+ (1 - our_model.gamma_s * dist_q_d * dist_q_d) * PS_new[middle_position + 1].aug_sum_dist_d_p_g_1_deg_0
						- 2 * our_model.gamma_s * dist_q_d * PS_new[middle_position + 1].aug_sum_dist_d_p_g_1_deg_1 - our_model.gamma_s * PS_new[middle_position + 1].aug_sum_dist_d_p_g_1_deg_2;

					density_g_2 = (1 - our_model.gamma_s * dist_q_c * dist_q_c) * PS_new[middle_position].aug_sum_dist_c_p_g_2_deg_0
						- 2 * our_model.gamma_s * dist_q_c * PS_new[middle_position].aug_sum_dist_c_p_g_2_deg_1 - our_model.gamma_s * PS_new[middle_position].aug_sum_dist_c_p_g_2_deg_2
						+ (1 - our_model.gamma_s * dist_q_d * dist_q_d) * PS_new[middle_position + 1].aug_sum_dist_d_p_g_2_deg_0
						- 2 * our_model.gamma_s * dist_q_d * PS_new[middle_position + 1].aug_sum_dist_d_p_g_2_deg_1 - our_model.gamma_s * PS_new[middle_position + 1].aug_sum_dist_d_p_g_2_deg_2;

					edge_KA_value = (1 - our_model.gamma_t * our_model.time_vec[our_model.cur_t_index] * our_model.time_vec[our_model.cur_t_index]) * density_g_0
						+ 2 * our_model.gamma_t * our_model.time_vec[our_model.cur_t_index] * density_g_1 - our_model.gamma_t * density_g_2;

					return edge_KA_value;
				}
				else
				{
					if (middle_position == -1)
					{
						density_g_0 = (1 - our_model.gamma_s * dist_q_d * dist_q_d) * PS_new[0].aug_sum_dist_d_p_g_0_deg_0
							- 2 * our_model.gamma_s * dist_q_d * PS_new[0].aug_sum_dist_d_p_g_0_deg_1 - our_model.gamma_s * PS_new[0].aug_sum_dist_d_p_g_0_deg_2;

						density_g_1 = (1 - our_model.gamma_s * dist_q_d * dist_q_d) * PS_new[0].aug_sum_dist_d_p_g_1_deg_0
							- 2 * our_model.gamma_s * dist_q_d * PS_new[0].aug_sum_dist_d_p_g_1_deg_1 - our_model.gamma_s * PS_new[0].aug_sum_dist_d_p_g_1_deg_2;

						density_g_2 = (1 - our_model.gamma_s * dist_q_d * dist_q_d) * PS_new[0].aug_sum_dist_d_p_g_2_deg_0
							- 2 * our_model.gamma_s * dist_q_d * PS_new[0].aug_sum_dist_d_p_g_2_deg_1 - our_model.gamma_s * PS_new[0].aug_sum_dist_d_p_g_2_deg_2;

						edge_KA_value = (1 - our_model.gamma_t * our_model.time_vec[our_model.cur_t_index] * our_model.time_vec[our_model.cur_t_index]) * density_g_0
							+ 2 * our_model.gamma_t * our_model.time_vec[our_model.cur_t_index] * density_g_1 - our_model.gamma_t * density_g_2;

						return edge_KA_value;
					}
					else
					{
						if (middle_position == PS_new.size() - 1)
						{
							density_g_0 = (1 - our_model.gamma_s * dist_q_c * dist_q_c) * PS_new[PS_new.size() - 1].aug_sum_dist_c_p_g_0_deg_0
								- 2 * our_model.gamma_s * dist_q_c * PS_new[PS_new.size() - 1].aug_sum_dist_c_p_g_0_deg_1 - our_model.gamma_s * PS_new[PS_new.size() - 1].aug_sum_dist_c_p_g_0_deg_2;

							density_g_1 = (1 - our_model.gamma_s * dist_q_c * dist_q_c) * PS_new[PS_new.size() - 1].aug_sum_dist_c_p_g_1_deg_0
								- 2 * our_model.gamma_s * dist_q_c * PS_new[PS_new.size() - 1].aug_sum_dist_c_p_g_1_deg_1 - our_model.gamma_s * PS_new[PS_new.size() - 1].aug_sum_dist_c_p_g_1_deg_2;

							density_g_2 = (1 - our_model.gamma_s * dist_q_c * dist_q_c) * PS_new[PS_new.size() - 1].aug_sum_dist_c_p_g_2_deg_0
								- 2 * our_model.gamma_s * dist_q_c * PS_new[PS_new.size() - 1].aug_sum_dist_c_p_g_2_deg_1 - our_model.gamma_s * PS_new[PS_new.size() - 1].aug_sum_dist_c_p_g_2_deg_2;

							edge_KA_value = (1 - our_model.gamma_t * our_model.time_vec[our_model.cur_t_index] * our_model.time_vec[our_model.cur_t_index]) * density_g_0
								+ 2 * our_model.gamma_t * our_model.time_vec[our_model.cur_t_index] * density_g_1 - our_model.gamma_t * density_g_2;

							return edge_KA_value;
						}
					}
				}
			}
		}

		return edge_KA_value;
	}

	return 0;
}

void NTKDV_TIC_phi(model& our_model, Lixel& cur_l)
{
	double edge_KA_value;
	int e_index;

	for (int e = 0; e < our_model.access_edge_list.size(); e++)
	{
		e_index = our_model.access_edge_list[e];
		for (int t = 0; t < our_model.T; t++)
		{
			our_model.cur_t_index = t;
			edge_KA_value = edge_KAF_TIC_phi(our_model, cur_l, e_index);
			cur_l.KDE_value_vec[our_model.cur_t_index] += edge_KA_value;
		}
	}
}

void TIC_phi_algorithm(model& our_model)
{
	int prev_edge_index = -1;
	int edge_index;

	init_TIC_phi(our_model);
	PS_new_distribution_phi(our_model);
	augment_preprocess_phi(our_model);
	bucket_construction_phi(our_model);

	for (int l = 0; l < (int)our_model.lixel_set.size(); l++)
	{
		Lixel& cur_l = our_model.lixel_set[l];
		edge_index = cur_l.edge_index;

		if (edge_index != prev_edge_index)
		{
			clear_scan_edge_list(our_model);
			prev_edge_index = edge_index;
			our_model.sel_node_index = our_model.edge_set[edge_index].n1;
			our_model.node_index_pi_a = our_model.sel_node_index;
			dijkstra(our_model, cur_l);
			copy_sp_info(our_model, true);

			our_model.sel_node_index = our_model.edge_set[edge_index].n2;
			our_model.node_index_pi_b = our_model.sel_node_index;
			dijkstra(our_model, cur_l);
			copy_sp_info(our_model, false);
		}

		NTKDV_TIC_phi(our_model, cur_l);
	}

	PS_new_phi_clear(our_model);
	bucket_clear_phi(our_model);
	our_model.max_inteval_vec.clear();
}

void augment_preprocess_phi(model& our_model)
{
	double incr_dist_c_p_g_0_deg_0, incr_dist_c_p_g_0_deg_1, incr_dist_c_p_g_0_deg_2;
	double incr_dist_c_p_g_1_deg_0, incr_dist_c_p_g_1_deg_1, incr_dist_c_p_g_1_deg_2;
	double incr_dist_c_p_g_2_deg_0, incr_dist_c_p_g_2_deg_1, incr_dist_c_p_g_2_deg_2;
	double incr_dist_d_p_g_0_deg_0, incr_dist_d_p_g_0_deg_1, incr_dist_d_p_g_0_deg_2;
	double incr_dist_d_p_g_1_deg_0, incr_dist_d_p_g_1_deg_1, incr_dist_d_p_g_1_deg_2;
	double incr_dist_d_p_g_2_deg_0, incr_dist_d_p_g_2_deg_1, incr_dist_d_p_g_2_deg_2;
	int size;

	for (int t = 0; t < our_model.T; t++)
	{
		for (int e = 0; e < our_model.m; e++)
		{
			//size = our_model.edge_set[e].PS_new_vec[t].size();
			incr_dist_c_p_g_0_deg_0 = 0; incr_dist_c_p_g_0_deg_1 = 0; incr_dist_c_p_g_0_deg_2 = 0;
			incr_dist_c_p_g_1_deg_0 = 0; incr_dist_c_p_g_1_deg_1 = 0; incr_dist_c_p_g_1_deg_2 = 0;
			incr_dist_c_p_g_2_deg_0 = 0; incr_dist_c_p_g_2_deg_1 = 0; incr_dist_c_p_g_2_deg_2 = 0;
			incr_dist_d_p_g_0_deg_0 = 0; incr_dist_d_p_g_0_deg_1 = 0; incr_dist_d_p_g_0_deg_2 = 0;
			incr_dist_d_p_g_1_deg_0 = 0; incr_dist_d_p_g_1_deg_1 = 0; incr_dist_d_p_g_1_deg_2 = 0;
			incr_dist_d_p_g_2_deg_0 = 0; incr_dist_d_p_g_2_deg_1 = 0; incr_dist_d_p_g_2_deg_2 = 0;

			vector<Point>& PS_new = our_model.edge_set[e].PS_new_vec[t];
			size = PS_new.size();

			for (int i = 0; i < size; i++)
			{
				our_model.edge_set[e].multi_dist_n1_vec[t][i] = PS_new[i].dist_n1;

				incr_dist_c_p_g_0_deg_0 += 1;
				incr_dist_c_p_g_0_deg_1 += PS_new[i].dist_n1;
				incr_dist_c_p_g_0_deg_2 += PS_new[i].dist_n1 * PS_new[i].dist_n1;
				incr_dist_c_p_g_1_deg_0 += PS_new[i].time;
				incr_dist_c_p_g_1_deg_1 += PS_new[i].time * PS_new[i].dist_n1;
				incr_dist_c_p_g_1_deg_2 += PS_new[i].time * PS_new[i].dist_n1 * PS_new[i].dist_n1;
				incr_dist_c_p_g_2_deg_0 += PS_new[i].time * PS_new[i].time;
				incr_dist_c_p_g_2_deg_1 += PS_new[i].time * PS_new[i].time * PS_new[i].dist_n1;
				incr_dist_c_p_g_2_deg_2 += PS_new[i].time * PS_new[i].time * PS_new[i].dist_n1 * PS_new[i].dist_n1;

				PS_new[i].aug_sum_dist_c_p_g_0_deg_0 = incr_dist_c_p_g_0_deg_0; PS_new[i].aug_sum_dist_c_p_g_0_deg_1 = incr_dist_c_p_g_0_deg_1; PS_new[i].aug_sum_dist_c_p_g_0_deg_2 = incr_dist_c_p_g_0_deg_2;
				PS_new[i].aug_sum_dist_c_p_g_1_deg_0 = incr_dist_c_p_g_1_deg_0; PS_new[i].aug_sum_dist_c_p_g_1_deg_1 = incr_dist_c_p_g_1_deg_1; PS_new[i].aug_sum_dist_c_p_g_1_deg_2 = incr_dist_c_p_g_1_deg_2;
				PS_new[i].aug_sum_dist_c_p_g_2_deg_0 = incr_dist_c_p_g_2_deg_0; PS_new[i].aug_sum_dist_c_p_g_2_deg_1 = incr_dist_c_p_g_2_deg_1; PS_new[i].aug_sum_dist_c_p_g_2_deg_2 = incr_dist_c_p_g_2_deg_2;

				incr_dist_d_p_g_0_deg_0 += 1;
				incr_dist_d_p_g_0_deg_1 += PS_new[size - 1 - i].dist_n2;
				incr_dist_d_p_g_0_deg_2 += PS_new[size - 1 - i].dist_n2 * PS_new[size - 1 - i].dist_n2;
				incr_dist_d_p_g_1_deg_0 += PS_new[size - 1 - i].time;
				incr_dist_d_p_g_1_deg_1 += PS_new[size - 1 - i].time * PS_new[size - 1 - i].dist_n2;
				incr_dist_d_p_g_1_deg_2 += PS_new[size - 1 - i].time * PS_new[size - 1 - i].dist_n2 * PS_new[size - 1 - i].dist_n2;
				incr_dist_d_p_g_2_deg_0 += PS_new[size - 1 - i].time * PS_new[size - 1 - i].time;
				incr_dist_d_p_g_2_deg_1 += PS_new[size - 1 - i].time * PS_new[size - 1 - i].time * PS_new[size - 1 - i].dist_n2;
				incr_dist_d_p_g_2_deg_2 += PS_new[size - 1 - i].time * PS_new[size - 1 - i].time * PS_new[size - 1 - i].dist_n2 * PS_new[size - 1 - i].dist_n2;

				PS_new[size - 1 - i].aug_sum_dist_d_p_g_0_deg_0 = incr_dist_d_p_g_0_deg_0; PS_new[size - 1 - i].aug_sum_dist_d_p_g_0_deg_1 = incr_dist_d_p_g_0_deg_1; PS_new[size - 1 - i].aug_sum_dist_d_p_g_0_deg_2 = incr_dist_d_p_g_0_deg_2;
				PS_new[size - 1 - i].aug_sum_dist_d_p_g_1_deg_0 = incr_dist_d_p_g_1_deg_0; PS_new[size - 1 - i].aug_sum_dist_d_p_g_1_deg_1 = incr_dist_d_p_g_1_deg_1; PS_new[size - 1 - i].aug_sum_dist_d_p_g_1_deg_2 = incr_dist_d_p_g_1_deg_2;
				PS_new[size - 1 - i].aug_sum_dist_d_p_g_2_deg_0 = incr_dist_d_p_g_2_deg_0; PS_new[size - 1 - i].aug_sum_dist_d_p_g_2_deg_1 = incr_dist_d_p_g_2_deg_1; PS_new[size - 1 - i].aug_sum_dist_d_p_g_2_deg_2 = incr_dist_d_p_g_2_deg_2;
			}
		}
	}
}

void bucket_construction_phi(model& our_model)
{
	Edge* edge_set = our_model.edge_set;
	Bucket temp_bucket;
	int num_buckets;
	int b_index;

	for (int e = 0; e < our_model.m; e++)
	{
		num_buckets = (int)ceil(edge_set[e].length / our_model.lixel_reg_length);

		for (int t = 0; t < our_model.T; t++)
		{
			for (int b = 0; b < num_buckets; b++)
				edge_set[e].multi_bucket_vec[t].push_back(temp_bucket);

			vector<Point>& PS_new = edge_set[e].PS_new_vec[t];
			for (int i = 0; i < PS_new.size(); i++)
			{
				b_index = (int)floor(PS_new[i].dist_n1 / our_model.lixel_reg_length);
				edge_set[e].multi_bucket_vec[t][b_index].dist_n1_vec.push_back(PS_new[i].dist_n1);
			}

			for (int b = 0; b < num_buckets; b++)
			{
				if (b == num_buckets - 1)
					edge_set[e].multi_bucket_vec[t][b].bucket_length = edge_set[e].length - (num_buckets - 1) * our_model.lixel_reg_length;
				else
					edge_set[e].multi_bucket_vec[t][b].bucket_length = our_model.lixel_reg_length;

				if (b == 0)
					edge_set[e].multi_bucket_vec[t][b].num_accum_index = -1;
				else
					edge_set[e].multi_bucket_vec[t][b].num_accum_index = edge_set[e].multi_bucket_vec[t][b - 1].num_accum_index + edge_set[e].multi_bucket_vec[t][b - 1].dist_n1_vec.size();
			}
		}
	}
}

int bucket_find_position_phi(double value, model& our_model, Edge& e)
{
	int num_buckets = (int)e.multi_bucket_vec[our_model.cur_t_index].size();
	int bucket_index = max(min((int)floor(value / our_model.lixel_reg_length), num_buckets - 1), 0);
	vector<double>::iterator position_iterator = lower_bound(e.multi_bucket_vec[our_model.cur_t_index][bucket_index].dist_n1_vec.begin(), e.multi_bucket_vec[our_model.cur_t_index][bucket_index].dist_n1_vec.end(), value);
	int pos_index = e.multi_bucket_vec[our_model.cur_t_index][bucket_index].num_accum_index + (int)(position_iterator - e.multi_bucket_vec[our_model.cur_t_index][bucket_index].dist_n1_vec.begin());

	return pos_index;
}

void bucket_clear_phi(model& our_model)
{
	Edge* edge_set = our_model.edge_set;
	int num_buckets;

	for (int e = 0; e < our_model.m; e++)
	{
		for (int t = 0; t < our_model.T; t++)
		{
			num_buckets = edge_set[e].multi_bucket_vec[t].size();
			for (int b = 0; b < num_buckets; b++)
				edge_set[e].multi_bucket_vec[t][b].dist_n1_vec.clear();

			edge_set[e].multi_bucket_vec[t].clear();
		}

		edge_set[e].multi_bucket_vec.clear();
	}
}