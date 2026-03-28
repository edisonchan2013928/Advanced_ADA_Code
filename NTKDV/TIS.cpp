#include "TIS.h"

void merge_timestamps(model& our_model, vector<double>& t_bandwidth_small_vec, vector<double>& t_bandwidth_large_vec)
{
	int s_index = 0;
	int l_index = 0;

	while (s_index < our_model.T || l_index < our_model.T)
	{
		if (s_index == our_model.T)
		{
			our_model.timestamps_vec.push_back(t_bandwidth_large_vec[l_index]);
			l_index++;
			continue;
		}

		if (t_bandwidth_small_vec[s_index] < t_bandwidth_large_vec[l_index])
		{
			our_model.timestamps_vec.push_back(t_bandwidth_small_vec[s_index]);
			s_index++;
		}
		else
		{
			our_model.timestamps_vec.push_back(t_bandwidth_large_vec[l_index]);
			l_index++;
		}
	}
}

void construct_intervals(model& our_model)
{
	vector<double> t_bandwidth_small_vec;
	vector<double> t_bandwidth_large_vec;
	//vector<double> timestamps_vec;
	Lixel_den_g lixel_den_g;
	time_interval t_interval;

	for (int t = 0; t < our_model.T; t++)
	{
		t_bandwidth_small_vec.push_back(our_model.time_vec[t] - our_model.bandwidth_t);
		t_bandwidth_large_vec.push_back(our_model.time_vec[t] + our_model.bandwidth_t);
	}

	merge_timestamps(our_model, t_bandwidth_small_vec, t_bandwidth_large_vec);

	for (int ii = 0; ii < our_model.timestamps_vec.size() - 1; ii++)
	{
		t_interval.start_time = our_model.timestamps_vec[ii];
		t_interval.end_time = our_model.timestamps_vec[ii + 1];
		
		for (int l = 0; l < our_model.lixel_set.size(); l++)
		{
			lixel_den_g.dist_n1 = our_model.lixel_set[l].dist_n1;
			lixel_den_g.dist_n2 = our_model.lixel_set[l].dist_n2;
			lixel_den_g.edge_index = our_model.lixel_set[l].edge_index;
			lixel_den_g.density_g_0 = 0;
			lixel_den_g.density_g_1 = 0;
			lixel_den_g.density_g_2 = 0;
			t_interval.lixel_g_set.push_back(lixel_den_g);
		}

		our_model.interval_vec.push_back(t_interval);
		t_interval.lixel_g_set.clear();
	}

	//construct the slide line
	for (int l = 0; l < our_model.lixel_set.size(); l++)
	{
		lixel_den_g.dist_n1 = our_model.lixel_set[l].dist_n1;
		lixel_den_g.dist_n2 = our_model.lixel_set[l].dist_n2;
		lixel_den_g.edge_index = our_model.lixel_set[l].edge_index;
		lixel_den_g.density_g_0 = 0;
		lixel_den_g.density_g_1 = 0;
		lixel_den_g.density_g_2 = 0;
		our_model.slide_line.push_back(lixel_den_g);
	}

	our_model.num_intervals = our_model.interval_vec.size();
}

void init_TIS(model& our_model)
{
	vector<Point> PS_vec;
	vector<double> dist_n1_vec;
	vector<Bucket> bucket_vec;

	for (int e = 0; e < our_model.m; e++)
	{
		for (int i = 0; i < our_model.num_intervals; i++)
		{
			our_model.edge_set[e].PS_new_vec.push_back(PS_vec);
			our_model.edge_set[e].multi_dist_n1_vec.push_back(dist_n1_vec);
			our_model.edge_set[e].multi_bucket_vec.push_back(bucket_vec);
		}
	}
}

void PS_new_distribution_sharing(model& our_model)
{
	Point p;
	int t_index;
	int i_index;

	for (int e_index = 0; e_index < our_model.m; e_index++)
	{
		Edge& e = our_model.edge_set[e_index];

		for (int p_index = 0; p_index < e.PS.size(); p_index++)
		{
			vector<double>::iterator position_iterator = lower_bound(our_model.timestamps_vec.begin(), our_model.timestamps_vec.end(), e.PS[p_index].time);
			t_index = (int)(position_iterator - our_model.timestamps_vec.begin());

			if (t_index == 0 || t_index == our_model.timestamps_vec.size())
				continue;

			i_index = t_index - 1;

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

			e.PS_new_vec[i_index].push_back(p);
			e.multi_dist_n1_vec[i_index].push_back(p.dist_n1);
			//e.PS_new_vec[t_index].push_back(p);
			//e.multi_dist_n1_vec[t_index].push_back(p.dist_n1);
		}
	}
}

void PS_new_sharing_clear(model& our_model)
{
	for (int e_index = 0; e_index < our_model.m; e_index++)
	{
		Edge& e = our_model.edge_set[e_index];
		for (int i_index = 0; i_index < our_model.num_intervals; i_index++)
		{
			e.PS_new_vec[i_index].clear();
			e.multi_dist_n1_vec[i_index].clear();
		}

		e.PS_new_vec.clear();
		e.multi_dist_n1_vec.clear();
	}

	our_model.timestamps_vec.clear();
}

void interval_clear(model& our_model)
{
	for (int i = 0; i < our_model.num_intervals; i++)
		our_model.interval_vec[i].lixel_g_set.clear();
	our_model.interval_vec.clear();

	our_model.slide_line.clear();
}

void edge_KAF_TIS(model& our_model, int lixel_index, int interval_index, int edge_index)
{
	int node_index_pi_c = our_model.edge_set[edge_index].n1;
	int node_index_pi_d = our_model.edge_set[edge_index].n2;
	int node_index_pi_a = our_model.node_index_pi_a;
	int node_index_pi_b = our_model.node_index_pi_b;
	int left_position, right_position, middle_position;
	double left_value, right_value;
	double dist_q_c, dist_q_d;
	//double density_g_0, density_g_1, density_g_2;
	Lixel& cur_l = our_model.lixel_set[lixel_index];
	Lixel_den_g& cur_l_g = our_model.interval_vec[our_model.cur_interval_index].lixel_g_set[lixel_index];

	Edge& e = our_model.edge_set[edge_index];
	//vector<Point>& PS_new = e.PS_new;
	vector<Point>& PS_new = e.PS_new_vec[our_model.cur_interval_index];
	if (PS_new.size() == 0)
		return;

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

		left_position = bucket_find_position_sharing(left_value, our_model, e);
		right_position = bucket_find_position_sharing(right_value, our_model, e);
		q_position = bucket_find_position_sharing(q_value, our_model, e);

		if (right_position == -1) //left, q and right positions are all -1
			return;

		if (q_position == -1) //left and q positions are all -1
		{
			a_P_right_g_0_deg_0 = PS_new[right_position].aug_sum_dist_c_p_g_0_deg_0; a_P_right_g_0_deg_1 = PS_new[right_position].aug_sum_dist_c_p_g_0_deg_1; a_P_right_g_0_deg_2 = PS_new[right_position].aug_sum_dist_c_p_g_0_deg_2;
			a_P_right_g_1_deg_0 = PS_new[right_position].aug_sum_dist_c_p_g_1_deg_0; a_P_right_g_1_deg_1 = PS_new[right_position].aug_sum_dist_c_p_g_1_deg_1; a_P_right_g_1_deg_2 = PS_new[right_position].aug_sum_dist_c_p_g_1_deg_2;
			a_P_right_g_2_deg_0 = PS_new[right_position].aug_sum_dist_c_p_g_2_deg_0; a_P_right_g_2_deg_1 = PS_new[right_position].aug_sum_dist_c_p_g_2_deg_1; a_P_right_g_2_deg_2 = PS_new[right_position].aug_sum_dist_c_p_g_2_deg_2;
			
			cur_l_g.density_g_0 += (1 - our_model.gamma_s * cur_l.dist_n1 * cur_l.dist_n1) * a_P_right_g_0_deg_0
				+ 2 * our_model.gamma_s * cur_l.dist_n1 * a_P_right_g_0_deg_1 - our_model.gamma_s * a_P_right_g_0_deg_2;
			cur_l_g.density_g_1 += (1 - our_model.gamma_s * cur_l.dist_n1 * cur_l.dist_n1) * a_P_right_g_1_deg_0
				+ 2 * our_model.gamma_s * cur_l.dist_n1 * a_P_right_g_1_deg_1 - our_model.gamma_s * a_P_right_g_1_deg_2;
			cur_l_g.density_g_2 += (1 - our_model.gamma_s * cur_l.dist_n1 * cur_l.dist_n1) * a_P_right_g_2_deg_0
				+ 2 * our_model.gamma_s * cur_l.dist_n1 * a_P_right_g_2_deg_1 - our_model.gamma_s * a_P_right_g_2_deg_2;

			return;
		}

		if (left_position == -1) //left position is -1
		{
			a_P_left_g_0_deg_0 = PS_new[q_position].aug_sum_dist_c_p_g_0_deg_0; a_P_left_g_0_deg_1 = PS_new[q_position].aug_sum_dist_c_p_g_0_deg_1; a_P_left_g_0_deg_2 = PS_new[q_position].aug_sum_dist_c_p_g_0_deg_2;
			a_P_left_g_1_deg_0 = PS_new[q_position].aug_sum_dist_c_p_g_1_deg_0; a_P_left_g_1_deg_1 = PS_new[q_position].aug_sum_dist_c_p_g_1_deg_1; a_P_left_g_1_deg_2 = PS_new[q_position].aug_sum_dist_c_p_g_1_deg_2;
			a_P_left_g_2_deg_0 = PS_new[q_position].aug_sum_dist_c_p_g_2_deg_0; a_P_left_g_2_deg_1 = PS_new[q_position].aug_sum_dist_c_p_g_2_deg_1; a_P_left_g_2_deg_2 = PS_new[q_position].aug_sum_dist_c_p_g_2_deg_2;

			a_P_right_g_0_deg_0 = PS_new[right_position].aug_sum_dist_c_p_g_0_deg_0 - PS_new[q_position].aug_sum_dist_c_p_g_0_deg_0; a_P_right_g_0_deg_1 = PS_new[right_position].aug_sum_dist_c_p_g_0_deg_1 - PS_new[q_position].aug_sum_dist_c_p_g_0_deg_1; a_P_right_g_0_deg_2 = PS_new[right_position].aug_sum_dist_c_p_g_0_deg_2 - PS_new[q_position].aug_sum_dist_c_p_g_0_deg_2;
			a_P_right_g_1_deg_0 = PS_new[right_position].aug_sum_dist_c_p_g_1_deg_0 - PS_new[q_position].aug_sum_dist_c_p_g_1_deg_0; a_P_right_g_1_deg_1 = PS_new[right_position].aug_sum_dist_c_p_g_1_deg_1 - PS_new[q_position].aug_sum_dist_c_p_g_1_deg_1; a_P_right_g_1_deg_2 = PS_new[right_position].aug_sum_dist_c_p_g_1_deg_2 - PS_new[q_position].aug_sum_dist_c_p_g_1_deg_2;
			a_P_right_g_2_deg_0 = PS_new[right_position].aug_sum_dist_c_p_g_2_deg_0 - PS_new[q_position].aug_sum_dist_c_p_g_2_deg_0; a_P_right_g_2_deg_1 = PS_new[right_position].aug_sum_dist_c_p_g_2_deg_1 - PS_new[q_position].aug_sum_dist_c_p_g_2_deg_1; a_P_right_g_2_deg_2 = PS_new[right_position].aug_sum_dist_c_p_g_2_deg_2 - PS_new[q_position].aug_sum_dist_c_p_g_2_deg_2;

			cur_l_g.density_g_0 += (1 - our_model.gamma_s * cur_l.dist_n1 * cur_l.dist_n1) * a_P_left_g_0_deg_0
				+ 2 * our_model.gamma_s * cur_l.dist_n1 * a_P_left_g_0_deg_1 - our_model.gamma_s * a_P_left_g_0_deg_2
				+ (1 - our_model.gamma_s * cur_l.dist_n1 * cur_l.dist_n1) * a_P_right_g_0_deg_0
				+ 2 * our_model.gamma_s * cur_l.dist_n1 * a_P_right_g_0_deg_1 - our_model.gamma_s * a_P_right_g_0_deg_2;

			cur_l_g.density_g_1 += (1 - our_model.gamma_s * cur_l.dist_n1 * cur_l.dist_n1) * a_P_left_g_1_deg_0
				+ 2 * our_model.gamma_s * cur_l.dist_n1 * a_P_left_g_1_deg_1 - our_model.gamma_s * a_P_left_g_1_deg_2
				+ (1 - our_model.gamma_s * cur_l.dist_n1 * cur_l.dist_n1) * a_P_right_g_1_deg_0
				+ 2 * our_model.gamma_s * cur_l.dist_n1 * a_P_right_g_1_deg_1 - our_model.gamma_s * a_P_right_g_1_deg_2;

			cur_l_g.density_g_2 += (1 - our_model.gamma_s * cur_l.dist_n1 * cur_l.dist_n1) * a_P_left_g_2_deg_0
				+ 2 * our_model.gamma_s * cur_l.dist_n1 * a_P_left_g_2_deg_1 - our_model.gamma_s * a_P_left_g_2_deg_2
				+ (1 - our_model.gamma_s * cur_l.dist_n1 * cur_l.dist_n1) * a_P_right_g_2_deg_0
				+ 2 * our_model.gamma_s * cur_l.dist_n1 * a_P_right_g_2_deg_1 - our_model.gamma_s * a_P_right_g_2_deg_2;

			return;
		}

		//left position >= 0
		a_P_left_g_0_deg_0 = PS_new[q_position].aug_sum_dist_c_p_g_0_deg_0 - PS_new[left_position].aug_sum_dist_c_p_g_0_deg_0; a_P_left_g_0_deg_1 = PS_new[q_position].aug_sum_dist_c_p_g_0_deg_1 - PS_new[left_position].aug_sum_dist_c_p_g_0_deg_1; a_P_left_g_0_deg_2 = PS_new[q_position].aug_sum_dist_c_p_g_0_deg_2 - PS_new[left_position].aug_sum_dist_c_p_g_0_deg_2;
		a_P_left_g_1_deg_0 = PS_new[q_position].aug_sum_dist_c_p_g_1_deg_0 - PS_new[left_position].aug_sum_dist_c_p_g_1_deg_0; a_P_left_g_1_deg_1 = PS_new[q_position].aug_sum_dist_c_p_g_1_deg_1 - PS_new[left_position].aug_sum_dist_c_p_g_1_deg_1; a_P_left_g_1_deg_2 = PS_new[q_position].aug_sum_dist_c_p_g_1_deg_2 - PS_new[left_position].aug_sum_dist_c_p_g_1_deg_2;
		a_P_left_g_2_deg_0 = PS_new[q_position].aug_sum_dist_c_p_g_2_deg_0 - PS_new[left_position].aug_sum_dist_c_p_g_2_deg_0; a_P_left_g_2_deg_1 = PS_new[q_position].aug_sum_dist_c_p_g_2_deg_1 - PS_new[left_position].aug_sum_dist_c_p_g_2_deg_1; a_P_left_g_2_deg_2 = PS_new[q_position].aug_sum_dist_c_p_g_2_deg_2 - PS_new[left_position].aug_sum_dist_c_p_g_2_deg_2;
		a_P_right_g_0_deg_0 = PS_new[right_position].aug_sum_dist_c_p_g_0_deg_0 - PS_new[q_position].aug_sum_dist_c_p_g_0_deg_0; a_P_right_g_0_deg_1 = PS_new[right_position].aug_sum_dist_c_p_g_0_deg_1 - PS_new[q_position].aug_sum_dist_c_p_g_0_deg_1; a_P_right_g_0_deg_2 = PS_new[right_position].aug_sum_dist_c_p_g_0_deg_2 - PS_new[q_position].aug_sum_dist_c_p_g_0_deg_2;
		a_P_right_g_1_deg_0 = PS_new[right_position].aug_sum_dist_c_p_g_1_deg_0 - PS_new[q_position].aug_sum_dist_c_p_g_1_deg_0; a_P_right_g_1_deg_1 = PS_new[right_position].aug_sum_dist_c_p_g_1_deg_1 - PS_new[q_position].aug_sum_dist_c_p_g_1_deg_1; a_P_right_g_1_deg_2 = PS_new[right_position].aug_sum_dist_c_p_g_1_deg_2 - PS_new[q_position].aug_sum_dist_c_p_g_1_deg_2;
		a_P_right_g_2_deg_0 = PS_new[right_position].aug_sum_dist_c_p_g_2_deg_0 - PS_new[q_position].aug_sum_dist_c_p_g_2_deg_0; a_P_right_g_2_deg_1 = PS_new[right_position].aug_sum_dist_c_p_g_2_deg_1 - PS_new[q_position].aug_sum_dist_c_p_g_2_deg_1; a_P_right_g_2_deg_2 = PS_new[right_position].aug_sum_dist_c_p_g_2_deg_2 - PS_new[q_position].aug_sum_dist_c_p_g_2_deg_2;

		cur_l_g.density_g_0 += (1 - our_model.gamma_s * cur_l.dist_n1 * cur_l.dist_n1) * a_P_left_g_0_deg_0
			+ 2 * our_model.gamma_s * cur_l.dist_n1 * a_P_left_g_0_deg_1 - our_model.gamma_s * a_P_left_g_0_deg_2
			+ (1 - our_model.gamma_s * cur_l.dist_n1 * cur_l.dist_n1) * a_P_right_g_0_deg_0
			+ 2 * our_model.gamma_s * cur_l.dist_n1 * a_P_right_g_0_deg_1 - our_model.gamma_s * a_P_right_g_0_deg_2;

		cur_l_g.density_g_1 += (1 - our_model.gamma_s * cur_l.dist_n1 * cur_l.dist_n1) * a_P_left_g_1_deg_0
			+ 2 * our_model.gamma_s * cur_l.dist_n1 * a_P_left_g_1_deg_1 - our_model.gamma_s * a_P_left_g_1_deg_2
			+ (1 - our_model.gamma_s * cur_l.dist_n1 * cur_l.dist_n1) * a_P_right_g_1_deg_0
			+ 2 * our_model.gamma_s * cur_l.dist_n1 * a_P_right_g_1_deg_1 - our_model.gamma_s * a_P_right_g_1_deg_2;

		cur_l_g.density_g_2 += (1 - our_model.gamma_s * cur_l.dist_n1 * cur_l.dist_n1) * a_P_left_g_2_deg_0
			+ 2 * our_model.gamma_s * cur_l.dist_n1 * a_P_left_g_2_deg_1 - our_model.gamma_s * a_P_left_g_2_deg_2
			+ (1 - our_model.gamma_s * cur_l.dist_n1 * cur_l.dist_n1) * a_P_right_g_2_deg_0
			+ 2 * our_model.gamma_s * cur_l.dist_n1 * a_P_right_g_2_deg_1 - our_model.gamma_s * a_P_right_g_2_deg_2;

		return;
	}

	//q is in different edge
	//compute the distances between q to c and q to d 
	dist_q_c = min(cur_l.dist_n1 + our_model.sp_node_vec_node_a[node_index_pi_c].cur_sp_value,
		cur_l.dist_n2 + our_model.sp_node_vec_node_b[node_index_pi_c].cur_sp_value);
	dist_q_d = min(cur_l.dist_n1 + our_model.sp_node_vec_node_a[node_index_pi_d].cur_sp_value,
		cur_l.dist_n2 + our_model.sp_node_vec_node_b[node_index_pi_d].cur_sp_value);

	//Case 1:
	if (dist_q_c > our_model.bandwidth_s && dist_q_d > our_model.bandwidth_s)
		return;

	//Case 2:
	if (dist_q_c <= our_model.bandwidth_s && dist_q_d > our_model.bandwidth_s)
	{
		left_value = our_model.bandwidth_s - dist_q_c;
		left_position = bucket_find_position_sharing(left_value, our_model, e);

		if (left_position == -1) //left position is -1
			return;
		else
		{
			cur_l_g.density_g_0 += (1 - our_model.gamma_s * dist_q_c * dist_q_c) * PS_new[left_position].aug_sum_dist_c_p_g_0_deg_0
				- 2 * our_model.gamma_s * dist_q_c * PS_new[left_position].aug_sum_dist_c_p_g_0_deg_1 - our_model.gamma_s * PS_new[left_position].aug_sum_dist_c_p_g_0_deg_2;
			cur_l_g.density_g_1 += (1 - our_model.gamma_s * dist_q_c * dist_q_c) * PS_new[left_position].aug_sum_dist_c_p_g_1_deg_0
				- 2 * our_model.gamma_s * dist_q_c * PS_new[left_position].aug_sum_dist_c_p_g_1_deg_1 - our_model.gamma_s * PS_new[left_position].aug_sum_dist_c_p_g_1_deg_2;
			cur_l_g.density_g_2 += (1 - our_model.gamma_s * dist_q_c * dist_q_c) * PS_new[left_position].aug_sum_dist_c_p_g_2_deg_0
				- 2 * our_model.gamma_s * dist_q_c * PS_new[left_position].aug_sum_dist_c_p_g_2_deg_1 - our_model.gamma_s * PS_new[left_position].aug_sum_dist_c_p_g_2_deg_2;

			return;
		}
	}

	//Case 3:
	if (dist_q_c > our_model.bandwidth_s && dist_q_d <= our_model.bandwidth_s)
	{
		right_value = our_model.bandwidth_s - dist_q_d;
		right_position = bucket_find_position_sharing(e.length - right_value, our_model, e) + 1;

		//if (right_position == PS.size() - 1)
		if (right_position == PS_new.size())
			return;
		else
		{
			cur_l_g.density_g_0 += (1 - our_model.gamma_s * dist_q_d * dist_q_d) * PS_new[right_position].aug_sum_dist_d_p_g_0_deg_0
				- 2 * our_model.gamma_s * dist_q_d * PS_new[right_position].aug_sum_dist_d_p_g_0_deg_1 - our_model.gamma_s * PS_new[right_position].aug_sum_dist_d_p_g_0_deg_2;
			cur_l_g.density_g_1 += (1 - our_model.gamma_s * dist_q_d * dist_q_d) * PS_new[right_position].aug_sum_dist_d_p_g_1_deg_0
				- 2 * our_model.gamma_s * dist_q_d * PS_new[right_position].aug_sum_dist_d_p_g_1_deg_1 - our_model.gamma_s * PS_new[right_position].aug_sum_dist_d_p_g_1_deg_2;
			cur_l_g.density_g_2 += (1 - our_model.gamma_s * dist_q_d * dist_q_d) * PS_new[right_position].aug_sum_dist_d_p_g_2_deg_0
				- 2 * our_model.gamma_s * dist_q_d * PS_new[right_position].aug_sum_dist_d_p_g_2_deg_1 - our_model.gamma_s * PS_new[right_position].aug_sum_dist_d_p_g_2_deg_2;

			return;
		}
	}

	//Case 4:
	if (dist_q_c <= our_model.bandwidth_s && dist_q_d <= our_model.bandwidth_s)
	{
		left_value = our_model.bandwidth_s - dist_q_c;
		left_position = bucket_find_position_sharing(left_value, our_model, e);

		right_value = our_model.bandwidth_s - dist_q_d;
		right_position = bucket_find_position_sharing(e.length - right_value, our_model, e) + 1;

		if (left_position == -1 || right_position == PS_new.size())
		{
			if (left_position == -1 && right_position == PS_new.size())
				return;
			if (right_position == PS_new.size())
			{
				cur_l_g.density_g_0 += (1 - our_model.gamma_s * dist_q_c * dist_q_c) * PS_new[left_position].aug_sum_dist_c_p_g_0_deg_0
					- 2 * our_model.gamma_s * dist_q_c * PS_new[left_position].aug_sum_dist_c_p_g_0_deg_1 - our_model.gamma_s * PS_new[left_position].aug_sum_dist_c_p_g_0_deg_2;
				cur_l_g.density_g_1 += (1 - our_model.gamma_s * dist_q_c * dist_q_c) * PS_new[left_position].aug_sum_dist_c_p_g_1_deg_0
					- 2 * our_model.gamma_s * dist_q_c * PS_new[left_position].aug_sum_dist_c_p_g_1_deg_1 - our_model.gamma_s * PS_new[left_position].aug_sum_dist_c_p_g_1_deg_2;
				cur_l_g.density_g_2 += (1 - our_model.gamma_s * dist_q_c * dist_q_c) * PS_new[left_position].aug_sum_dist_c_p_g_2_deg_0
					- 2 * our_model.gamma_s * dist_q_c * PS_new[left_position].aug_sum_dist_c_p_g_2_deg_1 - our_model.gamma_s * PS_new[left_position].aug_sum_dist_c_p_g_2_deg_2;

				return;
			}

			if (left_position == -1)
			{
				cur_l_g.density_g_0 += (1 - our_model.gamma_s * dist_q_d * dist_q_d) * PS_new[right_position].aug_sum_dist_d_p_g_0_deg_0
					- 2 * our_model.gamma_s * dist_q_d * PS_new[right_position].aug_sum_dist_d_p_g_0_deg_1 - our_model.gamma_s * PS_new[right_position].aug_sum_dist_d_p_g_0_deg_2;
				cur_l_g.density_g_1 += (1 - our_model.gamma_s * dist_q_d * dist_q_d) * PS_new[right_position].aug_sum_dist_d_p_g_1_deg_0
					- 2 * our_model.gamma_s * dist_q_d * PS_new[right_position].aug_sum_dist_d_p_g_1_deg_1 - our_model.gamma_s * PS_new[right_position].aug_sum_dist_d_p_g_1_deg_2;
				cur_l_g.density_g_2 += (1 - our_model.gamma_s * dist_q_d * dist_q_d) * PS_new[right_position].aug_sum_dist_d_p_g_2_deg_0
					- 2 * our_model.gamma_s * dist_q_d * PS_new[right_position].aug_sum_dist_d_p_g_2_deg_1 - our_model.gamma_s * PS_new[right_position].aug_sum_dist_d_p_g_2_deg_2;

				return;
			}
		}
		else
		{
			if (left_position < right_position)
			{
				cur_l_g.density_g_0 += (1 - our_model.gamma_s * dist_q_c * dist_q_c) * PS_new[left_position].aug_sum_dist_c_p_g_0_deg_0
					- 2 * our_model.gamma_s * dist_q_c * PS_new[left_position].aug_sum_dist_c_p_g_0_deg_1 - our_model.gamma_s * PS_new[left_position].aug_sum_dist_c_p_g_0_deg_2
					+ (1 - our_model.gamma_s * dist_q_d * dist_q_d) * PS_new[right_position].aug_sum_dist_d_p_g_0_deg_0
					- 2 * our_model.gamma_s * dist_q_d * PS_new[right_position].aug_sum_dist_d_p_g_0_deg_1 - our_model.gamma_s * PS_new[right_position].aug_sum_dist_d_p_g_0_deg_2;

				cur_l_g.density_g_1 += (1 - our_model.gamma_s * dist_q_c * dist_q_c) * PS_new[left_position].aug_sum_dist_c_p_g_1_deg_0
					- 2 * our_model.gamma_s * dist_q_c * PS_new[left_position].aug_sum_dist_c_p_g_1_deg_1 - our_model.gamma_s * PS_new[left_position].aug_sum_dist_c_p_g_1_deg_2
					+ (1 - our_model.gamma_s * dist_q_d * dist_q_d) * PS_new[right_position].aug_sum_dist_d_p_g_1_deg_0
					- 2 * our_model.gamma_s * dist_q_d * PS_new[right_position].aug_sum_dist_d_p_g_1_deg_1 - our_model.gamma_s * PS_new[right_position].aug_sum_dist_d_p_g_1_deg_2;

				cur_l_g.density_g_2 += (1 - our_model.gamma_s * dist_q_c * dist_q_c) * PS_new[left_position].aug_sum_dist_c_p_g_2_deg_0
					- 2 * our_model.gamma_s * dist_q_c * PS_new[left_position].aug_sum_dist_c_p_g_2_deg_1 - our_model.gamma_s * PS_new[left_position].aug_sum_dist_c_p_g_2_deg_2
					+ (1 - our_model.gamma_s * dist_q_d * dist_q_d) * PS_new[right_position].aug_sum_dist_d_p_g_2_deg_0
					- 2 * our_model.gamma_s * dist_q_d * PS_new[right_position].aug_sum_dist_d_p_g_2_deg_1 - our_model.gamma_s * PS_new[right_position].aug_sum_dist_d_p_g_2_deg_2;

				return;
			}
			else
			{
				middle_position = bucket_find_position_sharing((dist_q_d - dist_q_c + e.length) / 2.0, our_model, e);

				if (middle_position > -1 && middle_position < (int)PS_new.size() - 1)
				{
					cur_l_g.density_g_0 += (1 - our_model.gamma_s * dist_q_c * dist_q_c) * PS_new[middle_position].aug_sum_dist_c_p_g_0_deg_0
						- 2 * our_model.gamma_s * dist_q_c * PS_new[middle_position].aug_sum_dist_c_p_g_0_deg_1 - our_model.gamma_s * PS_new[middle_position].aug_sum_dist_c_p_g_0_deg_2
						+ (1 - our_model.gamma_s * dist_q_d * dist_q_d) * PS_new[middle_position + 1].aug_sum_dist_d_p_g_0_deg_0
						- 2 * our_model.gamma_s * dist_q_d * PS_new[middle_position + 1].aug_sum_dist_d_p_g_0_deg_1 - our_model.gamma_s * PS_new[middle_position + 1].aug_sum_dist_d_p_g_0_deg_2;

					cur_l_g.density_g_1 += (1 - our_model.gamma_s * dist_q_c * dist_q_c) * PS_new[middle_position].aug_sum_dist_c_p_g_1_deg_0
						- 2 * our_model.gamma_s * dist_q_c * PS_new[middle_position].aug_sum_dist_c_p_g_1_deg_1 - our_model.gamma_s * PS_new[middle_position].aug_sum_dist_c_p_g_1_deg_2
						+ (1 - our_model.gamma_s * dist_q_d * dist_q_d) * PS_new[middle_position + 1].aug_sum_dist_d_p_g_1_deg_0
						- 2 * our_model.gamma_s * dist_q_d * PS_new[middle_position + 1].aug_sum_dist_d_p_g_1_deg_1 - our_model.gamma_s * PS_new[middle_position + 1].aug_sum_dist_d_p_g_1_deg_2;

					cur_l_g.density_g_2 += (1 - our_model.gamma_s * dist_q_c * dist_q_c) * PS_new[middle_position].aug_sum_dist_c_p_g_2_deg_0
						- 2 * our_model.gamma_s * dist_q_c * PS_new[middle_position].aug_sum_dist_c_p_g_2_deg_1 - our_model.gamma_s * PS_new[middle_position].aug_sum_dist_c_p_g_2_deg_2
						+ (1 - our_model.gamma_s * dist_q_d * dist_q_d) * PS_new[middle_position + 1].aug_sum_dist_d_p_g_2_deg_0
						- 2 * our_model.gamma_s * dist_q_d * PS_new[middle_position + 1].aug_sum_dist_d_p_g_2_deg_1 - our_model.gamma_s * PS_new[middle_position + 1].aug_sum_dist_d_p_g_2_deg_2;

					return;
				}
				else
				{
					if (middle_position == -1)
					{
						cur_l_g.density_g_0 += (1 - our_model.gamma_s * dist_q_d * dist_q_d) * PS_new[0].aug_sum_dist_d_p_g_0_deg_0
							- 2 * our_model.gamma_s * dist_q_d * PS_new[0].aug_sum_dist_d_p_g_0_deg_1 - our_model.gamma_s * PS_new[0].aug_sum_dist_d_p_g_0_deg_2;

						cur_l_g.density_g_1 += (1 - our_model.gamma_s * dist_q_d * dist_q_d) * PS_new[0].aug_sum_dist_d_p_g_1_deg_0
							- 2 * our_model.gamma_s * dist_q_d * PS_new[0].aug_sum_dist_d_p_g_1_deg_1 - our_model.gamma_s * PS_new[0].aug_sum_dist_d_p_g_1_deg_2;

						cur_l_g.density_g_2 += (1 - our_model.gamma_s * dist_q_d * dist_q_d) * PS_new[0].aug_sum_dist_d_p_g_2_deg_0
							- 2 * our_model.gamma_s * dist_q_d * PS_new[0].aug_sum_dist_d_p_g_2_deg_1 - our_model.gamma_s * PS_new[0].aug_sum_dist_d_p_g_2_deg_2;

						return;
					}
					else
					{
						if (middle_position == PS_new.size() - 1)
						{
							cur_l_g.density_g_0 += (1 - our_model.gamma_s * dist_q_c * dist_q_c) * PS_new[PS_new.size() - 1].aug_sum_dist_c_p_g_0_deg_0
								- 2 * our_model.gamma_s * dist_q_c * PS_new[PS_new.size() - 1].aug_sum_dist_c_p_g_0_deg_1 - our_model.gamma_s * PS_new[PS_new.size() - 1].aug_sum_dist_c_p_g_0_deg_2;

							cur_l_g.density_g_1 += (1 - our_model.gamma_s * dist_q_c * dist_q_c) * PS_new[PS_new.size() - 1].aug_sum_dist_c_p_g_1_deg_0
								- 2 * our_model.gamma_s * dist_q_c * PS_new[PS_new.size() - 1].aug_sum_dist_c_p_g_1_deg_1 - our_model.gamma_s * PS_new[PS_new.size() - 1].aug_sum_dist_c_p_g_1_deg_2;

							cur_l_g.density_g_2 += (1 - our_model.gamma_s * dist_q_c * dist_q_c) * PS_new[PS_new.size() - 1].aug_sum_dist_c_p_g_2_deg_0
								- 2 * our_model.gamma_s * dist_q_c * PS_new[PS_new.size() - 1].aug_sum_dist_c_p_g_2_deg_1 - our_model.gamma_s * PS_new[PS_new.size() - 1].aug_sum_dist_c_p_g_2_deg_2;

							return;
						}
					}
				}
			}
		}

		return;
	}
}

void update_density_g(model& our_model)
{
	int prev_edge_index = -1;
	int edge_index;

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

		for (int e = 0; e < our_model.access_edge_list.size(); e++)
		{
			for (int ii = 0; ii < our_model.interval_vec.size(); ii++)
			{
				our_model.cur_interval_index = ii;
				edge_KAF_TIS(our_model, l, ii, our_model.access_edge_list[e]);
			}
		}
	}
}

void window_density(model& our_model)
{
	int left_interval_index = 0;
	int right_interval_index = 0;
	double start_win_time;
	double end_win_time;

	for (int t = 0; t < our_model.T; t++)
	{
		start_win_time = our_model.time_vec[t] - our_model.bandwidth_t;
		end_win_time = our_model.time_vec[t] + our_model.bandwidth_t;

		if (t == 0)
		{
			while (true)
			{
				for (int l = 0; l < our_model.lixel_set.size(); l++)
				{
					our_model.slide_line[l].density_g_0 += our_model.interval_vec[right_interval_index].lixel_g_set[l].density_g_0;
					our_model.slide_line[l].density_g_1 += our_model.interval_vec[right_interval_index].lixel_g_set[l].density_g_1;
					our_model.slide_line[l].density_g_2 += our_model.interval_vec[right_interval_index].lixel_g_set[l].density_g_2;
				}

				if (fabs(our_model.interval_vec[right_interval_index].end_time - end_win_time) <= eps)
				{
					for (int l = 0; l < our_model.lixel_set.size(); l++)
					{
						our_model.lixel_set[l].KDE_value_vec[t] = (1 - our_model.gamma_t * our_model.time_vec[t] * our_model.time_vec[t]) * our_model.slide_line[l].density_g_0
							+ 2 * our_model.gamma_t * our_model.time_vec[t] * our_model.slide_line[l].density_g_1 - our_model.gamma_t * our_model.slide_line[l].density_g_2;
					}

					right_interval_index++;
					//cout << left_interval_index << " " << right_interval_index << endl;

					break;
				}

				right_interval_index++;
			}

			continue;
		}

		while (true)
		{
			right_interval_index++;

			for (int l = 0; l < our_model.lixel_set.size(); l++)
			{
				our_model.slide_line[l].density_g_0 += our_model.interval_vec[right_interval_index].lixel_g_set[l].density_g_0;
				our_model.slide_line[l].density_g_1 += our_model.interval_vec[right_interval_index].lixel_g_set[l].density_g_1;
				our_model.slide_line[l].density_g_2 += our_model.interval_vec[right_interval_index].lixel_g_set[l].density_g_2;
			}

			if (fabs(our_model.interval_vec[right_interval_index].end_time - end_win_time) <= eps)
				break;

			/*for (int l = 0; l < our_model.lixel_set.size(); l++)
			{
				our_model.slide_line[l].density_g_0 += our_model.interval_vec[right_interval_index].lixel_g_set[l].density_g_0;
				our_model.slide_line[l].density_g_1 += our_model.interval_vec[right_interval_index].lixel_g_set[l].density_g_1;
				our_model.slide_line[l].density_g_2 += our_model.interval_vec[right_interval_index].lixel_g_set[l].density_g_2;
			}

			if (fabs(our_model.interval_vec[right_interval_index].end_time - end_win_time) <= eps)
			{
				right_interval_index++;
				break;
			}

			right_interval_index++;*/
		}

		while (true)
		{
			for (int l = 0; l < our_model.lixel_set.size(); l++)
			{
				our_model.slide_line[l].density_g_0 -= our_model.interval_vec[left_interval_index].lixel_g_set[l].density_g_0;
				our_model.slide_line[l].density_g_1 -= our_model.interval_vec[left_interval_index].lixel_g_set[l].density_g_1;
				our_model.slide_line[l].density_g_2 -= our_model.interval_vec[left_interval_index].lixel_g_set[l].density_g_2;
			}

			if (fabs(our_model.interval_vec[left_interval_index].end_time - start_win_time) <= eps)
			{
				left_interval_index++;
				break;
			}

			left_interval_index++;

			/*for (int l = 0; l < our_model.lixel_set.size(); l++)
			{
				our_model.slide_line[l].density_g_0 -= our_model.interval_vec[left_interval_index].lixel_g_set[l].density_g_0;
				our_model.slide_line[l].density_g_1 -= our_model.interval_vec[left_interval_index].lixel_g_set[l].density_g_1;
				our_model.slide_line[l].density_g_2 -= our_model.interval_vec[left_interval_index].lixel_g_set[l].density_g_2;
			}

			if (fabs(our_model.interval_vec[left_interval_index].start_time - start_win_time) <= eps)
			{
				//left_interval_index++;
				break;
			}

			left_interval_index++;*/
		}

		//cout << left_interval_index << " " << right_interval_index << endl;

		for (int l = 0; l < our_model.lixel_set.size(); l++)
		{
			our_model.lixel_set[l].KDE_value_vec[t] = (1 - our_model.gamma_t * our_model.time_vec[t] * our_model.time_vec[t]) * our_model.slide_line[l].density_g_0
				+ 2 * our_model.gamma_t * our_model.time_vec[t] * our_model.slide_line[l].density_g_1 - our_model.gamma_t * our_model.slide_line[l].density_g_2;
		}
	}
}

void TIS_algorithm(model& our_model)
{
	construct_intervals(our_model);
	init_TIS(our_model);
	PS_new_distribution_sharing(our_model);
	augment_preprocess_sharing(our_model);
	bucket_construction_sharing(our_model);

	update_density_g(our_model);
	window_density(our_model);

	PS_new_sharing_clear(our_model);
	bucket_clear_sharing(our_model);
	interval_clear(our_model);
}

void augment_preprocess_sharing(model& our_model)
{
	double incr_dist_c_p_g_0_deg_0, incr_dist_c_p_g_0_deg_1, incr_dist_c_p_g_0_deg_2;
	double incr_dist_c_p_g_1_deg_0, incr_dist_c_p_g_1_deg_1, incr_dist_c_p_g_1_deg_2;
	double incr_dist_c_p_g_2_deg_0, incr_dist_c_p_g_2_deg_1, incr_dist_c_p_g_2_deg_2;
	double incr_dist_d_p_g_0_deg_0, incr_dist_d_p_g_0_deg_1, incr_dist_d_p_g_0_deg_2;
	double incr_dist_d_p_g_1_deg_0, incr_dist_d_p_g_1_deg_1, incr_dist_d_p_g_1_deg_2;
	double incr_dist_d_p_g_2_deg_0, incr_dist_d_p_g_2_deg_1, incr_dist_d_p_g_2_deg_2;
	int size;

	for (int interval_index = 0; interval_index < our_model.interval_vec.size(); interval_index++)
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

			vector<Point>& PS_new = our_model.edge_set[e].PS_new_vec[interval_index];
			size = PS_new.size();

			for (int i = 0; i < size; i++)
			{
				our_model.edge_set[e].multi_dist_n1_vec[interval_index][i] = PS_new[i].dist_n1;

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

void bucket_construction_sharing(model& our_model)
{
	Edge* edge_set = our_model.edge_set;
	Bucket temp_bucket;
	int num_buckets;
	int b_index;

	for (int e = 0; e < our_model.m; e++)
	{
		num_buckets = (int)ceil(edge_set[e].length / our_model.lixel_reg_length);

		for (int interval_index = 0; interval_index < our_model.interval_vec.size(); interval_index++)
		{
			for (int b = 0; b < num_buckets; b++)
				edge_set[e].multi_bucket_vec[interval_index].push_back(temp_bucket);

			vector<Point>& PS_new = edge_set[e].PS_new_vec[interval_index];
			for (int i = 0; i < PS_new.size(); i++)
			{
				b_index = (int)floor(PS_new[i].dist_n1 / our_model.lixel_reg_length);
				edge_set[e].multi_bucket_vec[interval_index][b_index].dist_n1_vec.push_back(PS_new[i].dist_n1);
			}

			for (int b = 0; b < num_buckets; b++)
			{
				if (b == num_buckets - 1)
					edge_set[e].multi_bucket_vec[interval_index][b].bucket_length = edge_set[e].length - (num_buckets - 1) * our_model.lixel_reg_length;
				else
					edge_set[e].multi_bucket_vec[interval_index][b].bucket_length = our_model.lixel_reg_length;

				if (b == 0)
					edge_set[e].multi_bucket_vec[interval_index][b].num_accum_index = -1;
				else
					edge_set[e].multi_bucket_vec[interval_index][b].num_accum_index = edge_set[e].multi_bucket_vec[interval_index][b - 1].num_accum_index + edge_set[e].multi_bucket_vec[interval_index][b - 1].dist_n1_vec.size();
			}
		}
	}
}

int bucket_find_position_sharing(double value, model& our_model, Edge& e)
{
	int num_buckets = (int)e.multi_bucket_vec[our_model.cur_interval_index].size();
	int bucket_index = max(min((int)floor(value / our_model.lixel_reg_length), num_buckets - 1), 0);
	vector<double>::iterator position_iterator = lower_bound(e.multi_bucket_vec[our_model.cur_interval_index][bucket_index].dist_n1_vec.begin(), e.multi_bucket_vec[our_model.cur_interval_index][bucket_index].dist_n1_vec.end(), value);
	int pos_index = e.multi_bucket_vec[our_model.cur_interval_index][bucket_index].num_accum_index + (int)(position_iterator - e.multi_bucket_vec[our_model.cur_interval_index][bucket_index].dist_n1_vec.begin());

	return pos_index;
}

void bucket_clear_sharing(model& our_model)
{
	Edge* edge_set = our_model.edge_set;
	int num_buckets;

	for (int e = 0; e < our_model.m; e++)
	{
		for (int i = 0; i < our_model.num_intervals; i++)
		{
			num_buckets = edge_set[e].multi_bucket_vec[i].size();
			for (int b = 0; b < num_buckets; b++)
				edge_set[e].multi_bucket_vec[i][b].dist_n1_vec.clear();

			edge_set[e].multi_bucket_vec[i].clear();
		}

		edge_set[e].multi_bucket_vec.clear();
	}
}