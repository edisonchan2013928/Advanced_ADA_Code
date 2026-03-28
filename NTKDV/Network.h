#pragma once
#ifndef NETWORK_H
#define NETWORK_H

#include "Library.h"

const double inf = 999999999999999;
const double eps = 0.000001;

struct Point
{
	double dist_n1;
	double dist_n2;
	double time;
	int edge_index;

	//Used in ADA method
	double aug_sum_dist_c_p_g_0_deg_0; double aug_sum_dist_c_p_g_0_deg_1; double aug_sum_dist_c_p_g_0_deg_2;
	double aug_sum_dist_c_p_g_1_deg_0; double aug_sum_dist_c_p_g_1_deg_1; double aug_sum_dist_c_p_g_1_deg_2;
	double aug_sum_dist_c_p_g_2_deg_0; double aug_sum_dist_c_p_g_2_deg_1; double aug_sum_dist_c_p_g_2_deg_2;

	double aug_sum_dist_d_p_g_0_deg_0; double aug_sum_dist_d_p_g_0_deg_1; double aug_sum_dist_d_p_g_0_deg_2;
	double aug_sum_dist_d_p_g_1_deg_0; double aug_sum_dist_d_p_g_1_deg_1; double aug_sum_dist_d_p_g_1_deg_2;
	double aug_sum_dist_d_p_g_2_deg_0; double aug_sum_dist_d_p_g_2_deg_1; double aug_sum_dist_d_p_g_2_deg_2;

	//double**aug_sum_dist_c_p_g_deg;
};

struct Bucket
{
	int num_accum_index;
	double bucket_length;
	vector<double> dist_n1_vec;
};

struct Edge
{
	int n1;
	int n2;
	double length;
	vector<Point> PS;

	//ADA_{TIC}^{B+}
	vector<Point> PS_new;
	vector<double> dist_n1_vec;
	vector<Bucket> bucket_vec;

	//ADA_{TIC}^{B+, \phi}
	vector<vector<Point>> PS_new_vec;
	vector<vector<double>> multi_dist_n1_vec;
	vector<vector<Bucket>> multi_bucket_vec;
};

struct Lixel : Point
{
	vector<double> KDE_value_vec;
};

struct Lixel_den_g : Point //Used in TIS
{
	double density_g_0;
	double density_g_1;
	double density_g_2;
};

struct time_interval //Used in TIS
{
	double start_time;
	double end_time;
	vector<Lixel_den_g> lixel_g_set;
};

struct sp_node
{
	int node_index;
	double cur_sp_value;
	bool is_opt;
};

//This is the minimum heap for dijkstra's algorithm
struct comparePriority
{
	bool operator()(sp_node& n1, sp_node& n2)
	{
		return n1.cur_sp_value > n2.cur_sp_value;
	}
};

struct model
{
	int method;
	int n, m; //number of nodes and edges
	double lixel_reg_length;
	char* network_fileName;
	char* out_NTKDV_fileName;
	char* timestamp_fileName;
	Edge* edge_set;
	vector<vector<int> > Network;
	vector<Lixel> lixel_set;

	//NTKDV
	int T; //number of timestamps
	vector<double> time_vec; //timestamps
	int cur_t_index;

	//kernel functions
	int k_type; //k_type = 0: Gaussian kernel, k_type = 1: Triangular kernel 2: Epanechnikov kernel 3: Quartic kernel
	double bandwidth_s;
	double bandwidth_t;
	double gamma_s;
	double gamma_t;

	//Used in ADA_{TIC, phi}
	vector<double> max_inteval_vec;

	//Used in ADA_{TIS}^{B+}
	vector<double> timestamps_vec;
	vector<time_interval> interval_vec;
	int num_intervals;
	int cur_interval_index;
	vector<Lixel_den_g> slide_line;

	//Used in Dijkstra's algorithm
	//Lixel cur_l;
	vector<sp_node> sp_node_vec;

	//Used in method = 2
	int sel_node_index;
	int node_index_pi_a;
	int node_index_pi_b;
	vector<sp_node> sp_node_vec_node_a;
	vector<sp_node> sp_node_vec_node_b;

	//Avoid scanning some edges that are far away from the bandwidth b
	vector<int> access_edge_list;
	bool* is_scan_edge_list;
};

#endif