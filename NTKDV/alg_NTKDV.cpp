#include "alg_NTKDV.h"

void alg_NTKDV::init_parameters(int argc, char** argv)
{
	our_model.network_fileName = argv[1];
	our_model.out_NTKDV_fileName = argv[2];
	our_model.method = atoi(argv[3]);
	our_model.lixel_reg_length = atoi(argv[4]);
	our_model.k_type = atoi(argv[5]);
	our_model.bandwidth_s = atof(argv[6]);
	our_model.bandwidth_t = atof(argv[7]);

	//debug
	/*our_model.network_fileName = (char*)"../../../Datasets/Extension_Datasets/NTKDV_Datasets/Testing/Testing_network";
	our_model.out_NTKDV_fileName = (char*)"./Results/Testing_network_m2_l2_k2_bs15_bt3_T8";
	our_model.method = 2;
	our_model.lixel_reg_length = 2;
	our_model.k_type = 2;
	our_model.bandwidth_s = 15;
	our_model.bandwidth_t = 3;
	our_model.gamma_s = 1;
	our_model.gamma_t = 1;*/

	//Epanechnikov and Quartic kernels
	if (our_model.k_type == 2 || our_model.k_type == 3)
	{
		our_model.gamma_s = 1 / (our_model.bandwidth_s * our_model.bandwidth_s);
		our_model.gamma_t = 1 / (our_model.bandwidth_t * our_model.bandwidth_t);
	}
		
}

void alg_NTKDV::load_network(char** argv)
{
	fstream network_file;
	int num_points;
	//int degree;
	//int edge_index;
	Point p;
	vector<int> tempVector;
	sp_node temp_sp_Node;

	our_model.network_fileName = argv[1];
	our_model.method = atoi(argv[3]);

	//debug
	//our_model.network_fileName = (char*)"../../../Datasets/Extension_Datasets/NTKDV_Datasets/Testing/Testing_network";
	//our_model.method = 2;

	network_file.open(our_model.network_fileName);
	if (network_file.is_open() == false)
	{
		cout << "Cannot open network file!" << endl;
		exit(0);
	}

	network_file >> our_model.n;
	network_file >> our_model.m;
	our_model.edge_set = new Edge[our_model.m];

	for (int e = 0; e < our_model.m; e++)
	{
		network_file >> our_model.edge_set[e].n1;
		network_file >> our_model.edge_set[e].n2;
		network_file >> our_model.edge_set[e].length;
		network_file >> num_points;

		for (int i = 0; i < num_points; i++)
		{
			p.edge_index = e;
			network_file >> p.dist_n1;
			network_file >> p.time;
			p.dist_n2 = our_model.edge_set[e].length - p.dist_n1;
			our_model.edge_set[e].PS.push_back(p);
			//our_model.edge_set[e].aug_dist_diff_vec.push_back(-inf);
			//our_model.edge_set[e].dist_n1_vec.push_back(-1);
			//our_model.edge_set[e].dist_n2_vec.push_back(-1);
		}
	}

	//Initialization
	for (int v = 0; v < our_model.n; v++)
		our_model.Network.push_back(tempVector);

	//Construct the road network
	for (int e = 0; e < our_model.m; e++)
	{
		our_model.Network[our_model.edge_set[e].n1].push_back(e);
		our_model.Network[our_model.edge_set[e].n2].push_back(e);
	}

	//init the Dijkstra's algorithm
	for (int v = 0; v < our_model.n; v++)
	{
		our_model.sp_node_vec.push_back(temp_sp_Node);

		if (our_model.method >= 2)
		{
			our_model.sp_node_vec_node_a.push_back(temp_sp_Node);
			our_model.sp_node_vec_node_b.push_back(temp_sp_Node);
		}
	}

	our_model.is_scan_edge_list = new bool[our_model.m];
	for (int e = 0; e < our_model.m; e++)
		our_model.is_scan_edge_list[e] = false;

	network_file.close();
}

void alg_NTKDV::NTKDV_algorithm()
{
	//Used in method = 2
	int prev_edge_index = -1;
	int num_edge_index = 0;
	int edge_index;

	double run_time;

	//Different algorithms
	auto start_s = chrono::high_resolution_clock::now();

	if (our_model.method == 1) //better version of NKDV baseline method
	{
		for (int l = 0; l < (int)our_model.lixel_set.size(); l++)
		{
			clear_scan_edge_list(our_model);
			Lixel& cur_l = our_model.lixel_set[l];
			dijkstra(our_model, cur_l);
			NTKDV_basic(our_model, cur_l);
		}
	}

	if (our_model.method == 2) //Save the number of calls in dijstra algorithm
	{
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

			NTKDV_basic(our_model, cur_l);
		}
	}

	if (our_model.method == 3) //TIC
	{
		if (check_validity(our_model) == true)
			TIC_phi_algorithm(our_model);
		else
			TIC_algorithm(our_model);
	}

	if (our_model.method == 4) //TIS
	{
		if (check_validity(our_model) == true)
			TIC_phi_algorithm(our_model);
		else
			TIS_algorithm(our_model);
	}

	auto end_s = chrono::high_resolution_clock::now();

	run_time = (chrono::duration_cast<chrono::nanoseconds>(end_s - start_s).count()) / 1000000000.0;

	std::cout << "method " << our_model.method << ":" << run_time << endl;
}

void alg_NTKDV::NTKDV_compute_file(int argc, char** argv)
{
	init_parameters(argc, argv);
	obtain_timestamps_min_max_range(argc, argv, our_model);
	obtain_lixel_set(our_model);
	NTKDV_algorithm();
	output_Visual();
}

string alg_NTKDV::NTKDV_compute(int argc, char** argv)
{
	init_parameters(argc, argv);
	obtain_timestamps_user_input(argc, argv, our_model);
	//obtain_timestamps_min_max_range(argc, argv, our_model);
	obtain_lixel_set(our_model);
	NTKDV_algorithm();

	return output_Visual_String();
}

void alg_NTKDV::output_Visual()
{
	int edge_index;
	int num_lixels;
	double dist_n1, dist_n2;
	double KDE_value;
	fstream out_NTKDV_file;

	out_NTKDV_file.open(our_model.out_NTKDV_fileName, ios::in | ios::out | ios::trunc);
	if (out_NTKDV_file.is_open() == false)
	{
		cout << "Cannot open output file!" << endl;
		exit(0);
	}

	num_lixels = (int)our_model.lixel_set.size();
	out_NTKDV_file << "Number of lixels: " << num_lixels << endl;
	out_NTKDV_file << "Number of timestamps: " << our_model.T << endl;
	for (int t = 0; t < our_model.T; t++)
	{
		out_NTKDV_file << "timestamp value:" << our_model.time_vec[t] << endl;
		for (int l = 0; l < num_lixels; l++)
		{
			edge_index = our_model.lixel_set[l].edge_index;
			dist_n1 = our_model.lixel_set[l].dist_n1;
			dist_n2 = our_model.lixel_set[l].dist_n2;
			KDE_value = our_model.lixel_set[l].KDE_value_vec[t];
			out_NTKDV_file << edge_index << " " << dist_n1 << " " << dist_n2 << " " << KDE_value << endl;
		}
	}

	out_NTKDV_file.close();
	clear_memory();
}

string alg_NTKDV::output_Visual_String()
{
	stringstream outString_ss;
	int edge_index;
	int num_lixels;
	double dist_n1, dist_n2;
	double KDE_value;

	num_lixels = (int)our_model.lixel_set.size();
	outString_ss << num_lixels << endl;
	outString_ss << our_model.T << endl;

	for (int t = 0; t < our_model.T; t++)
	{
		outString_ss << our_model.time_vec[t] << endl;
		for (int l = 0; l < num_lixels; l++)
		{
			edge_index = our_model.lixel_set[l].edge_index;
			dist_n1 = our_model.lixel_set[l].dist_n1;
			dist_n2 = our_model.lixel_set[l].dist_n2;
			KDE_value = our_model.lixel_set[l].KDE_value_vec[t];
			outString_ss << edge_index << " " << dist_n1 << " " << dist_n2 << " " << KDE_value << endl;
		}
	}

	return outString_ss.str();
}

void alg_NTKDV::clear_memory()
{
	for (int l = 0; l < our_model.lixel_set.size(); l++)
		our_model.lixel_set[l].KDE_value_vec.clear();

	our_model.lixel_set.clear();
	our_model.time_vec.clear();
}

void alg_NTKDV::clear_basic_memory()
{
	for (int e = 0; e < our_model.m; e++)
		our_model.edge_set[e].PS.clear();
	
	delete[] our_model.edge_set;

	for (int v = 0; v < our_model.n; v++)
		our_model.Network[v].clear();
	our_model.Network.clear();

	our_model.sp_node_vec.clear();

	if (our_model.method >= 2) //Those methods that need to use shortest path sharing
	{
		our_model.sp_node_vec_node_a.clear();
		our_model.sp_node_vec_node_b.clear();
	}
}

void obtain_timestamps_min_max_range(int argc, char** argv, model& our_model)
{
	our_model.T = atoi(argv[8]); //number of timestamps

	//debug
	//our_model.T = 8;

	double min_time = inf;
	double max_time = -inf;
	double cur_time, new_time, middle_time;
	double time_interval;

	for (int e = 0; e < our_model.m; e++)
	{
		vector<Point>& PS = our_model.edge_set[e].PS;

		for (int p = 0; p < PS.size(); p++)
		{
			if (PS[p].time < min_time)
				min_time = PS[p].time;
			if (PS[p].time > max_time)
				max_time = PS[p].time;
		}
	}

	time_interval = (max_time - min_time) / our_model.T;
	for (int t = 0; t < our_model.T; t++)
	{
		cur_time = min_time + t * time_interval;
		new_time = cur_time + time_interval;
		middle_time = (cur_time + new_time) / 2.0;

		our_model.time_vec.push_back(middle_time);
	}
}

void obtain_timestamps_user_input(int argc, char** argv, model& our_model)
{
	our_model.T = atoi(argv[8]); //number of timestamps

	for (int t = 0; t < our_model.T; t++)
		our_model.time_vec.push_back(atof(argv[9 + t]));
}

void obtain_timestamps_file(int argc, char** argv, model& our_model)
{
	fstream timestamp_file;
	double time;

	our_model.timestamp_fileName = argv[8];
	timestamp_file.open(our_model.timestamp_fileName);
	if (timestamp_file.is_open() == false)
	{
		cout << "Cannot open timestamp file!" << endl;
		exit(0);
	}

	timestamp_file >> our_model.T;
	for (int t = 0; t < our_model.T; t++)
	{
		timestamp_file >> time;
		our_model.time_vec.push_back(time);
	}

	timestamp_file.close();
}

void add_lixels_for_edge(int e, model& our_model)
{
	double cur_dist = 0;
	double middle_dist;
	double next_dist;
	double length = our_model.edge_set[e].length;
	Lixel l;

	while (cur_dist < length)
	{
		next_dist = cur_dist + our_model.lixel_reg_length;
		if (next_dist > length)
			next_dist = length;

		middle_dist = (cur_dist + next_dist) / 2.0;
		l.dist_n1 = middle_dist;
		l.dist_n2 = length - middle_dist;
		l.edge_index = e;

		for (int t = 0; t < our_model.T; t++)
			l.KDE_value_vec.push_back(0);
		
		our_model.lixel_set.push_back(l);

		cur_dist += our_model.lixel_reg_length;
	}
}

void obtain_lixel_set(model& our_model)
{
	for (int e = 0; e < our_model.m; e++)
		add_lixels_for_edge(e, our_model);
}