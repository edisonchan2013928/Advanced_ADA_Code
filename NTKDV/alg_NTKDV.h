#pragma once
#ifndef ALG_NTKDV_H
#define ALG_NTKDV_H

#include "shortest_path.h"
#include "KAF.h"
#include "TIC.h"
#include "TIC_phi.h"
#include "TIS.h"

class alg_NTKDV
{
public:
	void init_parameters(int argc, char** argv);
	void load_network(char** argv);
	void NTKDV_algorithm();
	void NTKDV_compute_file(int argc, char** argv);
	string NTKDV_compute(int argc, char** argv);
	void output_Visual(); //Used for outputting the file
	string output_Visual_String();
	void clear_memory();
	void clear_basic_memory();

private:
	model our_model;
};

void obtain_timestamps_min_max_range(int argc, char** argv, model& our_model);
void obtain_timestamps_user_input(int argc, char** argv, model& our_model);
void obtain_timestamps_file(int argc, char** argv, model& our_model);
void add_lixels_for_edge(int e, model& our_model);
void obtain_lixel_set(model& our_model);

#endif