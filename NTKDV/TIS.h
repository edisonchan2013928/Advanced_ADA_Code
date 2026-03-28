#pragma once
#ifndef TIS_H
#define TIS_H

#include "Network.h"
#include "shortest_path.h"

void construct_intervals(model& our_model);
void init_TIS(model& our_model);
void PS_new_distribution_sharing(model& our_model);
void PS_new_sharing_clear(model& our_model);
void interval_clear(model& our_model);

void update_density_g(model& our_model);
void window_density(model& our_model);
void TIS_algorithm(model& our_model);

void augment_preprocess_sharing(model& our_model);
void bucket_construction_sharing(model& our_model);
int bucket_find_position_sharing(double value, model& our_model, Edge& e);
void bucket_clear_sharing(model& our_model);

#endif