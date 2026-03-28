#pragma once
#ifndef TIC_PHI_H
#define TIC_PHI_H

#include "Network.h"
#include "shortest_path.h"

bool check_validity(model& our_model);
void init_TIC_phi(model& our_model);
void PS_new_distribution_phi(model& our_model);
void PS_new_phi_clear(model& our_model);
void TIC_phi_algorithm(model& our_model);

void augment_preprocess_phi(model& our_model);
void bucket_construction_phi(model& our_model);
int bucket_find_position_phi(double value, model& our_model, Edge& e);
void bucket_clear_phi(model& our_model);

#endif