#pragma once
#ifndef TIC_H
#define TIC_H

#include "Network.h"
#include "shortest_path.h"

void PS_new_construction(model& our_model);
void PS_new_clear(model& our_model);
double edge_KAF_TIC(model& our_model, Lixel& cur_l, int edge_index);
void TIC_algorithm(model& our_model);

//Used in aggregate distance augmentation (ADA_{TIC}^{B+})
void augment_preprocess(model& our_model);
void bucket_construction(model& our_model);
int bucket_find_position(double value, model& our_model, Edge& e);
void bucket_clear(model& our_model);

#endif