#pragma once
#ifndef KAF_H
#define KAF_H

#include "Network.h"

double s_kernel_value(model& our_model, double s_dist);
double t_kernel_value(model& our_model, double t_dist);
void NTKDV_basic(model& our_model, Lixel& cur_l);

#endif