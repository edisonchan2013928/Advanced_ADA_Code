#pragma once
#ifndef SHORTEST_PATH_H
#define SHORTEST_PATH_H

//#include "init.h"
#include "Network.h"

typedef priority_queue<sp_node, vector<sp_node>, comparePriority> PQ;

void init_dijkstra(model& our_model, PQ& pq, Lixel& cur_l);
void dijkstra(model& our_model, Lixel& cur_l);
void copy_sp_info(model& our_model, bool isFirst);
void clear_scan_edge_list(model& our_model);

#endif