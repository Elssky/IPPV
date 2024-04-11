#pragma once

#ifndef KCLISTCORE_H
#define KCLISTCORE_H



#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <vector>
#include <unordered_map>

std::vector<std::vector<int>>get_kcliques(std::vector<std::vector<int>> adj, std::vector<std::vector<int>>& v2cliques, unsigned k_max, std::vector<bool>selected);
std::vector<std::vector<int>>get_kloops(std::vector<std::vector<int>> adj, std::vector<std::vector<int>>& v2cliques, unsigned k_max, std::vector<bool>selected);
std::vector<std::vector<int>>get_doubletriangle(std::vector<std::vector<int>> adj, std::vector<std::vector<int>>& v2cliques, unsigned k_max, std::vector<bool>selected);
std::vector<std::vector<int>>get_tritriangle(std::vector<std::vector<int>> adj, std::vector<std::vector<int>>& v2cliques, unsigned k_max, std::vector<bool>selected);
std::vector<std::vector<int>>get_kstar(std::vector<std::vector<int>> adj, std::vector<std::vector<int>>& v2cliques, unsigned k_max, std::vector<bool>selected);
std::vector<std::vector<int>>get_ctriangle(std::vector<std::vector<int>> adj, std::vector<std::vector<int>>& v2cliques, unsigned k_max, std::vector<bool>selected);


#endif