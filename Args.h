#ifndef LhCDSCVX_ARGS_H
#define LhCDSCVX_ARGS_H

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include "getopt.h"

#define USAGE_TXT							   \
    "usage: \n"                                \
    "\t[-g input directed graph file]\n"       \
    "\t[-t repeated iterations]\n"             \
    "\t[-h h-cliques]\n"                       \
    "\t[-k topk]\n"                             \
    "\t[-p pattern]\n"                          \

class Args {
public:
    char *address{};
    char ds_address[50];
    int NT = 100;
    int topk = 5;
    int h = 3;
    int p = 0; // p=0 hcliques, p=1 4loops, p=2, doubletriangle, p=3, tritriangle, p=4, kstar;
    ~Args();
    void usage(char *msg, int exit_status);
    void parse_args(int argc, char *argv[]);
};


#endif //LhCDSCVX_ARGS_H
