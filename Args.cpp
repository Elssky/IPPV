#include "Args.h"
#include <string>
#include <string.h>

Args::~Args() {
    delete address;
}

void Args::usage(char *msg, int exit_status) {
    fprintf(exit_status == 0 ? stdout : stderr, "%s", USAGE_TXT);

    if (msg) {
        fprintf(exit_status == 0 ? stdout : stderr, "\n%s\n", msg);
    }
    exit(exit_status);
}

void Args::parse_args(int argc, char **argv) {
    int c;
    opterr = 0;
    if (argc < 2) {
        usage(nullptr, 0);
    }

    while ((c = getopt(argc, argv, "g:t:h:k:p:")) != -1) {
        switch (c) {
            case 'g': {
                printf("%s\n", optarg);
                address = _strdup(optarg);
                break;
            }

            case 't': {
                sscanf(optarg, "%d", &NT);
                break;
            }
            case 'h': {
                sscanf(optarg, "%d", &h);
                break;
            }
            case 'k': {
                sscanf(optarg, "%d", &topk);
                break;
            }
            case 'p': {
                sscanf(optarg, "%d", &p);
            }


            default:
                break;
        }
    }

    std::string data_name(address+10, strlen(address)-14);

    sprintf(ds_address, "output/%s-%d.txt", data_name.c_str(), topk);

    printf("%s\n", ds_address);
}
