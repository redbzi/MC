#include <iostream>

#include "print.h"


int main(int argc, char *argv[]) {
    while (*++argv) {
        if ((*argv)[0] == '-') {
            switch ((*argv)[1]) {
                case '1': {
                    print_results_GBM();
                    print_results_Heston();
                    break;
                }
                case '2': {
                    print_graphs_GBM();
                    print_graphs_Heston();
                    break;
                }
                case '3': {
                    print_graphs_both();
                    break;
                }
                default: {
                    std::cerr << "Option does not exist." << std::endl;
                    exit(-1);
                }
            }
        } else {
            std::cerr << "Syntax error. Options must be preceded by \"-\"." << std::endl;
            exit(-1);
        }
    }

    return 0;
}
