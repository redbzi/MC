#include "print.h"

int main(int argc, char *argv[]) {
    print_results_GBM();
    print_results_Heston();

    print_graphs_GBM();
    print_graphs_Heston();
    print_graphs_both();
    return 0;
}