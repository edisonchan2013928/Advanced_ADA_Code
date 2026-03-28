#include "alg_NTKDV.h"

int main(int argc, char** argv)
{
	alg_NTKDV algorithm;
	algorithm.load_network(argv);
	algorithm.NTKDV_compute_file(argc, argv);
	//algorithm.NTKDV_compute(argc, argv);
	algorithm.clear_basic_memory();
}