// TODO: indirect_taxes(): more branch? APL: if txva1 is not 0, do more branch.
// TODO: uppsats

#include "moses.h"
using namespace moses;

int main() {
	Moses moses;
	const int simulation_length = moses.start_simulation();
	const int verbose = 2; 
	/* --- verbose ---
	 * 0: minimal info
	 * 1: minimal initialisation info & end of simulation info
	 * 2: extended initialisation info & end of simulation info
	 * 3: maximum info
	 */
	
	moses.initialise(verbose);
	moses.simulate(simulation_length, verbose);

	return 0;
}
