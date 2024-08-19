#include "moses.h"
#include "data.h"


int main(int argc, char ** argv) {

	// TODO
	//if (argc <= 1) return 1;
	//char * path = argv[1];

	static FILE * out = stdout;
	const char * path = "data/bin";

	std::vector<Firm>      firms;
	std::vector<Household> households;
	External               externals[NEX];
	Market                 markets[NMKT];
	Bank                   bank;
	Government             government;
	Labour                 labour;


	fprintf(out, "[ Initialising MOSES ]\n");
	if (Moses::initialise(firms, households, externals, markets)) {
		fprintf(out, "[ Failed MOSES initialisation ]\n");
		return 1;
	}
	fprintf(out, "[ Successfully initialised MOSES ]\n");


	fprintf(out, "[ Starting MOSES simulation... ]\n");
	int simulation_length = 30;
	//Moses::dump_initialization(out, simulation_length);
	if (Moses::simulate(out, simulation_length, firms, households, externals, markets, bank, government, labour)) {
		fprintf(out, "[ Failed MOSES simulation ]\n");
		return 2;
	}
	fprintf(out, "[ Finished MOSES simulation successfully ]\n");



	return 0;
}
