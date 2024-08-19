#include "data.h"


void read_var(FILE * in, const char * format, void * var, const int n) {
	void * tmp = var;
	for (int i = 0; i < n; i++) {
		if (format == "%i") var = (void *) ((int *) tmp + i);
		else var = (void *) ((double *) tmp + i);
		fscanf(in, format, var); // need to increment void pointer somehow but cannot do that lol!

//		if (format == "%i") printf("%i ", *(int *)var);
//		else printf("%lf ", *(double *)var);
	}
}

void read_file(FILE * in, const int n_vars, void * vars[], const size_t var_sizes[], const int var_types[]) {
	//printf("new session:\n");
	for (int var = 0; var < n_vars; var++) {
		if (var_types[var]) {
//			printf("double %i:  ", var);
			read_var(in, "%lf", vars[var], var_sizes[var] / sizeof(double));
		} else {
//			printf("int %i:  ", var);
			read_var(in, "%i", vars[var], var_sizes[var] / sizeof(int));
		}
	}
}
