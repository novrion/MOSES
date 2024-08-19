#include <stdio.h>
#include <vector>
#include <utility>

int a, b[3];
double c, d;

void * static_vars[] = {
	& a,
	& b,
	& c,
	& d
};

size_t static_var_sizes[] {
	sizeof(a),
	sizeof(b),
	sizeof(c),
	sizeof(d)
};

int static_var_types[] = {
	0,
	0,
	1,
	1
};

int id;
double ok[5];

void * vars[] = {
	& id,
	& ok
};
size_t var_sizes[] = {
	sizeof(id),
	sizeof(ok)
};
int var_types[] = {
	0,
	1
};


int n_static_vars = 4;
int n_instances = 1;
int n_vars = 2;




void read_var(FILE * in, const char * format, void * var, const int n) { // TODO FILE ** in ?
	void * tmp = var;
	for (int i = 0; i < n; i++) {
		if (format == "%i") var = (void *) ((int *) tmp + i);
		else var = (void *) ((double *) tmp + i);
		fscanf(in, format, var); // need to increment void pointer somehow but cannot do that lol!
	}
}

void read_file(FILE * in, const int n_static_vars, void * static_vars[], size_t static_var_sizes[], int static_var_types[], const int n_vars, void * vars[], size_t var_sizes[], int var_types[], const int n_instances) {
	
	for (int i = 0; i < n_static_vars; i++) {
		read_var(in, static_var_types[i] ? "%lf" : "%i", static_vars[i], (int) (static_var_types[i] ? static_var_sizes[i] / sizeof(double) : static_var_sizes[i] / sizeof(int)));
	}

	if (!n_vars) return;

	char newline;
	fscanf(in, "%c", &newline);

	for (int instance = 0; instance < n_instances; instance++) {
		for (int var = 0; var < n_vars; var++) {
			read_var(in, var_types[var] ? "%lf" : "%i", vars[var], (int) (var_types[var] ? var_sizes[var] / sizeof(double) : var_sizes[var] / sizeof(int)));
		}

		fscanf(in, "%c", &newline);
	}
}





//void get_var(FILE * in, FILE * out, const char * format, size_t size) {
//	double db;
//	int nt;
//	fscanf(in, format, size == sizeof(int) ? &nt : &db);
//	fwrite(size == sizeof(int) ? &nt : &db, size, 1, out);
//}

void get_int(FILE * in, FILE * out) {
	int var;
	fscanf(in, "%i", &var);
	fwrite(&var, sizeof(int), 1, out);
}
void get_double(FILE * in, FILE * out) {
	double var;
	fscanf(in, "%lf", &var);
	fwrite(&var, sizeof(double), 1, out);
}

void compile_bin(FILE * in, FILE * out, const std::vector<std::pair<int, int>> instructions) {
	for (auto & instruction : instructions) {
		if (instruction.first == -1) {
			char newline;
			for (int i = 0; i < instruction.second; i++) {
				fscanf(in, "%c", &newline);
			}
		} else if (instruction.first == 0) {
			int var;
			for (int i = 0; i < instruction.second; i++) {
				fscanf(in, "%i", &var);
				fwrite(&var, sizeof(int), 1, out);
			}
		} else {
			double var;
			for (int i = 0; i < instruction.second; i++) {
				fscanf(in, "%lf", &var);
				fwrite(&var, sizeof(double), 1, out);
			}
		}
	}
}


//	for (int type = 0; type < types.size(); type++) {
//		for (int i = 0; i < types[type]; i++) {
//			type ? get_double(in, out) : get_int(in, out);
//		}
//	}
//}

















int main() {

	FILE * in = fopen("testing", "r");
	FILE * out = fopen("out", "w");

	read_file(in, n_static_vars, static_vars, static_var_sizes, static_var_types, n_vars, vars, var_sizes, var_types, n_instances);

	printf("%i [ %i %i %i ] %lf %lf\n", a, b[0], b[1], b[2], c, d);
	for (int i = 0; i < 5; i++) {
		printf("%lf ", ok[i]);
	}
	printf("\nid: %i\n", id);

	const std::vector<std::pair<int, int>> instructions = { {0, 2}, {1, 2} };
	compile_bin(in, out, instructions);



	return 0;
}


