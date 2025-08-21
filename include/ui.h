#pragma once

namespace moses {

class Moses;

class ui {
public:
	static int start_simulation();
	static void print_simulation(const Moses& moses);
	static void print_year(const Moses& moses);
private:
	static int get_simulation_length();
};

} // namespace moses
