#pragma once
#include "constants.h"
#include <string>
#include <vector>
#include <array>
#include <memory>
#include <stdexcept>
#include <numeric>

namespace moses {
	
class Bank;
class Firm;
class Government;
class Household;
class Labour;
class Market;
class External;

// Stores help variables to labour_search()
struct LabourSearchData {
    std::vector<double> chl;
    std::vector<double> ww;
    std::vector<double> ll;
};

class Moses {
	friend class Bank;
	friend class Firm;
	friend class Government;
	friend class Household;
	friend class Labour;
	friend class Market;
	friend class External;
	friend class ui;
private:

	/* ===== MOSES actors ===== */
	std::unique_ptr<Bank> bank;
    std::unique_ptr<Government> gov;
    std::unique_ptr<Household> hh;
    std::unique_ptr<Labour> lab;
    std::vector<std::unique_ptr<Firm>> firms;
    std::vector<std::unique_ptr<Market>> markets;
    std::vector<std::unique_ptr<External>> externals;


	/* ===== Parameters ===== */
	int niter;      // (csv)
	int marketiter; // (csv)


	/* ===== Time ===== */
    int year{0};
    int quarter{0};
    int relative_quarter{0};


	/* ===== Yearly Variables ===== */

	// National Accounting
	std::array<double, NGNPCUR> gnpcur;
	std::array<double, NGNPFIX> gnpfix;
	double export_;
	double import_;


	// Cumulative
	std::array<double, NGNPCUR> cumgnpcur;
	std::array<double, NGNPFIX> cumgnpfix;
	double cumimport;
	double cumexport;


	/* ===== Quarterly Variables ===== */

	// National Accounting
	double qimport;
	double qexport;

	// Firm Aggregate Interest Receipts
	double qintf;
	double qintk2;

public:
	Moses();
	~Moses();

	int start_simulation();
	void initialise(const int verbose = 0);
	void simulate(const int simulation_length, const int verbose = 0);

private:

	// UI
	void print_initialisation(const int verbose = 0);
	void print_simulation();

	// Initialisation
	void initialise_extra();

	// Main Simulation
	void simulate_year();
	void simulate_quarter();

	// Yearly init
	void yearly_init();

	// Quarterly Expectations
	void quarterly_exp();

	// Production Planning
	void luupdate();
	void prodfront();
	void prodplan();

	// Labour Market
	void handle_successful_attack(LabourSearchData& data, size_t attacker, size_t target);
	LabourSearchData labour_search_input();
	void confront(LabourSearchData& data);
	void labour_search_output(LabourSearchData& data);
	void labour_search();
	void labour_update();
	void indalabour();
	void labour_market();

	// Export Market
	void export_market();

	// Domestic Market
	std::vector<double> compute_trial_prices();
	void market_entrance();
	void compute_external_production();
	void compute_total_buying();
	void compute_buying(const std::vector<double>& qsp, const std::vector<double>& pt);
	void market_confront(std::vector<double>& qsp, std::vector<double>& pt);
	void compute_imports(const std::vector<double>& pt);
	void domestic_result(const std::vector<double>& pt);
	void external_sectors(const std::vector<double>& pt);
	void indirect_taxes(const std::vector<double>& qsp);
	void domestic_market();

	// Inventory System
	void firm_sto();
	void sto_system();

	// Investment & Financing
	void invfin();

	// Monetary Sector
	void invfin_adjustments();
	void monetary_sector();

	// National Accounting
	void national_accounting();

	// IO
	void print_year();

	// Helper Methods
	void nullify_firm(const size_t index);
	void check_market_health(const bool verbose = false);
	size_t choose_labour_target(const LabourSearchData& data);
	std::vector<std::vector<double>> invert_matrix(std::vector<std::vector<double>> m);

	double qdpk() const;

	// Total quarterly wage sum (firm)
	double qwsf() const;

	// Total wage sum (firm) - txw adjusted
	double wstx() const;

	// Sum across a vector
	inline double sum(std::vector<double>& vec) const;

	// Sum across an array
	template<std::size_t N>
	inline double sum(const std::array<double, N>& ar) const;

	// Sum across a vector per market
	inline double sum(const std::vector<double>& vec, const int market_id) const;

	// Sum one variable across firms
	template<typename T>
	inline double sum(T Firm::*member) const;

	// Sum one variable across markets
	template<typename T>
	inline double sum(T Market::*member) const;

	// Sum one variable across firms in a sector
	template<typename T>
	inline double sum(T Firm::*member, const int market_id) const;

	// Weighted average sum across firms
	template<typename WeightT, typename ValueT>
	inline double sum_avg(WeightT Firm::*weight, ValueT Firm::*value) const;

	// Weighted average sum across firms per market
	template<typename WeightT, typename ValueT>
	inline double sum_avg(WeightT Firm::*weight, ValueT Firm::*value, const int market_id) const;

	// Sum the multiplication across two vectors
	inline double sum_mult(const std::vector<double>& v1, const std::vector<double>& v2) const;

	// Sum the multiplication across an array and a vector
	template<std::size_t N>
	inline double sum_mult(const std::array<double, N>& ar, const std::vector<double> vec) const;
	
	// Sum the multiplication between two variables across firms
	template<typename T1, typename T2>
	inline double sum_mult(T1 Firm::*member1, T2 Firm::*member2, const int market_id) const;

	// Sum the multiplication between a variable across firms per market, and a specific variable
	template<typename T1>
	inline double sum_mult(T1 Firm::*member, const double var, const int market_id) const;

	// Sum the ratio between two variables across an array and a vector
	template<std::size_t N>
	inline double sum_div(const std::array<double, N>& ar, const std::vector<double>& vec) const;

	// Sum the ratio between two variables across two arrays
	template<std::size_t N1, std::size_t N2>
	inline double sum_div(const std::array<double, N1>& ar1, const std::array<double, N2>& ar2) const;

	// Sum qbuy from a specific sector
	double sum_qbuy(const int mkt_id) const;

	// Sum qdeschbw over firms, but only take positive values
	double sum_qdeschbw_positive() const;
};

} // namespace moses
