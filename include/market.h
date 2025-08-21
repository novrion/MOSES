#pragma once
#include "constants.h"
#include "sector.h"
#include <vector>
#include <array>

namespace moses {

class Moses;
class Sector;
class Params;

class Market : public Sector {
	friend class Moses;
	friend class External;
	friend class Firm;
	friend class ui;
private:

	/* ===== Market Wide Parameters */

	const double maxdp; // (csv)
	

	/* ===== Market Specific Parameters ===== */

	// Prices
	const std::array<double, MAXT> _qdpfor; // (csv)

	// Technology
	const double qdmtec; // (csv)
    
	// Timing
	const double tminv; // (csv)
    const double tmx;   // (csv)

	// Subsidies
	const double rsubs; // (csv)


	/* ===== Yearly Variables ===== */
   
	// Technology 
	double mtec; // (csv)

	// Inventory
	double tstocurf; // (csv)
	double tstocurm; // (csv)


	/* ===== Quarterly Variables ===== */
    
	// Prices
	double qprelpdom;

	// Sales
	double qtsudom;

	// Inventory
    double qchtsto;
    double qwaste;
	double qchtstocurf;
	double qchtstocurm;

public:
	Market() = default;
	~Market() = default;

    Market(Moses& moses_ref, const Params& params,
			int id_, double qpdom_, double qpfor_, 
            double mtec_, double qdmtec_, double tmx_, 
            double rsubs_, double imp_, double tmimp_, 
            double tminv_, double pref_,
			double tstocurf_, double tstocurm_,
            const std::array<double, NSEC>& io_,
            const std::array<double, MAXT>& qdpfor);

private:
	void print_initialisation(const bool print_params = false, const bool verbose = false);

	double qdpfor() const;

	void quarterly_exp(Government& g) override;
    void market_entrance();
    double compute_trial_price(Government& g);
    void compute_buying();
    void price_adjust(std::vector<double>& pt);
    void compute_imports(const std::vector<double>& pt);
    void domestic_result(const std::vector<double>& pt) override;
    void firm_sto(Government& g,
			const std::vector<double>& limsto, 
			const std::vector<double>& upper, 
			const std::vector<double>& lower);
	void national_accounting(std::array<double, NGNPCUR>& qgnpcur, std::array<double, NGNPFIX>& qgnpfix) override;
};

} // namespace moses
