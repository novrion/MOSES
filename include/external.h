#pragma once
#include "constants.h"
#include "sector.h"
#include <vector>
#include <array>

namespace moses {

class Moses;
class Market;
class Government;

class External : public Sector {
	friend class Moses;
	friend class ui;
private:

	/* ===== Parameters ===== */

	// Prices
    const std::array<double, MAXT> _qdp; // (csv)
	
	// Export
    const double x; // (csv)


	/* ===== Variables ===== */

	// Core
	double qq;
	double qva;

public:
	External() = default;
	~External() = default;

    External(Moses& moses_ref,
			int id_, double qpdom_, double qpfor_, 
            double imp_, double tmimp_, double x_, double pref_,
            const std::array<double, NSEC>& io_,
            const std::array<double, MAXT>& qdp);

private:
	void print_initialisation(const bool verbose = false);

	double qdp() const;
	double qpfor_() const;
	
    void quarterly_exp(Government& g) override;
    double compute_trial_price();
    void compute_buying();
    void compute_imports(Government& g, const std::vector<double>& pt);
    void domestic_result(const std::vector<double>& pt) override;
    void external_sectors(Government& g, Household& hh, const std::vector<double>& pt);
	void national_accounting(std::array<double, NGNPCUR>& qgnpcur, std::array<double, NGNPFIX>& qgnpfix) override;

    // Helper Methods
    double qimpurch(Government& g, const std::vector<double>& pt) const;
};

} // namespace moses
