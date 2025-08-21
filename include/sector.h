#pragma once
#include "constants.h"
#include <vector>
#include <array>

namespace moses {
	
class Moses;
class Bank;
class Government;
class Household;

class Sector {
	friend class Moses;
	friend class Firm;
protected:
		    
	// Back References
	Moses& moses;	

	/* ===== Parameters ===== */
	
	// Identification
	const int id; // (csv)

	// Timing
    const double tmimp; // (csv)

	// Input-Output
    const std::array<double, NSEC> io; // (csv)

	// Price Reference
	const double pref; // (csv)


	/* ===== Variables ===== */

	// Prices
    double qpdom; // (csv)
    double qdpdom;
    double qpfor; // (csv)
    double qexppim;

	// Buying
    std::array<double, NSEC> qbuy;
    double qtbuy;
    double qtbuydom;
    double qtbuyfor;

	// Import
    double imp; // (csv)


    Sector(Moses& moses_ref,
			int id_, double qpdom_, double qpfor_, 
            double imp_, double tmimp_,
			double pref_,
            const std::array<double, NSEC>& io_);

	Sector() = default;
    ~Sector() = default;

    virtual void quarterly_exp(Government& g) {};
    virtual void domestic_result(const std::vector<double>& pt) {};
	virtual void national_accounting(std::array<double, NGNPCUR>& qgnpcur, std::array<double, NGNPFIX>& qgnpfix) {};

    void indirect_taxes(Bank& bank, Government& g, Household& hh, const std::vector<double>& qsp);
};

} // namespace moses
