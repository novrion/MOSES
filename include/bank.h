#pragma once
#include "constants.h"
#include <vector>
#include <array>

namespace moses {

class Moses;
class Government;
class Household;

class Bank {
	friend class Moses;
	friend class Firm;
	friend class Government;
	friend class Household;
	friend class Sector;
	friend class ui;
private:
    
	// Back References
	const Moses& moses;


	/* ===== Parameters ===== */
	
	// Interest Rates
    const std::array<double, MAXT> _ridepfor; // (csv)
    const std::array<double, MAXT> _ribwfor;  // (csv)
    const double mb;        // (csv)
    const double minri;     // (csv)
    const double maxri;     // (csv)
    const double maxqchri;  // (csv)
    const double maxridiff; // (csv)

	// Liquidity Restriction
	const double rfund1; // (csv)
	const double rfund2; // (csv)

	// Timing
    const double tmfass; // (csv)
    const double tmfd;   // (csv)

	// Market Adjustment
	const double kappa1; // (csv)
    const double kappa2; // (csv)
    const double lamda1; // (csv)
    const double lamda2; // (csv)


	/* ===== Variables ===== */

	// Interest Rates   
    double ri;    // (csv)
    double qchri; // (csv)
	double rik2;
	double rif;
    double rih;
	double ribwg;
	double ridepg;
	double ribwgfor;
	double ridepgfor;

	// Balance Sheet
    double liqb;    // (csv)
   	double qchliqb;
    double liqbfor; // (csv)
    double nw;      // (csv)
    double fass;    // (csv)
    double fd;      // (csv)
   
	// Bank Operations
	double money;
   	std::array<double, NSEC> qbuy{0};

	// Credit Market
   	double qsupfund;
   	double qtchbw;
   	double qtchinv;
   	double qtchk2;
    
public:
	Bank() = default;
	~Bank() = default;

    Bank(const Moses& moses_ref,
		const std::array<double, MAXT>& ridepfor,
		const std::array<double, MAXT>& ribwfor,
		double mb, double rfund1, double rfund2,
        double liqb, double liqbfor, double nw,
        double fass, double tmfass, double fd, double tmfd,
        double ri, double qchri, double maxri, double maxqchri,
        double minri, double maxridiff,
        double kappa1, double kappa2, double lamda1, double lamda2);

private:
	void initialise_extra(Government& g, Household& hh);
	void print_initialisation(const bool verbose = false);

	double ridepfor() const;
	double ribwfor() const;

	void compute_buying(Government& g, const std::vector<double>& pt);

	void bank_transactions(Government& g, Household& hh);
   	void credit_market(Government& g, Household& hh);
   	void bank_update(Government& g, Household& hh);
};

} // namespace moses
