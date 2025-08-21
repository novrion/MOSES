#pragma once
#include "constants.h"
#include <vector>
#include <array>

namespace moses {

class Moses;
class Labour;
class Household;
class Bank;

class Government {
	friend class Moses;
	friend class Bank;
	friend class External;
	friend class Firm;
	friend class Household;
	friend class Labour;
	friend class Market;
	friend class Sector;
private:

	// Back References
	const Moses& moses;


	/* ===== Parameters ===== */

	// Taxes
    const std::array<double, MAXT> _qchtxva1; // (csv)
    const std::array<double, MAXT> _qchtxva2; // (csv)
    const std::array<double, MAXY> _txw; // (csv)
    const std::array<double, MAXY> _txwg; // (csv)
    const std::array<double, MAXY> _txc; // (csv)
    const std::array<double, MAXY> _txi1; // (csv)
    const double txva1; // (csv)
    const double txva2; // (csv)

	// Investment Distribution
    const std::array<double, NSEC> omegag; // (csv)
    const std::array<double, NSEC> omegabld; // (csv)
    const std::array<double, NSEC> omegain; // (csv)
    const std::array<double, NSEC> omega; // (csv)
    const std::array<double, NSEC> gkoff; // (csv)

	// References
    const double wgref; // (csv)

	// Subsidies
	const double rsubscash; // (csv)


	/* ===== Yearly Variables ===== */
	
	// Core
	double sp;
	double inv;
	double bw;
	double bwfor;
	double surplus;
	double mprint;
	double trans;
	double dep;
	double depfor;
	std::array<double, NSEC> purch{0};

	// Position
    double pos;    // (csv)
    double posfor; // (csv)

	// Labour & Wages
    double l;     // (csv)
	double w;     // (csv)
    double ws;    // (csv)

	// Taxes	
	double wtax;
	double itax;
	double vatax;
	double ctax;
	
	// Interest Receipts
	double intg;

	// Subsidies
	double subs;

	// Cumulative
	double cumwtax;
	double cumitax;
	double cumvatax;
	double cumctax;
	double cuminv;
	double cummprint;
	double cumws;
	double cuml;
	double cumtrans;
	double cumint;
	double cumsubs;
	std::array<double, NSEC> cumpurch{0};


	/* ===== Quarterly Variables ===== */

	// Core
	double qsp;
	double qws;
	double qsurplus;
	double qmprint;
	double qtrans;
	std::array<double, NSEC> qpurch{0};
	std::array<double, NSEC> qbuy{0};

	// Labour & Wages
    double qrealchl; // (csv)
    double qw;       // (csv)

	// Position
	double qchpos;
	double qchposfor; // (csv)

	// Taxes
    double qttax; // (csv)
	double qitax;
	double qwtax;
	std::array<double, NSEC> qvatax{0};
	double qvataximp;
	double qctax;

	// Interest Receipts
	double qint;
	double qintfor;
	
	// Investments
    double qinv;   // (csv)
    double qinvbld; // (csv)
    double qinvin;  // (csv)

	// Subsidies
	double qsubs;
	double qsubsdom;
	double qsubsfor;
	double qsubscash;

public:
	Government() = default;
	~Government() = default;

    Government(const Moses& moses_ref,
        const std::array<double, MAXT>& qchtxva1,
        const std::array<double, MAXT>& qchtxva2,
        const std::array<double, MAXY>& txw,
        const std::array<double, MAXY>& txwg,
        const std::array<double, MAXY>& txc,
        const std::array<double, MAXY>& txi1,
        const std::array<double, NSEC>& omegag,
        const std::array<double, NSEC>& omegabld,
        const std::array<double, NSEC>& omegain,
        const std::array<double, NSEC>& omega,
        const std::array<double, NSEC>& gkoff,
        double txva2, double txva1, double qrealchl,
        double l, double qw, double qttax, double wgref,
        double qinv, double qinvbld, double qinvin,
        double rsubscash, double pos, double posfor,
        double w, double ws, double qchposfor);

private:
	void initialise_extra(Bank& bank);
	void print_initialisation(const bool verbose = false);

	double qchtxva1() const;
	double qchtxva2() const;
	double txw() const;
	double txwg() const;
	double txc() const;
	double txi1() const;

	void yearly_init();
	void glabour(Labour& lab);
	void compute_consumption(const std::vector<double>& pt);
	void compute_buying(const std::vector<double>& pt);
	void accounting(Bank& bank, Household& hh);
};

} // namespace moses
