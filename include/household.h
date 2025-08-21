#pragma once
#include "constants.h"
#include <vector>
#include <array>

namespace moses {

class Moses;
class Government;
class Labour;
class Bank;

class Household {
	friend class Moses;
	friend class Bank;
	friend class Government;
	friend class Sector;
	friend class External;
	friend class ui;
private:

	// Back References
	const Moses& moses;


	/* ===== Parameters ===== */

	// Core
    const double nh; // (csv)

	// Durable Depreciation
    const double rhodur; // (csv)

	// Transfer Payment Ratio
    const double rtrans; // (csv)

	// Expenditure
	const double alfa3; // (csv)
    const double alfa4; // (csv)
	const std::array<double, NEXPH> beta1;  // (csv)
    const std::array<double, NEXPH> beta2;  // (csv)
	const std::array<double, NEXPH> beta3;  // (csv)
	const std::array<double, NEXPH> smooth; // (csv)


	/* ===== Yearly Variables ===== */

	// Assets
    double wh;     // (csv)
    double whra;   // (csv)
    double stodur; // (csv)

	// Consumption
    std::array<double, NSEC> cva; // (csv)


	/* ===== Quarterly Variables ===== */

	// Income
	double qti;
	double qdi;
    double qinpay; // (csv)
    double qtdiv;  // (csv)

	// Prices
    double qcpi;  // (csv)
    double qdcpi; // (csv)
    std::array<double, NSEC> qp; // (csv)

	// Consumption
    std::array<double, NSEC> qc;  // (csv)
	std::vector<double> qsp;
	std::array<double, NSEC> qbuy{0};

	// Saving
	double qsav;
    double qsavhreq; // (csv)
	double qspsavreq;

	// Interest Receipts
	double qint;

public:
	Household() = default;
	~Household() = default;

    Household(const Moses& moses_ref,
        double qinpay, double qtdiv, double nh,
        double qsavhreq, double rtrans,
        double wh, double whra, double stodur, double rhodur,
        double qcpi, double qdcpi, double alfa3, double alfa4,
        const std::array<double, NEXPH>& beta1,
        const std::array<double, NEXPH>& beta2,
        const std::array<double, NEXPH>& beta3,
        const std::array<double, NSEC>& qp,
        const std::array<double, NSEC>& cva,
        const std::array<double, NSEC>& qc,
        const std::array<double, NEXPH>& smooth);

private:
	void print_initialisation(const bool verbose = false);

	void init(Bank& bank, Government& g, Labour& lab);
    std::vector<double> compute_essential_spending(Bank& bank, Labour& lab, const std::vector<double>& pt, const double qchdcpi);
    std::vector<double> compute_spending(const std::vector<double>& qspe, const double qprelcpi);
    std::vector<double> compute_expenditures(const std::vector<double>& pt);
    void compute_buying(const std::vector<double>& qsp, const std::vector<double>& pt);
    void household_update(std::vector<double>& qsp, const std::vector<double>& pt);
};

} // namespace moses
