#pragma once
#include "constants.h"
#include <vector>
#include <array>

namespace moses {

struct Params;
struct LabourSearchData;
class Market;
class Labour;
class Bank;
class Government;
class Moses;

class Firm {
	friend class Moses;
	friend class Bank;
	friend class Government;
	friend class Household;
	friend class Labour;
	friend class Market;
	friend class ui;
private:

	// Back Reference
    const Moses& moses;


	/* ===== Firm Wide Parameters ===== */
	const double e1;	  // (csv)
	const double e2;	  // (csv)
	const double eps;	  // (csv)
	const double r;		  // (csv)
	const double smp;	  // (csv)
	const double sms;	  // (csv)
	const double smt;	  // (csv)
	const double smw;	  // (csv)
	const double fip;	  // (csv)
	const double fis;	  // (csv)
	const double fiw;	  // (csv)
	const double expxdp;  // (csv)
	const double expxds;  // (csv)
	const double expxdw;  // (csv)
	const double loss;	  // (csv)
	const double rho;	  // (csv)
	const double resdown; // (csv)
	const double resmax;  // (csv)
	const double tmsto;	  // (csv)
	const double tmimsto; // (csv)
	const double beta;	  // (csv)
	const double imbeta;  // (csv)
	const double wtix;	  // (csv)
	const double iota;	  // (csv)
	const double rhobook; // (csv)
	const double rtd;	  // (csv)
	const double alfabw;  // (csv)
	const double betabw;  // (csv)
	const double redchbw; // (csv)
	const double elinv;   // (csv)
	const double utref;	  // (csv)


	/* ===== Firm Specific Parameters ===== */
	
	// Identification
	const int id;        // (csv)
	const int market_id; // (csv)

	// Inventory
	const double small;   // (csv)
	const double imsmall; // (csv)
	const double big;     // (csv)
	const double imbig;   // (csv)

	// Production & Planning
	const double share; // (csv)

	// Investment & Capital Structure
	const double rw; // (csv)


	/* ===== Yearly Variables ===== */

	// Core
	double m;   // (csv)
	double chm; // (csv)
	double p;   // (csv)
	double dp;  // (csv)
	double s;   // (csv)
	double ds;  // (csv)
	double snet;
	double w;   // (csv)
	double dw;  // (csv)
	double q;   // (csv)
	double dq;
	double va;  // (csv)
	double dva;

	// Inventory    
	double sto; // (csv)
	double minsto;
    double maxsto;
    double optsto;
	std::array<double, NSEC> imsto; // (csv)
    std::array<double, NSEC> minimsto{};
    std::array<double, NSEC> maximsto{};
    std::array<double, NSEC> optimsto{};	

	// Export
	double x; // (csv)

	// Labour
	double l; // (csv)
	std::array<double, 3> aman; // (csv)

	// Expectations
	double expdp; // (csv)
	double expds; // (csv)
	double expdw; // (csv)

	// History
	double mhist; // (csv)
	double histdp, histdpdev, histdpdev2; // (csv)
	double histds, histdsdev, histdsdev2; // (csv)
	double histdw, histdwdev, histdwdev2; // (csv)

	// Target
	double targm;

	// Cumulative
	double cumq;
    double cumm;
    double cumsu;
    double cums;
    double cump;
    double cumws;
    double cuml;
    double cuminv;
    double cumva;
    double cumsnet;

	// Investment & Capital Structure
	double nw;
	double bw;     // (csv)
	double k1;     // (csv)
	double k1book; // (csv)
	double k2;     // (csv)
   	double k3;
	double k3finish;

	double inveff; // (csv)

    std::array<double, 3> inv_delay{};


	/* ===== Quarterly Variables ===== */

	// Core
	double qm;
	double qp;  // (csv)
	double qdp;
	double qs;  // (csv)
	double qsnet;
	double qds;
	double qw;  // (csv)
	double qdw;
	double qq;  // (csv)
	double qva; // (csv)
	double qdva;

	// Production Planning
	double qtop; // (csv)
	double tec;  // (csv)
	double res;  // (csv)

	double qplanq;
    double qplanqsave;
	double qplanl;

	double qoptsu;
	double qoptsudom;

	double qexpsu;
	double qexpp;
	double qexpdp;
    double qexps;
    double qexpw;	

	// Production
	double qsufor;
	double qsfor;
	double qsudom;
	double qsdom;
	double qsu;

	// Target
	double qtargm;

	// Labour
	double qchl;
	double qchw;

	// Material Inputs
    std::array<double, NSEC> qimq{};

	// Investment & Capital Structure
	double qinv;    // (csv)
	double qinvlag; // (csv)
	
    double qdeschbw;
    double qchbw;

	double qdeschk2;
	double qchk2;

    double qrev;
    double qrr;
    double qcash;

    double qdepr;
    double qdeprbook;

    double qtax;
    double qdiv;


	/* ===== Helper Variables ===== */

    int bad{0};

public:
	Firm() = default;
	~Firm() = default;

    Firm(const Moses& moses_ref, const Params& params,
         int id_, int market_id_,
         double m_, double mhist_,
         double expdp_, double expds_, double expdw_,
         double dp_, double ds_, double dw_,
         double histdp_, double histdpdev_, double histdpdev2_,
         double histds_, double histdsdev_, double histdsdev2_,
         double histdw_, double histdwdev_, double histdwdev2_,
         double qp_, double qs_, double qw_,
         double qq_, double l_, double qtop_,
         double tec_, double qinv_, double inveff_,
         double res_, double sto_, double s_, double p_,
         double small_, double big_, double share_,
         double chm_, double imsmall_, double imbig_,
         double x_, double qinvlag_, double qva_,
         double q_, double va_, double w_,
         double bw_, double k1_, double k2_,
         double k1book_, double rw_,
         const std::array<double, 3>& aman_,
         const std::array<double, NSEC>& imsto_);

private:
	// Initialisation
	void print_initialisation(const bool print_params = false, const bool verbose = false);
    void initialise_extra(const Market& market);
    
	void yearly_init();
    void yearly_exp();
    void yearly_targ(); 
    void yearly_update();

    void quarterly_exp(bool first_quarter);
    void quarterly_targ();
    void quarterly_cum();

    void luupdate(Labour& lab);
    void prodfront(const Market& market);
    void init_prodplan();
    bool target_search();
    
    void labour_search_input(LabourSearchData& data);
    void handle_labour_search_layoffs();
    void handle_labour_update_layoffs(Labour& lab);
    void update_labour_force_and_wage();
    
    void planqrevise(const Market& market);
    void adjust_foreign(const Market& market, Government& g);
    void reference_inventory_levels(const Market& market);
    void finalqpqsqm();
    
    void invfin(Bank& bank, Government& g);
    bool invfin_adjustments();
	
    // Helper
	double qexppnet() const;
	double sum_io_qpdom_txva2() const;
	double sum_qimq_qpdom_txva2() const;
	double k3imed() const;

    bool try_reduce_production(double qexppnet);
    bool try_increase_production(double qexppnet);
    bool plan_implies_recruitment(double qexppnet);
    bool try_reduce_slack(double q7, double qexppnet);
    void exit_successfully(double qexppnet);
    void handle_layoffs();

    double qfr(double l);
    double rfq(double q);
    bool sat(double q, double l, double qexppnet);
    void solve(double qexppnet);
    double curs() const;
    double qcurs() const;
    double curp() const;
    double check_bw(int bw_ratio);
    double delay();
};

} // namespace moses
