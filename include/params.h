#pragma once
#include "constants.h"

namespace moses {

struct Params {
	// Moses
    int niter;
    int marketiter;

    // Firm
    double e1;
    double e2;
    double eps;
    double r;
    double smp;
    double sms;
    double smt;
    double smw;
    double fip;
    double fis;
    double fiw;
    double expxdp;
    double expxds;
    double expxdw;
    double loss;
    double rho;
    double resdown;
    double resmax;
    double tmsto;
    double tmimsto;
    double beta;
    double imbeta;
    double wtix;
    double iota;
    double rhobook;
    double rtd;
    double alfabw;
    double betabw;
    double redchbw;
    double elinv;
    double utref;

    // Market
    double maxdp;

	Params(int niter_, int marketiter_,
	        double e1_, double e2_, double eps_,
	        double r_, double smp_, double sms_, double smt_, double smw_,
	        double fip_, double fis_, double fiw_,
	        double expxdp_, double expxds_, double expxdw_,
	        double loss_, double rho_, double resdown_, double resmax_,
	        double tmsto_, double tmimsto_,
	        double beta_, double imbeta_, double wtix_, double iota_,
	        double rhobook_, double rtd_,
	        double alfabw_, double betabw_, double redchbw_,
	        double elinv_, double utref_,
		    double maxdp_);
};

} // namespace moses
