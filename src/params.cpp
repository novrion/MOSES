#include "params.h"

namespace moses {

Params::Params(
		int niter_, int marketiter_,
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
		double maxdp_)
	: niter(niter_)
	, marketiter(marketiter_)
	, e1(e1_)
	, e2(e2_)
	, eps(eps_)
	, r(r_)
	, smp(smp_)
	, sms(sms_)
	, smt(smt_)
	, smw(smw_)
	, fip(fip_)
	, fis(fis_)
	, fiw(fiw_)
	, expxdp(expxdp_)
	, expxds(expxds_)
	, expxdw(expxdw_)
	, loss(loss_)
	, rho(rho_)
	, resdown(resdown_)
	, resmax(resmax_)
	, tmsto(tmsto_)
	, tmimsto(tmimsto_)
	, beta(beta_)
	, imbeta(imbeta_)
	, wtix(wtix_)
	, iota(iota_)
	, rhobook(rhobook_)
	, rtd(rtd_)
	, alfabw(alfabw_)
	, betabw(betabw_)
	, redchbw(redchbw_)
	, elinv(elinv_)
	, utref(utref_)
	, maxdp(maxdp_)
{}

} // namespace moses
