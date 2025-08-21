#include "firm.h"
#include "moses.h"
#include "params.h"
#include "market.h"
#include "labour.h"
#include "government.h"
#include "bank.h"
#include "external.h"
#include <math.h>
#include <algorithm>
#include <stdexcept>
#include <iostream>

namespace moses {

Firm::Firm(const Moses& moses_ref, const Params& params,
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
    	const std::array<double, NSEC>& imsto_)
	: moses(moses_ref)

	, e1(params.e1)
	, e2(params.e2)
	, eps(params.eps)
	, r(params.r)
	, smp(params.smp)
	, sms(params.sms)
	, smt(params.smt)
	, smw(params.smw)
	, fip(params.fip)
	, fis(params.fis)
	, fiw(params.fiw)
	, expxdp(params.expxdp)
	, expxds(params.expxds)
	, expxdw(params.expxdw)
	, loss(params.loss)
	, rho(params.rho)
	, resdown(params.resdown)
	, resmax(params.resmax)
	, tmsto(params.tmsto)
	, tmimsto(params.tmimsto)
	, beta(params.beta)
	, imbeta(params.imbeta)
	, wtix(params.wtix)
	, iota(params.iota)
	, rhobook(params.rhobook)
	, rtd(params.rtd)
	, alfabw(params.alfabw)
	, betabw(params.betabw)
	, redchbw(params.redchbw)
	, elinv(params.elinv)
	, utref(params.utref)

	, id(id_)
	, market_id(market_id_)
	, m(m_)
	, mhist(mhist_)
	, expdp(expdp_)
	, expds(expds_)
	, expdw(expdw_)
	, dp(dp_)
	, ds(ds_)
	, dw(dw_)
	, histdp(histdp_)
	, histdpdev(histdpdev_)
	, histdpdev2(histdpdev2_)
	, histds(histds_)
	, histdsdev(histdsdev_)
	, histdsdev2(histdsdev2_)
	, histdw(histdw_)
	, histdwdev(histdwdev_)
	, histdwdev2(histdwdev2_)
	, qp(qp_)
	, qs(qs_)
	, qw(qw_)
	, qq(qq_)
	, l(l_)
	, qtop(qtop_)
	, tec(tec_)
	, qinv(qinv_)
	, inveff(inveff_)
	, res(res_)
	, sto(sto_)
	, s(s_)
	, p(p_)
	, small(small_)
	, big(big_)
	, share(share_)
	, chm(chm_)
	, imsmall(imsmall_)
	, imbig(imbig_)
	, x(x_)
	, qinvlag(qinvlag_)
	, qva(qva_)
	, q(q_)
	, va(va_)
	, w(w_)
	, bw(bw_)
	, k1(k1_)
	, k2(k2_)
	, k1book(k1book_)
	, rw(rw_)
	, aman(aman_)
	, imsto(imsto_)
{}

void Firm::print_initialisation(const bool print_params, const bool verbose) {
	std::cout << "========== Firm [" << id << "] ==========" << std::endl;
	std::cout << "market_id : " << market_id  << std::endl;
	std::cout << "m         : " << m          << std::endl;
	std::cout << "mhist     : " << mhist      << std::endl;
	std::cout << "expdp     : " << expdp      << std::endl;
	std::cout << "expds     : " << expds      << std::endl;
	std::cout << "expdw     : " << expdw      << std::endl;
	std::cout << "dp        : " << dp         << std::endl;
	std::cout << "ds        : " << ds         << std::endl;
	std::cout << "dw        : " << dw         << std::endl;
	std::cout << "histdp    : " << histdp     << std::endl;
	std::cout << "histdpdev : " << histdpdev  << std::endl;
	std::cout << "histdpdev2: " << histdpdev2 << std::endl;
	std::cout << "histds    : " << histds     << std::endl;
	std::cout << "histdsdev : " << histdsdev  << std::endl;
	std::cout << "histdsdev2: " << histdsdev2 << std::endl;
	std::cout << "histdw    : " << histdw     << std::endl;
	std::cout << "histdwdev : " << histdwdev  << std::endl;
	std::cout << "histdwdev2: " << histdwdev2 << std::endl;
	std::cout << "qp        : " << qp         << std::endl;
	std::cout << "qs        : " << qs         << std::endl;
	std::cout << "qw        : " << qw         << std::endl;
	std::cout << "qq        : " << qq         << std::endl;
	std::cout << "l         : " << l          << std::endl;
	std::cout << "qtop      : " << qtop       << std::endl;
	std::cout << "tec       : " << tec        << std::endl;
	std::cout << "qinv      : " << qinv       << std::endl;
	std::cout << "inveff    : " << inveff     << std::endl;
	std::cout << "res       : " << res        << std::endl;
	std::cout << "sto       : " << sto        << std::endl;
	std::cout << "s         : " << s          << std::endl;
	std::cout << "p         : " << p          << std::endl;
	std::cout << "small     : " << small      << std::endl;
	std::cout << "big       : " << big        << std::endl;
	std::cout << "share     : " << share      << std::endl;
	std::cout << "chm       : " << chm        << std::endl;
	std::cout << "imsmall   : " << imsmall    << std::endl;
	std::cout << "imbig     : " << imbig      << std::endl;
	std::cout << "x         : " << x          << std::endl;
	std::cout << "qinvlag   : " << qinvlag    << std::endl;
	std::cout << "qva       : " << qva        << std::endl;
	std::cout << "q         : " << q          << std::endl;
	std::cout << "va        : " << va         << std::endl;
	std::cout << "w         : " << w          << std::endl;
	std::cout << "bw        : " << bw         << std::endl;
	std::cout << "k1        : " << k1         << std::endl;
	std::cout << "k2        : " << k2         << std::endl;
	std::cout << "k1book    : " << k1book     << std::endl;
	std::cout << "rw        : " << rw         << std::endl;
	std::cout << std::endl;

	if (verbose) {
		std::cout << "aman: ";
		for (const auto& val : aman) std::cout << val << " ";
		std::cout << std::endl;
		std::cout << "imsto: ";
		for (const auto& val : imsto) std::cout << val << " ";
		std::cout << std::endl;
		std::cout << std::endl;
	}

	if (print_params) {
		std::cout << "========== Params ==========" << std::endl;
		std::cout << "e1	 : " << e1	    << std::endl;
		std::cout << "e2     : " << e2      << std::endl;
		std::cout << "eps    : " << eps     << std::endl;
		std::cout << "r      : " << r       << std::endl;
		std::cout << "smp    : " << smp     << std::endl;
		std::cout << "sms    : " << sms     << std::endl;
		std::cout << "smt    : " << smt     << std::endl;
		std::cout << "smw    : " << smw     << std::endl;
		std::cout << "fip    : " << fip     << std::endl;
		std::cout << "fis    : " << fis     << std::endl;
		std::cout << "fiw    : " << fiw     << std::endl;
		std::cout << "expxdp : " << expxdp  << std::endl;
		std::cout << "expxds : " << expxds  << std::endl;
		std::cout << "expxdw : " << expxdw  << std::endl;
		std::cout << "loss   : " << loss    << std::endl;
		std::cout << "rho    : " << rho     << std::endl;
		std::cout << "resdown: " << resdown << std::endl;
		std::cout << "resmax : " << resmax  << std::endl;
		std::cout << "tmsto  : " << tmsto   << std::endl;
		std::cout << "tmimsto: " << tmimsto << std::endl;
		std::cout << "beta   : " << beta    << std::endl;
		std::cout << "imbeta : " << imbeta  << std::endl;
		std::cout << "wtix   : " << wtix    << std::endl;
		std::cout << "iota   : " << iota    << std::endl;
		std::cout << "rhobook: " << rhobook << std::endl;
		std::cout << "rtd    : " << rtd     << std::endl;
		std::cout << "alfabw : " << alfabw  << std::endl;
		std::cout << "betabw : " << betabw  << std::endl;
		std::cout << "redchbw: " << redchbw << std::endl;
		std::cout << "elinv  : " << elinv   << std::endl;
		std::cout << "utref  : " << utref   << std::endl;
		std::cout << std::endl;
	}
}

double Firm::curs() const {
	const int j = moses.relative_quarter + 1; // 1-index
	return (1 / (4.0 + j))
		* (4 * s * std::pow(1 + histds/4.0, j + 1.5)
		+ j * cums * std::pow(1 + histds/4.0, (j - 1)/2.0));
}

double Firm::qcurs() const {
	const int j = moses.relative_quarter + 1; // 1-index
	return (1 / (4.0 + j))
		* (s * std::pow(1 + histds/4.0, j + 1.5)
		+ cums * std::pow(1 + histds/4.0, (j - 1)/2.0));
}

double Firm::curp() const {
	const int j = moses.relative_quarter + 1; // 1-index
	return (1 / (4.0 + j))
		* (4 * p * std::pow(1 + histdp/4.0, j + 1.5)
		+ j * cump * std::pow(1 + histdp/4.0, (j - 1)/2.0));
}

double Firm::qexppnet() const {
    double sum = 0.0;
	const auto& firm_market = moses.markets[market_id];
	for (const auto& mkt : moses.markets) {
		sum += firm_market->io[mkt->id] * mkt->qexppim;
	}
	for (const auto& ext : moses.externals) {
		sum += firm_market->io[ext->id] * ext->qexppim;
	}

    return qexpp - share * sum;
}

double Firm::sum_io_qpdom_txva2() const {
	double sum = 0.0;
	double txva2 = moses.gov->txva2;
	const auto& firm_market = moses.markets[market_id];
	for (const auto& mkt : moses.markets) {
		sum += firm_market->io[mkt->id] * mkt->qpdom * (1 - txva2);
	}
	for (const auto& ext : moses.externals) {
		sum += firm_market->io[ext->id] * ext->qpdom * (1 - txva2);
	}
	return sum;
}

double Firm::sum_qimq_qpdom_txva2() const {
	double sum = 0.0;
	double txva2 = moses.gov->txva2;
	for (const auto& mkt : moses.markets) {
		sum += qimq[mkt->id] * mkt->qpdom * (1 - txva2);
	}
	for (const auto& ext : moses.externals) {
		sum += qimq[ext->id] * ext->qpdom * (1 - txva2);
	}
	return sum;
}

double Firm::k3imed() const {
	double sum = 0.0;
	double txva2 = moses.gov->txva2;
	for (const auto& mkt : moses.markets) {
		sum += imsto[mkt->id] * mkt->qpdom * (1 - txva2);
	}
	for (const auto& ext : moses.externals) {
		sum += imsto[ext->id] * ext->qpdom * (1 - txva2);
	}
	return sum;
}

void Firm::initialise_extra(const Market& market) {
	bad = 0;

	yearly_init();
	reference_inventory_levels(market);
}

double Firm::qfr(const double l) {
	return wtix * (1 - res) * qtop * (1 - std::exp(-(tec / qtop) * l));
}

double Firm::rfq(const double q) {
	return (qtop / tec) * log((wtix * (1 - res) * qtop) / (wtix * (1 - res) * qtop - q));
}

bool Firm::sat(const double q, const double l, const double qexppnet) {
	double margin;
	if (l > 0) {
		margin = 1.001 - (l * (qexpw/4.0)) / (q * qexppnet);
	} else {
		margin = 1 - (qexpw/4.0) / (wtix * (1 - res) * tec * qexppnet);
	}

	if (margin >= qtargm) return true;
	else return false;
}

void Firm::solve(const double qexppnet) {
    double b = qexpw / ((1 - qtargm) * wtix * (1 - res) * tec * qexppnet * 4);
    if (b <= 0) {
		throw std::runtime_error("Entering solve() b <= 0");
	}

	const int MAX_ITERATIONS = 100;
	const double TOLERANCE = 0.001;
	double y = 1 / b;

	for (int i = 0; i < MAX_ITERATIONS; i++) {
		double f = b * y + std::exp(-y) - 1.0;
		double f_prime = b - std::exp(-y);

		if (std::abs(f_prime) < 1e-10) {
			throw std::runtime_error("Zero derivative in solve()");
		}
		
		double delta = f / f_prime;
        y -= delta;

		if (std::abs(delta) < TOLERANCE) {
			qplanl = y * (qtop / tec);
		    qplanq = qfr(qplanl);
			return;
		}
    }

	throw std::runtime_error("Failed to converge in solve()");
}

double Firm::check_bw(const int bw_ratio) {
	const double val = 1.0 - bw_ratio;
	if (val > 0.3) {
		return 1.0;
	} else if (val <= 0.1) {
		return 0.0;
	} else {
		return 0.0025 * std::pow(std::floor(100 * val + 0.5) - 10, 2.0);
	}
}

double Firm::delay() {
	const int tminv = moses.markets[market_id]->tminv;
	double delta = 1.0 / (1.0 + (4.0 * tminv) / 3.0);

	// Update each stage of the delay process
	inv_delay[2] += qinvlag;			   // Add new investment to stage 3
    inv_delay[1] += delta * inv_delay[2];  // Flow from stage 3 to 2
    inv_delay[2] -= delta * inv_delay[2];  // Reduce stage 3
    inv_delay[0] += delta * inv_delay[1];  // Flow from stage 2 to 1
    inv_delay[1] -= delta * inv_delay[1];  // Reduce stage 2
    
    // Calculate output investment
    double qinv = delta * inv_delay[0];
    
    // Update final stage
    inv_delay[0] *= (1.0 - delta);
    
    return qinv;
}

void Firm::yearly_init() {
	cumq    = 0;
	cumm    = 0;
	cumsu   = 0;
	cums    = 0;
	cump    = 0;
	cumws   = 0;
	cuml    = 0;
	cuminv  = 0;
	cumva   = 0;
	cumsnet = 0;
}

void Firm::yearly_exp() {
    histdp = smp * histdp + (1 - smp) * dp;
    histdpdev = smp * histdpdev + (1 - smp) * (dp - expdp);
    histdpdev2 = smp * histdpdev2 + (1 - smp) * pow(dp - expdp, 2);
    double expidp = histdp + e1 * histdpdev - e2 * sqrt(histdpdev2);
    expdp = (1 - r) * expidp + r * expxdp;

    histdw = smw * histdw + (1 - smw) * dw;
    histdwdev = smw * histdwdev + (1 - smw) * (dw - expdw);
    histdwdev2 = smw * histdwdev2 + (1 - smw) * pow(dw - expdw, 2);
    double expidw = histdw + e1 * histdwdev + e2 * sqrt(histdwdev2);
    expdw = (1 - r) * expidw + r * expxdw;

    histds = sms * histds + (1 - sms) * ds;
    histdsdev = sms * histdsdev + (1 - sms) * (ds - expds);
    histdsdev2 = sms * histdsdev2 + (1 - sms) * pow(ds - expds, 2);
    double expids = histds + e1 * histdsdev - e2 * sqrt(histdsdev2);
    expds = (1 - r) * expids + r * expxds;
}

void Firm::yearly_targ() {
	mhist = std::max(0.0, smt * mhist + (1 - smt) * m);
	targm = mhist * (1 + eps);
}

void Firm::yearly_update() {
	dq = cumq/q - 1;
	q *= (1 + dq);
	dp = cums/(cumsu * p) - 1;
	p *= (1 + dp);
	dw = cumws/(cuml * w) - 1;
	w *= (1 + dw);
	ds = cums/s - 1;
	s *= (1 + ds);
	dva = cumva/va - 1;
	va *= (1 + dva);
	snet = cumsnet;
	chm = cumm - m;
	m += chm;
}

void Firm::quarterly_exp(const bool first_quarter) {
	qexpdp = expdp/4.0;
	double qexpdw = expdw/4.0;
	double qexpds = expds/4.0;

	if (!first_quarter) {
		qexpdp += fip * (qdp - qexpdp);
		qexpdw += fiw * (qdw - qexpdw);
		qexpds += fis * (qds - qexpds);
	}

	qexpp = qp * (1 + qexpdp);
	qexpw = qw * (1 + qexpdw);
	qexps = qs * (1 + qexpds);
}

void Firm::quarterly_targ() {
	qtargm = targm;
}

void Firm::luupdate(Labour& lab) {
	l *= (1 - lab.ret);
	for (auto& val : aman) {
		val *= (1 - lab.ret);
	}
}

void Firm::prodfront(const Market& mkt) {
	qtop *= (1 - rho);

	double qchqtop1 = (1 - loss) * (qinv * inveff) / qp;
	double qchqtop2 = std::min(loss * ((qinv * inveff) / qp) * ((resmax - res) / resmax),
						((resmax - res) / (1 - resmax)) * (qtop + qchqtop1));
	double qchqtop = qchqtop1 + qchqtop2;

	res = (res * (qtop + qchqtop1) + qchqtop2) / (qtop + qchqtop);
	if (res < 0 || res > resmax) {
		throw std::runtime_error("res exceeded resmax or res < 0 in prodfront()");
	}

	tec = (qtop + qchqtop) / ((qtop / tec) + (qchqtop / mkt.mtec));
	qtop += qchqtop;
}

void Firm::init_prodplan() {
	qexpsu = qexps/qexpp;
	qplanq = std::max(0.0, qexpsu + (optsto - sto)/(4 * tmsto));
}

bool Firm::try_reduce_production(const double qexppnet) {

	// path 4
	// try reduce production to qplanq with corresponding decrease in labour
	if (sat(qplanq, rfq(qplanq), qexppnet)) {
		solve(qexppnet);
		exit_successfully(qexppnet);
		return true;
	}

	return try_reduce_slack(qplanq, qexppnet);
}

bool Firm::try_increase_production(const double qexppnet) {

	// path 2
	double q2 = std::min(qfr(l), qexpsu + maxsto - sto);
	if (sat(q2, l, qexppnet)) {
		qplanq = (l * qexpw/4.0) / ((1 - qtargm) * qexppnet);
		qplanl = l;
		exit_successfully(qexppnet);
		return true;
	}

	if (q2 == qfr(l)) {
		return try_reduce_production(qexppnet);
	}

	// path 3
	// try reduce employment, still producing q2
	if (sat(q2, rfq(q2), qexppnet)) {
		qplanq = q2;
		qplanl = ((1 - qtargm) * q2 * qexppnet) / (qexpw/4.0);
		exit_successfully(qexppnet);
		return true;
	}

	return try_reduce_production(qexppnet);
}

bool Firm::plan_implies_recruitment(const double qexppnet) {

	// path 6
	if (sat(qfr(l), l, qexppnet)) {
		solve(qexppnet);
		exit_successfully(qexppnet);
		return true;
	}

	return try_reduce_slack(qfr(l), qexppnet);
}

bool Firm::try_reduce_slack(const double q7, const double qexppnet) {

	// path 7
	// keep production at q7 and reduce slack (res)
	if (sat(q7, rfq((1 - res) / (1 - resdown * res) * q7), qexppnet)) {
		qplanq = q7;
		qplanl = ((1 - qtargm) * q7 * qexppnet) / (qexpw/4.0);
		res = 1 - (q7 * (1 - res)) / qfr(qplanl);
		exit_successfully(qexppnet);
		return true;
	}

	// reduce slack
	res *= resdown;

	// path 8
	if (sat(0, 0, qexppnet)) {
		solve(qexppnet);
		exit_successfully(qexppnet);
		return true;
	}

	return false;
}

void Firm::handle_layoffs() {
	double layoff = std::max(0.0, l - qplanl);
	aman[0] = std::min(layoff, aman[1]);
	aman[1] = std::min(layoff - aman[0], aman[2]);
	aman[2] = layoff - aman[0] - aman[1];
}

void Firm::exit_successfully(const double qexppnet) {
	handle_layoffs();

	if (qplanq < 0) {
		throw std::runtime_error("qplanq < 0 after target_search()");
	}
	if (qplanl < 0) {
		throw std::runtime_error("qplanl < 0 after target_search()");
	}
	if (!sat(qplanq, qplanl, qexppnet)) {
		throw std::runtime_error("exiting target_search(), qplanq & qplanl doesn't satisfy margin");
	}
}

bool Firm::target_search() {
	double qexppnet_ = qexppnet();
	if (qexppnet_ < 1e-6) {
		throw std::runtime_error("qexppnet <= 0 in target_search()");
	}

	// plan requires more than maximum production => recruitment
	if (qplanq >= qtop * (1 - res) * wtix) {
		return plan_implies_recruitment(qexppnet_);
	}

	// plan implies recruitment
	if (qplanq > qfr(l)) {

		// path 5
		if (sat(qplanq, rfq(qplanq), qexppnet_)) {
			qplanl = rfq(qplanq);
			exit_successfully(qexppnet_);
			return true;
		}

		return plan_implies_recruitment(qexppnet_);
	}

	// path 1
	// initial plan satisfies profit target?
	if (sat(qplanq, l, qexppnet_)) {
		qplanl = l;
		exit_successfully(qexppnet_);
		return true;
	}

	return try_increase_production(qexppnet_);
}

void Firm::labour_search_input(LabourSearchData& data) {
	if (chm > 0) {
		data.chl.push_back(qplanl - l);
	} else {
		data.chl.push_back(rfq(qplanq) - l);
	}

	data.ww.push_back(qw + iota * (qexpw - qw));
	data.ll.push_back(l);
}

void Firm::handle_labour_search_layoffs() {
	double exit = std::max(0.0, -qchl);

	if (exit > aman[0] + aman[1]) {
		aman[2] = std::max(0.0, aman[2] - (exit - aman[0] - aman[1]));
	}

	if (exit > aman[0]) {
		aman[1] = std::max(0.0, aman[1] - (exit - aman[0]));
	}

	if (exit > 0) {
		aman[0] = std::max(0.0, aman[0] - exit);
	}
}

void Firm::handle_labour_update_layoffs(Labour& lab) {
	double sack = std::min(aman[0], std::max(0.0, l + qchl - qplanl));
	qchl -= sack;
	aman[0] -= sack;
	lab.lu += sack;
}

void Firm::update_labour_force_and_wage() {
	l += qchl;
	qdw = qchw/qw;
	qw += qchw;
}

void Firm::planqrevise(const Market& market) {
	qplanq = std::min(qplanq, qfr(l));	
	double qdq = qplanq/qq - 1;
	qq *= (1 + qdq);

	if (qq <= 0) {
		throw std::runtime_error("qq <= 0. Division by 0 or negative production in planqrevise().");
	}	

	for (size_t i = 0; i < NSEC; i++) {
		qimq[i] = std::max(0.0, share * market.io[i] * qplanqsave + (optimsto[i] - imsto[i]) / (4 * tmimsto));

		if (qimq[i] < 0) {
			throw std::runtime_error("qimq[] < 0 in planqrevise()");
		}
	}

	for (size_t i = 0; i < NSEC; i++) {
		imsto[i] = std::min(maximsto[i], imsto[i] + qimq[i] - share * market.io[i] * qq);

		if (imsto[i] < -1e-6) {
			throw std::runtime_error("imsto[] < 0 in planqrevise()");
		}
		imsto[i] = std::max(0.0, imsto[i]); // account for precision of type double
	}

	qoptsu = std::max(0.0, qq * (qexpsu / (qexpsu + (optsto - sto) / (4 * tmsto))));
}

void Firm::adjust_foreign(const Market& market, Government& g) {
	if (market.qpdom * (1 - g.txva2) > market.qpfor) {
		x -= x * (1/(4*market.tmx)) * (market.qpdom * (1 - g.txva2) - market.qpfor)/market.qpfor;
	} else {
		x += (1 - x) * (1/(4*market.tmx)) * (market.qpfor - market.qpdom * (1 - g.txva2)) / (market.qpdom * (1 - g.txva2));
	}
	x = std::clamp(x, 0.0, 1.0);

	qsufor = x * qoptsu;
}

void Firm::reference_inventory_levels(const Market& market) {
	double curs_ = 4 * qcurs();
	double curp_ = curp();

	minsto = small * curs_/curp_;
	maxsto = big * curs_/curp_;
	optsto = minsto + beta * (maxsto - minsto);

	for (size_t i = 0; i < NSEC; i++) {
		minimsto[i] = imsmall * market.io[i] * share * curs_/curp_;
		maximsto[i] = imbig * market.io[i] * share * curs_/curp_;
		optimsto[i] = minimsto[i] + imbeta * (maximsto[i] - minimsto[i]);
	}
}

void Firm::finalqpqsqm() {
	qsu = qsufor + qsudom;
	qds = (qsfor + qsdom - qs)/qs;
	qs = qsfor + qsdom;
	if (qs <= 0) {
		throw std::runtime_error("qs <= 0 in finalqpqsqm()");
	}

	qdp = ((qs/qsu)-qp)/qp;
	qp = qs/qsu;
	qdva = qq * (qp - sum_io_qpdom_txva2()) / qva - 1;
	qva *= (1 + qdva);
	qsnet = qs - sum_qimq_qpdom_txva2();
	qm = 1 - (l * qw/4.0)/qsnet;
}

void Firm::quarterly_cum() {
	const int j = moses.relative_quarter + 1; // 1-index relative quarter

	cuminv += qinvlag;
	cumq += qq;
	cumva += qva;
	cums += qs;
	cumsu += qsu;
	cumsnet += qsnet;
	cumws += l * qw/4.0;
	cuml = ((j - 1) * cuml + l)/j;
	cumm = 1 - cumws/cumsnet;
	cump = cums/cumsu;
}

void Firm::invfin(Bank& bank, Government& g) {
	qrev = qm * qsnet + (k2 * bank.rik2/4.0) - (bw * bank.rif/4.0);

	double qdpk_ = moses.qdpk();
	k1 = qinv + k1 * (1 + qdpk_);
	qdepr = rho * k1;
	k1 -= qdepr;

	k1book += qinv;
	qdeprbook = std::max(0.0, std::min(qrev, rhobook * k1book));
	k1book -= qdeprbook;

	k3finish = sto * qp;
	k3 = k3imed() + k3finish;

	qtax = g.txc() * std::max(0.0, qrev - qdeprbook);
	qdiv = rtd * qtax;
	qcash = qrev - (qtax + qdiv) + qs * g.rsubscash;

	qrr = 4 * (qrev - qdepr)/(k1 + k2 + k3);
	qdeschk2 = (rw * 4 * qcurs()) - k2;

	qdeschbw = std::max(0.0, 1 - elinv * (utref - qq/(qtop * wtix * (1 - res)))) 
			   * bw * (alfabw + betabw * (qdpk_ + (qrr - bank.rif)/4.0));
	qdeschbw = std::min(qdeschbw, redchbw * bw);
	if (bw + qdeschbw <= 0) {
		throw std::runtime_error("bw + qdeschbw <= 0 in Firm::invfin()");
	}

	double a = k1 + k2 + qdeschk2 + k3;
	qdeschbw *= check_bw((bw + qdeschbw)/a);

	qinv = delay();
	inveff = qtop * qp/k1;
}

bool Firm::invfin_adjustments() {
	bw += qchbw;
	if (bw <= 0) {
		throw std::runtime_error("bw <= 0 in invfin_adjustments()");
	}
	
	qinvlag = std::max(0.0, qcash + qchbw - qdeschk2);
	qchk2 = qcash + qchbw - qinvlag;
	k2 += qchk2;

	nw = k1 + k2 + k3 - bw;
	if (nw < 0) {
		bad++;
	}
	if (bad >= 6) {
		return false;
	}

	return true;
}

} // namespace moses
