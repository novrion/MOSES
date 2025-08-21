#include "government.h"
#include "moses.h"
#include "moses_impl.h"
#include "labour.h"
#include "household.h"
#include "bank.h"
#include <iostream>

namespace moses {

Government::Government(const Moses& moses_ref,
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
        double w, double ws, double qchposfor)
	: moses(moses_ref)
	, _qchtxva1(qchtxva1)
    , _qchtxva2(qchtxva2)
    , _txw(txw)
    , _txwg(txwg)
    , _txc(txc)
    , _txi1(txi1)
    , omegag(omegag)
    , omegabld(omegabld)
    , omegain(omegain)
    , omega(omega)
    , gkoff(gkoff)
    , txva2(txva2)
    , txva1(txva1)
    , qrealchl(qrealchl)
    , l(l)
    , qw(qw)
    , qttax(qttax)
    , wgref(wgref)
    , qinv(qinv)
    , qinvbld(qinvbld)
    , qinvin(qinvin)
    , rsubscash(rsubscash)
    , pos(pos)
    , posfor(posfor)
    , w(w)
    , ws(ws)
    , qchposfor(qchposfor)
{}

void Government::initialise_extra(Bank& bank) {
	dep = std::max(0.0, pos);
	bw = std::max(0.0, -pos);
	depfor = std::max(0.0, posfor);
	bwfor = std::max(0.0, -posfor);
}

void Government::print_initialisation(const bool verbose) {
	std::cout << "========== Government Class ==========" << std::endl;
	std::cout << "txva2:     " << txva2 << std::endl;
	std::cout << "txva1:     " << txva1 << std::endl;
	std::cout << "qrealchl:  " << qrealchl << std::endl;
	std::cout << "l:         " << l << std::endl;
	std::cout << "qw:        " << qw << std::endl;
	std::cout << "qttax:     " << qttax << std::endl;
	std::cout << "wgref:     " << wgref << std::endl;
	std::cout << "qinv:     " << qinv << std::endl;
	std::cout << "qinvbld:   " << qinvbld << std::endl;
	std::cout << "qinvin:    " << qinvin << std::endl;
	std::cout << "rsubscash: " << rsubscash << std::endl;
	std::cout << "pos:       " << pos << std::endl;
	std::cout << "posfor:    " << posfor << std::endl;
	std::cout << "w:         " << w << std::endl;
	std::cout << "ws:        " << ws << std::endl;
	std::cout << "qchposfor: " << qchposfor << std::endl;
	std::cout << std::endl;

	if (verbose) {
		std::cout << "omegag: ";
		for (const auto& val : omegag) std::cout << val << " ";
		std::cout << std::endl;
		std::cout << "omegabld: ";
		for (const auto& val : omegabld) std::cout << val << " ";
		std::cout << std::endl;
		std::cout << "omegain: ";
		for (const auto& val : omegain) std::cout << val << " ";
		std::cout << std::endl;
		std::cout << "omega: ";
		for (const auto& val : omega) std::cout << val << " ";
		std::cout << std::endl;
		std::cout << "qchtxva2: ";
		for (const auto& val : _qchtxva2) std::cout << val << " ";
		std::cout << std::endl;
		std::cout << "txw: ";
		for (const auto& val : _txw) std::cout << val << " ";
		std::cout << std::endl;
		std::cout << "txwg: ";
		for (const auto& val : _txwg) std::cout << val << " ";
		std::cout << std::endl;
		std::cout << "txi1: ";
		for (const auto& val : _txi1) std::cout << val << " ";
		std::cout << std::endl;
		std::cout << "gkoff: ";
		for (const auto& val : gkoff) std::cout << val << " ";
		std::cout << std::endl;
		std::cout << "txc: ";
		for (const auto& val : _txc) std::cout << val << " ";
		std::cout << std::endl;
		std::cout << "qchtxva1: ";
		for (const auto& val : _qchtxva1) std::cout << val << " ";
		std::cout << std::endl;
		std::cout << std::endl;
	}
}

double Government::qchtxva1() const { return _qchtxva1[moses.quarter]; }
double Government::qchtxva2() const { return _qchtxva2[moses.quarter]; }
double Government::txw() const { return _txw[moses.year]; }
double Government::txwg() const { return _txwg[moses.year]; }
double Government::txc() const { return _txc[moses.year]; }
double Government::txi1() const { return _txi1[moses.year]; }

void Government::glabour(Labour& lab) {
	double qchl = std::min(lab.lu, l * lab.ret + qrealchl);
	l += qchl - lab.ret * l;
	lab.lu -= qchl;

	double qdw = lab.qdwind;
	qw *= (1 + qdw);
}

void Government::compute_consumption(const std::vector<double>& pt) {
	for (size_t i = 0; i < NSEC; i++) {
		qpurch[i] = gkoff[i] * l * qw/4.0 * pt[i] * (wgref/(100 * qw));
	}
}

void Government::compute_buying(const std::vector<double>& pt) {
	for (size_t i = 0; i < NSEC; i++) {
		qbuy[i] = qpurch[i]/pt[i];

		if (qbuy[i] < 0) {
			throw std::runtime_error("qbuy < 0 in Government::compute_buying()");
		}
	}
}

void Government::accounting(Bank& bank, Household& hh) {
	const int j = moses.relative_quarter + 1; // 1-index relative_quarter

	qint = (dep * bank.ridepg/4.0) - (bw * bank.ribwg/4.0);
	qintfor = (depfor * bank.ridepgfor/4.0) - (bwfor * bank.ribwgfor/4.0);
	
	qws = l * qw/4.0;
	qsubs = qsubsfor + qsubsdom + qsubscash;
	qsp = qws + qsubs + qtrans + moses.sum(qpurch);
	qttax = qwtax + qitax + moses.sum(qvatax) + qctax;
	qsurplus = qttax + qint + qintfor - qsp - qinv;
	qmprint = std::max(0.0, bank.money * 0.25 * moses.sum_avg(&Firm::s, &Firm::ds));

	posfor += qchposfor; // policy option
	depfor = std::max(0.0, pos);
  	bwfor = std::max(0.0, -posfor);

	qchpos = qsurplus + qmprint - qchposfor;
	pos += qchpos;
	dep = std::max(0.0, pos);
	bw = std::max(0.0, -pos);

	cumwtax += qwtax;
	cumitax += qitax;
	cumvatax += moses.sum(qvatax);
	cumctax += qctax;
	cumws += qws;
	cuml = (l + cuml * (j - 1))/j;
	cumsubs += qsubs;
	cumtrans += qtrans;
	for (size_t i = 0; i < NSEC; i++) {
		cumpurch[i] += qpurch[i];
	}
	cuminv += qinv;
	cumint += qint + qintfor;
	cummprint += qmprint;

	// last quarter
	if (moses.relative_quarter == 3) {
		wtax = cumwtax;
		itax = cumitax;
		vatax = cumvatax;
		ctax = cumctax;

		double dws = cumws/ws - 1;
		ws *= (1 + dws);

		double dw = cumws/(cuml * w) - 1;
		w *= (1 + dw);

		subs = cumsubs;
		trans = cumtrans;

		for (size_t i = 0; i < NSEC; i++) {
			purch[i] = cumpurch[i];
		}

		sp = ws + subs + trans + moses.sum(purch);
		inv = cuminv;
		intg = cumint;
		surplus = wtax + itax + vatax + ctax + intg - sp - inv;
		mprint = cummprint;
	}
}

void Government::yearly_init() {
	cumwtax   = 0.0;
	cumitax   = 0.0;
	cumvatax  = 0.0;
	cumctax   = 0.0;
	cumws     = 0.0;
	cuml      = 0.0;
	cuminv    = 0.0;
	cumpurch.fill(0.0);
	cumtrans  = 0.0;
	cumsubs   = 0.0;
	cummprint = 0.0;
	cumint    = 0.0;
}

} // namespace moses
