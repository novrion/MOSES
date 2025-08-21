#include "bank.h"
#include "moses.h"
#include "moses_impl.h"
#include "government.h"
#include "firm.h"
#include "household.h"
#include <iostream>
#include <algorithm>

namespace moses {

Bank::Bank(const Moses& moses_ref,
		const std::array<double, MAXT>& ridepfor,
		const std::array<double, MAXT>& ribwfor,
        double mb, double rfund1, double rfund2,
        double liqb, double liqbfor, double nw,
        double fass, double tmfass, double fd, double tmfd,
        double ri, double qchri, double maxri, double maxqchri,
        double minri, double maxridiff,
        double kappa1, double kappa2, double lamda1, double lamda2)
	: moses(moses_ref)
	, _ridepfor(ridepfor)
    , _ribwfor(ribwfor)
    , mb(mb)
    , rfund1(rfund1)
    , rfund2(rfund2)
    , liqb(liqb)
    , liqbfor(liqbfor)
    , nw(nw)
    , fass(fass)
    , tmfass(tmfass)
    , fd(fd)
    , tmfd(tmfd)
    , ri(ri)
    , qchri(qchri)
    , maxri(maxri)
    , maxqchri(maxqchri)
    , minri(minri)
    , maxridiff(maxridiff)
    , kappa1(kappa1)
    , kappa2(kappa2)
    , lamda1(lamda1)
    , lamda2(lamda2)
{}

void Bank::initialise_extra(Government& g, Household& hh) {
	rik2 = ri - mb;
	rif = ri;
	rih = ri - mb;
	ribwg = ri;
	ridepg = ri - mb;
	ribwgfor = ribwfor();
	ridepgfor = ridepfor();

	money = moses.sum(&Firm::k2) + g.dep + hh.wh * hh.nh;
}

void Bank::print_initialisation(const bool verbose) {
	std::cout << "========== Bank Class ==========" << std::endl;
	std::cout << "fass:      " << fass << std::endl;
	std::cout << "tmfass:    " << tmfass << std::endl;
	std::cout << "fd:        " << fd << std::endl;
	std::cout << "tmfd:      " << tmfd << std::endl;
	std::cout << "kappa1:    " << kappa1 << std::endl;
	std::cout << "kappa2:    " << kappa2 << std::endl;
	std::cout << "lamda1:    " << lamda1 << std::endl;
	std::cout << "lamda2:    " << lamda2 << std::endl;
	std::cout << "ri:        " << ri << std::endl;
	std::cout << "qchri:     " << qchri << std::endl;
	std::cout << "maxri:     " << maxri << std::endl;
	std::cout << "maxqchri:  " << maxqchri << std::endl;
	std::cout << "maxridiff: " << maxridiff << std::endl;
	std::cout << "minri:     " << minri << std::endl;
	std::cout << "mb:        " << mb << std::endl;
	std::cout << "rfund1:    " << rfund1 << std::endl;
	std::cout << "rfund2:    " << rfund2 << std::endl;
	std::cout << "liqb:      " << liqb << std::endl;
	std::cout << "liqbfor:   " << liqbfor << std::endl;
	std::cout << "nw:        " << nw << std::endl;
	std::cout << std::endl;

	if (verbose) {
		std::cout << "ridepfor: ";
		for (const auto& val : _ridepfor) std::cout << val << " ";
		std::cout << std::endl;
		std::cout << "ribwfor: ";
		for (const auto& val : _ribwfor)  std::cout << val << " ";
		std::cout << std::endl;
		std::cout << std::endl;
	}
}

double Bank::ridepfor() const { return _ridepfor[moses.quarter]; }
double Bank::ribwfor() const { return _ribwfor[moses.quarter]; }

void Bank::compute_buying(Government& g, const std::vector<double>& pt) {
	for (size_t i = 0; i < NSEC; i++) {
		double qinvtot = 0;
		qinvtot += g.omegag[i] * g.qinv;
		qinvtot += g.omegabld[i] * g.qinvbld;
		qinvtot += g.omegain[i] * g.qinvin;
		qinvtot += g.omega[i] * moses.sum(&Firm::qinvlag);

		qbuy[i] = qinvtot / (pt[i] * (1 - g.txva2)/(1 - g.txva1));

		if (qbuy[i] < 0) {
			throw std::runtime_error("qbuy < 0 in Bank::compute_buying()");
		}
	}
}

void Bank::bank_transactions(Government& g, Household& hh) {
	double rfpay = (ri - ribwfor())/maxridiff;
	rfpay = std::clamp(rfpay, -lamda2, lamda2);

	double qfasspay = (fass + moses.qexport)/(1 + 4 * tmfass * (1 - rfpay));
	double qchfass = moses.qexport - qfasspay;
	fass += qchfass;

	double qfdpay = (fd + moses.qimport)/(1 + 4 * tmfd * (1 + rfpay));
	double qchfd = moses.qimport - qfdpay;
	fd += qchfd;
	
	double qchliqbfor = qfasspay + g.qintfor - qfdpay - g.qchposfor;
	liqbfor += qchliqbfor;

	qchliqb = moses.qintf + hh.qsav + g.qchpos + moses.qimport + g.qchposfor 
			  - moses.qintk2 - hh.qint - g.qint - moses.qexport - g.qintfor;
}

void Bank::credit_market(Government& g, Household& hh) {
	double qdemfund = std::max(0.0, -g.qchpos) + moses.sum_qdeschbw_positive();

	double qsupfund1 = (liqb + qchliqb + moses.sum(&Firm::qdeschk2) - rfund1 
		* (moses.sum(&Firm::bw) + g.bw + std::max(0.0, -hh.wh * hh.nh)))
		/(1 + rfund1);
	double qsupfund2 = liqb + qchliqb + moses.sum(&Firm::qdeschk2) 
					   - rfund2 * (moses.sum(&Firm::k2) + moses.sum(&Firm::qdeschk2) + g.dep + std::max(0.0, hh.wh * hh.nh));
	qsupfund = std::max(0.0, std::min(qsupfund1, qsupfund2));

	hh.qsavhreq = std::min(kappa1 * std::max(0.0, hh.qdi * hh.nh), std::max(0.0, qdemfund - qsupfund));

	double sum_qdeschbw = moses.sum_qdeschbw_positive();
	double qredtbw = std::min(kappa2 * sum_qdeschbw, std::max(0.0, qdemfund - qsupfund - hh.qsavhreq));
	for (auto& firm : moses.firms) {
		firm->qchbw = firm->qdeschbw - qredtbw * std::max(0.0, firm->qdeschbw)/sum_qdeschbw;
	}

	qchri = lamda1 * (qdemfund - qsupfund)/(std::max(1.0, qdemfund));
	qchri = std::clamp(qchri, -maxqchri, maxqchri);

	ri += qchri;
	ri = std::clamp(ri, minri, maxri);

	rif = ri;
	rih = ri - mb;
	rik2 = ri - mb;
	
	ribwg = ri;
	ridepg = ri - mb;
	ribwgfor = ribwfor();
	ridepgfor = ridepfor();
}

void Bank::bank_update(Government& g, Household& hh) {
	qchliqb += qtchk2 - qtchbw;
	liqb += qchliqb;

	double assets = moses.sum(&Firm::bw) + fass + liqbfor + liqb;
	double liabilities = moses.sum(&Firm::k2) + hh.wh * hh.nh + g.pos + fd;
	nw = assets - liabilities;

	money = moses.sum(&Firm::k2) + g.dep + hh.wh * hh.nh;
}

} // namespace moses
