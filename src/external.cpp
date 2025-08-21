#include "external.h"
#include "moses.h"
#include "market.h"
#include "government.h"
#include "household.h"
#include <iostream>

namespace moses {

External::External(Moses& moses_ref,
		int id_, double qpdom_, double qpfor_, 
		double imp_, double tmimp_, double x_, double pref_,
        const std::array<double, NSEC>& io_,
        const std::array<double, MAXT>& qdp)
	: Sector(moses_ref, id_, qpdom_, qpfor_, imp_, tmimp_, pref_, io_)
    , x(x_)
    , _qdp(qdp)
{}

void External::print_initialisation(const bool verbose) {
	std::cout << "========== External [" << id << "] ==========" << std::endl;
	std::cout << "qpdom: " << qpdom << std::endl;
	std::cout << "qpfor: " << qpfor << std::endl;
	std::cout << "imp  : " << imp   << std::endl;
	std::cout << "tmimp: " << tmimp << std::endl;
	std::cout << "x    : " << x     << std::endl;
	std::cout << std::endl;

	if (verbose) {
		std::cout << "io: ";
		for (const auto& val : io) std::cout << val << " ";
		std::cout << std::endl;
		std::cout << "qdp: ";
		for (const auto& val : _qdp) std::cout << val << " ";
		std::cout << std::endl;
		std::cout << std::endl;
	}
}

double External::qdp() const { return _qdp[moses.quarter]; }

double External::qpfor_() const {return qpdom * (1 - moses.gov->txva2); }

void External::quarterly_exp(Government& g) {
	double qexpdpim = qdp() - g.qchtxva2();
	qexppim = (1 - g.txva2) * qpdom * (1 + qexpdpim);
}

double External::compute_trial_price() {
	return qpdom * (1 + qdp());
}

void External::compute_buying() {
	for (size_t i = 0; i < NSEC; i++) {
		qbuy[i] = io[i] * qq;

		if (qbuy[i] < 0) {
			throw std::runtime_error("qbuy < 0 in External::compute_buying()");
		}
	}
}

void External::compute_imports(Government& g, const std::vector<double>& pt) {
	qtbuyfor = imp * qtbuy;
	qtbuydom = qtbuy - qtbuyfor;

	if (qtbuyfor < 0) {
		throw std::runtime_error("qtbuyfor < 0 in compute_imports()");
	}

	moses.qimport += qtbuyfor * ((1 - g.txva2) * qpdom * (1 + qdp()));
}

void External::domestic_result(const std::vector<double>& pt) {
	qdpdom = pt[id]/qpdom - 1;
	qpdom *= (1 + qdpdom);
}

double External::qimpurch(Government& g, const std::vector<double>& pt) const {
	double ret = 0;
	for (const auto& mkt : moses.markets) {
		ret += qbuy[mkt->id] * mkt->qpdom * (1 - g.txva2);
	}
	for (const auto& ext : moses.externals) {
		ret += qbuy[ext->id] * ext->qpdom * (1 - g.txva2);
	}
	return ret;
}

void External::external_sectors(Government& g, Household& hh, const std::vector<double>& pt) {
	double qexportin = qq * x * qpfor_();
	moses.qexport += qexportin;

	double qsdom = qtbuydom * qpfor_();

	qva = qsdom + qexportin - qimpurch(g, pt);
	hh.qinpay += qva;
}

void External::national_accounting(std::array<double, NGNPCUR>& qgnpcur, std::array<double, NGNPFIX>& qgnpfix) {
	// Production in external sectors
	qgnpcur[id] = qva;
	qgnpfix[id] = pref * qq;
}

} // namespace moses
