#include "market.h"
#include "moses.h"
#include "moses_impl.h"
#include "params.h"
#include "firm.h"
#include "government.h"
#include <iostream>
#include <algorithm>

namespace moses {

Market::Market(Moses& moses_ref, const Params& params,
		int id_, double qpdom_, double qpfor_, 
		double mtec_, double qdmtec_, double tmx_, 
        double rsubs_, double imp_, double tmimp_, 
        double tminv_, double pref_,
		double tstocurf_, double tstocurm_,
        const std::array<double, NSEC>& io_,
        const std::array<double, MAXT>& qdpfor)
	: Sector(moses_ref, id_, qpdom_, qpfor_, imp_, tmimp_, pref_, io_)

	, maxdp(params.maxdp)

    , mtec(mtec_)
    , qdmtec(qdmtec_)
    , tminv(tminv_)
    , tmx(tmx_)
    , rsubs(rsubs_)
	, tstocurf(tstocurf_)
	, tstocurm(tstocurm_)
    , _qdpfor(qdpfor)
{}

void Market::print_initialisation(const bool print_params, const bool verbose) {
	std::cout << "========== Market [" << id << "] ==========" << std::endl;
	std::cout << "qpdom : " << qpdom  << std::endl;
	std::cout << "qpfor : " << qpfor  << std::endl;
	std::cout << "mtec  : " << mtec   << std::endl;
	std::cout << "qdmtec: " << qdmtec << std::endl;
	std::cout << "tmx   : " << tmx    << std::endl;
	std::cout << "rsubs : " << rsubs  << std::endl;
	std::cout << "imp   : " << imp    << std::endl;
	std::cout << "tmimp : " << tmimp  << std::endl;
	std::cout << "tminv : " << tminv  << std::endl;
	std::cout << std::endl;

	if (verbose) {
		std::cout << "io: ";
		for (const auto& val : io) std::cout << val << " ";
		std::cout << std::endl;
		std::cout << "qdpfor: ";
		for (const auto& val : _qdpfor) std::cout << val << " ";
		std::cout << std::endl;
		std::cout << std::endl;
	}

	if (print_params) {
		std::cout << "========== Params ==========" << std::endl;
		std::cout << "maxdp: " << maxdp << std::endl;
		std::cout << std::endl;
	}
}

double Market::qdpfor() const { return _qdpfor[moses.quarter]; }

void Market::quarterly_exp(Government& g) {
	double qexpdpim = moses.sum_avg(&Firm::qq, &Firm::qexpdp, id);
	qexppim = (1 - g.txva2) * qpdom * (1 + qexpdpim);
}

double Market::compute_trial_price(Government& g) {
	if (qpdom > qpfor) {
		imp += (1 - imp)/(4 * tmimp) * (qpdom * (1 - g.txva2) - qpfor)/qpfor;
	} else {
		imp -= imp/(4 * tmimp) * (qpfor - qpdom * (1 - g.txva2))/(qpdom * (1 - g.txva2));
	}
	imp = std::clamp(imp, 0.0, 1.0);

	return (1 - imp) * qprelpdom + imp * qpfor / (1 - g.txva2);
}

void Market::compute_buying() {
	for (size_t i = 0; i < NSEC; i++) {
		qbuy[i] = 0;
		for (auto& firm : moses.firms) {
			if (firm->market_id == id) {
				qbuy[i] += firm->qimq[i];
			}
		}

		if (qbuy[i] < 0) {
			throw std::runtime_error("qbuy < 0 in Market::compute_buying()");
		}
	}
}

void Market::price_adjust(std::vector<double>& pt) {
	double ch = ((1 - imp) * maxdp * pt[id])/(4 * (moses.marketiter - 1));
	if (qtbuy * (1 - imp) < moses.sum(&Firm::qoptsudom, id)) {
		pt[id] -= ch;
	} else {
		pt[id] += ch;
	}
}

void Market::compute_imports(const std::vector<double>& pt) {
	double qtbuyfor1 = imp * qtbuy;
	qtbuydom = qtbuy - qtbuyfor1;

	double qmaxtsudom = std::max(0.0, moses.sum(&Firm::qq, id) + moses.sum(&Firm::sto, id) - moses.sum(&Firm::minsto, id) - moses.sum(&Firm::qsufor, id));
	qtbuydom = std::min(qtbuydom, qmaxtsudom);

	double qtbuyfor2 = qtbuy - (qtbuydom + qtbuyfor1);
	qtbuyfor = qtbuyfor1 + qtbuyfor2;

	if (qtbuyfor < 0) {
		throw std::runtime_error("qtbuyfor < 0 in compute_imports()");
	}

	moses.qimport += qtbuyfor * qpfor;
}

void Market::domestic_result(const std::vector<double>& pt) {
	qdpdom = (pt[id] - imp * qpfor)/((1 - imp) * qprelpdom) * (qprelpdom/qpdom - 1);
	qpdom *= (1 + qdpdom);
	qtsudom = qtbuydom;
}

void Market::firm_sto(Government& g, const std::vector<double>& limsto, const std::vector<double>& upper, const std::vector<double>& lower) {
	double totchsto = std::min(moses.sum(&Firm::qq, id) - moses.sum(&Firm::qsufor, id) - qtsudom, moses.sum(upper, id) - moses.sum(&Firm::sto, id));
	qchtsto = totchsto;
	qwaste = moses.sum(&Firm::qq, id) - moses.sum(&Firm::qsufor, id) - qtsudom - totchsto;

	if (moses.sum(lower, id) > totchsto + moses.sum(&Firm::sto, id) + 1e-2) {
		throw std::runtime_error("sum(lower) > totchsto + sum(sto) [mkt] in Market::firm_sto()");
	}
	if (moses.sum(upper, id) < totchsto + moses.sum(&Firm::sto, id) - 1e-2) {
		throw std::runtime_error("sum(upper) < totchsto + sum(sto) [mkt] in Market::firm_sto()");
	}

	std::vector<double> qchsto(moses.firms.size(), 0);

	// adjust firms outside prespecified sto limits
	for (size_t i = 0; i < moses.firms.size(); i++) {
		auto& firm = moses.firms[i];
		if (firm->market_id != id) continue;

		double ch = 0.0;;
		if (firm->sto > upper[i]) {
			ch = upper[i] - firm->sto;
		} else if (firm->sto < lower[i]) {
			ch = lower[i] - firm->sto;
		}

		firm->sto += ch;
		qchsto[i] += ch;
		totchsto -= ch;
	}

	// Make sure these are here! DO NOT REDO THE SUM FOR EACH ITERATION!
	double upper_sto_diff = moses.sum(upper, id) - moses.sum(&Firm::sto, id);
	double lower_sto_diff = moses.sum(lower, id) - moses.sum(&Firm::sto, id);
		
	if (totchsto > 1e-6) {
		if (upper_sto_diff < 1e-6 && upper_sto_diff > -1e-6) {
			throw std::runtime_error("sum(upper - sto) = 0 in firm_sto()");
		}
	} else if (totchsto < -1e-6) {
		if (lower_sto_diff < 1e-6 && lower_sto_diff > -1e-6) {
			throw std::runtime_error("sum(lower - sto) = 0 in firm_sto()");
		}
	}

	for (size_t i = 0; i < moses.firms.size(); i++) {
		auto& firm = moses.firms[i];
		if (firm->market_id != id) continue;
	
		double ch = 0.0;
		if (totchsto > 1e-6) {
			ch = (upper[i] - firm->sto)/upper_sto_diff * totchsto;
		} else if (totchsto < -1e-6) {
			ch = (lower[i] - firm->sto)/lower_sto_diff * totchsto;
		}

		firm->sto += ch;
		qchsto[i] += ch;

		firm->qsudom = firm->qq - firm->qsufor - qchsto[i];
		if (firm->qsudom < -1e-6) {
			throw std::runtime_error("Negative qsudom in firm_sto().");
		}
		firm->qsudom = std::max(0.0, firm->qsudom); // account for precision error of type double

		firm->qsdom = firm->qsudom * (1 + rsubs) * qpdom * (1 - g.txva2);

		g.qsubsdom += firm->qsudom * rsubs * qpdom * (1 - g.txva2);
	}
}

void Market::national_accounting(std::array<double, NGNPCUR>& qgnpcur, std::array<double, NGNPFIX>& qgnpfix) {
	qchtstocurf = moses.sum(&Firm::k3finish, id) - tstocurf;
	tstocurf += qchtstocurf;

	qchtstocurm = moses.sum_mult(&Firm::sto, qpdom, id) - tstocurm;
	tstocurm += qchtstocurm;

	// Production in explicit model sectors
	qgnpcur[id] = moses.sum(&Firm::qsnet, id) + qchtstocurf;
	qgnpfix[id] = pref * (moses.sum(&Firm::qq, id) - qwaste);
}

} // namespace moses
