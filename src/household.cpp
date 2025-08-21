#include "household.h"
#include "moses.h"
#include "moses_impl.h"
#include "government.h"
#include "labour.h"
#include "firm.h"
#include "bank.h"
#include <iostream>
#include <algorithm>

namespace moses {

Household::Household(const Moses& moses_ref,
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
        const std::array<double, NEXPH>& smooth)
	: moses(moses_ref)
	, qinpay(qinpay)
    , qtdiv(qtdiv)
    , nh(nh)
    , qsavhreq(qsavhreq)
    , rtrans(rtrans)
    , wh(wh)
    , whra(whra)
    , stodur(stodur)
    , rhodur(rhodur)
    , qcpi(qcpi)
    , qdcpi(qdcpi)
    , alfa3(alfa3)
    , alfa4(alfa4)
    , beta1(beta1)
    , beta2(beta2)
    , beta3(beta3)
    , qp(qp)
    , cva(cva)
    , qc(qc)
    , smooth(smooth)
{}

void Household::print_initialisation(const bool verbose) {
	std::cout << "========== Household Class ==========" << std::endl;
	std::cout << "qinpay:   " << qinpay << std::endl;
	std::cout << "qtdiv:    " << qtdiv << std::endl;
	std::cout << "nh:       " << nh << std::endl;
	std::cout << "qsavhreq: " << qsavhreq << std::endl;
	std::cout << "rtrans:   " << rtrans << std::endl;
	std::cout << "wh:       " << wh << std::endl;
	std::cout << "whra:     " << whra << std::endl;
	std::cout << "stodur:   " << stodur << std::endl;
	std::cout << "rhodur:   " << rhodur << std::endl;
	std::cout << "qcpi:     " << qcpi << std::endl;
	std::cout << "qdcpi:    " << qdcpi << std::endl;
	std::cout << "alfa3:    " << alfa3 << std::endl;
	std::cout << "alfa4:    " << alfa4 << std::endl;
	std::cout << std::endl;

	if (verbose) {
		std::cout << "beta1: ";
		for (const auto& val : beta1) std::cout << val << " ";
		std::cout << std::endl;
		std::cout << "beta2: ";
		for (const auto& val : beta2) std::cout << val << " ";
		std::cout << std::endl;
		std::cout << "beta3: ";
		for (const auto& val : beta3) std::cout << val << " ";
		std::cout << std::endl;
		std::cout << "qp: ";
		for (const auto& val : qp) std::cout << val << " ";
		std::cout << std::endl;
		std::cout << "cva: ";
		for (const auto& val : cva) std::cout << val << " ";
		std::cout << std::endl;
		std::cout << "qc: ";
		for (const auto& val : qc) std::cout << val << " ";
		std::cout << std::endl;
		std::cout << "smooth: ";
		for (const auto& val : smooth) std::cout << val << " ";
		std::cout << std::endl;
		std::cout << std::endl;
	}
}

void Household::init(Bank& bank, Government& g, Labour& lab) {
	g.qtrans = (rtrans * g.qttax) + lab.lu * lab.rlu * moses.wstx() / (4 * moses.sum(&Firm::l));
	qint = bank.rih * wh/4.0 * nh;

	double qtws = g.l * g.qw/4.0 + qinpay + moses.qwsf();
	g.qwtax = (g.l * g.qw/4.0) * (g.txwg() / (1 + g.txwg())) + (qinpay + moses.qwsf()) * g.txw() / (1 + g.txw());
	
	qti = (qtws - g.qwtax) + qint + g.qtrans + qtdiv;

	g.qitax = qti * g.txi1();
	
	qdi = (qti - g.qitax) / nh;
	qspsavreq = qsavhreq / nh;
}

std::vector<double> Household::compute_essential_spending(Bank& bank, Labour& lab, const std::vector<double>& pt, const double qchdcpi) {
	std::vector<double> qspe(NEXPH);

	for (size_t i = 0; i < NSEC; i++) {
		if (i != DUR) {
			qspe[i] = cva[i] * pt[i];
		}
	}

	double swap = alfa3 * (bank.qchri/4.0 - qchdcpi) + alfa4 * lab.qchru;
	swap = std::clamp(swap, -0.05, 0.05);

	qspe[DUR] = std::max(0.0, (pt[DUR] * cva[DUR])/rhodur - (pt[DUR]/qp[DUR] * (1 - rhodur) * stodur) - qdi * swap);
	qspe[SAV] = (whra * qdi - wh) + qdi * swap;

	return qspe;
}

std::vector<double> Household::compute_spending(const std::vector<double>& qspe, const double qprelcpi) {
	double sum = moses.sum_mult(beta1, qspe);
	std::vector<double> qsp(NEXPH);
	for (size_t i = 0; i < NEXPH; i++) {
		qsp[i] = beta1[i] * qspe[i] + (beta2[i] + ((beta3[i] * qprelcpi)/(qdi - qspsavreq))) * (qdi - qspsavreq - sum);

		if (i != SAV) {
			qsp[i] = std::max(0.0, qsp[i]);
		}
	}
	qsp[SAV] += qspsavreq;

	return qsp;
}

std::vector<double> Household::compute_expenditures(const std::vector<double>& pt) {
	double qprelcpi = moses.sum(qc) / moses.sum_div(qc, pt);
	double qchdcpi = qprelcpi/qcpi - 1 - qdcpi;

	std::vector<double> qspe = compute_essential_spending(*moses.bank, *moses.lab, pt, qchdcpi);
	qsp = compute_spending(qspe, qprelcpi);
	return qsp;
}

void Household::compute_buying(const std::vector<double>& qsp, const std::vector<double>& pt) {
	for (size_t i = 0; i < NSEC; i++) {
		qbuy[i] = nh * qsp[i]/pt[i];

		if (qbuy[i] < 0) {
			throw std::runtime_error("qbuy < 0 in Household::compute_buying()");
		}
	}
}

void Household::household_update(std::vector<double>& qsp, const std::vector<double>& pt) {
	for (size_t i = 0; i < NSEC; i++) {
		if (i != DUR) {
			qc[i] = qsp[i];
		}
	}

	stodur = pt[DUR]/qp[DUR] * stodur + qsp[DUR];
	qc[DUR] = rhodur * stodur;
	stodur *= (1 - rhodur);

	qsp[SAV] = qdi - (moses.sum(qsp) - qsp[SAV]);
	wh += qsp[SAV];
	qsav = qsp[SAV] * nh;

	for (size_t i = 0; i < NSEC; i++) {
		cva[i] = smooth[i] * cva[i] + (1 - smooth[i]) * qc[i]/pt[i];
	}
	whra = smooth[SAV] * whra + (1 - smooth[SAV]) * wh/qdi;

	for (size_t i = 0; i < NSEC; i++) {
		qp[i] = pt[i];
	}
	double oldqcpi = qcpi;
	qcpi = moses.sum(qc)/moses.sum_div(qc, qp);
	qdcpi = (qcpi - oldqcpi)/oldqcpi;
}

} // namespace moses
