#include "sector.h"
#include "bank.h"
#include "government.h"
#include "household.h"

namespace moses {

Sector::Sector(Moses& moses_ref,
		int id_, double qpdom_, double qpfor_, 
		double imp_, double tmimp_,
		double pref_,
        const std::array<double, NSEC>& io_)
	: moses(moses_ref)
	, id(id_)
    , qpdom(qpdom_)
    , qpfor(qpfor_)
    , imp(imp_)
    , tmimp(tmimp_)
    , pref(pref_)
    , io(io_)
    , qexppim(0)
    , qdpdom(0)
    , qbuy()  // Initialize array to zeros
    , qtbuy(0)
    , qtbuyfor(0)
    , qtbuydom(0)
{}

void Sector::indirect_taxes(Bank& bank, Government& g, Household& hh, const std::vector<double>& qsp) {
	g.qvatax[id] = g.txva2 * (g.qpurch[id] + hh.nh * qsp[id]);
	if (g.txva1 < 0) {
		double more = (qpdom * (1 - g.txva2))/(1 - g.txva1) * bank.qbuy[id];
		g.qvatax[id] += g.txva1 * more;
	}

	g.qvataximp += qtbuyfor/qtbuy * g.qvatax[id];
}

} // namespace moses
