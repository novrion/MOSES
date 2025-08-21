#include "labour.h"
#include "moses.h"
#include "moses_impl.h"
#include "government.h"
#include "firm.h"
#include <iostream>

namespace moses {

Labour::Labour(const Moses& moses_ref,
		double lu, double entry, double ret,
		double qdwind, double theta, double ksi,
        double gamma, double ru, double rlu,
        double skrepa)
	: moses(moses_ref)
	, lu(lu)
    , entry(entry)
    , ret(ret)
    , qdwind(qdwind)
    , theta(theta)
    , ksi(ksi)
    , gamma(gamma)
    , ru(ru)
    , rlu(rlu)
    , skrepa(skrepa)
{}

void Labour::print_initialisation() {
	std::cout << "========== Labour Class ==========" << std::endl;
	std::cout << "lu:     " << lu << std::endl;
	std::cout << "entry:  " << entry << std::endl;
	std::cout << "ret:    " << ret << std::endl;
	std::cout << "qdwind: " << qdwind << std::endl;
	std::cout << "theta:  " << theta << std::endl;
	std::cout << "ksi:    " << ksi << std::endl;
	std::cout << "gamma:  " << gamma << std::endl;
	std::cout << "ru:     " << ru << std::endl;
	std::cout << "rlu:    " << rlu << std::endl;
	std::cout << "skrepa: " << skrepa << std::endl;
	std::cout << std::endl;	
}

void Labour::luupdate(Government& g) {
	lf = lu + g.l + moses.sum(&Firm::l);
	lu *= (1 - ret);
	lu += entry * lf;
}

} // namespace moses
