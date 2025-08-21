#include "moses.h"
#include "params.h"
#include "bank.h"
#include "government.h"
#include "household.h"
#include "labour.h"
#include "market.h"
#include "external.h"
#include "firm.h"
#include <sstream>
#include <fstream>
#include <stdexcept>

namespace moses::csv {

std::vector<std::string> split_line(const std::string& line) {
			std::vector<std::string> tokens;
			std::istringstream ss(line);
			std::string token;

			while (std::getline(ss, token, ',')) {
				tokens.push_back(token);
			}
			return tokens;
		}

std::vector<std::string> read_file(const std::string& path) {
	std::vector<std::string> lines;
	std::ifstream file(path);

	if (!file.is_open()) {
		throw std::runtime_error("Failed to open CSV file: " + path);
	}

	std::string line;
	std::getline(file, line); // remove csv header
	while (std::getline(file, line)) {
		if (!line.empty()) {
			lines.push_back(line);
		}
	}
	return lines;
}

Params parse_params(const std::string& line) {
	auto tokens = split_line(line);
	size_t t = 0;

	int niter      = std::stoi(tokens.at(t++));
	int marketiter = std::stoi(tokens.at(t++));

	double e1      = std::stod(tokens.at(t++));
	double e2      = std::stod(tokens.at(t++));
	double eps     = std::stod(tokens.at(t++));
	double r       = std::stod(tokens.at(t++));
	double smp     = std::stod(tokens.at(t++));
	double sms     = std::stod(tokens.at(t++));
	double smt     = std::stod(tokens.at(t++));
	double smw     = std::stod(tokens.at(t++));
	double fip     = std::stod(tokens.at(t++));
	double fis     = std::stod(tokens.at(t++));
	double fiw     = std::stod(tokens.at(t++));
	double expxdp  = std::stod(tokens.at(t++));
	double expxds  = std::stod(tokens.at(t++));
	double expxdw  = std::stod(tokens.at(t++));
	double loss    = std::stod(tokens.at(t++));
	double rho     = std::stod(tokens.at(t++));
	double resdown = std::stod(tokens.at(t++));
	double resmax  = std::stod(tokens.at(t++));
	double tmsto   = std::stod(tokens.at(t++));
	double tmimsto = std::stod(tokens.at(t++));
	double beta    = std::stod(tokens.at(t++));
	double imbeta  = std::stod(tokens.at(t++));
	double wtix    = std::stod(tokens.at(t++));
	double iota    = std::stod(tokens.at(t++));
	double rhobook = std::stod(tokens.at(t++));
	double rtd     = std::stod(tokens.at(t++));
	double alfabw  = std::stod(tokens.at(t++));
	double betabw  = std::stod(tokens.at(t++));
	double redchbw = std::stod(tokens.at(t++));
	double elinv   = std::stod(tokens.at(t++));
	double utref   = std::stod(tokens.at(t++));
	
	double maxdp = std::stod(tokens.at(t++));

	return Params(niter, marketiter,
			e1, e2, eps, r, smp, sms, smt, smw, fip, fis, fiw,
			expxdp, expxds, expxdw, loss, rho, resdown, resmax,
			tmsto, tmimsto, beta, imbeta, wtix, iota, rhobook, rtd,
			alfabw, betabw, redchbw, elinv, utref,
			maxdp);
}

std::unique_ptr<Bank> parse_bank(const Moses& moses_ref, const std::string& line) {
	auto tokens = split_line(line);
	size_t t = 0;

	double fass      = std::stod(tokens.at(t++));
	double tmfass    = std::stod(tokens.at(t++));
	double fd        = std::stod(tokens.at(t++));
	double tmfd      = std::stod(tokens.at(t++));
	double kappa1    = std::stod(tokens.at(t++));
	double kappa2    = std::stod(tokens.at(t++));
	double lamda1    = std::stod(tokens.at(t++));
	double lamda2    = std::stod(tokens.at(t++));
	double ri        = std::stod(tokens.at(t++));
	double qchri     = std::stod(tokens.at(t++));
	double maxri     = std::stod(tokens.at(t++));
	double maxqchri  = std::stod(tokens.at(t++));
	double maxridiff = std::stod(tokens.at(t++));
	double minri     = std::stod(tokens.at(t++));
	double mb        = std::stod(tokens.at(t++));
	double rfund1    = std::stod(tokens.at(t++));
	double rfund2    = std::stod(tokens.at(t++));
	double liqb      = std::stod(tokens.at(t++));
	double liqbfor   = std::stod(tokens.at(t++));
	double nw        = std::stod(tokens.at(t++));

	std::array<double, MAXT> ridepfor;
    std::array<double, MAXT> ribwfor;
	for (auto& val : ridepfor) val = std::stod(tokens.at(t++));
    for (auto& val : ribwfor)  val = std::stod(tokens.at(t++));

	return std::make_unique<Bank>(moses_ref, ridepfor, ribwfor, mb, rfund1, rfund2, 
           liqb, liqbfor, nw, fass, tmfass, fd, tmfd,
           ri, qchri, maxri, maxqchri, minri, maxridiff,
           kappa1, kappa2, lamda1, lamda2);
}

std::unique_ptr<Government> parse_government(const Moses& moses_ref, const std::string& line) {
	auto tokens = split_line(line);
    size_t t = 0;

    double txva2     = std::stod(tokens.at(t++));
    double txva1     = std::stod(tokens.at(t++));
    double qrealchl  = std::stod(tokens.at(t++));
    double l         = std::stod(tokens.at(t++));
    double qw        = std::stod(tokens.at(t++));
    double qttax     = std::stod(tokens.at(t++));
    double wgref     = std::stod(tokens.at(t++));
    double qinv     = std::stod(tokens.at(t++));
    double qinvbld   = std::stod(tokens.at(t++));
    double qinvin    = std::stod(tokens.at(t++));
    double rsubscash = std::stod(tokens.at(t++));
    double pos       = std::stod(tokens.at(t++));
    double posfor    = std::stod(tokens.at(t++));
    double w         = std::stod(tokens.at(t++));
    double ws        = std::stod(tokens.at(t++));
    double qchposfor = std::stod(tokens.at(t++));

    std::array<double, NSEC> omegag;
    std::array<double, NSEC> omegabld;
    std::array<double, NSEC> omegain;
    std::array<double, NSEC> omega;
    std::array<double, NSEC> gkoff;
    std::array<double, MAXT> qchtxva2;
    std::array<double, MAXY> txw;
    std::array<double, MAXY> txwg;
    std::array<double, MAXY> txi1;
    std::array<double, MAXY> txc;
    std::array<double, MAXT> qchtxva1;

    for (auto& val : omegag)   val = std::stod(tokens.at(t++));
    for (auto& val : omegabld) val = std::stod(tokens.at(t++));
    for (auto& val : omegain)  val = std::stod(tokens.at(t++));
    for (auto& val : omega)    val = std::stod(tokens.at(t++));
    for (auto& val : qchtxva2) val = std::stod(tokens.at(t++));
    for (auto& val : txw)      val = std::stod(tokens.at(t++));
    for (auto& val : txwg)     val = std::stod(tokens.at(t++));
    for (auto& val : txi1)     val = std::stod(tokens.at(t++));
    for (auto& val : gkoff)    val = std::stod(tokens.at(t++));
    for (auto& val : txc)      val = std::stod(tokens.at(t++));
    for (auto& val : qchtxva1) val = std::stod(tokens.at(t++));

	return std::make_unique<Government>(moses_ref, qchtxva1, qchtxva2, txw, txwg, txc, txi1,
                 omegag, omegabld, omegain, omega, gkoff,
                 txva2, txva1, qrealchl, l, qw, qttax, wgref,
                 qinv, qinvbld, qinvin, rsubscash, pos, posfor,
                 w, ws, qchposfor);
}

std::unique_ptr<Household> parse_household(const Moses& moses_ref, const std::string& line) {
	auto tokens = split_line(line);
    size_t t = 0;

    double qinpay   = std::stod(tokens.at(t++));
    double qtdiv    = std::stod(tokens.at(t++));
    double nh       = std::stod(tokens.at(t++));
    double qsavhreq = std::stod(tokens.at(t++));
    double rtrans   = std::stod(tokens.at(t++));
    double wh       = std::stod(tokens.at(t++));
    double whra     = std::stod(tokens.at(t++));
    double stodur   = std::stod(tokens.at(t++));
    double rhodur   = std::stod(tokens.at(t++));
    double qcpi     = std::stod(tokens.at(t++));
    double qdcpi    = std::stod(tokens.at(t++));
    double alfa3    = std::stod(tokens.at(t++));
    double alfa4    = std::stod(tokens.at(t++));

    std::array<double, NEXPH> beta1;
    std::array<double, NEXPH> beta2;
    std::array<double, NEXPH> beta3;
    std::array<double, NSEC> qp;
    std::array<double, NSEC> cva;
    std::array<double, NSEC> qc;
    std::array<double, NEXPH> smooth;

    for (auto& val : beta1)  val = std::stod(tokens.at(t++));
    for (auto& val : beta2)  val = std::stod(tokens.at(t++));
    for (auto& val : beta3)  val = std::stod(tokens.at(t++));
    for (auto& val : qp)    val = std::stod(tokens.at(t++));
    for (auto& val : cva)    val = std::stod(tokens.at(t++));
    for (auto& val : qc)     val = std::stod(tokens.at(t++));
    for (auto& val : smooth) val = std::stod(tokens.at(t++));

    return std::make_unique<Household>(moses_ref, qinpay, qtdiv, nh, qsavhreq, rtrans,
                wh, whra, stodur, rhodur, qcpi, qdcpi, 
                alfa3, alfa4, beta1, beta2, beta3,
                qp, cva, qc, smooth);
}

std::unique_ptr<Labour> parse_labour(const Moses& moses_ref, const std::string& line) {
	auto tokens = split_line(line);
    size_t t = 0;

    double lu     = std::stod(tokens.at(t++));
    double entry  = std::stod(tokens.at(t++));
    double ret    = std::stod(tokens.at(t++));
    double qdwind = std::stod(tokens.at(t++));
    double theta  = std::stod(tokens.at(t++));
    double ksi    = std::stod(tokens.at(t++));
    double gamma  = std::stod(tokens.at(t++));
    double ru     = std::stod(tokens.at(t++));
    double rlu    = std::stod(tokens.at(t++));
    double skrepa = std::stod(tokens.at(t++));

    return std::make_unique<Labour>(moses_ref, lu, entry, ret, qdwind, theta, ksi,
             gamma, ru, rlu, skrepa);
}

std::unique_ptr<External> parse_external(Moses& moses_ref, const std::string& line) {
	auto tokens = split_line(line);
	size_t t = 0;
	
	int    id    = std::stoi(tokens.at(t++));
	double qpdom = std::stod(tokens.at(t++));
	double qpfor = std::stod(tokens.at(t++));
	double imp   = std::stod(tokens.at(t++));
	double tmimp = std::stod(tokens.at(t++));
	double x     = std::stod(tokens.at(t++));
	double pref  = std::stod(tokens.at(t++));

	std::array<double, NSEC> io;
	for (auto& val : io) val = std::stod(tokens.at(t++));

	std::array<double, MAXT> qdp;
	for (auto& val : qdp) val = std::stod(tokens.at(t++));

	return std::make_unique<External>(moses_ref, id, qpdom, qpfor, imp, tmimp, x, pref, io, qdp);
}

std::unique_ptr<Market> parse_market(Moses& moses_ref, const Params& params, const std::string& line) {
	auto tokens = split_line(line);
	size_t t = 0;

	int    id       = std::stoi(tokens.at(t++));
	double qpdom    = std::stod(tokens.at(t++));
	double qpfor    = std::stod(tokens.at(t++));
	double mtec     = std::stod(tokens.at(t++));
	double qdmtec   = std::stod(tokens.at(t++));
	double tmx      = std::stod(tokens.at(t++));
	double rsubs    = std::stod(tokens.at(t++));
	double imp      = std::stod(tokens.at(t++));
	double tmimp    = std::stod(tokens.at(t++));
	double tminv    = std::stod(tokens.at(t++));
	double pref     = std::stod(tokens.at(t++));
	double tstocurf = std::stod(tokens.at(t++));
	double tstocurm = std::stod(tokens.at(t++));

	std::array<double, NSEC> io;
	for (auto& val : io) val = std::stod(tokens.at(t++));
	
	std::array<double, MAXT> qdpfor;
	for (auto& val : qdpfor) val = std::stod(tokens.at(t++));

	return std::make_unique<Market>(moses_ref, params, id, qpdom, qpfor, mtec, qdmtec, tmx, rsubs, imp, tmimp, tminv, pref, tstocurf, tstocurm, io, qdpfor);
}

std::unique_ptr<Firm> parse_firm(const Moses& moses_ref, const Params& params, const std::string& line) {
	auto tokens = split_line(line);
	size_t t = 0;

	int    id         = std::stoi(tokens.at(t++));
	int    market_id  = std::stoi(tokens.at(t++));
	double m          = std::stod(tokens.at(t++));
	double mhist      = std::stod(tokens.at(t++));
	double expdp      = std::stod(tokens.at(t++));
	double expds      = std::stod(tokens.at(t++));
	double expdw      = std::stod(tokens.at(t++));
	double dp         = std::stod(tokens.at(t++));
	double ds         = std::stod(tokens.at(t++));
	double dw         = std::stod(tokens.at(t++));
	double histdp     = std::stod(tokens.at(t++));
	double histdpdev  = std::stod(tokens.at(t++));
	double histdpdev2 = std::stod(tokens.at(t++));
	double histds     = std::stod(tokens.at(t++));
	double histdsdev  = std::stod(tokens.at(t++));
	double histdsdev2 = std::stod(tokens.at(t++));
	double histdw     = std::stod(tokens.at(t++));
	double histdwdev  = std::stod(tokens.at(t++));
	double histdwdev2 = std::stod(tokens.at(t++));
	double qp         = std::stod(tokens.at(t++));
	double qs         = std::stod(tokens.at(t++));
	double qw         = std::stod(tokens.at(t++));
	double qq         = std::stod(tokens.at(t++));
	double l          = std::stod(tokens.at(t++));
	double qtop       = std::stod(tokens.at(t++));
	double tec        = std::stod(tokens.at(t++));
	double qinv       = std::stod(tokens.at(t++));
	double inveff     = std::stod(tokens.at(t++));
	double res        = std::stod(tokens.at(t++));
	double sto        = std::stod(tokens.at(t++));
	double s          = std::stod(tokens.at(t++));
	double p          = std::stod(tokens.at(t++));
	double small      = std::stod(tokens.at(t++));
	double big        = std::stod(tokens.at(t++));
	double share      = std::stod(tokens.at(t++));
	double chm        = std::stod(tokens.at(t++));
	double imsmall    = std::stod(tokens.at(t++));
	double imbig      = std::stod(tokens.at(t++));
	double x          = std::stod(tokens.at(t++));
	double qinvlag    = std::stod(tokens.at(t++));
	double qva        = std::stod(tokens.at(t++));
	double q          = std::stod(tokens.at(t++));
	double va         = std::stod(tokens.at(t++));
	double w          = std::stod(tokens.at(t++));
	double bw         = std::stod(tokens.at(t++));
	double k1         = std::stod(tokens.at(t++));
	double k2         = std::stod(tokens.at(t++));
	double k1book     = std::stod(tokens.at(t++));
	double rw         = std::stod(tokens.at(t++));

	std::array<double, 3> aman;
	for (auto& val : aman) val = std::stod(tokens.at(t++));

	std::array<double, NSEC> imsto;
	for (auto& val : imsto) val = std::stod(tokens.at(t++));

	return std::make_unique<Firm>(moses_ref, params, id, market_id, m, mhist, expdp, expds, expdw, dp, ds, dw, histdp, histdpdev, histdpdev2, histds, histdsdev, histdsdev2, histdw, histdwdev, histdwdev2, qp, qs, qw, qq, l, qtop, tec, qinv, inveff, res, sto, s, p, small, big, share, chm, imsmall, imbig, x, qinvlag, qva, q, va, w, bw, k1, k2, k1book, rw, aman, imsto);
}

} //namespace moses::csv
