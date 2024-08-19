#include "moses.h"



void read_double_array(FILE * in, double list[], const int n) {

	for (int i = 0; i < n; i++) {
		fscanf(in, "%lf", &list[i]);
	}
}

void read_aman(FILE * in, double list[][3], const int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < 3; j++) {
			fscanf(in, "%lf", &list[i][j]);
		}
		char newline; fscanf(in, "%c", &newline);
	}
}


int Moses::initialise(std::vector<Firm> & firms, std::vector<Household> & households, External externals[NEX], Market markets[NMKT]) {

	FILE * BANK       = fopen(PATH_BANK,       "r");
	FILE * EXTERNAL   = fopen(PATH_EXTERNAL,   "r");
	FILE * FIRM       = fopen(PATH_FIRM,       "r");
	FILE * GOVERNMENT = fopen(PATH_GOVERNMENT, "r");
	FILE * HOUSEHOLD  = fopen(PATH_HOUSEHOLD,  "r");
	FILE * LABOUR     = fopen(PATH_LABOUR,     "r");
	FILE * MARKET     = fopen(PATH_MARKET,     "r");
	FILE * MISC       = fopen(PATH_MISC,       "r");
	FILE * MOSES      = fopen(PATH_MOSES,      "r");

	if (BANK       == NULL) return 1;
	if (EXTERNAL   == NULL) return 2;
	if (FIRM       == NULL) return 3;
	if (GOVERNMENT == NULL) return 4;
	if (HOUSEHOLD  == NULL) return 5;
	if (LABOUR     == NULL) return 6;
	if (MARKET     == NULL) return 7;
	if (MISC       == NULL) return 8;
	if (MOSES      == NULL) return 9;

	char newline;

	// ----- Order -----
	// misc
	// bank
	// external
	// firm
	// government
	// hosuehold
	// labour
	// market
	// moses
	




	// Misc
	int n_firms, n_households;
	void * misc_static[] = {
		 & n_firms,
		 & n_households
	};

	constexpr size_t misc_static_size[] = {
		sizeof(n_firms),
		sizeof(n_households)
	};

	constexpr int misc_static_type[] = {
		0,
		0
	};

	read_file(MISC, sizeof(misc_static) / sizeof(void *), misc_static, misc_static_size, misc_static_type);





	// Bank
	void * bank_static[] = {
		& Bank::fass,
		& Bank::tmfass,
		& Bank::fd,
		& Bank::tmfd,
		& Bank::kappa1,
		& Bank::kappa2,
		& Bank::lamda1,
		& Bank::lamda2,	
		& Bank::ri,
		& Bank::qchri,
		& Bank::maxqchri,
		& Bank::maxri,
		& Bank::maxridiff,
		& Bank::minri,
		& Bank::mb,
		& Bank::rfund1,
		& Bank::rfund2,
		& Bank::liqb,
		& Bank::liqbfor,
		& Bank::nwb
	};

	constexpr size_t bank_static_size[] = {
		sizeof(Bank::fass),
		sizeof(Bank::tmfass),
		sizeof(Bank::fd),
		sizeof(Bank::tmfd),
		sizeof(Bank::kappa1),
		sizeof(Bank::kappa2),
		sizeof(Bank::lamda1),
		sizeof(Bank::lamda2),	
		sizeof(Bank::ri),
		sizeof(Bank::qchri),
		sizeof(Bank::maxqchri),
		sizeof(Bank::maxri),
		sizeof(Bank::maxridiff),
		sizeof(Bank::minri),
		sizeof(Bank::mb),
		sizeof(Bank::rfund1),
		sizeof(Bank::rfund2),
		sizeof(Bank::liqb),
		sizeof(Bank::liqbfor),
		sizeof(Bank::nwb)
	};

	constexpr int bank_static_type[] = {
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1
	};

	read_file(BANK, sizeof(bank_static) / sizeof(void *), bank_static, bank_static_size, bank_static_type);





	// External
	for (int external = 0; external < NEX; external++) {
		External ex;

		void * var_external[] = {
			& ex.qpdom,
			ex.io, // MKT + NEX array
			ex.qdp // MAXT array
		};

		constexpr size_t var_external_size[] = {
			sizeof(ex.qpdom),
			sizeof(ex.io), // MKT + NEX array
			sizeof(ex.qdp) // MAXT array
		};

		constexpr int var_external_type[] = {
			1,
			1,
			1
		};

		read_file(EXTERNAL, sizeof(var_external) / sizeof(void *), var_external, var_external_size, var_external_type);
		externals[external] = ex;
	}





	// Firms
	void * firm_static[] = {
		& Firm::e1,
		& Firm::e2,
		& Firm::eps,
		& Firm::r,
		& Firm::smp,
		& Firm::sms,
		& Firm::smt,
		& Firm::smw,
		& Firm::fip,
		& Firm::fis,
		& Firm::fiw,
		& Firm::expxdp,
		& Firm::expxds,
		& Firm::expxdw,
		& Firm::loss,
		& Firm::rho,
		& Firm::resmax,
		& Firm::tmsto,
		& Firm::beta,
		& Firm::wtix,
		& Firm::resdown
	};

	constexpr size_t firm_static_size[] = {
		sizeof(Firm::e1),
		sizeof(Firm::e2),
		sizeof(Firm::eps),
		sizeof(Firm::r),
		sizeof(Firm::smp),
		sizeof(Firm::sms),
		sizeof(Firm::smt),
		sizeof(Firm::smw),
		sizeof(Firm::fip),
		sizeof(Firm::fis),
		sizeof(Firm::fiw),
		sizeof(Firm::expxdp),
		sizeof(Firm::expxds),
		sizeof(Firm::expxdw),
		sizeof(Firm::loss),
		sizeof(Firm::rho),
		sizeof(Firm::resmax),
		sizeof(Firm::tmsto),
		sizeof(Firm::beta),
		sizeof(Firm::wtix),
		sizeof(Firm::resdown)
	};

	constexpr int firm_static_type[] = {
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1
	};

	read_file(FIRM, sizeof(firm_static) / sizeof(void *), firm_static, firm_static_size, firm_static_type);


//	double sto[225];
//	double s[225];
//	double p[225];
//	double small[225];
//	double big[225];
//	double share[225];
//	FILE * SHARE = fopen("data/r_share", "r");
//	read_double_array(SHARE, share, 225);
//
//	FILE * STO = fopen("data/r_sto", "r");
//	FILE * S = fopen("data/r_s", "r");
//	FILE * P = fopen("data/r_p", "r");
//	FILE * SMALL = fopen("data/r_small", "r");
//	FILE * BIG = fopen("data/r_big", "r");
//
//	read_double_array(STO, sto, 225);
//	read_double_array(S, s, 225);
//	read_double_array(P, p, 225);
//	read_double_array(SMALL, small, 225);
//	read_double_array(BIG, big, 225);
//
//	printf("\n\n\n\n\n");

	for (int firm = 0; firm < n_firms; firm++) {
		Firm f;

		void * var_firm[] = {
			& f.id,
			& f.mkt,
			& f.m,
			& f.mhist,
			& f.expdp,
			& f.expds,
			& f.expdw,
			& f.dp,
			& f.ds,
			& f.dw,
			& f.histdp,
			& f.histdpdev,
			& f.histdpdev2,
			& f.histds,
			& f.histdsdev,
			& f.histdsdev2,
			& f.histdw,
			& f.histdwdev,
			& f.histdwdev2,
			& f.qp,
			& f.qs,
			& f.qw,
			& f.qq,
			& f.l,
			f.aman, // 3 array
			& f.qtop,
			& f.tec,
			& f.qinv,
			& f.inveff,
			& f.res,
			& f.sto,
			& f.s,
			& f.p,
			& f.small,
			& f.big,
			& f.share
		};

		constexpr size_t var_firm_size[] = {
			sizeof(f.id),
			sizeof(f.mkt),
			sizeof(f.m),
			sizeof(f.mhist),
			sizeof(f.expdp),
			sizeof(f.expds),
			sizeof(f.expdw),
			sizeof(f.dp),
			sizeof(f.ds),
			sizeof(f.dw),
			sizeof(f.histdp),
			sizeof(f.histdpdev),
			sizeof(f.histdpdev2),
			sizeof(f.histds),
			sizeof(f.histdsdev),
			sizeof(f.histdsdev2),
			sizeof(f.histdw),
			sizeof(f.histdwdev),
			sizeof(f.histdwdev2),
			sizeof(f.qp),
			sizeof(f.qs),
			sizeof(f.qw),
			sizeof(f.qq),
			sizeof(f.l),
			sizeof(f.aman),
			sizeof(f.qtop),
			sizeof(f.tec),
			sizeof(f.qinv),
			sizeof(f.inveff),
			sizeof(f.res),
			sizeof(f.sto),
			sizeof(f.s),
			sizeof(f.p),
			sizeof(f.small),
			sizeof(f.big),
			sizeof(f.share)
		};

		constexpr int var_firm_type[] = {
			0,
			0,
			1,
			1,
			1,
			1,
			1,
			1,
			1,
			1,
			1,
			1,
			1,
			1,
			1,
			1,
			1,
			1,
			1,
			1,
			1,
			1,
			1,
			1,
			1,
			1,
			1,
			1,
			1,
			1,
			1,
			1,
			1,
			1,
			1,
			1
		};



		fscanf(FIRM, "%c", &newline);
		read_file(FIRM, sizeof(var_firm) / sizeof(void *), var_firm, var_firm_size, var_firm_type);
		firms.push_back(f);

//		printf("%lf %lf %lf %lf %lf\n", sto[firm], s[firm], p[firm], small[firm], big[firm]);
	}
	




	// Government
	void * government_static[] = {
		& Government::txva2,
		& Government::l,
		Government::qchtxva2 // MAXT array
	};

	constexpr size_t government_static_size[] = {
		sizeof(Government::txva2),
		sizeof(Government::l),
		sizeof(Government::qchtxva2) // MAXT array
	};

	constexpr int government_static_type[] = {
		1,
		1,
		1
	};

	read_file(GOVERNMENT, sizeof(government_static) / sizeof(void *), government_static, government_static_size, government_static_type);



	// Labour
	void * labour_static[] = {
		& Labour::lu,
		& Labour::entry,
		& Labour::ret
	};

	constexpr size_t labour_static_size[] = {
		sizeof(Labour::lu),
		sizeof(Labour::entry),
		sizeof(Labour::ret)

	};

	constexpr int labour_static_type[] = {
		1,
		1,
		1
	};

	read_file(LABOUR, sizeof(labour_static) / sizeof(void *), labour_static, labour_static_size, labour_static_type);





	// Market
	for (int market = 0; market < NMKT; market++) {
		Market m;

		void * var_market[] = {
			& m.qpdom,
			& m.mtec,
			& m.qdmtec,
			m.io // NMKT + NEX array
		};
	
		constexpr size_t var_market_size[] = {
			sizeof(m.qpdom),
			sizeof(m.mtec),
			sizeof(m.qdmtec),
			sizeof(m.io) // NMKT + NEX array
		};

		constexpr int var_market_type[] = {
			1,
			1,
			1,
			1
		};

		read_file(MARKET, sizeof(var_market) / sizeof(void *), var_market, var_market_size, var_market_type);
		markets[market] = m;
	}





	// Moses
	void * moses_static[] = {
		& Moses::niter,
		& Moses::marketiter
	};

	constexpr size_t moses_static_size[] = {
		sizeof(Moses::niter),
		sizeof(Moses::marketiter)
	};

	constexpr int moses_static_type[] = {
		0,
		0
	};

	read_file(MOSES, sizeof(moses_static) / sizeof(void *), moses_static, moses_static_size, moses_static_type);





	initialise_extra(firms, households, externals, markets);


	return 0;
}

inline void Firm::initialise_extra() {

	expids = histds + e1 * histdsdev - e2 * sqrt(histdsdev2);
	expidp = histdp + e1 * histdpdev - e2 * sqrt(histdpdev2);

	curs = 4.0 * (s * pow(1.0 + expids / 4.0, j + 1.5) + cums * pow(1.0 + histds / 4.0, (j - 1.0) / 2.0)) / (4.0 + j);
	curp = (4.0 * p * pow(1.0 + expidp / 4.0, j + 1.5) + cump * pow(1.0 + histdp / 4.0, (j - 1.0) / 2.0)) / (4.0 + j);

	minsto = small * curs / curp;
	maxsto = big * curs / curp;
	optsto = minsto + beta * (maxsto - minsto);
}

void Moses::initialise_extra(std::vector<Firm> & firms, std::vector<Household> & households, External externals[NEX], Market markets[NMKT]) {

	for (auto & firm : firms) {
		firm.initialise_extra();
	}
}


































inline void nullify_firm(std::vector<Firm> & firms, const int firm_index) {
	firms.erase(firms.begin() + firm_index);
}


inline double Firm::qfr(const double l) {
	return wtix * (1.0 - res) * qtop * (1.0 - pow(e, -(tec / qtop) * l));
}

inline double Firm::rfq(const double q) {
	return (qtop / tec) * log((wtix * (1.0 - res) * qtop) / (wtix * (1.0 - res) * qtop - q));
}


inline bool Firm::sat(const double q, const double l) {
	double margin;
	if (l > 0.0) {
		margin = 1.0 - (l * (qexpw / 4.0)) / (q * qexppnet);
	} else {
		margin = 1.0 - (qexpw / 4.0) / (wtix * (1.0 - res) * tec * qexppnet);
	}

	if (margin >= qtargm) return true;
	else return false;
}

inline bool Firm::solve() {

    double b = qexpw / ((1.0 - qtargm) * wtix * (1.0 - res) * tec * qexppnet * 4.0);
    if (b <= 0.0) {
		fprintf(stderr, "[ERROR] b < 0 when entering solve. Exit program.\n");
		fflush(stderr);
		return false;
	} else if (b < 1.0) {
	//	fprintf(stderr, "[ERROR] 0 < b < 1 when entering solve.\n");
	//	fflush(stderr);
	}

	double y = 1.0 / b;
    double h = (b * y + pow(e, -y) - 1.0) / (b - pow(e, -y));

	while (h > 0.001) {
        h = (b * y + pow(e, -y) - 1.0) / (b - pow(e, -y));
        y -= h;
    }

    qplanl = y * (qtop / tec);
    qplanq = qfr(qplanl);


	return true;
}







inline void Firm::quarterly_exp(const bool first_quarter) {
	qexpdp = expdp / 4.0;
	qexpdw = expdw / 4.0;
	qexpds = expds / 4.0;

	if (!first_quarter) {
		qexpdp += fip * (qdp - qexpdp);
		qexpdw += fiw * (qdw - qexpdw);
		qexpds += fis * (qds - qexpds);
	}

	qexpp = qp * (1.0 + qexpdp);
	qexpw = qw * (1.0 + qexpdw);
	qexps = qs * (1.0 + qexpds);
}

inline void Moses::quarterly_exp(const bool first_quarter, std::vector<Firm> & firms, Market markets[NMKT], External externals[NEX]) {
	double sum_qq_qexpdp[NMKT] = {}, sum_qq[NMKT] = {};

	for (auto & firm : firms) {
		firm.quarterly_exp(first_quarter);
		sum_qq_qexpdp[firm.mkt] += firm.qq * firm.qexpdp;
		sum_qq[firm.mkt] += firm.qq;
	}

	for (int i = 0; i < NMKT; i++) {
		markets[i].qexpdpim = sum_qq_qexpdp[i] / sum_qq[i];
		markets[i].qexppim = (1.0 - Government::txva2) * markets[i].qpdom * (1.0 + markets[i].qexpdpim);
	}

	for (int i = 0; i < NEX; i++) {
		externals[i].qexpdpim = externals[i].qdp[t] - Government::qchtxva2[t];
		externals[i].qexppim = (1.0 - Government::txva2) * externals[i].qpdom * (1.0 + externals[i].qexpdpim);
	}
}

inline void Firm::quarterly_targ() {
	qtargm = targm;
}

inline void Moses::luupdate(std::vector<Firm> & firms) {
	double sum_l = 0;
	for (auto & firm : firms) {
		sum_l += firm.l;
		firm.l *= (1.0 - Labour::ret);
		for (int i = 0; i < 3; i++) {
			firm.aman[i] *= (1.0 - Labour::ret);
		}
	}

	Labour::lf = Labour::lu + Government::l + sum_l;
	Labour::lu *= (1.0 - Labour::ret);
	Labour::lu += Labour::entry * Labour::lf;
}

inline void Firm::prodfront(Market markets[NMKT]) {
	qtop *= (1.0 - rho);
	qchqtop1 = (1.0 - loss) * (qinv * inveff) / qp;
	qchqtop2 = std::min(loss * ((qinv * inveff) / qp) * ((resmax - res) / resmax),
						((resmax - res) / (1.0 - resmax)) * (qtop + qchqtop1));
	qchqtop = qchqtop1 + qchqtop2;

	res = (res * (qtop + qchqtop1) + qchqtop2) / (qtop + qchqtop);
	tec = (qtop + qchqtop) / ((qtop / tec) + (qchqtop / markets[mkt].mtec));

	qtop += qchqtop;
}


inline void Moses::prodfront(std::vector<Firm> & firms, Market markets[NMKT]) {
	for (int market = 0; market < NMKT; market++) {
		markets[market].mtec *= (1.0 + markets[market].qdmtec);
	}

	for (auto & firm : firms) {
		firm.prodfront(markets);
	}
}

inline void Firm::initprodplan() {
	qexpsu = qexps / qexpp;
	qplanq = std::max(0.0, qexpsu + (optsto - sto) / (4.0 * tmsto));
}

inline int Firm::target_search(External externals[NEX], Market markets[NMKT]) {

	double sum_io = 0;
	for (int market = 0; market < NMKT; market++) {
		sum_io += markets[mkt].io[market] * markets[market].qexppim;
	}
	for (int external = 0; external < NEX; external++) {
		sum_io += markets[mkt].io[NMKT + external] * externals[external].qexppim;
	}
	qexppnet = qexpp - share * sum_io;


	bool path[12] = {};
	double q2 = 0, q7 = 0;

	// 4.3.1
	if (qplanq >= qtop * (1.0 - res) * wtix) {
		path[7] = true;
	} else if (qplanq > qfr(l)) {
		path[6] = true;
	} else {
		path[2] = true;
	}

	// 4.3.2
	if (path[2]) {
		if (sat(qplanq, l)) {
			qplanl = l;
			path[11] = true;
		} else {
			path[3] = true;
		}
	}

	// 4.3.3
	if (path[3]) {
		//q2 = std::min(qfr(l), qexpsu + maxsto - sto);
		q2 = std::min(qplanq, qfr(l)); // bypass 4.3.3 with this
		if (sat(q2, l)) {
			qplanq = (l * (qexpw / 4.0)) / ((1.0 - qtargm) * qexppnet);
			qplanl = l;
			path[11] = true;
		} else if (q2 == qfr(l)) {
			path[5] = true;
		} else {
			path[4] = true;
		}
	}

	// 4.3.4
	if (path[4]) {
		if (sat(q2, rfq(q2))) {
			qplanq = q2;
			qplanl = ((1.0 - qtargm) * q2 * qexppnet) / (qexpw / 4.0);
			path[11] = true;
		} else {
			path[5] = true;
		}
	}

	// 4.3.5
	if (path[5]) {
		if (sat(qplanq, rfq(qplanq))) {
			solve();
			path[11] = true;
		} else {
			q7 = qplanq;
			path[8] = true;
		}
	}

	// 4.3.6
	if (path[6]) {
		if (sat(qplanq, rfq(qplanq))) {
			qplanl = rfq(qplanq);
			path[11] = true;
		} else {
			path[7] = true;
		}
	}

	// 4.3.7
	if (path[7]) {
		if (sat(qfr(l), l)) {
			solve();
			path[11] = true;
		} else {
			q7 = qfr(l);
			path[8] = true;
		}
	}

	// 4.3.8
	if (path[8]) {
		if (sat(q7, rfq(((1.0 - res) / (1.0 - resdown * res)) * q7))) {
			qplanq = q7;
			qplanl = ((1.0 - qtargm) * q7 * qexppnet) / (qexpw / 4.0);
			res = 1.0 - ((q7 * (1.0 - res)) / qfr(qplanl));
			path[11] = true;
		} else {
			res *= resdown;
			path[9] = true;
		}
	}

	// 4.3.9
	if (path[9]) {
		if (sat(0.0, 0.0)) {
			solve();
			path[11] = true;
		} else {
			path[10] = true;
		}
	}

	// 4.3.10
	if (path[10]) {
		Labour::lu += l;
		return 0;
	}

	// 4.3.11
	if (path[11]) {
		double layoff = std::max(l - qplanl, 0.0);
		aman[0] = std::min(layoff, aman[1]);
		aman[1] = std::min(layoff - aman[0], aman[2]);
		aman[2] = layoff - aman[0] - aman[1];
	}


	return 1;
}


inline void Moses::prodplan(std::vector<Firm> & firms, External externals[NEX], Market markets[NMKT]) {
	luupdate(firms);
	prodfront(firms, markets);

	for (int firm_index = 0; firm_index < firms.size(); firm_index++) {
		firms[firm_index].initprodplan();
		if (!firms[firm_index].target_search(externals, markets)) {
			if ((int) firms.size() <= 1) {
				fprintf(stderr, "[ERROR] attempting to nullify last firm. Terminate program\n");
				fflush(stderr);
			}
			nullify_firm(firms, firm_index);
			--firm_index;
		}
	}

	// 4.4
	for (auto & firm : firms) {
		firm.qplanqsave = firm.qplanq;
	}
}





inline int Moses::simulate_quarter(std::vector<Firm> & firms, std::vector<Household> & households, External externals[NEX], Market markets[NMKT], Bank & bank, Government & government, Labour & Labour) {

	if (t >= MAXT) {
		printf("[ ERROR ] Exceeded MAXT.\n");
		return 1;
	}

	quarterly_exp(j == 0, firms, markets, externals);
	for (auto & firm : firms) {
		firm.quarterly_targ();
	}

	prodplan(firms, externals, markets);

	labour_market(firms, externals, markets);
	

	t += 1;


	return 0;
}





inline void Firm::yearly_init() {
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

inline void Government::yearly_init() {
	cumwtax   = 0;
	cumitax   = 0;
	cumvatax  = 0;
	cumctax   = 0;
	cumws     = 0;
	cuml      = 0;
	cuminv    = 0;
	cumpurch  = 0;
	cumtrans  = 0;
	cumsubs   = 0;
	cummprint = 0;
	cumint    = 0;

	cumgnpcur = 0;
	cumgnpfix = 0;
	cumexport = 0;
	cumimport = 0;
}

inline void Firm::yearly_exp() {
	histdp     = smp * histdp     + (1.0 - smp) * dp;
	histdpdev  = smp * histdpdev  + (1.0 - smp) * (dp - expdp);
	histdpdev2 = smp * histdpdev2 + (1.0 - smp) * ((dp - expdp) * (dp - expdp));
	expidp = histdp + e1 * histdpdev - e2 * sqrt(histdpdev2);
	expdp = (1.0 - r) * expidp + r * expxdp;

	histdw     = smw * histdw     + (1.0 - smw) * dw;
	histdwdev  = smw * histdwdev  + (1.0 - smw) * (dw - expdw);
	histdwdev2 = smw * histdwdev2 + (1.0 - smw) * ((dw - expdw) * (dw - expdw));
	expidw = histdw + e1 * histdwdev - e2 * sqrt(histdwdev2);
	expdw = (1.0 - r) * expidw + r * expxdw;

	histds     = sms * histds     + (1.0 - sms) * ds;
	histdsdev  = sms * histdsdev  + (1.0 - sms) * (ds - expds);
	histdsdev2 = sms * histdsdev2 + (1.0 - sms) * ((ds - expds) * (ds - expds));
	expids = histds + e1 * histdsdev - e2 * sqrt(histdsdev2);
	expds = (1.0 - r) * expids + r * expxds;
}

inline void Firm::yearly_targ() {
	mhist = smt * mhist + (1.0 - smt) * m;
	targm = mhist * (1.0 + eps);
}





inline int Moses::simulate_year(std::vector<Firm> & firms, std::vector<Household> & households, External externals[NEX], Market markets[NMKT], Bank & bank, Government & government, Labour & labour) {

	Government::yearly_init();
	for (auto & firm : firms) {
		firm.yearly_init();
		firm.yearly_exp();
		firm.yearly_targ();
	}

	for (j = 0; j < 4; j++) {
		if (simulate_quarter(firms, households, externals, markets, bank, government, labour)) {
			printf("[ ERROR ] Failed to simulate quarter.\n");
			return 1;
		}
	}


	return 0;
}




int Moses::simulate(FILE * out, const int simulation_length, std::vector<Firm> & firms, std::vector<Household> & households, External externals[NEX], Market markets[NMKT], Bank & bank, Government & government, Labour & labour) {

	for (int year = 0; year < simulation_length; year++) {
		if (simulate_year(firms, households, externals, markets, bank, government, labour)) {
			printf("[ ERROR ] Failed to simulate year.\n");
			return 1;
		}
	}


	return 0;
}
