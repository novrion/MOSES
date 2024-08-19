// TODO: // <->? means possible class change
// TODO: qinv for one firm is negative ?? (-1517834...)
// TODO: double check target search

#ifndef MOSES_H
#define MOSES_H


#include "data.h"

#include <stdio.h>
#include <algorithm>
#include <vector>
#include <cmath>

#define PATH_BANK       "data/bank"
#define PATH_EXTERNAL   "data/external"
#define PATH_FIRM       "data/firm"
#define PATH_GOVERNMENT "data/government"
#define PATH_HOUSEHOLD  "data/household"
#define PATH_LABOUR     "data/labour"
#define PATH_MARKET     "data/market"
#define PATH_MISC       "data/misc"
#define PATH_MOSES      "data/moses"

#define NMKT 4 // number of explicitly producing sectors
#define NEX  6 // number of external sectors

#define MAXT 120 // maximum number of time steps

constexpr size_t size_int    = 4;
constexpr size_t size_double = 8;

constexpr double e = 2.718281828459;

class Moses;
class Bank;
class Firm;
class Government;
class Household;
class Labour;
class Sector;
class Market;
class External;


class Moses {

	public:

	static inline int niter, marketiter;
	static inline int t = 0;
	static inline int j = 0;


		static int initialise(std::vector<Firm> & firms, std::vector<Household> & households, External externals[NEX], Market markets[NMKT]);
		static void initialise_extra(std::vector<Firm> & firms, std::vector<Household> & households, External externals[NEX], Market markets[NMKT]);

		// Run MOSES model
		static int simulate(FILE * out, const int simulation_length, std::vector<Firm> & firms, std::vector<Household> & households, External externals[NEX], Market markets[NMKT], Bank & bank, Government & government, Labour & labour);

		// MOSES model blocks
		static inline int simulate_year(std::vector<Firm> & firms, std::vector<Household> & households, External externals[NEX], Market markets[NMKT], Bank & bank, Government & government, Labour & Labour);
		static inline int simulate_quarter(std::vector<Firm> & firms, std::vector<Household> & households, External externals[NEX], Market markets[NMKT], Bank & bank, Government & government, Labour & Labour);

		static inline void quarterly_exp(const bool first_quarter, std::vector<Firm> & firms, Market mkt[NMKT], External in[NEX]);
		static inline void prodplan(std::vector<Firm> & firms, External externals[NEX], Market markets[NMKT]);
		static inline void luupdate(std::vector<Firm> & firms);
		static inline void prodfront(std::vector<Firm> & firms, Market markets[NMKT]);
};

class Bank : public Moses {

	public:

		static inline double fass, tmfass, fd, tmfd; //const
		static inline double kappa1, kappa2, lamda1, lamda2; //const
		static inline double ri, qchri, maxqchri, maxri, maxridiff, minri, mb; //const
		static inline double rfund1, rfund2; //const

		static inline double liqb, liqbfor;
		static inline double nwb;
};

class Firm : public Moses {

	public:

		static inline double e1, e2, eps, r, smp, sms, smt, smw, fip, fiw, fis; //const
		static inline double expxdp, expxdw, expxds; //const
		static inline double loss, rho, resdown, resmax;
		static inline double tmsto, beta;
		static inline double wtix;

		int id; //const
		int mkt; //const

		double cumq, cumm, cumsu, cums, cumws, cuml, cuminv, cumva, cumsnet, cump;
		double s, p, m, targm, mhist;
		double expdp, expdw, expds, expidp, expidw, expids;
		double dp, dw, ds;
		double histdp, histdpdev, histdpdev2, histdw, histdwdev, histdwdev2, histds, histdsdev, histdsdev2;

		double l, aman[3];
		double curs, curp;
		double sto, minsto, optsto, maxsto, small, big;
		double share, qexppnet;

		double qp, qdp, qexpp, qexpdp,
		       qw, qdw, qexpw, qexpdw,
			   qs, qds, qexps, qexpds,
			   qq, qplanq, qplanqsave, qplanl,
			   qexpsu, qtargm;

		double qtop, qchqtop, qchqtop1, qchqtop2;
		double tec, res;
		double qinv, inveff;

		// yearly functions
		inline void yearly_init();
		inline void yearly_exp();
		inline void yearly_targ();

		// quarterly functions
		inline void quarterly_exp(const bool first_quarter);
		inline void quarterly_targ();
		inline void prodfront(Market markets[NMKT]);
		inline void initprodplan();
		inline int target_search(External externals[NEX], Market markets[NMKT]);


		// helper functions
		inline void initialise_extra();
		inline double qfr(const double l);
		inline double rfq(const double q);
		inline bool sat(double q, double l);
		inline bool solve();
};

class Government : public Moses {

	public:

		static inline double l;

		static inline double qchtxva2[MAXT]; //const
		static inline double txva2; //const

		static inline double cumwtax, cumitax, cumvatax, cumctax, cumws, cuml, cuminv,
							 cumpurch, cumtrans, cumsubs, cummprint, cumint;

		static inline double cumgnpcur, cumgnpfix, cumexport, cumimport; // <->?

		// yearly functions
		static inline void yearly_init();
};

class Household : public Moses {
};

class Labour : public Moses {

	public:

		static inline double lu, lf;
		static inline double entry, ret;
};

class Sector : public Moses {

	public:

		double io[NMKT + NEX];
		double qexppim, qexpdpim;
		double qpdom; // <->?
};

class Market : public Sector {

	public:

		//add
		double mtec, qdmtec;
};

class External : public Sector {

	public:

		double qdp[MAXT]; //qdpin
};


#endif // MOSES_H
