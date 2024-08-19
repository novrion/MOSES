#pragma once

#include <stdio.h>
#include <algorithm>
#include <cmath>
#include <time.h>
#include <chrono>


// Constants
const double e = 2.718281828459;

constexpr int DUR  = 0;
constexpr int NDUR = 1;
constexpr int SAV  = 2;
constexpr int iZ   = 3;

constexpr int MAXNTOT = 1000;
constexpr int MAXNM = 2;
constexpr int MAXNMH = 4;

// Parameters
const int YEARS = 50;
const int NM = 2;
const int NMH = 4;
const int NITER = 5, MARKETITER = 3;
const double RI = 0.144, CHRI = -0.005;

// Model wide variables
extern unsigned xrandom;
extern int NTOT;
extern int NTOTDUR;
extern int NTOTNDUR;


class Market {

  public:

  Market() {}


  double MTEC, QDMTEC;

  double BETA;

  double EXPXDP, EXPXDS, EXPXDW;

  double IMP, TMIMP;
  double TMX;

  double QTSP;
  double QTBUY;

  double QPRELPDOM;
  double QPDOM, QPFOR;
  double QDPDOM, QDPFOR;

  double PT, QPH;

  double QCHTSTO;
  double QMAXTSUDOM;
  double REDUCE;
};

class Household {

  public:

  Household() {}


  static inline double ALFA1[MAXNMH], ALFA2[MAXNMH], ALFA3, ALFA4;
  static inline double BETA1[MAXNMH], BETA2[MAXNMH], BETA3[MAXNMH];
  static inline double SMOOTH[MAXNMH];
  static inline double RHODUR;
  

  double NH;

  double CVA[MAXNMH];
  double CVE[MAXNMH];

  double QDI;
  double STODUR;
  double WH, WHRA;

  double OLDQCPI, QPRELCPI, QDCPI, QCPI, CHDCPI;
  double SWAP;

  double QSP[MAXNMH], QSPE[MAXNMH];
  double QC[MAXNMH];

  double QSAVH;
};

class Firm {

  public:

  Firm() {}


  static inline double ALFABW, BETABW;
  static inline double RHO;
  static inline double SMP, SMW, SMS, SMT;
  static inline double FIP, FIS, FIW;
  static inline double KSIFAIL, KSISUCC;
  static inline double TMSTO;
  static inline double IOTA;
  static inline double E1, E2, R;
  static inline double EPS;
  static inline double RESDOWN, RESMAX;
  static inline double NRS;
  static inline double LOSS;


  int iF;
  int iM;

  // Year
  double MHIST;

  double Q, P, W, S, M;
  double DQ, DP, DS, DW, CHM;
  double CUMQ, CUMM, CUMSU, CUMS, CUMWS, CUML;

  double TEC;

  double TARGM;

  double EXPIDP, EXPIDS, EXPIDW;
  double EXPDP, EXPDS, EXPDW;

  // Quarter
  double QQ, QP, QS, QW, QM, QSU;
  double QDQ, QDP, QDS, QDW;

  double QSDOM, QSFOR;
  double QSUDOM, QSUFOR;
  double QOPTSU, QOPTSUDOM;

  double L, RES, QTOP;
  double QCHL, QCHW;
  double QCHQTOP, QCHQTOP1, QCHQTOP2;
  double AMAN[3];
  double CHL, LL, WW;

  double STO, MINSTO, MAXSTO, OPTSTO;
  double QCHSTO;
  double SMALL, BIG;
  double X;

  double QPLANQ, QPLANL;

  double QTARGM;

  double QEXPP, QEXPS, QEXPW, QEXPSU;
  double QEXPDP, QEXPDS, QEXPDW;

  // Investment
  double INVEFF;
  double QINV, QINVLAG;

  double QRR;

  double NW;
  double BW, CHBW, QCHBW;
  double RW;

  double K1, K2;
  double QCHK2;
  double QCHS;
};

class Labour {

  public:

  Labour() {}


  double LF;
  double LU;
  double RET;
  double ENTRY;

  double RU, CHRU;

  double QDWIND;

  double GAMMA;
  double THETA;

  double LL;
};

class Service {

  public:

  Service() {}


  double QTSP, QTBUY;
  double REDUCE;

  double QQZ, QPZ, QSZ, QWZ, QMZ;
  double TARGMZ, QTARGMZ;

  double PT, QPH;
  double QPRELPZ;

  double IMP;

  double LZ, QCHLZ;

  double TECZ, QDTECZ;
};

class Government {

  public:

  Government() {}


  double LG, QWG;
  double QCHLG;
  double REALCHLG;
};


// Helper functions
double QFR(Firm& F, const double L);
double RFQ(Firm& F, const double Q);
bool SAT(Firm& F, const double Q, const double L);
bool SOLVE(Firm& F);
unsigned XorShift32(unsigned& x);
int CHOOSE(Firm F[MAXNTOT], Labour& L, const int kThisFirm, bool& is_firm);
void NULLIFY_FIRM(Firm F[MAXNTOT], const int firm);
