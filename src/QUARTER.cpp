#include "QUARTER.h"



inline void QUARTERLY_EXP(Firm& F, const bool kFirstQ) {
  F.QEXPDP = F.EXPDP / 4.0;
  F.QEXPDW = F.EXPDW / 4.0;
  F.QEXPDS = F.EXPDS / 4.0;

  if (!kFirstQ) {
    F.QEXPDP += Firm::FIP * (F.QDP - F.QEXPDP);
    F.QEXPDW += Firm::FIW * (F.QDW - F.QEXPDW);
    F.QEXPDS += Firm::FIS * (F.QDS - F.QEXPDS);
  }

  F.QEXPP = F.QP * (1.0 + F.QEXPDP);
  F.QEXPW = F.QW * (1.0 + F.QEXPDW);
  F.QEXPS = F.QS * (1.0 + F.QEXPDS);
}

inline void QUARTERLY_TARG(Firm& F) {
  F.QTARGM = F.TARGM + (Firm::NRS - 1.0) / (5.0 - Firm::NRS) * (F.TARGM - F.CUMM);
  //F.QTARGM = F.TARGM;
}





#pragma region PRODPLAN
inline void LUUPDATE(Firm F[MAXNTOT], Labour& L, Service& Z, Government& G) {
  double SUM_L = 0;
  for (int I = 0; I < NTOT; ++I) {
    SUM_L += F[I].L;
    F[I].L *= (1.0 - L.RET);

    for (int i = 0; i < 3; ++i) {
      F[I].AMAN[i] *= (1.0 - L.RET);
    }
  }

  L.LF = L.LU + Z.LZ + G.LG + SUM_L;
  L.LU *= (1.0 - L.RET);
  L.LU += L.ENTRY * L.LF;
}

inline void PRODFRONT(Firm F[MAXNTOT], Market M[MAXNM]) {
  for (int I = 0; I < NM; ++I) {
    M[I].MTEC *= (1.0 + M[I].QDMTEC);
  }

  for (int I = 0; I < NTOT; ++I) {
    F[I].QTOP *= (1.0 - Firm::RHO);

    F[I].QCHQTOP1 = (1.0 - Firm::LOSS) * (F[I].QINV * F[I].INVEFF) / F[I].QP;
    F[I].QCHQTOP2 = std::min(Firm::LOSS * ((F[I].QINV * F[I].INVEFF) / F[I].QP) * ((Firm::RESMAX - F[I].RES) / Firm::RESMAX), ((Firm::RESMAX - F[I].RES) / (1.0 - Firm::RESMAX)) * (F[I].QTOP + F[I].QCHQTOP1));
    F[I].QCHQTOP = F[I].QCHQTOP1 + F[I].QCHQTOP2;

    F[I].RES = (F[I].RES * (F[I].QTOP + F[I].QCHQTOP1) + F[I].QCHQTOP2) / (F[I].QTOP + F[I].QCHQTOP);

    F[I].TEC = (F[I].QTOP + F[I].QCHQTOP) / ((F[I].QTOP / F[I].TEC) + (F[I].QCHQTOP / M[F[I].iM].MTEC));
    F[I].QTOP += F[I].QCHQTOP;
  }
}

inline void INITPRODPLAN(Firm& F) {
  F.QEXPSU = F.QEXPS / F.QEXPP;
  F.QPLANQ = std::max(0.0, F.QEXPSU + (F.OPTSTO - F.STO) / (4.0 * Firm::TMSTO));
}

inline bool TARGSEARCH(Firm& F, Labour& L) {
  bool path[11] = {};
  double Q2, Q7;

  if (F.QPLANQ >= F.QTOP * (1.0 - F.RES)) path[6] = true;
  else if (F.QPLANQ > QFR(F, F.L)) path[5] = true;
  else path[1] = true;

  if (path[1]) {
    if (SAT(F, F.QPLANQ, F.L)) {
      F.QPLANL = F.L;
      path[10] = true;
    }

    else path[2] = true;
  }

  if (path[2]) {
    Q2 = std::min(QFR(F, F.L), F.QEXPSU + F.MAXSTO - F.STO);

    if (SAT(F, Q2, F.L)) {
      F.QPLANQ = (F.L * (F.QEXPW / 4.0)) / ((1.0 - F.QTARGM) * F.QEXPP);
      F.QPLANL = F.L;
      path[10] = true;
    }

    else if (Q2 == QFR(F, F.L)) path[4] = true;

    else path[3] = true;
  }

  if (path[3]) {
    if (SAT(F, Q2, RFQ(F, Q2))) {
      F.QPLANQ = Q2;
      F.QPLANL = ((1.0 - F.QTARGM) * Q2 * F.QEXPP) / (F.QEXPW / 4.0);
      path[10] = true;
    }

    else path[4] = true;
  }

  if (path[4]) {
    if (SAT(F, F.QPLANQ, RFQ(F, F.QPLANQ))) {
      SOLVE(F);
      path[10] = true;
    }

    else {
      Q7 = F.QPLANQ;
      path[7] = true;
    }
  }

  if (path[5]) {
    if (SAT(F, F.QPLANQ, RFQ(F, F.QPLANQ))) {
      F.QPLANL = RFQ(F, F.QPLANQ);
      path[10] = true;
    }

    else path[6] = true;
  }

  if (path[6]) {
    if (SAT(F, QFR(F, F.L), F.L)) {
      SOLVE(F);
      path[10] = true;
    }

    else {
      Q7 = QFR(F, F.L);
      path[7] = true;
    }
  }

  if (path[7]) {
    if (SAT(F, Q7, RFQ(F, ((1.0 - F.RES) / (1.0 - Firm::RESDOWN * F.RES)) * Q7))) {
      F.QPLANQ = Q7;
      F.QPLANL = ((1.0 - F.QTARGM) * Q7 * F.QEXPP) / (F.QEXPW / 4.0);
      F.RES = 1.0 - (Q7 * (1.0 - F.RES)) / (QFR(F, F.QPLANL));
      path[10] = true;
    }

    else {
      F.RES = Firm::RESDOWN * F.RES;
      path[8] = true;
    }
  }

  if (path[8]) {
    if (SAT(F, 0, 0)) {
      SOLVE(F);
      path[10] = true;
    }

    else path[9] = true;
  }

  if (path[9]) {
    L.LU += F.L;
    return false;
  }

  if (path[10]) {
    double LAYOFF = std::max(F.L - F.QPLANL, 0.0);
    F.AMAN[0] = std::min(LAYOFF, F.AMAN[1]);
    F.AMAN[1] = std::min(LAYOFF - F.AMAN[0], F.AMAN[2]);
    F.AMAN[2] = LAYOFF - F.AMAN[0] - F.AMAN[1];
  }

  return true;
}


inline void PRODPLAN(Firm F[MAXNTOT], Market M[MAXNM], Labour& L, Service& Z, Government& G) {

  // Update Labour Unemployment
  LUUPDATE(F, L, Z, G);

  // Determine Change In Production Frontier
  PRODFRONT(F, M);

  for (int I = 0; I < NTOT; ++I) {

    // Initialize Quarterly Production Plan
    INITPRODPLAN(F[I]);

    // Search For Target Satisfaction
    if (!TARGSEARCH(F[I], L)) {
      NULLIFY_FIRM(F, I);
      --I;
    }
  }
}
#pragma endregion PRODPLAN





inline void EXPORT(Firm F[MAXNTOT], Market M[MAXNM]) {
  for (int I = 0; I < NTOT; ++I) {
    double QPDOM = M[F[I].iM].QPDOM;
    double QPFOR = M[F[I].iM].QPFOR;
    double TMX = M[F[I].iM].TMX;

    if (QPDOM >= QPFOR) {
      F[I].X -= F[I].X * (1.0 / (4.0 * TMX)) * ((QPDOM - QPFOR) / QPFOR);
      F[I].X = std::max(0.0, F[I].X);
    }

    else {
      F[I].X += (1.0 - F[I].X) * (1.0 / (4.0 * TMX)) * ((QPFOR - QPDOM) / QPDOM);
      F[I].X = std::min(1.0, F[I].X);
    }

    F[I].QSUFOR = F[I].X * F[I].QOPTSU;
  }

  for (int I = 0; I < NM; ++I) {
    M[I].QPFOR *= (1.0 + M[I].QDPFOR);
  }

  for (int I = 0; I < NTOT; ++I) {
    F[I].QSFOR = F[I].QSUFOR * M[F[I].iM].QPFOR;
  }
}



#pragma region STOCK_SYSTEM
inline void FIRMSTO(Firm F[MAXNTOT], Market M[MAXNM]) {
  double SUM[MAXNM][2];
  for (int I = 0; I < NTOT; ++I) {
    F[I].QCHSTO = 0;

    if (F[I].STO > F[I].MAXSTO) {
      double CH = F[I].STO - F[I].MAXSTO;
      M[F[I].iM].QCHTSTO += CH;
      F[I].STO = F[I].MAXSTO;
      F[I].QCHSTO += CH;
    }

    else if (F[I].STO < F[I].MINSTO) {
      double CH = F[I].STO - F[I].MINSTO;
      M[F[I].iM].QCHTSTO += CH;
      F[I].STO = F[I].MINSTO;
      F[I].QCHSTO += CH;
    }

    SUM[F[I].iM][0] += F[I].MAXSTO - F[I].STO;
    SUM[F[I].iM][1] += F[I].MINSTO - F[I].STO;
  }

  for (int I = 0; I < NTOT; ++I) {
    if (M[F[I].iM].QCHTSTO > 0) {
      double CH = (F[I].MAXSTO - F[I].STO) / SUM[F[I].iM][0];
      F[I].STO += CH;
      F[I].QCHSTO += CH;
    }

    else {
      double CH = (F[I].MINSTO - F[I].STO) / SUM[F[I].iM][1];
      F[I].STO += CH;
      F[I].QCHSTO += CH;
    }
  }

  for (int I = 0; I < NTOT; ++I) {
    F[I].QSUDOM = F[I].QQ - F[I].QSUFOR - F[I].QCHSTO;
    F[I].QSDOM = F[I].QSUDOM * M[F[I].iM].QPDOM;
  }
}


inline void STOSYSTEM(Firm F[MAXNTOT], Market M[MAXNM]) {

  FIRMSTO(F, M);

  for (int I = 0; I < NTOT; ++I) {
    F[I].MINSTO = F[I].SMALL * (4.0 * (F[I].QS / F[I].QP));
    F[I].MAXSTO = F[I].BIG * (4.0 * (F[I].QS / F[I].QP));
    F[I].OPTSTO = F[I].MINSTO + M[F[I].iM].BETA * (F[I].MAXSTO - F[I].MINSTO);
  }
}
#pragma endregion STOCK_SYSTEM





inline void FINALQPQSQM(Firm F[MAXNTOT], Service& Z) {
  for (int I = 0; I < NTOT; ++I) {
    F[I].QSU = F[I].QSUFOR + F[I].QSUDOM;

    F[I].QDS = (F[I].QSFOR + F[I].QSDOM) / F[I].QS - 1.0;
    F[I].QS = F[I].QSFOR + F[I].QSDOM;

    F[I].QDP = (F[I].QS / F[I].QSU) / F[I].QP - 1.0;
    F[I].QP = F[I].QS / F[I].QSU;

    if (!F[I].QS) F[I].QM = 0;
    else F[I].QM = 1.0 - (F[I].L * (F[I].QW / 4.0)) / F[I].QS;
  }

  if (!Z.QSZ) Z.QMZ = 0;
  else Z.QMZ = 1.0 - (Z.LZ * (Z.QWZ / 4.0)) / Z.QSZ;
}

inline void QUARTERLY_CUM(Firm& F) {
  F.CUMQ += F.QQ;
  F.CUMS += F.QS;
  F.CUMSU += F.QSU;
  F.CUMWS += F.L * (F.QW / 4.0);
  F.CUML = ((Firm::NRS - 1.0) * F.CUML + F.L) / Firm::NRS;
  F.CUMM = 1.0 - F.CUMWS / F.CUMS;
}


inline void QUARTERLY_RESULT(Firm F[MAXNTOT], Service& Z) {

  FINALQPQSQM(F, Z);

  for (int I = 0; I < NTOT; ++I) {
    QUARTERLY_CUM(F[I]);
  }
}





inline void INVFIN(Firm& F, Market M[MAXNM]) {
  F.K1 = (1.0 - Firm::RHO + M[DUR].QDPDOM) * F.K1 + F.QINV * (1.0 - Firm::RHO);
  F.QRR = 4.0 * ((F.QM * F.QS - Firm::RHO * F.K1) / (F.K1 + F.K2 + F.STO * F.QP));

  F.QCHS = F.QS * (F.QDS / (1.0 + F.QDS));
  F.QCHK2 = F.RW * 4.0 * F.QCHS;
  F.K2 += F.QCHK2;

  F.QCHBW = F.BW * (Firm::ALFABW + Firm::BETABW * (F.QRR / 4.0 + M[DUR].QDPDOM - RI / 4.0));
  F.BW += F.QCHBW;
  F.NW = F.K1 + F.K2 + F.STO * F.QP - F.BW;

  F.QINV = F.QINVLAG;
  F.QINVLAG = std::max(0.0, F.QM * F.QS - F.QCHK2 + F.QCHBW - (RI / 4.0) * F.BW);

  F.INVEFF = (F.QTOP * F.QP) / F.K1;
}





void QUARTER(Firm F[MAXNTOT], Market M[MAXNM], Household& H, Labour& L, Service& Z, Government& G) {

  for (int q = 0; q < 4; ++q) {
    for (int I = 0; I < NTOT; ++I) {
      QUARTERLY_EXP(F[I], q == 0);
      QUARTERLY_TARG(F[I]);
    }

    PRODPLAN(F, M, L, Z, G);

    LABOUR_MARKET(F, L, Z, G);

    EXPORT(F, M);

    DOMESTIC_MARKET(F, M, H, L, Z, G);

    STOSYSTEM(F, M);

    QUARTERLY_RESULT(F, Z);

    for (int I = 0; I < NTOT; ++I) {
      INVFIN(F[I], M);
    }
  }
}
