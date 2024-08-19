#include "LABOR.h"


inline void LABOUR_SEARCH_INPUT(Firm F[MAXNTOT], Labour& L) {
  for (int I = 0; I < NTOT; ++I) {
    F[I].CHL = F[I].QPLANL - F[I].L;
    F[I].WW = F[I].QW + Firm::IOTA * (F[I].QEXPW - F[I].QW);
    F[I].LL = F[I].L;
  }

  L.LL = L.LU;
}

inline void CONFRONT(Firm F[MAXNTOT], Labour& L) {
  std::sort(F, F + NTOT, [](const Firm& a, const Firm& b) {return (a.CHL / a.L) > (b.CHL / b.L); });

  for (int NIT = 0; NIT < NITER; ++NIT) {
    for (int I = 0; I < NTOT; ++I) {
      if (F[I].CHL <= 0) continue;

      bool is_firm;
      int II = CHOOSE(F, L, I, is_firm);

      if (is_firm) {
        if (F[I].WW >= F[II].WW * (1.0 + L.GAMMA)) {
          F[II].WW += Firm::KSISUCC * (F[I].WW - F[II].WW);

          double CHLNOW = std::min(L.THETA * F[II].LL, F[I].CHL);
          F[I].LL += CHLNOW;
          F[I].CHL -= CHLNOW;
          F[II].LL -= CHLNOW;
          F[II].CHL += CHLNOW;
        }

        else {
          F[I].WW += Firm::KSIFAIL * (F[II].WW * (1.0 + L.GAMMA) - F[I].WW);
          continue;
        }
      }

      else {
        double CHLNOW = std::min(L.THETA * L.LL, F[I].CHL);
        F[I].LL += CHLNOW;
        F[I].CHL -= CHLNOW;
        L.LL -= CHLNOW;
      }
    }
  }
}

inline void LABOUR_SEARCH_OUTPUT(Firm F[MAXNTOT], Labour& L) {
  for (int I = 0; I < NTOT; ++I) {
    F[I].QCHL = F[I].LL - F[I].L;
    F[I].QCHW = F[I].WW - F[I].QW;

    double EXIT = std::max(0.0, -F[I].QCHL);

    if (EXIT > F[I].AMAN[0] + F[I].AMAN[1]) {
      F[I].AMAN[2] -= (EXIT - F[I].AMAN[0] - F[I].AMAN[1]);
      F[I].AMAN[2] = std::max(0.0, F[I].AMAN[2]);
    }

    if (EXIT > F[I].AMAN[0]) {
      F[I].AMAN[1] -= (EXIT - F[I].AMAN[0]);
      F[I].AMAN[1] = std::max(0.0, F[I].AMAN[1]);
    }

    if (EXIT > 0) {
      F[I].AMAN[0] -= EXIT;
      F[I].AMAN[0] = std::max(0.0, F[I].AMAN[0]);
    }
  }

  L.LU = L.LL;
}


inline void LABOUR_SEARCH(Firm F[MAXNTOT], Labour& L) {

  LABOUR_SEARCH_INPUT(F, L);

  CONFRONT(F, L);

  LABOUR_SEARCH_OUTPUT(F, L);
}

inline void LABOUR_UPDATE(Firm F[MAXNTOT], Labour& L, Service& Z, Government& G) {
  double SUM[5] = {};
  for (int I = 0; I < NTOT; ++I) {
    double SACK = std::min(F[I].AMAN[0], std::max(0.0, F[I].L + F[I].QCHL - F[I].QPLANL));
    F[I].QCHL -= SACK;
    F[I].AMAN[0] -= SACK;

    SUM[0] += SACK;
    SUM[1] += F[I].L * F[I].QW;
    SUM[2] += F[I].L;
    SUM[3] += (F[I].L + F[I].QCHL) * (F[I].QW + F[I].QCHW);
    SUM[4] += F[I].L + F[I].QCHL;

    F[I].L += F[I].QCHL;

    F[I].QDW = F[I].QCHW / F[I].QW;
    F[I].QW += F[I].QCHW;
  }

  L.LU += SUM[0];
  double OLDQW = SUM[1] / SUM[2];
  double NEWQW = SUM[3] / SUM[4];
  L.QDWIND = NEWQW / OLDQW - 1.0;
  L.CHRU = L.LU / (L.LU + Z.LZ + G.LG + SUM[2]) - L.RU;
  L.RU += L.CHRU;
}

inline void PLANQREVISE(Firm& F) {
  F.QPLANQ = std::min(F.QPLANQ, QFR(F, F.L));

  if (F.QQ) {
    F.QDQ = F.QPLANQ / F.QQ - 1.0;
    F.QQ *= (1.0 + F.QDQ);
  }

  else {
	fprintf(stderr, "[ERROR] QQ = 0, avoided division by 0\n");
    F.QQ = F.QPLANQ;
  }

  F.QOPTSU = std::max(0.0, F.QEXPSU * (F.QQ / (F.QEXPSU + ((F.OPTSTO - F.STO) / (4.0 * Firm::TMSTO)))));
}


inline void ZLABOUR(Labour& L, Service& Z) {
  Z.TECZ *= (1.0 + Z.QDTECZ);

  Z.QCHLZ = std::min(L.LU, ((Z.QMZ - Z.QTARGMZ) * Z.QPZ * Z.TECZ * Z.LZ) / (Z.QWZ / 4.0) + L.RET * Z.LZ);
  Z.LZ = std::max(0.0, Z.LZ + Z.QCHLZ - L.RET * Z.LZ);
  L.LU -= Z.QCHLZ;

  Z.QWZ *= (1.0 + L.QDWIND);
  Z.QQZ = Z.TECZ * Z.LZ;

  Z.QPRELPZ = Z.QPZ * (1.0 + L.QDWIND - Z.QDTECZ);
}

inline void GLABOUR(Labour& L, Government& G) {
  G.QCHLG = std::min(L.LU, G.LG * L.RET + G.REALCHLG);
  G.LG += G.QCHLG - L.RET * G.LG;
  L.LU -= G.QCHLG;

  G.QWG *= (1.0 + L.QDWIND);
}

inline void INDALABOUR(Firm F[MAXNTOT], Labour& L, Service& Z, Government& G) {

  LABOUR_SEARCH(F, L);

  LABOUR_UPDATE(F, L, Z, G);

  for (int I = 0; I < NTOT; ++I) {
    PLANQREVISE(F[I]);
  }
}


void LABOUR_MARKET(Firm F[MAXNTOT], Labour& L, Service& Z, Government& G) {

  ZLABOUR(L, Z);

  GLABOUR(L, G);

  INDALABOUR(F, L, Z, G);
}
