#include "YEAR.h"


inline void YEARLY_INIT(Firm& F) {
  F.CUMQ = 0;
  F.CUMM = 0;
  F.CUMSU = 0;
  F.CUMS = 0;
  F.CUMWS = 0;
  F.CUML = 0;
}

inline void YEARLY_EXP(Firm& F, Market M[MAXNM]) {
  F.EXPIDP = Firm::SMP * F.EXPIDP + (1.0 - Firm::SMP) * (F.DP + Firm::E1 * (F.DP - F.EXPDP) - Firm::E2 * ((F.DP - F.EXPDP) * (F.DP - F.EXPDP)));
  F.EXPDP = (1.0 - Firm::R) * F.EXPIDP + Firm::R * M[F.iM].EXPXDP;

  F.EXPIDW = Firm::SMW * F.EXPIDW + (1.0 - Firm::SMW) * (F.DW + Firm::E1 * (F.DW - F.EXPDW) - Firm::E2 * ((F.DW - F.EXPDW) * (F.DW - F.EXPDW)));
  F.EXPDW = (1.0 - Firm::R) * F.EXPIDW + Firm::R * M[F.iM].EXPXDW;

  F.EXPIDS = Firm::SMS * F.EXPIDS + (1.0 - Firm::SMS) * (F.DS + Firm::E1 * (F.DS - F.EXPDS) - Firm::E2 * ((F.DS - F.EXPDS) * (F.DS - F.EXPDS)));
  F.EXPDS = (1.0 - Firm::R) * F.EXPIDS + Firm::R * M[F.iM].EXPXDS;

  // Exogenous:
  // EXPXDP, EXPXDW, EXPXDS
}

inline void YEARLY_TARG(Firm& F) {
  F.MHIST = Firm::SMT * F.MHIST + (1.0 - Firm::SMT) * F.M;
  F.TARGM = F.MHIST * (1.0 + Firm::EPS);
}

inline void YEARLY_UPDATE(Firm& F) {
  F.DQ = F.CUMQ / F.Q - 1.0;
  F.Q *= (1.0 + F.DQ);

  F.DP = (F.CUMS / F.CUMSU) / F.P - 1.0;
  F.P *= (1.0 + F.DP);

  F.DW = (F.CUMWS / F.CUML) / F.W - 1.0;
  F.W *= (1.0 + F.DW);

  F.DS = F.CUMS / F.S - 1.0;
  F.S *= (1.0 + F.DS);

  F.CHM = F.CUMM - F.M;
  F.M += F.CHM;
}


void YEAR(Firm F[MAXNTOT], Market M[MAXNM], Household& H, Labour& L, Service& Z, Government& G) {

  // Yearly Initialization, Expectation, and Target
  for (int I = 0; I < NTOT; ++I) {
    YEARLY_INIT(F[I]);
    YEARLY_EXP(F[I], M);
    YEARLY_TARG(F[I]);
  }

  // Quarterly Execution
  QUARTER(F, M, H, L, Z, G);

  // Yearly Update
  for (int I = 0; I < NTOT; ++I) {
    YEARLY_UPDATE(F[I]);
  }
}
