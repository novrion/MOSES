#include "DOM.h"


inline void COMPUTE_SPENDING(Household& H, Market M[MAXNM], Labour& L, Service& Z) {
    double SUM1[2] = {};
    for (int j = 0; j < NM; ++j) {
      SUM1[0] += H.QC[j];
      SUM1[1] += H.QC[j] / M[j].PT;
    }
    SUM1[0] += H.QC[iZ];
    SUM1[1] += H.QC[iZ] / Z.PT;

    H.QPRELCPI = SUM1[0] / SUM1[1];
    H.CHDCPI = H.QPRELCPI / H.QCPI - 1.0 - H.QDCPI;

    H.QSPE[NDUR] = H.CVE[NDUR] * M[NDUR].PT;
    H.QSPE[iZ] = H.CVE[iZ] * Z.PT;

    H.SWAP = Household::ALFA3 * (CHRI / 4.0 - H.CHDCPI) + Household::ALFA4 * L.CHRU;

    H.QSPE[DUR] = ((M[DUR].PT * H.CVE[DUR]) / Household::RHODUR) - (M[DUR].PT / M[DUR].QPH) * H.STODUR - H.QDI * H.SWAP;
    H.QSPE[SAV] = (H.WHRA * H.QDI - H.WH) + H.QDI * H.SWAP;

	double SUM2 = 0;
	for (int i = 0; i < NMH; ++i) { SUM2 += Household::BETA1[i] * H.QSPE[i]; }

    H.QSP[DUR]  = std::max(0.0, Household::BETA1[DUR]  * H.QSPE[DUR]  + (Household::BETA2[DUR]  + Household::BETA3[DUR]  / (H.QDI / H.QPRELCPI)) * (H.QDI - SUM2));
    H.QSP[NDUR] = std::max(0.0, Household::BETA1[NDUR] * H.QSPE[NDUR] + (Household::BETA2[NDUR] + Household::BETA3[NDUR] / (H.QDI / H.QPRELCPI)) * (H.QDI - SUM2));
    H.QSP[SAV]  =               Household::BETA1[SAV]  * H.QSPE[SAV]  + (Household::BETA2[SAV]  + Household::BETA3[SAV]  / (H.QDI / H.QPRELCPI)) * (H.QDI - SUM2);
    H.QSP[iZ] = std::max(0.0,   Household::BETA1[iZ]   * H.QSPE[iZ]   + (Household::BETA2[iZ]   + Household::BETA3[iZ]   / (H.QDI / H.QPRELCPI)) * (H.QDI - SUM2));
}

inline void COMPUTE_BUYING(Firm F[MAXNTOT], Market M[MAXNM], Household& H, Service& Z) {
  double SUM2 = 0;
  for (int I = 0; I < NTOT; ++I) { SUM2 += F[I].QINVLAG; }

  for (int I = 0; I < NM; ++I) {
    M[I].QTSP = H.QSP[I] * H.NH;
    if (I == DUR) { M[I].QTSP += SUM2; }

    M[I].QTBUY = (1.0 - M[I].IMP) * (M[I].QTSP / M[I].PT);
  }

  Z.QTSP = H.QSP[iZ] * H.NH;
  Z.QTBUY = Z.QTSP / Z.PT;
}

inline void PRICE_ADJUST(Firm F[MAXNTOT], Market M[MAXNM]) {
  double SUM[MAXNM] = {};
  for (int I = 0; I < NTOT; ++I) { SUM[F[I].iM] += F[I].QOPTSUDOM; }

  for (int I = 0; I < NM; ++I) {
    if (M[I].QTBUY < SUM[I]) {
      M[I].PT -= (abs(M[I].EXPXDP) * M[I].PT) / (4.0 * (MARKETITER - 1.0));
    }

    else {
      M[I].PT += (abs(M[I].EXPXDP) * M[I].PT) / (4.0 * (MARKETITER - 1.0));
    }
  }
}


inline void MARKET_ENTRANCE(Firm F[MAXNTOT], Market M[MAXNM]) {
  double SUM1[MAXNM] = {};
  double SUM2[MAXNM] = {};
  for (int I = 0; I < NTOT; ++I) {
    F[I].QOPTSUDOM = (1.0 - F[I].X) * F[I].QOPTSU;
    SUM1[F[I].iM] += F[I].QOPTSUDOM * (F[I].QEXPP / F[I].QP);
    SUM2[F[I].iM] += F[I].QOPTSUDOM;
  }

  for (int I = 0; I < NM; ++I) {
    M[I].QPRELPDOM = M[I].QPDOM * (SUM1[I] / SUM2[I]);
  }
}

inline void HOUSEHOLD_INIT(Firm F[MAXNTOT], Market M[MAXNM], Household& H, Service& Z, Government& G) {
  double SUM = Z.QMZ * Z.QSZ + Z.LZ * (Z.QWZ / 4.0) + G.LG * (G.QWG / 4.0);
  for (int I = 0; I < NTOT; ++I) { SUM += F[I].L * (F[I].QW / 4.0); }
  H.QDI = SUM / H.NH + H.WH * (RI / 4.0);

  for (int I = 0; I < NMH; ++I) {
    if (I == SAV) continue;
	H.CVE[I] = Household::ALFA1[I] + Household::ALFA2[I] * H.CVA[I];
  }
}

inline void MARKET_CONFRONT(Firm F[MAXNTOT], Market M[MAXNM], Household& H, Labour& L, Service& Z) {
  for (int I = 0; I < NM; ++I) {
    if (M[I].QPDOM >= M[I].QPFOR) {
      M[I].IMP += ((1.0 - M[I].IMP) / (4.0 * M[I].TMIMP)) * ((M[I].QPDOM - M[I].QPFOR) / M[I].QPFOR);
      M[I].IMP = std::min(1.0, M[I].IMP);
    }
    else {
      M[I].IMP -= (M[I].IMP / (4.0 * M[I].TMIMP)) * ((M[I].QPFOR - M[I].QPDOM) / M[I].QPDOM);
      M[I].IMP = std::max(0.0, M[I].IMP);
    }

    M[I].PT = M[I].QPRELPDOM;
  }

  Z.PT = Z.QPRELPZ;

  for (int I = 0; I < MARKETITER; ++I) {
    COMPUTE_SPENDING(H, M, L, Z);
    COMPUTE_BUYING(F, M, H, Z);
    if (I < MARKETITER - 1.0) PRICE_ADJUST(F, M);
  }
}

inline void MINSTO_ADJUST(Firm F[MAXNTOT], Market M[MAXNM], Household& H, Service& Z) {
  double SUM[MAXNM] = {};
  for (int I = 0; I < NTOT; ++I) { SUM[F[I].iM] += F[I].QQ + (F[I].STO - F[I].MINSTO) - F[I].QSUFOR; }

  for (int I = 0; I < NM; ++I) {
    M[I].QMAXTSUDOM = std::max(0.0, SUM[I]);

    if (!M[I].QTBUY) {
	  fprintf(stderr, "[ERROR] Forced REDUCE = 1 to avoid division by zero\n");
	  M[I].REDUCE = 1.0;
	}
    else M[I].REDUCE = std::min(1.0, M[I].QMAXTSUDOM / M[I].QTBUY);
    H.QSP[I] *= M[I].REDUCE;
	M[I].QTBUY *= M[I].REDUCE;
  }


  if (!Z.QTBUY) {
    fprintf(stderr, "[ERROR] Forced REDUCE = 1 to avoid division by zero\n");
    Z.REDUCE = 1.0;
  }
  else Z.REDUCE = std::min(1.0, Z.QQZ / Z.QTBUY);
  H.QSP[iZ] *= Z.REDUCE;
  Z.QTBUY *= Z.REDUCE;


  for (int I = 0; I < NTOT; ++I) { F[I].QINVLAG *= M[DUR].REDUCE; }
}

inline void DOMESTIC_RESULT(Firm F[MAXNTOT], Market M[MAXNM], Service& Z) {
  double SUM[MAXNM][2];
  for (int I = 0; I < NTOT; ++I) {
    SUM[F[I].iM][0] += F[I].MAXSTO - F[I].STO;
    SUM[F[I].iM][1] += F[I].QQ - F[I].QSUFOR;
  }

  for (int I = 0; I < NM; ++I) {
    M[I].QDPDOM = M[I].PT / M[I].QPDOM - 1.0;
    M[I].QPDOM = M[I].PT;
    M[I].QCHTSTO = std::min(SUM[I][0], SUM[I][1] - M[I].QTBUY);
  }

  Z.QPZ = Z.PT;
  Z.QSZ = Z.QTBUY * Z.QPZ;
}

inline void HOUSEHOLD_UPDATE(Household& H, Market M[MAXNM], Service& Z) {
  H.QC[NDUR] = H.QSP[NDUR];
  H.QC[iZ] = H.QSP[iZ];

  H.STODUR = (M[DUR].PT / M[DUR].QPH) * H.STODUR + H.QSP[DUR];
  H.QC[DUR] = Household::RHODUR * H.STODUR;
  H.STODUR = (1.0 - Household::RHODUR) * H.STODUR;

  H.QSAVH = H.QDI - (H.QSP[DUR] + H.QSP[NDUR] + H.QSP[iZ]);
  H.QSP[SAV] = H.QSAVH;
  H.WH += H.QSAVH;

  H.CVA[DUR]  = Household::SMOOTH[DUR]  * H.CVA[DUR]  + (1.0 - Household::SMOOTH[DUR])  * (H.QC[DUR]  / M[DUR].PT);
  H.CVA[NDUR] = Household::SMOOTH[NDUR] * H.CVA[NDUR] + (1.0 - Household::SMOOTH[NDUR]) * (H.QC[NDUR] / M[NDUR].PT);
  H.CVA[iZ]   = Household::SMOOTH[iZ]   * H.CVA[iZ]   + (1.0 - Household::SMOOTH[iZ])   * (H.QC[iZ]   / Z.PT);

  H.WHRA = Household::SMOOTH[SAV] * H.WHRA + (1.0 - Household::SMOOTH[SAV]) * (H.WH / H.QDI);

  for (int j = 0; j < NM; ++j) { M[j].QPH = M[j].PT; }
  Z.QPH = Z.PT;

  H.OLDQCPI = H.QCPI;
  H.QCPI = (H.QC[DUR] + H.QC[NDUR] + H.QC[iZ]) / ((H.QC[DUR] + H.QC[NDUR] + H.QC[iZ]) / (M[DUR].QPH + M[NDUR].QPH + Z.QPH));
  H.QDCPI = (H.QCPI - H.OLDQCPI) / H.OLDQCPI;
}


void DOMESTIC_MARKET(Firm F[MAXNTOT], Market M[MAXNM], Household& H, Labour& L, Service& Z, Government& G) {

  MARKET_ENTRANCE(F, M);

  HOUSEHOLD_INIT(F, M, H, Z, G);

  MARKET_CONFRONT(F, M, H, L, Z);

  MINSTO_ADJUST(F, M, H, Z);

  DOMESTIC_RESULT(F, M, Z);

  HOUSEHOLD_UPDATE(H, M, Z);
}
