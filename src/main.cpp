#include "YEAR.h"

// Model wide variables (declared extern in MOSES.h)
unsigned xrandom = 878726348;
int NTOT = 256;
int NTOTDUR = 0;
int NTOTNDUR = 0;


inline void INPUT(Firm F[MAXNTOT], Market M[MAXNM], Household& H, Labour& L, Service& Z, Government& G) {

  // Firms
  scanf("%lf", &Firm::ALFABW);
  scanf("%lf", &Firm::BETABW);
  scanf("%lf", &Firm::RHO);
  scanf("%lf", &Firm::SMP);
  scanf("%lf", &Firm::SMW);
  scanf("%lf", &Firm::SMS);
  scanf("%lf", &Firm::SMT);
  scanf("%lf", &Firm::FIP);
  scanf("%lf", &Firm::FIS);
  scanf("%lf", &Firm::FIW);
  scanf("%lf", &Firm::KSIFAIL);
  scanf("%lf", &Firm::KSISUCC);
  scanf("%lf", &Firm::TMSTO);
  scanf("%lf", &Firm::IOTA);
  scanf("%lf", &Firm::E1);
  scanf("%lf", &Firm::E2);
  scanf("%lf", &Firm::R);
  scanf("%lf", &Firm::EPS);
  scanf("%lf", &Firm::RESDOWN);
  scanf("%lf", &Firm::RESMAX);
  scanf("%lf", &Firm::NRS);
  scanf("%lf", &Firm::LOSS);
  for (int i = 0; i < NTOT; ++i) { scanf("%lf %lf %lf", &F[i].AMAN[0], &F[i].AMAN[1], &F[i].AMAN[2]); }
  for (int i = 0; i < NTOT; ++i) { scanf("%lf", &F[i].BIG); }
  for (int i = 0; i < NTOT; ++i) { scanf("%lf", &F[i].BW); }
  for (int i = 0; i < NTOT; ++i) { scanf("%lf", &F[i].DP); }
  for (int i = 0; i < NTOT; ++i) { scanf("%lf", &F[i].DS); }
  for (int i = 0; i < NTOT; ++i) { scanf("%lf", &F[i].DW); }
  for (int i = 0; i < NTOT; ++i) { scanf("%lf", &F[i].EXPDP); }
  for (int i = 0; i < NTOT; ++i) { scanf("%lf", &F[i].EXPDS); }
  for (int i = 0; i < NTOT; ++i) { scanf("%lf", &F[i].EXPDW); } 
  for (int i = 0; i < NTOT; ++i) { scanf("%lf", &F[i].EXPIDP); }
  for (int i = 0; i < NTOT; ++i) { scanf("%lf", &F[i].EXPIDS); }
  for (int i = 0; i < NTOT; ++i) { scanf("%lf", &F[i].EXPIDW); }
  for (int i = 0; i < NTOT; ++i) { scanf("%i",  &F[i].iM); }
  for (int i = 0; i < NTOT; ++i) { scanf("%lf", &F[i].INVEFF); }
  for (int i = 0; i < NTOT; ++i) { scanf("%lf", &F[i].K1); }
  for (int i = 0; i < NTOT; ++i) { scanf("%lf", &F[i].K2); }
  for (int i = 0; i < NTOT; ++i) { scanf("%lf", &F[i].L); }
  for (int i = 0; i < NTOT; ++i) { scanf("%lf", &F[i].M); }
  for (int i = 0; i < NTOT; ++i) { scanf("%lf", &F[i].MHIST); }
  for (int i = 0; i < NTOT; ++i) { scanf("%lf", &F[i].P); }
  for (int i = 0; i < NTOT; ++i) { scanf("%lf", &F[i].Q); }
  for (int i = 0; i < NTOT; ++i) { scanf("%lf", &F[i].QINV); }
  for (int i = 0; i < NTOT; ++i) { scanf("%lf", &F[i].QINVLAG); }
  for (int i = 0; i < NTOT; ++i) { scanf("%lf", &F[i].QP); }
  for (int i = 0; i < NTOT; ++i) { scanf("%lf", &F[i].QQ); }
  for (int i = 0; i < NTOT; ++i) { scanf("%lf", &F[i].QS); }
  for (int i = 0; i < NTOT; ++i) { scanf("%lf", &F[i].QTOP); }
  for (int i = 0; i < NTOT; ++i) { scanf("%lf", &F[i].QW); }
  for (int i = 0; i < NTOT; ++i) { scanf("%lf", &F[i].RES); }
  for (int i = 0; i < NTOT; ++i) { scanf("%lf", &F[i].RW); }
  for (int i = 0; i < NTOT; ++i) { scanf("%lf", &F[i].S); }
  for (int i = 0; i < NTOT; ++i) { scanf("%lf", &F[i].SMALL); }
  for (int i = 0; i < NTOT; ++i) { scanf("%lf", &F[i].STO); }
  for (int i = 0; i < NTOT; ++i) { scanf("%lf", &F[i].TEC); }
  for (int i = 0; i < NTOT; ++i) { scanf("%lf", &F[i].W); }
  for (int i = 0; i < NTOT; ++i) { scanf("%lf", &F[i].X); }

  // Government
  scanf("%lf", &G.LG);
  scanf("%lf", &G.QWG);
  scanf("%lf", &G.REALCHLG);

  // Household
  for (int j = 0; j < NMH; ++j) { scanf("%lf", &Household::ALFA1[j]); }
  for (int j = 0; j < NMH; ++j) { scanf("%lf", &Household::ALFA2[j]); }
  scanf("%lf", &Household::ALFA3);
  scanf("%lf", &Household::ALFA4);
  for (int j = 0; j < NMH; ++j) { scanf("%lf", &Household::BETA1[j]); }
  for (int j = 0; j < NMH; ++j) { scanf("%lf", &Household::BETA2[j]); }
  for (int j = 0; j < NMH; ++j) { scanf("%lf", &Household::BETA3[j]); }
  for (int j = 0; j < NMH; ++j) { scanf("%lf", &Household::SMOOTH[j]); }
  scanf("%lf", &Household::RHODUR);
  for (int j = 0; j < NMH; ++j) { scanf("%lf", &H.CVA[j]); }
  for (int j = 0; j < NMH; ++j) { scanf("%lf", &H.QC[j]); }
  scanf("%lf", &H.NH);
  scanf("%lf", &H.QCPI);
  scanf("%lf", &H.QDCPI);
  scanf("%lf", &H.STODUR);
  scanf("%lf", &H.WH);
  scanf("%lf", &H.WHRA);

  scanf("%lf", &L.ENTRY);
  scanf("%lf", &L.GAMMA);
  scanf("%lf", &L.LU);
  scanf("%lf", &L.RET);
  scanf("%lf", &L.RU);
  scanf("%lf", &L.THETA);

  // Markets
  for (int i = 0; i < NM; ++i) { scanf("%lf", &M[i].BETA); }
  for (int i = 0; i < NM; ++i) { scanf("%lf", &M[i].EXPXDP); }
  for (int i = 0; i < NM; ++i) { scanf("%lf", &M[i].EXPXDS); }
  for (int i = 0; i < NM; ++i) { scanf("%lf", &M[i].EXPXDW); }
  for (int i = 0; i < NM; ++i) { scanf("%lf", &M[i].IMP); }
  for (int i = 0; i < NM; ++i) { scanf("%lf", &M[i].MTEC); }
  for (int i = 0; i < NM; ++i) { scanf("%lf", &M[i].QDMTEC); }
  for (int i = 0; i < NM; ++i) { scanf("%lf", &M[i].QDPDOM); }
  for (int i = 0; i < NM; ++i) { scanf("%lf", &M[i].QDPFOR); }
  for (int i = 0; i < NM; ++i) { scanf("%lf", &M[i].QPDOM); }
  for (int i = 0; i < NM; ++i) { scanf("%lf", &M[i].QPFOR); }
  for (int i = 0; i < NM; ++i) { scanf("%lf", &M[i].QPH); }
  for (int i = 0; i < NM; ++i) { scanf("%lf", &M[i].TMIMP); }
  for (int i = 0; i < NM; ++i) { scanf("%lf", &M[i].TMX); }

  // Service
  scanf("%lf", &Z.LZ);
  scanf("%lf", &Z.QDTECZ);
  scanf("%lf", &Z.QMZ);
  scanf("%lf", &Z.QPH);
  scanf("%lf", &Z.QPZ);
  scanf("%lf", &Z.QSZ);
  scanf("%lf", &Z.QTARGMZ);
  scanf("%lf", &Z.QWZ);
  scanf("%lf", &Z.TECZ);

  // Remove excess data
  for (int i = 0; i < NTOT; ++i) {
    if (F[i].iM == 1 || F[i].iM == 2) { F[i] = F[--NTOT]; --i; }
  }
  for (int i = 0; i < NTOT; ++i) {
    F[i].iM -= 3;
  }

  // Count durables and non durables firms
  for (int i = 0; i < NTOT; ++i) {
    F[i].iF = i;
    if (F[i].iM == DUR) ++NTOTDUR;
    else ++NTOTNDUR;
  }

  // Initialise missing variables
  for (int I = 0; I < NTOT; ++I) {
    F[I].MINSTO = F[I].SMALL * (4 * (F[I].QS / F[I].QP));
    F[I].MAXSTO = F[I].BIG   * (4 * (F[I].QS / F[I].QP));
    F[I].OPTSTO = F[I].MINSTO + M[F[I].iM].BETA * (F[I].MAXSTO - F[I].MINSTO);
  }
}

inline void PrintYearlyResult(Firm F[MAXNTOT], Government& G, Household& H, Labour& L, Market M[MAXNM], Service& Z, const int kYear) {
  printf("---------------------------------------- Yearly Result ----------------------------------------\n");
  printf("YEAR: %i\n",   kYear);
  printf("\n");
  printf("NTOT:   %i\n",   NTOT);
  printf("(DUR):  %i\n",  NTOTDUR);
  printf("(NDUR): %i\n", NTOTNDUR);
  printf("\n");
  printf("RI: %f\n",     RI);

  printf("\n------------------------------ Labour ------------------------------\n");
  printf("LF: %f\n",   L.LF);
  printf("LU: %f\n",   L.LU);
  printf("RU: %f\n",   L.RU);
  printf("CHRU: %f\n", L.CHRU);

  printf("\n------------------------------ Service ------------------------------\n");
  printf("LZ: %f\n",    Z.LZ);
  printf("PT: %f\n",    Z.PT);
  printf("\n");
  printf("QTBUY: %f\n", Z.QTBUY);
  printf("QPH: %f\n",   Z.QPH);
  printf("QQZ: %f\n",   Z.QQZ);
  printf("QMZ: %f\n",   Z.QMZ);
  printf("QSZ: %f\n",   Z.QSZ);

  printf("\n------------------------------ Government ------------------------------\n");
  printf("LG: %f\n",  G.LG);
  printf("QWG: %f\n", G.QWG);

  printf("\n------------------------------ Household ------------------------------\n");
  printf("QDI: %lf\n",        H.QDI);
  printf("STODUR: %lf\n",     H.STODUR);
  printf("WH: %lf\n",         H.WH);
  printf("\n");
  printf("QSP(DUR): %lf\n",   H.QSP[DUR]);
  printf("QSP(NDUR): %lf\n",  H.QSP[NDUR]);
  printf("QSP(SAV): %lf\n",   H.QSP[SAV]);
  printf("QSP(Z): %lf\n",     H.QSP[iZ]);

  printf("\n------------------------------ Market ------------------------------\n");
  printf("QTSP: ");  for (int i = 0; i < NM; ++i) { printf("%f          ", M[i].QTSP); }
  printf("\n");
  printf("QTBUY: "); for (int i = 0; i < NM; ++i) { printf("%f          ", M[i].QTBUY); }
  printf("\n");
  printf("QPH: ");   for (int i = 0; i < NM; ++i) { printf("%f          ", M[i].QPH); }
  printf("\n");
  printf("QPDOM: "); for (int i = 0; i < NM; ++i) { printf("%f          ", M[i].QPDOM); }
  printf("\n");
  printf("QPFOR: "); for (int i = 0; i < NM; ++i) { printf("%f          ", M[i].QPFOR); }
  printf("\n");
  printf("IMP: ");   for (int i = 0; i < NM; ++i) { printf("%f          ", M[i].IMP); }

  printf("\n\n\n\n");
}



int main() {

  // Guessed data: NRS, QDPFOR, EXPIDP, EXPIDW, EXPIDS, (all of SERVICE except Z.QPH), ALFA1, ALFA2

  Firm F[MAXNTOT];
  Market M[MAXNM];
  Household H;
  Labour L;
  Service Z;
  Government G;

  INPUT(F, M, H, L, Z, G);

  for (int i = 0; i < YEARS; ++i) {
    if (NTOT < 2 || !NTOTDUR || !NTOTNDUR) break;
    YEAR(F, M, H, L, Z, G);
    PrintYearlyResult(F, G, H, L, M, Z, i);
  }

  return 0;
}
