#include "MOSES.h"



double QFR(Firm& F, const double L) {
  return (1.0 - F.RES) * F.QTOP * (1.0 - pow(e, -(F.TEC / F.QTOP) * L));
}

double RFQ(Firm& F, const double Q) {
  return (F.QTOP / F.TEC) * log(((1.0 - F.RES) * F.QTOP) / ((1.0 - F.RES) * F.QTOP - Q));
}

bool SAT(Firm& F, const double Q, const double L) {
  double MARGIN;
  if (L > 0) MARGIN = 1.0 - (L * (F.QEXPW / 4.0)) / (Q * F.QEXPP);
  else MARGIN = 1.0 - (F.QEXPW / 4.0) / ((1.0 - F.RES) * F.TEC * F.QEXPP);

  if (MARGIN >= F.QTARGM) return true;
  return false;
}

bool SOLVE(Firm& F) {
  double B = F.QEXPW / ((1.0 - F.QTARGM) * (1.0 - F.RES) * F.TEC * F.QEXPP * 4.0);
  if (B <= 0) return false;

  double Y = 1.0 / B;
  double h = (B * Y + pow(e, -Y) - 1.0) / (B - pow(e, -Y));

  while (h > 0.001) {
    h = (B * Y + pow(e, -Y) - 1.0) / (B - pow(e, -Y));
    Y -= h;
  }

  F.QPLANL = Y * (F.QTOP / F.TEC);
  F.QPLANQ = QFR(F, F.QPLANL);

  return true;
}

// Random number generator from seed x
unsigned XorShift32(unsigned& x) {
  x ^= x << 13;
  x ^= x >> 17;
  x ^= x << 5;
  return x;
}

int CHOOSE(Firm F[MAXNTOT], Labour& L, const int kThisFirm, bool& is_firm) {
  double SUMLABOUR = 0;
  for (int I = 0; I < NTOT; ++I) { SUMLABOUR += F[I].LL; }
  SUMLABOUR += L.LL;

  int random = XorShift32(xrandom) % ((int)SUMLABOUR + 1);
  int cum_sum = 0;
  for (int i = 0; i < NTOT; ++i) {
    cum_sum += F[i].LL;
    if (cum_sum > random) {
      if (i == kThisFirm && i == NTOT - 1) break;
      if (i == kThisFirm) ++i;
      is_firm = true;
      return i;
    }
  }

  is_firm = false;
  return NTOT;
}

void NULLIFY_FIRM(Firm F[MAXNTOT], const int firm) {
  if (F[firm].iM == DUR) --NTOTDUR;
  else --NTOTNDUR;

  F[firm] = F[--NTOT];
}
