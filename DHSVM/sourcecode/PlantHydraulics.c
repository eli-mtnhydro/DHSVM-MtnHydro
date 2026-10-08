/*****************************************************************************
  PlantHydraulics.c

  Xylem hydraulics for DHSVM-MtnHydro: vulnerability curves, the Kirchhoff
  supply function from root collar to canopy, and the hydraulic failure
  threshold.

  GRAVITY (Phase-VI refactor).  The hydrostatic drop rho g h between the
  root collar and the canopy is applied at the entry to the xylem stage,
  Pxylem,0 = Hcollar - Pgrav, so the same convention serves the supply
  function used by ProfitMax/Medlyn-PHS and the SOX predawn potential
  ([ELL20] Notes S1.2, Psi_pd = Psi_r - h g rho 1e-6).  Before this the SOX
  path carried gravity and the Kirchhoff path did not.

  WHAT LEFT THIS FILE IN PHASE III
  -------------------------------------------------------------------------
  * The duplicate FvCB kernel (HydArrhenius, HydPeaked, HydKinetics,
    HydAssimilationAtGc).  Deleted.  Photosynthesis.c is now the single
    source, and its PhotoAssimilationAtGc() closes the diffusion balance on
    NET assimilation, which the deleted version did not.
  * The leaf energy balance.  Moved to Photosynthesis.c in Phase II.
  * The entire layered rhizosphere (HydSetRhizosphere, HydRhizosphereK,
    HydLayerK, HydLayerUptake, HydSupplyFunctionLayered) and the ROOTZONE
    struct that went with it.  Replaced by RootHydraulics.c.

  The three departures that machinery carried are retired, not patched:
  emergent layered uptake is now Couvreur/Vanderborght; the root radial
  conductance that suppressed winner-take-all is unnecessary because a
  SUF-weighted mean cannot be dominated by one layer; and the soil-plant
  conductance is monotone again, so the hydraulic cost can go back to
  Sperry's definition instead of the xylem-only substitute.

  SOURCES
  -------------------------------------------------------------------------
  [SPE17] Sperry JS et al. (2017) Plant Cell Environ 40:816-830.
          doi:10.1111/pce.12852.  Supply function, Ecrit, profit maximization.
  [SL15]  Sperry JS, Love DM (2015) New Phytol 207:14-27.
          The Kirchhoff transform used to integrate the supply function.
  [ELL20] Eller CB et al. (2020) New Phytol 226:1622-1637. Eqn 2: the
          sigmoidal vulnerability curve K = 1/[1 + (psi/psi50)^a], also used
          by CLM5-PHS (Kennedy et al. 2019).
  [LEI25] Leitner D et al. (2025) HESS 29:1759-1782.  Eqn 20 gives the
          analytic collar potential this file consumes.

  STRUCTURE OF THE SUPPLY FUNCTION
  -------------------------------------------------------------------------
  With the soil-to-collar stage analytic ([LEI25] Eqn 20),

      Hcollar(E) = Heff - E / Krs

  the xylem stage is a clean one-dimensional Kirchhoff integral ([SL15]):

      F(psi) = integral of k(psi') dpsi' from psi to 0
      E      = F(psi_leaf) - F(Hcollar)

  so for a given E the leaf potential follows EXPLICITLY:

      psi_leaf = F^-1( E + F(Heff - E/Krs) )

  No iteration, no nested solve, and the curve is monotone by construction
  because k > 0.  That is the whole reason for splitting the root stage out.

  No DHSVM dependencies.
*****************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "planthydraulics.h"

/* ========================================================================= */
/* 1  Vulnerability curves                                                   */
/* ========================================================================= */

/* Weibull, [SPE17]:  k/kmax = exp(-(|psi|/b)^c).  P50 = b (ln2)^(1/c). */
float HydWeibull(float Psi, float B, float C)
{
  double u;
  if (B <= 0.0f || C <= 0.0f) return 1.0f;
  if (Psi >= 0.0f) return 1.0f;
  u = fabs((double)Psi) / (double)B;
  return (float)exp(-pow(u, (double)C));
}

/* Sigmoidal, [ELL20] Eqn 2 / CLM5-PHS:  k/kmax = 1/[1 + (psi/psi50)^a].
   Cheaper than Weibull and its parameters are the measured traits. */
float HydSigmoidal(float Psi, float Psi50, float A)
{
  double u;
  if (Psi50 >= 0.0f || A <= 0.0f) return 1.0f;
  if (Psi >= 0.0f) return 1.0f;
  u = pow(fabs((double)Psi) / fabs((double)Psi50), (double)A);
  return (float)(1.0 / (1.0 + u));
}

float HydVulnerability(float Psi, const HYDCURVE *V)
{
  if (V == NULL) return 1.0f;
  if (V->Form == HYD_WEIBULL)
    return HydWeibull(Psi, V->P1, V->P2);
  return HydSigmoidal(Psi, V->P1, V->P2);
}

float HydP50ToWeibullB(float P50, float C)
{
  if (P50 >= 0.0f || C <= 0.0f) return 2.0f;
  return (float)(fabs((double)P50) / pow(log(2.0), 1.0 / (double)C));
}

void HydCurveWeibullFromP50(HYDCURVE *V, float P50, float C)
{
  if (V == NULL) return;
  V->Form = HYD_WEIBULL;
  V->P1 = HydP50ToWeibullB(P50, C);
  V->P2 = C;
}

void HydCurveSigmoidal(HYDCURVE *V, float P50, float A)
{
  if (V == NULL) return;
  V->Form = HYD_SIGMOIDAL;
  V->P1 = P50;
  V->P2 = A;
}

/* Pressure at which conductance falls to the form's critical fraction of
   maximum (HYD_KCRIT_FRAC for Weibull, HYD_KCRIT_FRAC_SIGMOIDAL for the
   long-tailed sigmoidal).  Bisection on a monotone function: bounded, and
   evaluated once per class at startup, not per timestep. */
float HydCriticalPressure(const HYDCURVE *V)
{
  float Lo = 0.0f, Hi = -200.0f, Mid = 0.0f;
  float Frac = (V != NULL && V->Form == HYD_SIGMOIDAL)
               ? (float)HYD_KCRIT_FRAC_SIGMOIDAL : (float)HYD_KCRIT_FRAC;
  int i;

  for (i = 0; i < 80; i++) {
    Mid = 0.5f * (Lo + Hi);
    if (HydVulnerability(Mid, V) > Frac) Lo = Mid;
    else                                 Hi = Mid;
    if ((float)fabs((double)(Hi - Lo)) < 1.0e-5f) break;
  }
  return 0.5f * (Lo + Hi);
}

/* ========================================================================= */
/* 2  Kirchhoff cumulant                                                     */
/* ========================================================================= */

/* F(psi) = integral of k(psi') dpsi' from psi up to 0, tabulated on a
   uniform grid from 0 down to Pcrit and integrated with the trapezoid rule.
   [SL15].  Depends only on the vulnerability curve and kmax, so it is built
   once per vegetation class at startup, NOT per pixel per timestep -- which
   is the other half of why Phase III is cheap. */
void HydBuildCumulant(HYDXYLEM *X, const HYDCURVE *V, float KxMax)
{
  float dP, kPrev, kCur;
  int i;

  if (X == NULL || V == NULL) return;
  memset(X, 0, sizeof(*X));

  X->Curve = *V;
  X->KxMax = (KxMax > 0.0f) ? KxMax : 1.0f;
  X->Pcrit = HydCriticalPressure(V);
  X->N = HYD_NCUMULANT;

  dP = X->Pcrit / (float)(X->N - 1);      /* negative step                  */
  X->dP = dP;

  X->P[0] = 0.0f;
  X->F[0] = 0.0f;
  kPrev = X->KxMax * HydVulnerability(0.0f, V);

  for (i = 1; i < X->N; i++) {
    X->P[i] = dP * (float)i;
    kCur = X->KxMax * HydVulnerability(X->P[i], V);
    X->F[i] = X->F[i - 1] + 0.5f * (kPrev + kCur) * (float)fabs((double)dP);
    kPrev = kCur;
  }

  X->Fcrit = X->F[X->N - 1];
}

/* Hydrostatic drop for the class's canopy height (MPa). */
void HydSetCanopyHeight(HYDXYLEM *X, float HeightM)
{
  if (X == NULL) return;
  X->Pgrav = (HeightM > 0.0f) ? HeightM * 9810.0f * 1.0e-6f : 0.0f;
}

/* F at an arbitrary pressure, by linear interpolation on the grid. */
float HydCumulant(const HYDXYLEM *X, float Psi)
{
  float t;
  int k;

  if (X == NULL || X->N < 2) return 0.0f;
  if (Psi >= 0.0f) return 0.0f;
  if (Psi <= X->Pcrit) return X->Fcrit;

  t = Psi / X->dP;                        /* both negative -> positive index */
  k = (int)t;
  if (k < 0) k = 0;
  if (k > X->N - 2) k = X->N - 2;
  t -= (float)k;

  return X->F[k] + t * (X->F[k + 1] - X->F[k]);
}

/* Inverse: the pressure at which the cumulant reaches Ftarget.  Monotone,
   so a direct index search plus one linear interpolation suffices. */
float HydCumulantInverse(const HYDXYLEM *X, float Ftarget)
{
  int Lo, Hi, Mid;
  float t, dF;

  if (X == NULL || X->N < 2) return 0.0f;
  if (Ftarget <= 0.0f) return 0.0f;
  if (Ftarget >= X->Fcrit) return X->Pcrit;

  Lo = 0; Hi = X->N - 1;
  while (Hi - Lo > 1) {
    Mid = (Lo + Hi) / 2;
    if (X->F[Mid] <= Ftarget) Lo = Mid;
    else                      Hi = Mid;
  }

  dF = X->F[Hi] - X->F[Lo];
  t = ((float)fabs((double)dF) > HYD_TINY) ? (Ftarget - X->F[Lo]) / dF : 0.0f;

  return X->P[Lo] + t * (X->P[Hi] - X->P[Lo]);
}

/* ========================================================================= */
/* 3  Supply function                                                        */
/* ========================================================================= */

/* Leaf pressure sustaining a transpiration E, given the effective soil
   potential and the root system conductance.  EXPLICIT -- no iteration.

       Hcollar  = Heff - E/Krs                         [LEI25] Eqn 20
       Pxylem0  = Hcollar - Pgrav                      gravity, see header
       psi_leaf = F^-1( E + F(Pxylem0) )               [SL15]

   Returns 0 and sets *Failed when E exceeds what the xylem can carry. */
float HydLeafPressure(const HYDXYLEM *X, float Heff, float Krs, float E,
                      float *CollarOut, int *Failed)
{
  float Hcollar, Pxylem0, Ftarget;

  if (Failed != NULL) *Failed = 0;
  if (CollarOut != NULL) *CollarOut = Heff;

  if (X == NULL || Krs <= HYD_TINY) {
    if (Failed != NULL) *Failed = 1;
    return 0.0f;
  }

  Hcollar = Heff - E / Krs;
  if (CollarOut != NULL) *CollarOut = Hcollar;

  Pxylem0 = Hcollar - X->Pgrav;
  if (Pxylem0 <= X->Pcrit) {
    if (Failed != NULL) *Failed = 1;
    return X->Pcrit;
  }

  Ftarget = E + HydCumulant(X, Pxylem0);

  if (Ftarget >= X->Fcrit) {
    if (Failed != NULL) *Failed = 1;
    return X->Pcrit;
  }

  return HydCumulantInverse(X, Ftarget);
}

/* Transpiration at hydraulic failure: the largest E for which a leaf
   pressure above Pcrit still exists.  The residual
       g(E) = Fcrit - E - F(Heff - E/Krs - Pgrav)
   decreases monotonically in E, so bisection is safe and needs ~30 cheap
   evaluations, once per pixel per timestep. */
float HydEcrit(const HYDXYLEM *X, float Heff, float Krs)
{
  float Lo, Hi, Mid, g, Htop;
  int i;

  if (X == NULL || Krs <= HYD_TINY) return 0.0f;
  Htop = Heff - X->Pgrav;                  /* xylem entry at zero flow       */
  if (Htop <= X->Pcrit) return 0.0f;

  Lo = 0.0f;
  Hi = Krs * (Htop - X->Pcrit);            /* xylem entry hits Pcrit here    */
  if (Hi <= 0.0f) return 0.0f;

  for (i = 0; i < 40; i++) {
    Mid = 0.5f * (Lo + Hi);
    g = X->Fcrit - Mid - HydCumulant(X, Htop - Mid / Krs);
    if (g > 0.0f) Lo = Mid;
    else          Hi = Mid;
    if ((Hi - Lo) < 1.0e-6f * (Hi + 1.0f)) break;
  }
  return 0.5f * (Lo + Hi);
}

/* Build the discretized supply curve the profit-maximization scans.
   Points are laid out in E, and each leaf pressure follows explicitly, so
   the whole curve costs N interpolations rather than N nested solves. */
void HydBuildSupply(const HYDXYLEM *X, float Heff, float Krs,
                    HYDSUPPLY *S)
{
  float Ecrit, dE;
  int i, Failed;

  if (S == NULL) return;
  memset(S, 0, sizeof(*S));

  if (X == NULL || Krs <= HYD_TINY) return;

  Ecrit = HydEcrit(X, Heff, Krs);
  S->Ecrit = Ecrit;
  S->Heff  = Heff;
  S->Krs   = Krs;
  S->Pcrit = X->Pcrit;

  if (Ecrit <= HYD_TINY) { S->N = 0; return; }

  S->N = HYD_NSUPPLY;
  dE = Ecrit / (float)(S->N - 1);

  for (i = 0; i < S->N; i++) {
    S->E[i] = dE * (float)i;
    S->Pleaf[i] = HydLeafPressure(X, Heff, Krs, S->E[i],
                                  &S->Pcollar[i], &Failed);
    S->Kc[i] = X->KxMax * HydVulnerability(S->Pleaf[i], &X->Curve);
  }

  S->Kcmax = S->Kc[0];
  S->Kcrit = S->Kc[S->N - 1];
}

/* NOTE.  HydKmaxFromVcmax() used to live here.  It has moved to
   StomatalScheme.c as StomKmaxFromVcmax(), because the Sperry coordination
   is DEFINED through the profit-maximization optimum and cannot be evaluated
   without it.  A version that avoided the optimizer was tried and was wrong:
   its test statistic was invariant under kmax, so the bisection never moved
   off its lower bracket and every class got the same (tiny) conductance. */
