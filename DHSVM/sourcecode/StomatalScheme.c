/*****************************************************************************
  StomatalScheme.c

  Medlyn, Sperry profit maximization, Wang ProfitMax2, and SOX (numerical and
  semi-analytical).  Sources and transcribed equations are in
  stomatalscheme.h.

  No DHSVM dependencies.
*****************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "stomatalscheme.h"
#include "roothydraulics.h"

/* mol H2O m-2 s-1  ->  mmol */
#define STOM_MOL_TO_MMOL   1000.0f

/* ========================================================================= */
/* Helpers                                                                   */
/* ========================================================================= */

void StomInit(STOMSOLUTION *S, const STOMMET *M)
{
  if (S == NULL) return;
  memset(S, 0, sizeof(*S));
  if (M != NULL) {
    S->Ci = M->Ca;
    S->Tleaf = M->Tair;
    S->VpdLeaf = M->VpdAir;
  }
  S->Failed = 1;
}

/* Departure S2: [ELL20] treat the whole soil-plant path as one compartment.
   DHSVM resolves two stages, so r_p,min is their series resistance. */
float StomRpMin(float Krs, float KxMax)
{
  float Ktot;

  if (Krs <= STOM_TINY || KxMax <= STOM_TINY) return 0.0f;
  Ktot = 1.0f / (1.0f / Krs + 1.0f / KxMax);   /* mmol/m2/s/MPa             */
  if (Ktot <= STOM_TINY) return 0.0f;

  /* mmol -> mol, and invert to a resistance */
  return STOM_MOL_TO_MMOL / Ktot;
}

/* [ELL20] Eqn 6:  r_p = r_p,min / K(Psi_pd). */
float StomSoxRp(float RpMin, float PsiPd, const HYDCURVE *V)
{
  float K = HydVulnerability(PsiPd, V);
  if (K < 1.0e-6f) K = 1.0e-6f;
  return RpMin / K;
}

/* [ELL20] Notes S2.9:
       (dK/dPsi_m)(1/K) = [K(Psi_pd) - K((Psi_pd+Psi_50)/2)]
                          / [Psi_pd - (Psi_pd+Psi_50)/2] / K(Psi_pd)

   The gradient is evaluated toward Psi_50 -- the steepest point of the
   vulnerability curve -- deliberately, so it stays positive even in the flat
   part of the curve.  That is what stops SOX from opening stomata without
   limit when the cost is near zero, the defect Buckley (2017) identified in
   Wolf et al. (2016) and Sperry et al. (2017). */
float StomSoxKGradient(float PsiPd, const HYDCURVE *V)
{
  float Psi50, PsiHalf, Kpd, Khalf, dPsi, Grad;

  if (V == NULL) return 0.0f;

  /* Psi_50 from whichever curve form is in use. */
  if (V->Form == HYD_WEIBULL)
    Psi50 = -(float)((double)V->P1 * pow(log(2.0), 1.0 / (double)V->P2));
  else
    Psi50 = V->P1;

  PsiHalf = 0.5f * (PsiPd + Psi50);
  dPsi = PsiPd - PsiHalf;
  if ((float)fabs((double)dPsi) < 1.0e-9f) return 0.0f;

  Kpd   = HydVulnerability(PsiPd,   V);
  Khalf = HydVulnerability(PsiHalf, V);
  if (Kpd < 1.0e-9f) Kpd = 1.0e-9f;

  Grad = (Kpd - Khalf) / dPsi / Kpd;
  return (Grad > 0.0f) ? Grad : 0.0f;
}

/* Leaf-to-air VPD as a mole fraction, which is what [ELL20] Eqn 5 wants. */
static float StomVpdMoleFraction(float VpdLeafKpa, float PressPa)
{
  float Pk = PressPa * 0.001f;
  if (Pk < STOM_TINY) Pk = 101.325f;
  return VpdLeafKpa / Pk;
}

/* The photosynthetic traits the kernel needs, from the scheme trait block. */
static void StomPhotoTrait(const STOMTRAIT *Tr, PHOTOTRAIT *P)
{
  PhotoTraitDefaults(P, Tr->Vcmax25);
  if (Tr->JmaxRatio > 0.0f) P->JmaxRatio = Tr->JmaxRatio;
  if (Tr->Rd25Ratio > 0.0f) P->Rd25Ratio = Tr->Rd25Ratio;
  if (Tr->LeafWidth > 0.0f) P->LeafWidth = Tr->LeafWidth;
}

/* Total leaf conductance to CO2 from the stomatal CO2 conductance and the
   boundary layer conductance to water. */
static float StomTotalGc(float GsCO2, float GbH2O)
{
  if (GsCO2 <= STOM_TINY) return 0.0f;
  if (GbH2O <= STOM_TINY) return GsCO2;
  return 1.0f / (1.0f / GsCO2 + (float)PHOTO_H2O_CO2_BL / GbH2O);
}

/* ========================================================================= */
/* Medlyn  [MED11]                                                           */
/* ========================================================================= */

void StomMedlyn(const STOMTRAIT *Tr, const STOMMET *M, float Beta,
                STOMSOLUTION *S)
{
  PHOTOTRAIT P;
  float An, GsH2O, Ci;

  StomInit(S, M);
  if (Tr == NULL || M == NULL || S == NULL) return;

  StomPhotoTrait(Tr, &P);
  PhotoLeafFluxTrait(&P, Tr->G1, Tr->G0, M->ParAbs, M->Tair, M->VpdAir,
                     M->Ca, M->Press, Tr->Gb, Beta, Tr->Dormancy,
                     &An, &GsH2O, &Ci);

  if (Tr->Gmax > 0.0f && GsH2O > Tr->Gmax) GsH2O = Tr->Gmax;

  S->Gw = GsH2O;
  S->Gs = GsH2O / (float)PHOTO_H2O_CO2_STOM;
  S->An = An;
  S->Ci = Ci;
  S->E  = GsH2O * StomVpdMoleFraction(M->VpdAir, M->Press) * STOM_MOL_TO_MMOL;
  S->Tleaf = M->Tair;
  S->VpdLeaf = M->VpdAir;
  S->Failed = 0;
}

/* ========================================================================= */
/* Profit maximization  [SPE17] and ProfitMax2 [WAN20]                       */
/* ========================================================================= */

/* K = dE/dPleaf along the supply curve.  [SPE17]'s definition, restored in
   Phase III (departure S3): with Krs/SUF the series conductance is monotone,
   so this no longer needs the xylem-only substitute. */
static void StomSupplyConductance(const HYDSUPPLY *Sup, float *K)
{
  int i;
  float dE, dP;

  if (Sup == NULL || Sup->N < 2) return;

  for (i = 0; i < Sup->N; i++) {
    int a = (i == 0) ? 0 : i - 1;
    int b = (i == Sup->N - 1) ? Sup->N - 1 : i + 1;
    dE = Sup->E[b] - Sup->E[a];
    dP = Sup->Pleaf[a] - Sup->Pleaf[b];      /* positive: Pleaf falls        */
    K[i] = (dP > 1.0e-9f) ? dE / dP : 0.0f;
  }
}

void StomProfitMax(const STOMTRAIT *Tr, const STOMMET *M,
                   const HYDSUPPLY *Sup, int Variant, STOMSOLUTION *S)
{
  static float K[HYD_MAXSUPPLY];
  static float Astore[HYD_MAXSUPPLY];
  static float Gwstore[HYD_MAXSUPPLY];
  static float TLstore[HYD_MAXSUPPLY];
  static float DLstore[HYD_MAXSUPPLY];
  static float Cistore[HYD_MAXSUPPLY];

  static int Done[HYD_MAXSUPPLY];

  PHOTOKIN Kin;
  PHOTOTRAIT P;
  float Kmax, Kcrit, Amax, Best, Obj;
  int i, iBest = 0, Stride, NEval = 0;

  StomInit(S, M);
  if (Tr == NULL || M == NULL || Sup == NULL || S == NULL) return;
  if (Sup->N < 3 || Sup->Ecrit <= STOM_TINY) return;

  StomPhotoTrait(Tr, &P);

  S->Ecrit = Sup->Ecrit;

  StomSupplyConductance(Sup, K);
  Kmax  = K[0];
  Kcrit = K[Sup->N - 1];
  if (Kmax <= Kcrit + STOM_TINY) return;

  /* --- lazy point evaluation ------------------------------------------ */
  /* Evaluating every supply point costs one leaf energy balance and one
     Newton-refined A(Gc) each -- about 1000 FvCB evaluations for a 120-point
     curve, which is why a full scan runs ~170x Medlyn.  The profit curve is
     unimodal ([SPE17] Fig. 1b), so a coarse pass locates the optimum and
     supplies Amax, and only a neighbourhood is refined.  The coarse pass
     gives Amax to within the stride, and assimilation rises monotonically
     with conductance apart from a weak leaf-temperature effect, so that is
     accurate enough for the normalization. */
  memset(Done, 0, sizeof(int) * (size_t)Sup->N);

# define STOM_EVAL(ii)                                                        \
  do {                                                                        \
    if (!Done[ii]) {                                                          \
      float Em_ = Sup->E[ii] / STOM_MOL_TO_MMOL;                              \
      float TL_, DL_, Dm_, Gw_, Gc_, An_, Ci_;                                \
      PhotoLeafEnergyBalanceW(Em_, M->Tair, M->VpdAir, M->Rabs, M->Wind,      \
                              M->Press, P.LeafWidth, &TL_, &DL_);             \
      Dm_ = StomVpdMoleFraction(DL_, M->Press);                               \
      if (Dm_ < 1.0e-6f) Dm_ = 1.0e-6f;                                       \
      Gw_ = Em_ / Dm_;                                                        \
      if (Tr->Gmax > 0.0f && Gw_ > Tr->Gmax) Gw_ = Tr->Gmax;                     \
      Gc_ = StomTotalGc(Gw_ / (float)PHOTO_H2O_CO2_STOM, Tr->Gb);              \
      PhotoKineticsTrait(&P, TL_, M->ParAbs, Tr->Dormancy, 1.0f, M->Press, &Kin); \
      An_ = PhotoAssimilationAtGc(&Kin, Gc_, M->Ca, &Ci_);                    \
      Astore[ii] = An_;  Gwstore[ii] = Gw_;                                   \
      TLstore[ii] = TL_; DLstore[ii] = DL_; Cistore[ii] = Ci_;                \
      Done[ii] = 1;  NEval++;                                                 \
    }                                                                         \
  } while (0)

  Stride = (Sup->N > 24) ? (Sup->N / 12) : 1;

  /* --- coarse pass: Amax and the neighbourhood of the optimum ---------- */
  Amax = 0.0f;
  for (i = 0; i < Sup->N; i += Stride) {
    STOM_EVAL(i);
    if (Astore[i] > Amax) Amax = Astore[i];
  }
  i = Sup->N - 1; STOM_EVAL(i);
  if (Astore[i] > Amax) Amax = Astore[i];

  if (Amax <= STOM_TINY) {
    STOM_EVAL(0);
    S->Gs = STOM_MIN_GS; S->Gw = STOM_MIN_GS * (float)PHOTO_H2O_CO2_STOM;
    S->An = Astore[0]; S->Ci = M->Ca; S->E = 0.0f;
    S->PsiLeaf = Sup->Pleaf[0]; S->PsiCollar = Sup->Pcollar[0];
    S->Tleaf = TLstore[0]; S->VpdLeaf = DLstore[0];
    S->Iterations = NEval;
    S->Failed = 0;
    return;
  }

  Best = -1.0e30f;
  for (i = 0; i < Sup->N; i += Stride) {
    Obj = (Variant == STOM_PROFITMAX2)
          ? Astore[i] * (Sup->Ecrit - Sup->E[i]) / Sup->Ecrit
          : Astore[i] / Amax - (Kmax - K[i]) / (Kmax - Kcrit);
    if (Obj > Best) { Best = Obj; iBest = i; }
  }

  /* --- refine within one stride either side ---------------------------- */
  {
    int Lo = iBest - Stride, Hi = iBest + Stride;
    if (Lo < 0) Lo = 0;
    if (Hi > Sup->N - 1) Hi = Sup->N - 1;
    Best = -1.0e30f;
    for (i = Lo; i <= Hi; i++) {
      STOM_EVAL(i);
      Obj = (Variant == STOM_PROFITMAX2)
            ? Astore[i] * (Sup->Ecrit - Sup->E[i]) / Sup->Ecrit
            : Astore[i] / Amax - (Kmax - K[i]) / (Kmax - Kcrit);
      if (Obj > Best) { Best = Obj; iBest = i; }
    }
  }
# undef STOM_EVAL

  S->Gw        = Gwstore[iBest];
  S->Gs        = S->Gw / (float)PHOTO_H2O_CO2_STOM;
  S->An        = Astore[iBest];
  S->Ci        = Cistore[iBest];
  S->E         = Sup->E[iBest];
  S->PsiLeaf   = Sup->Pleaf[iBest];
  S->PsiCollar = Sup->Pcollar[iBest];
  S->Tleaf     = TLstore[iBest];
  S->VpdLeaf   = DLstore[iBest];
  S->Objective = Best;
  S->Iterations = NEval;
  S->Failed    = 0;
}

/* ========================================================================= */
/* SOX  [ELL18] / [ELL20]                                                    */
/* ========================================================================= */

/* Psi_m and the resulting normalized conductance at a trial gs.
   [ELL20] S1.9:  Psi_m = Psi_pd - r_p 1.6 gs D / 2                          */
static float StomSoxPsiM(float PsiPd, float Rp, float GsCO2, float Dmol)
{
  return PsiPd - Rp * (float)PHOTO_H2O_CO2_STOM * GsCO2 * Dmol * 0.5f;
}

void StomSOX(const STOMTRAIT *Tr, const STOMMET *M, const HYDXYLEM *X,
             float Heff, float Krs, int Solver, STOMSOLUTION *S)
{
  PHOTOKIN Kin;
  PHOTOTRAIT P;
  float PsiPd, RpMin, Rp, Dmol, DL, TL;
  float dAdCi, CiCol, ACol, Xi, Gs, GwH2O, Gc, An, Ci;
  float PsiM, Kn, Emol;
  int i;

  StomInit(S, M);
  if (Tr == NULL || M == NULL || X == NULL || S == NULL) return;
  if (Krs <= STOM_TINY) return;

  StomPhotoTrait(Tr, &P);

  /* --- predawn potential, [ELL20] S1.2 ---------------------------------
     Psi_pd = Psi_r - h g rho 1e-6, with Psi_r the mean root-zone soil
     potential.  Departure S2: Heff (SUF-weighted) is that mean, better
     defined than an arithmetic average.  The gravity term is the xylem
     object's Pgrav, shared with the Kirchhoff supply function. */
  PsiPd = Heff - X->Pgrav;
  S->PsiPd = PsiPd;

  if (PsiPd <= X->Pcrit) return;          /* soil already past failure      */

  RpMin = (Tr->RpMin > 0.0f) ? Tr->RpMin : StomRpMin(Krs, X->KxMax);
  if (RpMin <= STOM_TINY) return;

  Rp = StomSoxRp(RpMin, PsiPd, &X->Curve);        /* Eqn 6                  */

  S->Ecrit = HydEcrit(X, Heff, Krs);

  /* Leaf temperature depends on E which depends on gs; [ELL18] iterate the
     energy balance.  Start from air temperature and refine once the optimum
     is known (loop below). */
  TL = M->Tair;
  DL = M->VpdAir;

  for (i = 0; i < 3; i++) {
    Dmol = StomVpdMoleFraction(DL, M->Press);
    if (Dmol < 1.0e-6f) Dmol = 1.0e-6f;

    PhotoKineticsTrait(&P, TL, M->ParAbs, Tr->Dormancy, 1.0f, M->Press, &Kin);

    if (Solver == STOM_SOX_ANALYTICAL) {
      /* ---------------- [ELL20] Eqns 4-5 --------------------------- */
      float KGrad = StomSoxKGradient(PsiPd, &X->Curve);   /* Notes S2.9    */
      float Denom;

      dAdCi = PhotoDAdCi(&Kin, M->Ca, &CiCol, &ACol);     /* Eqn S2.1      */

      if (dAdCi <= STOM_TINY) {
        /* Eqn S2.2: no carbon benefit from opening further. */
        Gs = (M->Ca - CiCol > 1.0f) ? ACol / (M->Ca - CiCol) : STOM_MIN_GS;
      }
      else {
        Denom = KGrad * Rp * (float)PHOTO_H2O_CO2_STOM * Dmol;
        if (Denom <= STOM_TINY) {
          Gs = (Tr->Gmax > 0.0f)
               ? Tr->Gmax / (float)PHOTO_H2O_CO2_STOM : 1.0f;
        }
        else {
          Xi = 2.0f / Denom;                              /* Eqn 5         */
          /* Eqn 4 = positive root of S1.17b */
          Gs = 0.5f * dAdCi *
               ((float)sqrt(1.0 + 4.0 * (double)Xi / (double)dAdCi) - 1.0f);
        }
      }
    }
    else {
      /* ---------------- [ELL18] numerical iteration ------------------
         Scan ci, evaluate the objective A K, take the maximum, then
         recover gs from the diffusion equation.  This is the reference
         solution: [JON22] treat it as truth. */
      float Best = -1.0e30f, BestGs = STOM_MIN_GS;
      int j;

      for (j = 1; j < STOM_NSCAN; j++) {
        float CiT = M->Ca * (float)j / (float)STOM_NSCAN;
        float AgT = PhotoAssimilationAtCi(&Kin, CiT, NULL, NULL);
        float AnT = AgT - Kin.Rd;
        float GsT, PsiMT, KT, ObjT;

        if (AnT <= 0.0f) continue;
        GsT = AnT / (M->Ca - CiT);                 /* Fick, S1.14a         */
        if (GsT <= STOM_TINY) continue;

        PsiMT = StomSoxPsiM(PsiPd, Rp, GsT, Dmol); /* S1.9                 */
        KT = HydVulnerability(PsiMT, &X->Curve);   /* Eqn 2                */
        ObjT = AnT * KT;                           /* Eqn 1                */

        if (ObjT > Best) { Best = ObjT; BestGs = GsT; }
      }
      Gs = BestGs;
    }

    if (Gs < STOM_MIN_GS) Gs = STOM_MIN_GS;
    GwH2O = Gs * (float)PHOTO_H2O_CO2_STOM;
    if (Tr->Gmax > 0.0f && GwH2O > Tr->Gmax) {
      GwH2O = Tr->Gmax;
      Gs = GwH2O / (float)PHOTO_H2O_CO2_STOM;
    }

    /* Refine the leaf energy balance at this conductance. */
    Emol = GwH2O * Dmol;
    PhotoLeafEnergyBalanceW(Emol, M->Tair, M->VpdAir, M->Rabs, M->Wind,
                            M->Press, P.LeafWidth, &TL, &DL);
    S->Iterations = i + 1;
  }

  /* The loop above ends by UPDATING the leaf temperature, so the gs it last
     computed belongs to the previous temperature.  Returning both would give
     a solution that is one step out of step with itself, and Eqn 4 would not
     satisfy its own quadratic at the reported (gs, D_leaf).  One more pass at
     the converged temperature closes it. */
  Dmol = StomVpdMoleFraction(DL, M->Press);
  if (Dmol < 1.0e-6f) Dmol = 1.0e-6f;
  PhotoKineticsTrait(&P, TL, M->ParAbs, Tr->Dormancy, 1.0f, M->Press, &Kin);

  if (Solver == STOM_SOX_ANALYTICAL) {
    float KGrad = StomSoxKGradient(PsiPd, &X->Curve);
    float Denom;
    dAdCi = PhotoDAdCi(&Kin, M->Ca, &CiCol, &ACol);
    if (dAdCi <= STOM_TINY) {
      Gs = (M->Ca - CiCol > 1.0f) ? ACol / (M->Ca - CiCol) : STOM_MIN_GS;
    }
    else {
      Denom = KGrad * Rp * (float)PHOTO_H2O_CO2_STOM * Dmol;
      if (Denom <= STOM_TINY) {
        Gs = (Tr->Gmax > 0.0f) ? Tr->Gmax / (float)PHOTO_H2O_CO2_STOM : 1.0f;
      }
      else {
        Xi = 2.0f / Denom;
        Gs = 0.5f * dAdCi *
             ((float)sqrt(1.0 + 4.0 * (double)Xi / (double)dAdCi) - 1.0f);
      }
    }
    if (Gs < STOM_MIN_GS) Gs = STOM_MIN_GS;
    GwH2O = Gs * (float)PHOTO_H2O_CO2_STOM;
    if (Tr->Gmax > 0.0f && GwH2O > Tr->Gmax) {
      GwH2O = Tr->Gmax;
      Gs = GwH2O / (float)PHOTO_H2O_CO2_STOM;
    }
  }
  Emol  = GwH2O * Dmol;
  PsiM  = StomSoxPsiM(PsiPd, Rp, Gs, Dmol);
  Kn    = HydVulnerability(PsiM, &X->Curve);

  PhotoKineticsTrait(&P, TL, M->ParAbs, Tr->Dormancy, 1.0f, M->Press, &Kin);
  Gc = StomTotalGc(Gs, Tr->Gb);
  An = PhotoAssimilationAtGc(&Kin, Gc, M->Ca, &Ci);

  S->Gs        = Gs;
  S->Gw        = GwH2O;
  S->An        = An;
  S->Ci        = Ci;
  S->E         = Emol * STOM_MOL_TO_MMOL;
  S->Tleaf     = TL;
  S->VpdLeaf   = DL;
  S->PsiLeaf   = 2.0f * PsiM - PsiPd;      /* S1.1b: Psi_c = 2 Psi_m - Psi_pd */
  S->PsiCollar = Heff - S->E / Krs;        /* [LEI25] Eqn 20                */
  S->Objective = An * Kn;
  S->Failed    = 0;
}

/* ========================================================================= */
/* Medlyn with plant hydraulic stress -- CLM5-PHS                            */
/* ========================================================================= */

/* [MED11] stomata driven by [SPE17]/[ELL20] hydraulics.  Kennedy et al.
   (2019) implement this in CLM5 as PHS: leaf water potential, not soil
   moisture, sets the stress factor, and the stress factor attenuates Vcmax
   (their Eqn 2; PHOTO_BETA_MODE 0).

   The coupling is circular -- gs sets E, E sets psi_leaf, psi_leaf sets beta,
   beta sets gs -- so it is iterated to a fixed point.  CLM5 Newton-solves the
   same system; a damped fixed point is enough here because beta is bounded in
   [0,1] and the map is a contraction over that range.

   beta = K(psi_leaf) / K(0), the fractional loss of conductance, which is the
   same attenuation the vulnerability curve already defines.  No new
   parameter, and it reduces to beta = 1 in saturated soil.

   Known approximation: E = gw D_air here (no leaf energy balance), so the
   transpiration Penman-Monteith later computes from gw differs by the
   radiation term.  CLM5 solves the two jointly. */
void StomMedlynPHS(const STOMTRAIT *Tr, const STOMMET *M, const HYDXYLEM *X,
                   float Heff, float Krs, STOMSOLUTION *S)
{
  PHOTOTRAIT P;
  float Beta = 1.0f, BetaNew, An, GwH2O, Ci, Dmol, Emol, Emmol;
  float PsiLeaf = 0.0f, PsiCollar = Heff, Ecrit;
  int Failed = 0, it;

  StomInit(S, M);
  if (Tr == NULL || M == NULL || X == NULL || S == NULL) return;
  if (Krs <= STOM_TINY) return;

  StomPhotoTrait(Tr, &P);

  Ecrit = HydEcrit(X, Heff, Krs);
  S->Ecrit = Ecrit;
  if (Ecrit <= STOM_TINY) {
    /* Soil past hydraulic failure: stomata shut, but report a valid state. */
    S->Gs = STOM_MIN_GS;
    S->Gw = STOM_MIN_GS * (float)PHOTO_H2O_CO2_STOM;
    S->PsiLeaf = X->Pcrit; S->PsiCollar = Heff;
    S->Failed = 0;
    return;
  }

  Dmol = StomVpdMoleFraction(M->VpdAir, M->Press);
  if (Dmol < 1.0e-6f) Dmol = 1.0e-6f;

  for (it = 0; it < 8; it++) {
    PhotoLeafFluxTrait(&P, Tr->G1, Tr->G0, M->ParAbs, M->Tair,
                       M->VpdAir, M->Ca, M->Press, Tr->Gb, Beta, Tr->Dormancy,
                       &An, &GwH2O, &Ci);

    if (Tr->Gmax > 0.0f && GwH2O > Tr->Gmax) GwH2O = Tr->Gmax;

    Emol  = GwH2O * Dmol;                       /* mol H2O/m2/s            */
    Emmol = Emol * STOM_MOL_TO_MMOL;

    /* The supply function cannot deliver more than Ecrit.  Clip there and
       let beta close the stomata rather than reporting an impossible flux. */
    if (Emmol > Ecrit) Emmol = Ecrit;

    PsiLeaf = HydLeafPressure(X, Heff, Krs, Emmol, &PsiCollar, &Failed);
    if (Failed) PsiLeaf = X->Pcrit;

    BetaNew = HydVulnerability(PsiLeaf, &X->Curve);
    if (BetaNew < 0.0f) BetaNew = 0.0f;
    if (BetaNew > 1.0f) BetaNew = 1.0f;

    if ((float)fabs((double)(BetaNew - Beta)) < 1.0e-4f) {
      Beta = BetaNew;
      break;
    }
    Beta = 0.5f * (Beta + BetaNew);             /* damping                  */
    S->Iterations = it + 1;
  }

  /* Final evaluation at the converged beta, so the reported gs, An and
     psi_leaf are mutually consistent rather than one step apart. */
  PhotoLeafFluxTrait(&P, Tr->G1, Tr->G0, M->ParAbs, M->Tair, M->VpdAir,
                     M->Ca, M->Press, Tr->Gb, Beta, Tr->Dormancy,
                     &An, &GwH2O, &Ci);
  if (Tr->Gmax > 0.0f && GwH2O > Tr->Gmax) GwH2O = Tr->Gmax;

  Emol  = GwH2O * Dmol;
  Emmol = Emol * STOM_MOL_TO_MMOL;
  if (Emmol > Ecrit) {
    Emmol = Ecrit;
    Emol  = Emmol / STOM_MOL_TO_MMOL;
    GwH2O = Emol / Dmol;
  }
  PsiLeaf = HydLeafPressure(X, Heff, Krs, Emmol, &PsiCollar, &Failed);

  S->Gw        = GwH2O;
  S->Gs        = GwH2O / (float)PHOTO_H2O_CO2_STOM;
  S->An        = An;
  S->Ci        = Ci;
  S->E         = Emmol;
  S->PsiLeaf   = PsiLeaf;
  S->PsiCollar = PsiCollar;
  S->Tleaf     = M->Tair;
  S->VpdLeaf   = M->VpdAir;
  S->Objective = Beta;
  S->Failed    = 0;
}

/* ========================================================================= */
/* Coordination: KxMax from Vcmax25  ([SPE17])                               */
/* ========================================================================= */

void StomReferencePoint(const PHOTOTRAIT *P, const HYDXYLEM *X, float Krs,
                        float Ca, float Press, STOMSOLUTION *S, float *Ecrit)
{
  HYDSUPPLY Sup;
  STOMTRAIT Tr;
  STOMMET M;

  if (Ecrit != NULL) *Ecrit = 0.0f;
  if (S == NULL) return;
  memset(S, 0, sizeof(*S));
  S->Failed = 1;
  if (P == NULL || X == NULL || Krs <= 0.0f) return;
  if (Ca <= 0.0f) Ca = (float)PHOTO_CA;

  /* [SPE17] reference conditions: bright, warm, moderate VPD, wet soil.  The
     leaf's absorbed radiation follows departure D3 for PAR = 2000. */
  memset(&Tr, 0, sizeof(Tr));
  Tr.Vcmax25 = P->Vcmax25; Tr.G1 = 4.0f; Tr.G0 = 0.0f;
  Tr.Dormancy = 1.0f; Tr.Gb = (float)PHOTO_GB; Tr.Gmax = 0.0f;
  Tr.JmaxRatio = P->JmaxRatio; Tr.Rd25Ratio = P->Rd25Ratio;
  Tr.LeafWidth = P->LeafWidth;

  memset(&M, 0, sizeof(M));
  M.Tair = 25.0f; M.VpdAir = 1.0f;
  M.Press = (Press > 0.0f) ? Press : (float)PHOTO_P0;
  M.Ca = Ca; M.ParAbs = 2000.0f;
  M.Rabs = PhotoLeafRabs(M.ParAbs, M.Tair); M.Wind = 2.0f;

  HydBuildSupply(X, 0.0f, Krs, &Sup);
  if (Ecrit != NULL) *Ecrit = Sup.Ecrit;
  if (Sup.N < 3 || Sup.Ecrit <= STOM_TINY) return;

  StomProfitMax(&Tr, &M, &Sup, STOM_PROFITMAX, S);
}

float StomKmaxFromVcmax(const PHOTOTRAIT *P, float P50, float VulnShape,
                        int VulnForm, float CiCaTarget, float CanopyHeight,
                        float Ca, float Press)
{
  HYDCURVE V;
  HYDXYLEM X;
  STOMSOLUTION S;
  float Lo = 0.01f, Hi = 500.0f, Mid = 1.0f, Krs, CiCa, Ecrit;
  int i;

  if (P == NULL || P->Vcmax25 <= 0.0f || P50 >= 0.0f) return -1.0f;
  if (VulnShape <= 0.0f) VulnShape = (float)HYD_WEIBULL_C_DEFAULT;
  if (CiCaTarget <= 0.0f || CiCaTarget >= 1.0f)
    CiCaTarget = (float)HYD_CICA_TARGET;
  if (Ca <= 0.0f) Ca = (float)PHOTO_CA;

  if (VulnForm == HYD_SIGMOIDAL) HydCurveSigmoidal(&V, P50, VulnShape);
  else                           HydCurveWeibullFromP50(&V, P50, VulnShape);

  /* Ci/Ca rises monotonically with kmax: more hydraulic capacity means the
     optimum sits at higher conductance and a smaller CO2 drawdown. */
  for (i = 0; i < 45; i++) {
    Mid = 0.5f * (Lo + Hi);
    HydBuildCumulant(&X, &V, Mid);
    HydSetCanopyHeight(&X, CanopyHeight);
    Krs = RootKrsFromXylem(Mid);

    StomReferencePoint(P, &X, Krs, Ca, Press, &S, &Ecrit);
    if (S.Failed || S.Gs <= STOM_MIN_GS) { Lo = Mid; continue; }

    CiCa = S.Ci / Ca;
    if (CiCa < CiCaTarget) Lo = Mid;
    else                   Hi = Mid;

    if (Hi - Lo < 1.0e-4f * Hi) break;
  }

  Mid = 0.5f * (Lo + Hi);

  /* No coordination point inside the bracket: refuse rather than return a
     bracket bound, which is exactly the silent failure this replaces. */
  if (Mid <= 1.05f * 0.01f || Mid >= 0.95f * 500.0f) return -1.0f;

  return Mid;
}
