/*****************************************************************************
  Photosynthesis.c

  Leaf photosynthesis, temperature kinetics and leaf energy balance for
  DHSVM-MtnHydro.  Rebuilt in Phase II from the sources listed in
  photosynthesis.h and docs/provenance.md.

  This file has no DHSVM dependencies.  It includes math.h, stdlib.h and its
  own header, and nothing else.  Keeping it that way is what lets the Phase 0
  harness link it directly; `make lint` in harness/ checks it.

  Layout:
    1  units and vapour pressure
    2  temperature response functions
    3  kinetics assembly
    4  FvCB:  A(Ci)
    5  FvCB:  A(Gc)   -- analytic, net-basis
    6  SOX support: Ci at co-limitation, dA/dCi
    7  Medlyn coupled solve
    8  leaf energy balance
    9  two-leaf canopy scaling
   10  canopy conductance
*****************************************************************************/

#include <math.h>
#include <stdlib.h>
#include "photosynthesis.h"


/* ========================================================================= */
/* 1  Units and vapour pressure                                              */
/* ========================================================================= */

float PhotoMolarToVelocity(float Tair, float Press)
{
  /* Ideal gas: 1 mol/m2/s of conductance corresponds to RT/P m/s. */
  if (Press < PHOTO_TINY)
    Press = (float)PHOTO_P0;
  return (float)(PHOTO_RGAS * (PHOTO_TFRZ + (double)Tair) / (double)Press);
}

/* Tetens/Murray form, kPa.  Standard; used only for the energy balance and
   for the harness, not inside the carbon kernel. */
float PhotoSatVaporPressure(float Tc)
{
  return 0.61078f * (float)exp(17.269 * (double)Tc / (237.3 + (double)Tc));
}

float PhotoSatVaporSlope(float Tc)
{
  float es = PhotoSatVaporPressure(Tc);
  return 4098.0f * es / ((Tc + 237.3f) * (Tc + 237.3f));
}

/* ========================================================================= */
/* 2  Temperature response functions                                         */
/* ========================================================================= */

/* Arrhenius, normalized to 1 at 25 degC.  [BER] */
static float PhotoArrhenius(float Tleaf, float Ea)
{
  double Tk = PHOTO_TFRZ + (double)Tleaf;
  return (float)exp((double)Ea / PHOTO_RGAS * (1.0 / 298.15 - 1.0 / Tk));
}

/* Peaked Arrhenius, normalized to 1 at 25 degC.  [MED02] functional form,
   [KUM19] parameters.

     f(T) = exp[Ea/R (1/298.15 - 1/T)]
            * (1 + exp((298.15 dS - Hd)/(R 298.15)))
            / (1 + exp((T dS - Hd)/(R T)))

   This is term-for-term the plantecophys TVcmax/TJmax expression (V1/V2 and
   J1/J2 in its vignette), reimplemented from [MED02] rather than transcribed:
   plantecophys is GPL and its code is used here only as a numerical oracle. */
static float PhotoPeaked(float Tleaf, float Ea, float Hd, float Ds)
{
  double Tk = PHOTO_TFRZ + (double)Tleaf;
  double Num, Den;

  Num = 1.0 + exp((298.15 * (double)Ds - (double)Hd) / (PHOTO_RGAS * 298.15));
  Den = 1.0 + exp((Tk * (double)Ds - (double)Hd) / (PHOTO_RGAS * Tk));

  if (Den < PHOTO_TINY)
    return 0.0f;

  return PhotoArrhenius(Tleaf, Ea) * (float)(Num / Den);
}

/* ========================================================================= */
/* 3  Kinetics assembly                                                      */
/* ========================================================================= */

void PhotoTraitDefaults(PHOTOTRAIT *P, float Vcmax25)
{
  if (P == NULL) return;
  P->Vcmax25   = Vcmax25;
  P->JmaxRatio = (float)PHOTO_JMAXRATIO;
  P->Rd25Ratio = (float)PHOTO_RD25RATIO;
  P->LeafWidth = (float)PHOTO_LEAF_WIDTH;
}

void PhotoKineticsTrait(const PHOTOTRAIT *P, float Tleaf, float ParAbs,
                        float Dormancy, float Beta, float Press, PHOTOKIN *K)
{
  double Aq, Bq, Cq, Disc;
  float O2, Vcmax25, JmaxRatio, Rd25Ratio, Stress;

  if (K == NULL || P == NULL)
    return;
  if (Press < PHOTO_TINY)
    Press = (float)PHOTO_P0;
  if (Dormancy < 0.0f) Dormancy = 0.0f;
  if (Dormancy > 1.0f) Dormancy = 1.0f;
  if (Beta < 0.0f) Beta = 0.0f;
  if (Beta > 1.0f) Beta = 1.0f;
  if (ParAbs < 0.0f) ParAbs = 0.0f;

  Vcmax25   = (P->Vcmax25   > 0.0f) ? P->Vcmax25   : 0.0f;
  JmaxRatio = (P->JmaxRatio > 0.0f) ? P->JmaxRatio : (float)PHOTO_JMAXRATIO;
  Rd25Ratio = (P->Rd25Ratio > 0.0f) ? P->Rd25Ratio : (float)PHOTO_RD25RATIO;

  K->Press = Press;

  /* --- Rubisco kinetics, [BER] ----------------------------------------- */
  K->Kc        = (float)PHOTO_KC25        * PhotoArrhenius(Tleaf, (float)PHOTO_EA_KC);
  K->Ko        = (float)PHOTO_KO25        * PhotoArrhenius(Tleaf, (float)PHOTO_EA_KO);
  K->GammaStar = (float)PHOTO_GAMMASTAR25 * PhotoArrhenius(Tleaf, (float)PHOTO_EA_GAMMASTAR);

  O2 = (float)PHOTO_O2_MOLFRAC * Press;          /* Pa */
  K->Km = K->Kc * (1.0f + O2 / (K->Ko + (float)PHOTO_TINY));

  /* Mole-fraction equivalents, for the Ci-in-umol/mol formulation. */
  K->KmMol    = K->Km        / Press * 1.0e6f;
  K->GStarMol = K->GammaStar / Press * 1.0e6f;

  /* --- Water stress ------------------------------------------------------
     PHOTO_BETA_MODE 0: beta attenuates Vcmax (and, through it, the export
     limit Ae), exactly as CLM4.5 BTRAN and CLM5-PHS [KEN19] Eqn 2.  Jmax and
     Rd are not stressed, as in CLM.  Mode 1 leaves the kinetics alone and the
     Medlyn solve applies beta to its slope [DEK15]. */
  Stress = (PHOTO_BETA_MODE == 0) ? Beta : 1.0f;

  /* --- Capacities, [MED02] form with [KUM19] values ---------------------- */
  K->Vcmax = Vcmax25 * Dormancy * Stress *
             PhotoPeaked(Tleaf, (float)PHOTO_EA_VCMAX,
                         (float)PHOTO_HD_VCMAX, (float)PHOTO_DS_VCMAX);

  K->Jmax  = Vcmax25 * Dormancy * JmaxRatio *
             PhotoPeaked(Tleaf, (float)PHOTO_EA_JMAX,
                         (float)PHOTO_HD_JMAX, (float)PHOTO_DS_JMAX);

  /* Dark respiration is NOT scaled by dormancy or stress: a dormant leaf
     still respires.  [COL91] fraction of Vcmax25, [BER] temperature
     response. */
  K->Rd = Vcmax25 * Rd25Ratio *
          PhotoArrhenius(Tleaf, (float)PHOTO_EA_RD);

  /* Export/product limitation, [COL91]; Ci-independent. */
  K->Ae = (float)PHOTO_KEXPORT * K->Vcmax;

  /* --- Electron transport: non-rectangular hyperbola --------------------- */
  /*   theta J^2 - (alpha I + Jmax) J + alpha I Jmax = 0, smaller root       */
  if (K->Jmax <= 0.0f || ParAbs <= 0.0f) {
    K->J = 0.0f;
    return;
  }

  Aq = (double)PHOTO_THETA;
  Bq = -((double)PHOTO_ALPHA * (double)ParAbs + (double)K->Jmax);
  Cq = (double)PHOTO_ALPHA * (double)ParAbs * (double)K->Jmax;

  Disc = Bq * Bq - 4.0 * Aq * Cq;
  if (Disc < 0.0) Disc = 0.0;

  K->J = (float)((-Bq - sqrt(Disc)) / (2.0 * Aq));
  if (K->J < 0.0f) K->J = 0.0f;
}

void PhotoKinetics(float Tleaf, float ParAbs, float Vcmax25, float Dormancy,
                   float Press, PHOTOKIN *K)
{
  PHOTOTRAIT P;
  PhotoTraitDefaults(&P, Vcmax25);
  PhotoKineticsTrait(&P, Tleaf, ParAbs, Dormancy, 1.0f, Press, K);
}

/* ========================================================================= */
/* 4  FvCB:  A(Ci)                                                           */
/* ========================================================================= */

/* Smoothed minimum of two rates.  [COL91] quadratic form:
     theta X^2 - (A + B) X + A B = 0,  smaller root.  */
static float PhotoColimit(float A, float B, float Theta)
{
  double Sum, Disc, X;

  if (A <= 0.0f || B <= 0.0f)
    return (A < B) ? A : B;
  if (Theta <= 0.0f || Theta >= 1.0f)
    return (A < B) ? A : B;

  Sum  = (double)A + (double)B;
  Disc = Sum * Sum - 4.0 * (double)Theta * (double)A * (double)B;
  if (Disc < 0.0) Disc = 0.0;

  X = (Sum - sqrt(Disc)) / (2.0 * (double)Theta);

  if (X > (double)A) X = (double)A;
  if (X > (double)B) X = (double)B;
  if (X < 0.0) X = 0.0;

  return (float)X;
}

/* Gross assimilation and its analytic slope dAg/dCi at the same Ci.

   Each branch W = a(Ci - G*)/(Ci + b) has
       dW/dCi = a (b + G*) / (Ci + b)^2
   and the [COL91] smoothing  theta X^2 - (Ac+Aj) X + Ac Aj = 0  gives, by
   implicit differentiation of F(X, Ac, Aj) = 0,
       dX/dAc = (X - Aj) / (2 theta X - Ac - Aj)
       dX/dAj = (X - Ac) / (2 theta X - Ac - Aj)
   so  dAg/dCi = dX/dAc dAc/dCi + dX/dAj dAj/dCi.

   Used to give PhotoAssimilationAtGc() a Newton step, and available to the
   stomatal schemes that want a tangent rather than a chord. */
void PhotoAgAndSlope(const PHOTOKIN *K, float Ci,
                     float *AgOut, float *SlopeOut)
{
  float Ac, Aj, Ag, dAc, dAj, Denom, dXdAc, dXdAj, X1, dX1;

  Ag = PhotoAssimilationAtCi(K, Ci, &Ac, &Aj);
  if (AgOut != NULL) *AgOut = Ag;
  if (SlopeOut == NULL) return;
  *SlopeOut = 0.0f;

  if (Ci < 0.0f) Ci = 0.0f;

  dAc = (K->Vcmax > 0.0f && Ci > K->GStarMol)
        ? K->Vcmax * (K->KmMol + K->GStarMol) /
          ((Ci + K->KmMol) * (Ci + K->KmMol) + (float)PHOTO_TINY)
        : 0.0f;

  dAj = (K->J > 0.0f && Ci > K->GStarMol)
        ? (0.25f * K->J) * (2.0f * K->GStarMol + K->GStarMol) /
          ((Ci + 2.0f * K->GStarMol) * (Ci + 2.0f * K->GStarMol)
           + (float)PHOTO_TINY)
        : 0.0f;

  /* Stage 1: X1 = colimit(Ac, Aj) */
  X1 = PhotoColimit(Ac, Aj, (float)PHOTO_COLIM1);

  if (Ac <= 0.0f || Aj <= 0.0f) {
    dX1 = (Ac < Aj) ? dAc : dAj;
  }
  else {
    Denom = 2.0f * (float)PHOTO_COLIM1 * X1 - Ac - Aj;
    if ((float)fabs((double)Denom) < 1.0e-8f) {
      dX1 = (dAc < dAj) ? dAc : dAj;
    }
    else {
      dXdAc = (X1 - Aj) / Denom;
      dXdAj = (X1 - Ac) / Denom;
      dX1 = dXdAc * dAc + dXdAj * dAj;
    }
  }

  /* Stage 2: Ag = colimit(X1, Ae).  Ae does not depend on Ci, so only the
     dAg/dX1 factor survives. */
  if (K->Ae <= 0.0f || X1 <= 0.0f) {
    *SlopeOut = dX1;
  }
  else {
    Denom = 2.0f * (float)PHOTO_COLIM2 * Ag - X1 - K->Ae;
    if ((float)fabs((double)Denom) < 1.0e-8f)
      *SlopeOut = (X1 < K->Ae) ? dX1 : 0.0f;
    else
      *SlopeOut = ((Ag - K->Ae) / Denom) * dX1;
  }

  if (*SlopeOut < 0.0f) *SlopeOut = 0.0f;
}

float PhotoAssimilationAtCi(const PHOTOKIN *K, float Ci,
                            float *AcOut, float *AjOut)
{
  float Ac, Aj, Ag;

  if (AcOut != NULL) *AcOut = 0.0f;
  if (AjOut != NULL) *AjOut = 0.0f;
  if (K == NULL) return 0.0f;

  if (Ci < 0.0f) Ci = 0.0f;

  /* Rubisco-limited, [FvCB]:  Ac = Vcmax (Ci - G*) / (Ci + Km) */
  Ac = (K->Vcmax > 0.0f)
       ? K->Vcmax * (Ci - K->GStarMol) / (Ci + K->KmMol + (float)PHOTO_TINY)
       : 0.0f;
  if (Ac < 0.0f) Ac = 0.0f;

  /* RuBP-limited, [FvCB]:  Aj = J (Ci - G*) / (4 Ci + 8 G*) */
  Aj = (K->J > 0.0f)
       ? K->J * (Ci - K->GStarMol) / (4.0f * Ci + 8.0f * K->GStarMol + (float)PHOTO_TINY)
       : 0.0f;
  if (Aj < 0.0f) Aj = 0.0f;

  /* Two-stage [COL91] smoothing: (Ac, Aj) then that with the export limit. */
  Ag = PhotoColimit(Ac, Aj, (float)PHOTO_COLIM1);
  if (K->Ae > 0.0f)
    Ag = PhotoColimit(Ag, K->Ae, (float)PHOTO_COLIM2);

  if (AcOut != NULL) *AcOut = Ac;
  if (AjOut != NULL) *AjOut = Aj;

  return Ag;
}

/* ========================================================================= */
/* 5  FvCB:  A(Gc)  -- analytic, net-basis                                   */
/* ========================================================================= */

/* Solve one limiting branch of the form
       W(Ci) = a (Ci - G*) / (Ci + b)
   simultaneously with the NET diffusion balance
       An = Gc (Ca - Ci) = W(Ci) - Rd
   which rearranges to the quadratic
       Gc Ci^2 + [(a - Rd) - Gc (Ca - b)] Ci - [Gc Ca b + a G* + b Rd] = 0
   The constant term is negative for all physical inputs, so the discriminant
   is positive and the physically meaningful root is the one with +sqrt.

   Returns the branch's net assimilation; CiOut receives its Ci. */
static float PhotoBranchAtGc(float a, float b, float GStar, float Rd,
                             float Gc, float Ca, float *CiOut)
{
  double Aq, Bq, Cq, Disc, Ci;

  if (CiOut != NULL) *CiOut = Ca;

  if (a <= 0.0 || Gc <= PHOTO_TINY)
    return -Rd;

  Aq = (double)Gc;
  Bq = ((double)a - (double)Rd) - (double)Gc * ((double)Ca - (double)b);
  Cq = -((double)Gc * (double)Ca * (double)b
         + (double)a * (double)GStar
         + (double)b * (double)Rd);

  Disc = Bq * Bq - 4.0 * Aq * Cq;
  if (Disc < 0.0) Disc = 0.0;

  Ci = (-Bq + sqrt(Disc)) / (2.0 * Aq);

  if (Ci < 0.0) Ci = 0.0;
  if (Ci > (double)Ca) Ci = (double)Ca;

  if (CiOut != NULL) *CiOut = (float)Ci;

  return (float)((double)Gc * ((double)Ca - Ci));
}

float PhotoAssimilationAtGc(const PHOTOKIN *K, float Gc, float Ca,
                            float *CiOut)
{
  float AnC, AnJ, AgC, AgJ, Ag, An, Ci;

  if (CiOut != NULL) *CiOut = Ca;
  if (K == NULL) return 0.0f;

  if (Gc <= (float)PHOTO_TINY || (K->Vcmax <= 0.0f && K->J <= 0.0f)) {
    if (CiOut != NULL) *CiOut = Ca;
    return -K->Rd;
  }

  /* Each branch solved at its own Ci, then the GROSS rates are co-limited and
     Ci is recovered from the net balance.  Co-limiting the rates rather than
     the Ci values is the standard treatment (CLM5, and [ELL20] Notes S2 do
     the same); it is an approximation only where the two branches are within
     the smoothing width of each other. */
  AnC = PhotoBranchAtGc(K->Vcmax, K->KmMol, K->GStarMol, K->Rd, Gc, Ca, NULL);

  if (K->J > 0.0f) {
    /* Aj = J(Ci - G*)/(4Ci + 8G*) = (J/4)(Ci - G*)/(Ci + 2G*)  */
    AnJ = PhotoBranchAtGc(0.25f * K->J, 2.0f * K->GStarMol, K->GStarMol,
                          K->Rd, Gc, Ca, NULL);
  }
  else {
    AnJ = -K->Rd;
  }

  AgC = AnC + K->Rd;
  AgJ = AnJ + K->Rd;
  if (AgC < 0.0f) AgC = 0.0f;
  if (AgJ < 0.0f) AgJ = 0.0f;

  Ag = PhotoColimit(AgC, AgJ, (float)PHOTO_COLIM1);
  if (K->Ae > 0.0f)
    Ag = PhotoColimit(Ag, K->Ae, (float)PHOTO_COLIM2);
  An = Ag - K->Rd;
  Ci = Ca - An / Gc;

  /* The two branch quadratics were each solved at their OWN Ci, so smoothing
     their rates gives a pair (An, Ci) that does not sit on the co-limited
     A(Ci) curve.  [ELL20] Notes S2-S3 accept that inconsistency; here it can
     reach 4 umol/m2/s at low conductance, which matters because ProfitMax
     evaluates this function tens of times per leaf and SOX differentiates it.

     The seed above is refined with Newton steps on
         f(Ci) = Gc (Ca - Ci) - [Ag(Ci) - Rd]
         f'(Ci) = -Gc - dAg/dCi
     using the analytic slope.  Convergence is quadratic and 4 steps are
     ample; each step costs one A(Ci) evaluation, so the whole call stays far
     cheaper than a bracketed solve.  DEPARTURE: this refinement is ours, not
     [ELL20]'s -- but it moves the answer ONTO the published A(Ci) curve
     rather than away from it, and P1 in the harness measures the residual. */
  {
    int it;
    for (it = 0; it < PHOTO_AGC_NEWTON; it++) {
      float AgL, Slope, Fx, Fp, Step;

      PhotoAgAndSlope(K, Ci, &AgL, &Slope);
      Fx = Gc * (Ca - Ci) - (AgL - K->Rd);
      Fp = -Gc - Slope;
      if ((float)fabs((double)Fp) < (float)PHOTO_TINY) break;

      Step = Fx / Fp;
      Ci -= Step;

      if (Ci < 0.0f) Ci = 0.0f;
      if (Ci > 2.0f * Ca) Ci = 2.0f * Ca;
      if ((float)fabs((double)Step) < (float)PHOTO_CI_TOL * Ca) break;
    }

    Ag = PhotoAssimilationAtCi(K, Ci, NULL, NULL);
    An = Ag - K->Rd;
  }

  /* Close the diffusion balance exactly.  An = Gc (Ca - Ci) must hold to
     machine precision because the stomatal schemes rely on it.  Note Ci > Ca
     is CORRECT for a respiring leaf: the CO2 flux runs outward. */
  Ci = Ca - An / Gc;
  if (Ci < 0.0f) {
    Ci = 0.0f;
    An = Gc * Ca;
  }
  else if (Ci > 2.0f * Ca) {
    Ci = 2.0f * Ca;
    An = Gc * (Ca - Ci);
  }

  if (CiOut != NULL) *CiOut = Ci;

  return An;
}

/* =========================================================================== */
/* 6  SOX support: Ci at co-limitation, and dA/dCi                           */
/* =========================================================================== */

float PhotoCiColimit(const PHOTOKIN *K, float Ca)
{
  float Wcol, Denom, CiCol;

  if (K == NULL) return Ca;

  /* [ELL20] Eqn S2.5: Wcol is the co-limited rate of the branches that do
     NOT depend strongly on Ci -- the export limit Ae and the light limit.
     In FvCB the light branch asymptotes to J/4 as Ci -> inf, which is the
     FvCB counterpart of Collatz's ci-independent Wl. */
  Wcol = PhotoColimit(K->Ae, 0.25f * K->J, (float)PHOTO_COLIM2);

  if (K->Vcmax <= 0.0f || Wcol <= 0.0f)
    return Ca;

  /* J/4 >= Vcmax: the branches never cross, the leaf is Rubisco-limited at
     every Ci.  Return the cap so the chord becomes a local gradient. */
  Denom = K->Vcmax - Wcol;
  if (Denom <= (float)PHOTO_TINY)
    return (float)PHOTO_CICOL_CAP * Ca;

  CiCol = (K->Vcmax * K->GStarMol + K->KmMol * Wcol) / Denom;

  /* ci,col ABOVE Ca is normal and must not be clamped: in bright light the
     Rubisco and light branches cross beyond ambient CO2, and [ELL20] Eqn S2.1
     is then a chord over [Ca, ci,col] rather than [ci,col, Ca].  Clamping to
     Ca collapses the chord and drives dA/dCi to zero in full sun, which is
     the opposite of the intended behaviour.
     The upper bound is a numerical guard only: as Vcmax - J/4 -> 0 the
     crossing runs to infinity, and a chord over [Ca, PHOTO_CICOL_CAP*Ca] is
     still a sensible local gradient.  DEPARTURE: the cap is ours. */
  if (CiCol < 0.0f) CiCol = 0.0f;
  if (CiCol > (float)PHOTO_CICOL_CAP * Ca)
    CiCol = (float)PHOTO_CICOL_CAP * Ca;

  return CiCol;
}

/* [ELL20] Eqn S2.1:
       dA/dci = [A(ca) - A(ci,col)] / (ca - ci,col)
   evaluated on NET assimilation (Rd cancels in the difference, so gross or
   net give the same gradient -- net is used for consistency with the rest of
   this file).

   Returns 0 when ci,col has reached Ca, i.e. the leaf gains nothing from
   opening further.  [ELL20] Eqn S2.2 makes this the signal for SOX to stop at
   gs(ci,col), which is how SOX avoids the unregulated-aperture problem
   Buckley (2017) identified in Wolf et al. (2016) and Sperry et al. (2017).
   That is the same failure mode HYD_ABSOLUTE_GAIN was added to patch; see
   docs/provenance.md sec. 3. */
float PhotoDAdCi(const PHOTOKIN *K, float Ca, float *CiColOut, float *AColOut)
{
  float CiCol, ACol, ACa, Denom;

  if (CiColOut != NULL) *CiColOut = Ca;
  if (AColOut  != NULL) *AColOut  = 0.0f;
  if (K == NULL) return 0.0f;

  CiCol = PhotoCiColimit(K, Ca);
  ACol  = PhotoAssimilationAtCi(K, CiCol, NULL, NULL) - K->Rd;
  ACa   = PhotoAssimilationAtCi(K, Ca,    NULL, NULL) - K->Rd;

  if (CiColOut != NULL) *CiColOut = CiCol;
  if (AColOut  != NULL) *AColOut  = ACol;

  Denom = Ca - CiCol;
  if ((float)fabs((double)Denom) <= (float)PHOTO_TINY)
    return 0.0f;

  /* Works for ci,col on either side of Ca: both the numerator and the
     denominator change sign together, so the chord stays positive. */
  return (ACa - ACol) / Denom;
}

/* ========================================================================= */
/* 7  Medlyn coupled solve                                                   */
/* ========================================================================= */

/* Residual of the coupled system at a trial Ci.
       Ag  = FvCB(Ci);  An = Ag - Rd
       Cs  = Ca - 1.37 An / Gb
       gs  = g0 + Beta * 1.6 (1 + g1/sqrt(D)) An / Cs      [MED11]
       Ci' = Cs - 1.6 An / gs
   Returns Ci' - Ci. */
static float PhotoCiResidual(float Ci, const PHOTOKIN *K, float G1, float G0,
                             float Vpd, float Ca, float Gb, float Beta,
                             float *AnOut, float *GsOut)
{
  float Ag, An, Cs, Gs, CiNew, Slope;

  Ag = PhotoAssimilationAtCi(K, Ci, NULL, NULL);
  An = Ag - K->Rd;

  Cs = Ca - (float)PHOTO_H2O_CO2_BL * An / (Gb + (float)PHOTO_TINY);
  if (Cs < (float)PHOTO_TINY) Cs = (float)PHOTO_TINY;

  if (Vpd < (float)PHOTO_MIN_VPD) Vpd = (float)PHOTO_MIN_VPD;

  /* [MED11].  With PHOTO_BETA_MODE 1 beta multiplies the slope term only,
     leaving g0 intact (the supply-side form of [DEK15]).  With mode 0 the
     stress is already in the kinetics ([KEN19] Eqn 2) and Beta is unused
     here. */
  Slope = ((PHOTO_BETA_MODE == 1) ? Beta : 1.0f) * (float)PHOTO_H2O_CO2_STOM *
          (1.0f + G1 / (float)sqrt((double)Vpd));

  if (An > 0.0f)
    Gs = G0 + Slope * An / Cs;
  else
    Gs = G0;          /* respiring leaf: the optimality form does not apply */

  if (Gs < (float)PHOTO_MIN_GS) Gs = (float)PHOTO_MIN_GS;
  if (Gs > (float)PHOTO_MAX_GS) Gs = (float)PHOTO_MAX_GS;

  CiNew = Cs - (float)PHOTO_H2O_CO2_STOM * An / Gs;

  /* Bound the returned Ci.  Near the compensation point An goes negative
     while gs sits on its floor, so the raw CiNew blows up to ~1e6 and the
     residual becomes wildly asymmetric across the bracket -- which starves
     any false-position method and forces it to degenerate to bisection.
     Clamping does not move the root (there CiNew == Ci, well inside the
     range); it only tames the residual away from it. */
  if (CiNew < 0.0f)          CiNew = 0.0f;
  else if (CiNew > 2.0f * Ca) CiNew = 2.0f * Ca;

  if (AnOut != NULL) *AnOut = An;
  if (GsOut != NULL) *GsOut = Gs;

  return CiNew - Ci;
}

void PhotoLeafFlux(float Vcmax25, float G1, float G0, float ParAbs,
                   float Tleaf, float Vpd, float Ca, float Press, float Gb,
                   float Beta, float Dormancy,
                   float *An, float *Gs, float *Ci)
{
  PHOTOTRAIT P;
  PhotoTraitDefaults(&P, Vcmax25);
  PhotoLeafFluxTrait(&P, G1, G0, ParAbs, Tleaf, Vpd, Ca, Press, Gb, Beta,
                     Dormancy, An, Gs, Ci);
}

void PhotoLeafFluxTrait(const PHOTOTRAIT *P, float G1, float G0, float ParAbs,
                        float Tleaf, float Vpd, float Ca, float Press,
                        float Gb, float Beta, float Dormancy,
                        float *An, float *Gs, float *Ci)
{
  PHOTOKIN K;
  float Lo, Hi, Mid, FLo, FHi, FMid;
  float AnLoc = 0.0f, GsLoc = (float)PHOTO_MIN_GS;
  int i;

  if (An != NULL) *An = 0.0f;
  if (Gs != NULL) *Gs = (float)PHOTO_MIN_GS;
  if (Ci != NULL) *Ci = Ca;
  if (P == NULL) return;

  if (Beta < 0.0f) Beta = 0.0f;
  if (Beta > 1.0f) Beta = 1.0f;

  PhotoKineticsTrait(P, Tleaf, ParAbs, Dormancy, Beta, Press, &K);

  /* No capacity at all: the leaf respires at closed stomata. */
  if (K.Vcmax <= 0.0f && K.J <= 0.0f) {
    if (An != NULL) *An = -K.Rd;
    if (Gs != NULL) *Gs = (G0 > (float)PHOTO_MIN_GS) ? G0 : (float)PHOTO_MIN_GS;
    if (Ci != NULL) *Ci = Ca;
    return;
  }

  /* Bracket Ci on (G*, Ca].  The residual is monotone decreasing in Ci over
     this interval for the coupled system, so a bracketed method is safe.

     Illinois (modified false position) rather than plain bisection: it keeps
     the bracket, so it cannot diverge on the kink the hard-minimum
     co-limitation puts in A(Ci), but converges superlinearly instead of
     linearly.  Plain bisection needed ~23 iterations for this tolerance and
     made PhotoLeafFlux nearly twice as expensive as the pre-Phase-II solver;
     Illinois gets there in 6-8. */
  Lo = K.GStarMol + 1.0e-4f;
  Hi = Ca;

  FLo = PhotoCiResidual(Lo, &K, G1, G0, Vpd, Ca, Gb, Beta, &AnLoc, &GsLoc);
  FHi = PhotoCiResidual(Hi, &K, G1, G0, Vpd, Ca, Gb, Beta, NULL, NULL);

  if (FLo <= 0.0f) {
    /* Even at the compensation point the diffusive supply cannot keep up:
       stomata are at their floor. */
    Mid = Lo;
    PhotoCiResidual(Mid, &K, G1, G0, Vpd, Ca, Gb, Beta, &AnLoc, &GsLoc);
  }
  else if (FHi >= 0.0f) {
    Mid = Hi;
    PhotoCiResidual(Mid, &K, G1, G0, Vpd, Ca, Gb, Beta, &AnLoc, &GsLoc);
  }
  else {
    int Side = 0;
    Mid = Lo;
    for (i = 0; i < PHOTO_MAX_ITER; i++) {
      float Denom = FLo - FHi;
      if ((float)fabs((double)Denom) < (float)PHOTO_TINY) break;

      Mid = (Lo * (-FHi) + Hi * FLo) / Denom;   /* false position */
      if (Mid <= Lo || Mid >= Hi)               /* guard against stalling */
        Mid = 0.5f * (Lo + Hi);

      FMid = PhotoCiResidual(Mid, &K, G1, G0, Vpd, Ca, Gb, Beta,
                             &AnLoc, &GsLoc);

      if (FMid > 0.0f) {
        Lo = Mid; FLo = FMid;
        if (Side == +1) FHi *= 0.5f;            /* Illinois relaxation */
        Side = +1;
      }
      else {
        Hi = Mid; FHi = FMid;
        if (Side == -1) FLo *= 0.5f;
        Side = -1;
      }

      if ((Hi - Lo) < (float)PHOTO_CI_TOL * Ca ||
          (float)fabs((double)FMid) < (float)PHOTO_CI_TOL * Ca)
        break;
    }
    /* AnLoc/GsLoc already correspond to the last Mid evaluated. */
  }

  if (An != NULL) *An = AnLoc;
  if (Gs != NULL) *Gs = GsLoc;
  if (Ci != NULL) *Ci = Mid;
}

/* ========================================================================= */
/* 8  Leaf energy balance                                                    */
/* ========================================================================= */

/* [CN98] Table 7.6:  gHa = 0.135 sqrt(u/d),  d = 0.72 * leaf width,
   with the 1.4 factor for outdoor turbulence.  mol/m2/s. */
float PhotoBoundaryLayerH(float Wind, float LeafWidth)
{
  float d;

  if (Wind < (float)PHOTO_WIND_MIN) Wind = (float)PHOTO_WIND_MIN;
  if (LeafWidth <= 0.0f) LeafWidth = (float)PHOTO_LEAF_WIDTH;

  d = (float)PHOTO_D_FACTOR * LeafWidth;

  return (float)PHOTO_BL_TURB * (float)PHOTO_BL_COEF *
         (float)sqrt((double)(Wind / d));
}

/* Absorbed radiation of one leaf, W/m2 leaf.  Departure D3.

   Total absorbed shortwave is the absorbed PAR (which the two-leaf partition
   already gives per unit leaf area) plus the NIR a leaf absorbs alongside it,
   PHOTO_NIR_PER_PAR of the PAR energy.  The longwave term is the isothermal
   one -- what the leaf would absorb from surroundings at air temperature --
   so that Rabs - eps sigma Tl^4 reduces to the [CN98] ch. 14 isothermal net
   radiation SW - eps sigma (Tl^4 - Ta^4).

   Never pass the canopy net radiation here: it is per m2 GROUND and already
   net of longwave emission, and the leaf balance subtracts emission itself. */
float PhotoLeafRabs(float ParAbsLeaf, float Tair)
{
  float Tk = Tair + (float)PHOTO_TFRZ;
  float Sw;

  if (ParAbsLeaf < 0.0f) ParAbsLeaf = 0.0f;
  Sw = ParAbsLeaf / (float)PHOTO_PAR_PER_WATT * (1.0f + (float)PHOTO_NIR_PER_PAR);

  return Sw + (float)PHOTO_EMISSIVITY * (float)PHOTO_STEFAN * Tk * Tk * Tk * Tk;
}

/* Residual of the leaf energy budget (W/m2 leaf):

     Rabs - eps sigma Tl^4 - N cp gHa (Tl - Ta) - lambda E = 0

   with Rabs from PhotoLeafRabs(), gHa the [CN98] Table 7.6 one-sided flat
   plate conductance and N = PHOTO_LEAF_SIDES the number of faces that
   exchange sensible heat.  Solved rather than linearized, so the answer does
   not depend on recalling a particular textbook rearrangement, and so the
   balance can be checked afterwards by evaluating this function at the
   returned Tleaf.  Iterating the leaf energy balance is what [ELL20] and
   Sperry et al. (2017) do. */
float PhotoEnergyResidualW(float Tleaf, float E, float Tair, float Rabs,
                           float Wind, float Press, float LeafWidth)
{
  float gHa, Lambda, Tk, LwOut, Sensible, Latent;

  (void)Press;

  gHa = PhotoBoundaryLayerH(Wind, LeafWidth);

  /* Latent heat of vaporization, J/mol (18.015 g/mol) */
  Lambda = (2.501e6f - 2370.0f * Tleaf) * 0.018015f;

  Tk = Tleaf + (float)PHOTO_TFRZ;
  LwOut    = (float)PHOTO_EMISSIVITY * (float)PHOTO_STEFAN * Tk * Tk * Tk * Tk;
  Sensible = (float)PHOTO_LEAF_SIDES * (float)PHOTO_CP_MOLAR * gHa * (Tleaf - Tair);
  Latent   = Lambda * E;

  return Rabs - LwOut - Sensible - Latent;
}

float PhotoEnergyResidual(float Tleaf, float E, float Tair, float Vpd,
                          float Rabs, float Wind, float Press)
{
  (void)Vpd;
  return PhotoEnergyResidualW(Tleaf, E, Tair, Rabs, Wind, Press,
                              (float)PHOTO_LEAF_WIDTH);
}

void PhotoLeafEnergyBalanceW(float E, float Tair, float Vpd, float Rabs,
                             float Wind, float Press, float LeafWidth,
                             float *Tleaf, float *VpdLeaf)
{
  float Lo, Hi, Mid, FLo, FMid, es, ea;
  int i;

  if (Tleaf   != NULL) *Tleaf   = Tair;
  if (VpdLeaf != NULL) *VpdLeaf = Vpd;

  if (E < 0.0f) E = 0.0f;
  if (LeafWidth <= 0.0f) LeafWidth = (float)PHOTO_LEAF_WIDTH;

  Lo = Tair - (float)PHOTO_TLEAF_SPAN;
  Hi = Tair + (float)PHOTO_TLEAF_SPAN;

  FLo  = PhotoEnergyResidualW(Lo, E, Tair, Rabs, Wind, Press, LeafWidth);
  FMid = PhotoEnergyResidualW(Hi, E, Tair, Rabs, Wind, Press, LeafWidth);

  /* The residual falls monotonically with Tleaf (both the longwave and the
     sensible terms increase), so a sign change inside the bracket is
     guaranteed unless the forcing is extreme; if it is, clamp. */
  if (FLo <= 0.0f) {
    Mid = Lo;
  }
  else if (FMid >= 0.0f) {
    Mid = Hi;
  }
  else {
    for (i = 0; i < PHOTO_MAX_ITER; i++) {
      Mid = 0.5f * (Lo + Hi);
      FMid = PhotoEnergyResidualW(Mid, E, Tair, Rabs, Wind, Press, LeafWidth);
      if (FMid > 0.0f) Lo = Mid;
      else             Hi = Mid;
      if ((Hi - Lo) < (float)PHOTO_TLEAF_TOL)
        break;
    }
    Mid = 0.5f * (Lo + Hi);
  }

  /* Leaf-to-air VPD: saturation at leaf temperature minus actual air vapour
     pressure.  The air vapour pressure is recovered from the AIR VPD, so the
     leaf temperature effect is applied exactly once. */
  es = PhotoSatVaporPressure(Mid);
  ea = PhotoSatVaporPressure(Tair) - Vpd;
  if (ea < 0.0f) ea = 0.0f;

  if (Tleaf != NULL) *Tleaf = Mid;
  if (VpdLeaf != NULL) {
    float D = es - ea;
    *VpdLeaf = (D < (float)PHOTO_MIN_VPD) ? (float)PHOTO_MIN_VPD : D;
  }
}

void PhotoLeafEnergyBalance(float E, float Tair, float Vpd, float Rabs,
                            float Wind, float Press,
                            float *Tleaf, float *VpdLeaf)
{
  PhotoLeafEnergyBalanceW(E, Tair, Vpd, Rabs, Wind, Press,
                          (float)PHOTO_LEAF_WIDTH, Tleaf, VpdLeaf);
}

/* ========================================================================= */
/* 9  Two-leaf canopy scaling  -- [DPF97]                                    */
/* ========================================================================= */

void PhotoTwoLeafPartition(float Lai, float SinAlt, float ParBeam,
                           float ParDiff, float Vcmax25,
                           float *LaiSun, float *LaiSha,
                           float *ParSun, float *ParSha,
                           float *VcmaxSun, float *VcmaxSha)
{
  float Kb, Kn, LSun, LSha, FDiffSun;
  float VcTot, VcSunInt, VcShaInt;
  float ParSunTot, ParShaTot;

  Kn = (float)PHOTO_KN_NITROGEN;

  if (Lai < (float)PHOTO_TINY) {
    *LaiSun = 0.0f;   *LaiSha = 0.0f;
    *ParSun = 0.0f;   *ParSha = 0.0f;
    *VcmaxSun = 0.0f; *VcmaxSha = 0.0f;
    return;
  }

  /* Canopy carboxylation capacity under an exponential nitrogen profile. */
  VcTot = Vcmax25 * (1.0f - (float)exp(-Kn * Lai)) / Kn;

  /* Night: no beam, all leaf area behaves as shaded. */
  if (SinAlt <= 0.0f) {
    LSun = 0.0f;
    LSha = Lai;

    *LaiSun = LSun;
    *LaiSha = LSha;
    *ParSun = 0.0f;
    *ParSha = (ParDiff > 0.0f)
              ? (ParDiff * (float)PHOTO_PAR_PER_WATT / LSha) : 0.0f;
    *VcmaxSun = 0.0f;
    *VcmaxSha = VcTot / LSha;
    return;
  }

  if (SinAlt < (float)PHOTO_MIN_SINALT)
    SinAlt = (float)PHOTO_MIN_SINALT;

  Kb = (float)PHOTO_G_SPHERICAL / SinAlt;

  LSun = (1.0f - (float)exp(-Kb * Lai)) / Kb;
  if (LSun > Lai) LSun = Lai;
  LSha = Lai - LSun;

  /* Diffuse PAR is absorbed with extinction kd through the canopy, and the
     sunlit fraction at depth L is exp(-kb L).  The sunlit share of the
     absorbed diffuse flux is therefore  [DPF97]

         fDiffSun = kd/(kd+kb) * (1 - exp(-(kd+kb) L)) / (1 - exp(-kd L)),

     larger than LSun/Lai because sunlit leaves sit where diffuse light is
     strongest.  (Departure D4: the pre-refactor code used LSun/Lai.) */
  {
    float Kd = (float)PHOTO_KD_DIFFUSE;
    float DenD = 1.0f - (float)exp(-Kd * Lai);
    if (DenD > (float)PHOTO_TINY)
      FDiffSun = Kd / (Kd + Kb) * (1.0f - (float)exp(-(Kd + Kb) * Lai)) / DenD;
    else
      FDiffSun = (Lai > 0.0f) ? LSun / Lai : 0.0f;
    if (FDiffSun < 0.0f) FDiffSun = 0.0f;
    if (FDiffSun > 1.0f) FDiffSun = 1.0f;
  }

  /* Absorbed PAR, W/m2 ground -> umol/m2 leaf/s (departure D2). */
  if (LSun > (float)PHOTO_TINY) {
    ParSunTot = ParBeam + ParDiff * FDiffSun;
    *ParSun = ParSunTot * (float)PHOTO_PAR_PER_WATT / LSun;
  }
  else {
    *ParSun = 0.0f;
  }

  if (LSha > (float)PHOTO_TINY) {
    ParShaTot = ParDiff * (1.0f - FDiffSun);
    *ParSha = ParShaTot * (float)PHOTO_PAR_PER_WATT / LSha;
  }
  else {
    *ParSha = 0.0f;
  }

  VcSunInt = Vcmax25 * (1.0f - (float)exp(-(Kn + Kb) * Lai)) / (Kn + Kb);
  VcShaInt = VcTot - VcSunInt;
  if (VcShaInt < 0.0f) VcShaInt = 0.0f;

  *VcmaxSun = (LSun > (float)PHOTO_TINY) ? (VcSunInt / LSun) : 0.0f;
  *VcmaxSha = (LSha > (float)PHOTO_TINY) ? (VcShaInt / LSha) : 0.0f;

  *LaiSun = LSun;
  *LaiSha = LSha;
}

/* ========================================================================= */
/* 10  Canopy conductance                                                    */
/* ========================================================================= */

float PhotoCanopyConductance(float Vcmax25, float G1, float G0, float Lai,
                             float SinAlt, float ParBeam, float ParDiff,
                             float Tair, float Vpd, float Ca, float Press,
                             float Gb, float Beta, float Dormancy,
                             float *AnCanopy)
{
  float LaiSun, LaiSha, ParSun, ParSha, VcmaxSun, VcmaxSha;
  float AnSun = 0.0f, GsSun = 0.0f, CiSun = 0.0f;
  float AnSha = 0.0f, GsSha = 0.0f, CiSha = 0.0f;
  float GsCanopyMol, AnTot;

  if (Lai < (float)PHOTO_TINY) {
    if (AnCanopy != NULL) *AnCanopy = 0.0f;
    return 0.0f;
  }

  PhotoTwoLeafPartition(Lai, SinAlt, ParBeam, ParDiff, Vcmax25,
    &LaiSun, &LaiSha, &ParSun, &ParSha, &VcmaxSun, &VcmaxSha);

  /* Each big leaf keeps the class water-use strategy but carries its own
     photosynthetic capacity from the canopy nitrogen profile. */
  if (LaiSun > (float)PHOTO_TINY)
    PhotoLeafFlux(VcmaxSun, G1, G0, ParSun, Tair, Vpd, Ca, Press, Gb,
                  Beta, Dormancy, &AnSun, &GsSun, &CiSun);

  if (LaiSha > (float)PHOTO_TINY)
    PhotoLeafFlux(VcmaxSha, G1, G0, ParSha, Tair, Vpd, Ca, Press, Gb,
                  Beta, Dormancy, &AnSha, &GsSha, &CiSha);

  /* Conductances add in parallel across leaf area. */
  GsCanopyMol = GsSun * LaiSun + GsSha * LaiSha;
  AnTot       = AnSun * LaiSun + AnSha * LaiSha;

  if (AnCanopy != NULL) *AnCanopy = AnTot;

  return GsCanopyMol * PhotoMolarToVelocity(Tair, Press);
}
