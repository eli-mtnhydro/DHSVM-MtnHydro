
#include <math.h>
#include <stdlib.h>
#include "photosynthesis.h"

static float ArrheniusFactor(float Tleaf, float Ea);
static float PeakedFactor(float Tleaf, float Ea, float Hd, float Ds);
static float ElectronTransport(float ParAbs, float Jmax);
static float CiResidual(float Ci, float Vcmax25, float G1, float Vcmax, float Jmax,
  float Rd, float GammaStar, float Km, float Ca, float Press, float Gb,
  float Slope, float G0, float *AnOut, float *GsOut);

/*****************************************************************************
 PhotoMolarToVelocity()
 Conductance in mol/m2/s to m/s
 *****************************************************************************/
float PhotoMolarToVelocity(float Tair, float Press)
{
  float TairK;

  TairK = Tair + PHOTO_TFRZ;
  if (Press < PHOTO_TINY)
    Press = PHOTO_P0;

  return (float)(PHOTO_RGAS * TairK / Press);
}

/*****************************************************************************
 ArrheniusFactor()
 Normalized to 1.0 at 25 degC.  Used for Kc, Ko, Gamma* and Rd
 *****************************************************************************/
static float ArrheniusFactor(float Tleaf, float Ea)
{
  float TleafK, Tref;

  TleafK = Tleaf + PHOTO_TFRZ;
  Tref = 25.0f + PHOTO_TFRZ;

  return (float)exp(Ea * (TleafK - Tref) / (Tref * PHOTO_RGAS * TleafK));
}

/*****************************************************************************
 PeakedFactor()
 Arrhenius with high-temperature deactivation, normalized to 1.0 at 25 degC
 Used for Vcmax and Jmax
 *****************************************************************************/
static float PeakedFactor(float Tleaf, float Ea, float Hd, float Ds)
{
  float TleafK, Tref, Arrh, NumTerm, DenTerm;

  TleafK = Tleaf + PHOTO_TFRZ;
  Tref = 25.0f + PHOTO_TFRZ;

  Arrh = ArrheniusFactor(Tleaf, Ea);
  NumTerm = 1.0f + (float)exp((Tref * Ds - Hd) / (Tref * PHOTO_RGAS));
  DenTerm = 1.0f + (float)exp((TleafK * Ds - Hd) / (TleafK * PHOTO_RGAS));

  return Arrh * NumTerm / DenTerm;
}

/*****************************************************************************
 ElectronTransport()
 Smaller root of the non-rectangular hyperbola
   theta*J^2 - (I2 + Jmax)*J + I2*Jmax = 0
 where I2 = alpha * absorbed PAR
 *****************************************************************************/
static float ElectronTransport(float ParAbs, float Jmax)
{
  float I2, A, B, C, Disc;

  if (ParAbs <= 0.0f || Jmax <= 0.0f)
    return 0.0f;

  I2 = PHOTO_ALPHA * ParAbs;

  A = PHOTO_THETA;
  B = -(I2 + Jmax);
  C = I2 * Jmax;

  Disc = B * B - 4.0f * A * C;
  if (Disc < 0.0f)
    Disc = 0.0f;

  return (float)((-B - sqrt(Disc)) / (2.0f * A));
}

/*****************************************************************************
 CiResidual()
 Given a trial Ci (umol/mol), evaluate the biochemical demand, then the
 stomatal and diffusive supply, and return
   residual = Ci_implied_by_supply - Ci_trial
 A root of this function is the simultaneous solution of the three coupled
 relations.  The function is evaluated in partial pressures internally but
 takes and returns mole fractions, because the Medlyn and diffusion
 relations are naturally mole-fraction quantities.
 *****************************************************************************/
static float CiResidual(float Ci, float Vcmax25, float G1, float Vcmax, float Jmax,
  float Rd, float GammaStar, float Km, float Ca, float Press, float Gb,
  float Slope, float G0, float *AnOut, float *GsOut)
{
  float CiPa, GammaStarMf;
  float Wc, Wj, We, Ag, An;
  float Cs, Gs, CiNew;

  /* Partial pressure of the trial Ci (Pa) */
  CiPa = Ci * 1.0e-6f * Press;

  /* --- Biochemical demand (Farquhar et al. 1980) ------------------------- */
  /* Rubisco-limited */
  Wc = Vcmax * (CiPa - GammaStar) / (CiPa + Km);

  /* RuBP-regeneration (light) limited.  ElectronTransport() has already
     folded absorbed PAR into Jmax-limited J upstream; here Jmax carries J. */
  Wj = Jmax * (CiPa - GammaStar) / (4.0f * (CiPa + 2.0f * GammaStar));

  /* Product/export limited.  Noah-MP's WE = 0.5*Vcmax for C3. */
  We = 0.5f * Vcmax;

  Ag = Wc;
  if (Wj < Ag)
    Ag = Wj;
  if (We < Ag)
    Ag = We;
  if (Ag < 0.0f)
    Ag = 0.0f;

  An = Ag - Rd;

  /* --- Leaf surface CO2 (umol/mol) --------------------------------------- */
  Cs = Ca - PHOTO_H2O_CO2_BL * An / (Gb + PHOTO_TINY);
  if (Cs < PHOTO_TINY)
    Cs = PHOTO_TINY;

  /* --- Medlyn stomatal conductance (mol H2O/m2/s) ------------------------ */
  /* Slope has already absorbed 1.6*(1 + G1/sqrt(D)) and the Beta factor.
     When An <= 0 (respiring leaf: night, or deep shade) the optimality
     relation has no meaning and conductance collapses to the residual term. */
  if (An > 0.0f)
    Gs = G0 + Slope * An / Cs;
  else
    Gs = G0;

  if (Gs < PHOTO_MIN_GS)
    Gs = PHOTO_MIN_GS;
  if (Gs > PHOTO_MAX_GS)
    Gs = PHOTO_MAX_GS;

  /* --- Diffusive supply: Ci implied by this An and Gs -------------------- */
  CiNew = Cs - PHOTO_H2O_CO2_STOM * An / Gs;

  GammaStarMf = GammaStar / (1.0e-6f * Press);
  if (CiNew < 0.0f)
    CiNew = 0.0f;
  if (CiNew > Ca && An > 0.0f)
    CiNew = Ca;

  if (AnOut != NULL)
    *AnOut = An;
  if (GsOut != NULL)
    *GsOut = Gs;

  /* GammaStarMf is computed for clarity/debugging of the bracket bounds */
  (void)GammaStarMf;

  return CiNew - Ci;
}

/*****************************************************************************
 PhotoLeafFlux()
 Solve the coupled leaf system at fixed leaf temperature
*****************************************************************************/
void PhotoLeafFlux(float Vcmax25, float G1, float ParAbs, float Tleaf, float Vpd,
  float Ca, float Press, float Gb, float Beta, float Dormancy,
  float *An, float *Gs, float *Ci)
{
  float Vcmax, Jmax, J, Rd;
  float GammaStar, Kc, Ko, Km, O2Pa;
  float Slope, G0;
  float CiLo, CiHi, CiMid;
  float FLo, FHi;
  float AnLoc, GsLoc;
  float A, B, C, D, E, Fa, Fb, Fc, P, Q, R, S, Tol1, Xm;
  int Iter;

  AnLoc = 0.0f;
  GsLoc = 0.0f;

  /* --- Guard rails on the inputs ----------------------------------------- */
  if (Press < PHOTO_TINY)
    Press = PHOTO_P0;
  if (Gb < PHOTO_TINY)
    Gb = PHOTO_TINY;
  if (Beta < 0.0f)
    Beta = 0.0f;
  if (Beta > 1.0f)
    Beta = 1.0f;
  if (ParAbs < 0.0f)
    ParAbs = 0.0f;

  /* Medlyn is undefined at D = 0 and the 1/sqrt(D) term becomes very large at
     small D.  Floor at 0.05 kPa, which is the low end of the range the model
     was fitted over (De Kauwe et al. 2015 Fig. 1). */
  if (Vpd < 0.05f)
    Vpd = 0.05f;

  G0 = PHOTO_G0 * Dormancy; /* Temporary approximation heuristic for soil freezing */
  if (G0 < PHOTO_MIN_GS)
    G0 = PHOTO_MIN_GS;

  /* --- Temperature-adjusted kinetics ------------------------------------- */
  Kc = PHOTO_KC25 * ArrheniusFactor(Tleaf, PHOTO_EA_KC);
  Ko = PHOTO_KO25 * ArrheniusFactor(Tleaf, PHOTO_EA_KO);
  GammaStar = PHOTO_GAMMASTAR25 * ArrheniusFactor(Tleaf, PHOTO_EA_GAMMASTAR);

  /* O2 partial pressure scales with total pressure: this is the elevation
     effect that the mole-fraction formulation would lose. */
  O2Pa = PHOTO_O2_MOLFRAC * Press;
  Km = Kc * (1.0f + O2Pa / Ko);

  /* Temperature dormancy only scales carboxylation, not respiration */
  Vcmax = Vcmax25 * Dormancy * PeakedFactor(Tleaf, PHOTO_EA_VCMAX, PHOTO_HD_VCMAX, PHOTO_DS_VCMAX);
  Jmax  = Vcmax25 * Dormancy * PHOTO_JMAXRATIO * PeakedFactor(Tleaf, PHOTO_EA_JMAX, PHOTO_HD_JMAX, PHOTO_DS_JMAX);
  Rd = Vcmax25 * PHOTO_RD25RATIO * ArrheniusFactor(Tleaf, PHOTO_EA_RD);

  if (Vcmax < 0.0f)
    Vcmax = 0.0f;
  if (Jmax < 0.0f)
    Jmax = 0.0f;

  J = ElectronTransport(ParAbs, Jmax);

  /* --- Medlyn slope, with soil moisture stress applied ------------------- */
  /* m = 1.6 * (1 + G1/sqrt(D)); Beta scales the whole slope so that
     gs -> G0 as Beta -> 0.  Scaling G1 alone would leave a residual
     1.6*An/Cs term and no clean drought floor. */
  Slope = Beta * PHOTO_H2O_CO2_STOM * (1.0f + G1 / (float)sqrt(Vpd));

  /* --- Night / non-assimilating shortcut --------------------------------- */
  /* At zero light the leaf respires, Medlyn has no meaning, and gs = G0.
     This reproduces the behavior the Jarvis scheme gets from RsMax. */
  if (J <= 0.0f || Vcmax <= 0.0f || Slope <= 0.0f) {
    AnLoc = -Rd;
    GsLoc = G0;
    if (An != NULL)
      *An = AnLoc;
    if (Gs != NULL)
      *Gs = GsLoc;
    if (Ci != NULL) {
      /* A respiring leaf has Ci > Ca.  With a very small G0 this expression
         diverges, so bound it: the value is diagnostic only at night. */
      CiLo = Ca - PHOTO_H2O_CO2_BL * AnLoc / Gb
        - PHOTO_H2O_CO2_STOM * AnLoc / GsLoc;
      if (CiLo > 2.0f * Ca)
        CiLo = 2.0f * Ca;
      if (CiLo < 0.0f)
        CiLo = 0.0f;
      *Ci = CiLo;
    }
    return;
  }

  /* --- Bracket the root -------------------------------------------------- */
  /* Lower bound: Ci at the compensation point, where An <= 0, gs = G0, and
     the implied Ci exceeds the trial value (residual > 0).
     Upper bound: Ci = Ca, where assimilation is maximal and the drawdown
     pulls the implied Ci below Ca (residual < 0).
     A sign change across this interval is guaranteed. */
  CiLo = GammaStar / (1.0e-6f * Press);
  CiHi = Ca;
  if (CiLo >= CiHi) {
    CiLo = 0.0f;
    CiHi = Ca;
  }

  FLo = CiResidual(CiLo, Vcmax25, G1, Vcmax, J, Rd, GammaStar, Km, Ca, Press, Gb,
    Slope, G0, &AnLoc, &GsLoc);
  FHi = CiResidual(CiHi, Vcmax25, G1, Vcmax, J, Rd, GammaStar, Km, Ca, Press, Gb,
    Slope, G0, &AnLoc, &GsLoc);

  if (FLo * FHi > 0.0f) {
    /* No sign change: the leaf is not assimilating over the whole interval.
       Fall back to the residual-conductance state rather than iterating. */
    CiMid = (FLo < 0.0f) ? CiHi : CiLo;
    (void)CiResidual(CiMid, Vcmax25, G1, Vcmax, J, Rd, GammaStar, Km, Ca, Press, Gb,
      Slope, G0, &AnLoc, &GsLoc);
    if (An != NULL)
      *An = AnLoc;
    if (Gs != NULL)
      *Gs = GsLoc;
    if (Ci != NULL)
      *Ci = CiMid;
    return;
  }

  /* --- Brent's method ---------------------------------------------------- */
  A = CiLo;
  B = CiHi;
  Fa = FLo;
  Fb = FHi;
  C = B;
  Fc = Fb;
  D = B - A;
  E = D;

  for (Iter = 0; Iter < PHOTO_MAX_ITER; Iter++) {

    if ((Fb > 0.0f && Fc > 0.0f) || (Fb < 0.0f && Fc < 0.0f)) {
      C = A;
      Fc = Fa;
      D = B - A;
      E = D;
    }

    if (fabs(Fc) < fabs(Fb)) {
      A = B;  B = C;  C = A;
      Fa = Fb;  Fb = Fc;  Fc = Fa;
    }

    Tol1 = 2.0f * 1.0e-7f * (float)fabs(B) + 0.5f * PHOTO_CI_TOL;
    Xm = 0.5f * (C - B);

    if (fabs(Xm) <= Tol1 || Fb == 0.0f)
      break;

    if (fabs(E) >= Tol1 && fabs(Fa) > fabs(Fb)) {
      S = Fb / Fa;
      if (A == C) {
        P = 2.0f * Xm * S;
        Q = 1.0f - S;
      }
      else {
        Q = Fa / Fc;
        R = Fb / Fc;
        P = S * (2.0f * Xm * Q * (Q - R) - (B - A) * (R - 1.0f));
        Q = (Q - 1.0f) * (R - 1.0f) * (S - 1.0f);
      }
      if (P > 0.0f)
        Q = -Q;
      P = (float)fabs(P);

      if (2.0f * P < ((3.0f * Xm * Q - (float)fabs(Tol1 * Q)) <
          (float)fabs(E * Q) ? (3.0f * Xm * Q - (float)fabs(Tol1 * Q))
          : (float)fabs(E * Q))) {
        E = D;
        D = P / Q;
      }
      else {
        D = Xm;
        E = D;
      }
    }
    else {
      D = Xm;
      E = D;
    }

    A = B;
    Fa = Fb;

    if (fabs(D) > Tol1)
      B += D;
    else
      B += (Xm > 0.0f ? Tol1 : -Tol1);

    Fb = CiResidual(B, Vcmax25, G1, Vcmax, J, Rd, GammaStar, Km, Ca, Press, Gb,
      Slope, G0, &AnLoc, &GsLoc);
  }

  /* Re-evaluate at the converged Ci so An and Gs are exactly consistent
     with the returned Ci rather than with the last trial point. */
  (void)CiResidual(B, Vcmax25, G1, Vcmax, J, Rd, GammaStar, Km, Ca, Press, Gb,
    Slope, G0, &AnLoc, &GsLoc);

  if (An != NULL)
    *An = AnLoc;
  if (Gs != NULL)
    *Gs = GsLoc;
  if (Ci != NULL)
    *Ci = B;
}

/*****************************************************************************
 PhotoTwoLeafPartition()
 Sunlit/shaded split after de Pury & Farquhar (1997)
 The sunlit fraction shrinks as the sun drops and as LAI grows:
   Lsun = (1 - exp(-kb*L)) / kb,  kb = G / sin(solar altitude)
 Photosynthetic capacity is distributed through the canopy following a
 nitrogen profile Vcmax(l) = Vcmax0*exp(-kn*l), integrated separately over
 the sunlit and shaded fractions.  Sunlit leaves therefore carry both more
 light AND more capacity.
 PAR partitioning here is simplified relative to de Pury & Farquhar: the
 direct beam is assigned entirely to sunlit leaves and the diffuse flux is
 split in proportion to leaf area, rather than using separate kb', kd'
 scattering coefficients.  This conserves absorbed energy exactly, which
 keeps the kernel consistent with whatever DHSVM's radiation scheme already
 absorbed, and it captures the dominant sunlit/shaded asymmetry.
 Could upgrade to the full scattering treatment here eventually.
 *****************************************************************************/
void PhotoTwoLeafPartition(float Lai, float SinAlt, float ParBeam, float ParDiff,
  float Vcmax25, float *LaiSun, float *LaiSha, float *ParSun, float *ParSha,
  float *VcmaxSun, float *VcmaxSha)
{
  float Kb, Kn, LSun, LSha;
  float VcTot, VcSunInt, VcShaInt;
  float ParSunTot, ParShaTot;

  Kn = PHOTO_KN_NITROGEN;

  if (Lai < PHOTO_TINY) {
    *LaiSun = 0.0f;   *LaiSha = 0.0f;
    *ParSun = 0.0f;   *ParSha = 0.0f;
    *VcmaxSun = 0.0f; *VcmaxSha = 0.0f;
    return;
  }

  /* Total canopy carboxylation capacity under the N profile */
  VcTot = Vcmax25 * (1.0f - (float)exp(-Kn * Lai)) / Kn;

  /* --- Night: no beam, all leaf area is "shaded" ------------------------- */
  if (SinAlt <= 0.0f) {
    LSun = 0.0f;
    LSha = Lai;

    *LaiSun = LSun;
    *LaiSha = LSha;
    *ParSun = 0.0f;
    *ParSha = (ParDiff > 0.0f) ?
      (ParDiff * PHOTO_PAR_PER_WATT / LSha) : 0.0f;
    *VcmaxSun = 0.0f;
    *VcmaxSha = VcTot / LSha;
    return;
  }

  /* --- Sunlit / shaded leaf area ----------------------------------------- */
  if (SinAlt < PHOTO_MIN_SINALT)
    SinAlt = PHOTO_MIN_SINALT;

  Kb = PHOTO_G_SPHERICAL / SinAlt;

  LSun = (1.0f - (float)exp(-Kb * Lai)) / Kb;
  if (LSun > Lai)
    LSun = Lai;
  LSha = Lai - LSun;

  /* --- Absorbed PAR per fraction (W/m2 ground -> umol/m2/s leaf) ---------- */
  if (LSun > PHOTO_TINY) {
    ParSunTot = ParBeam + ParDiff * (LSun / Lai);
    *ParSun = ParSunTot * PHOTO_PAR_PER_WATT / LSun;
  }
  else {
    *ParSun = 0.0f;
  }

  if (LSha > PHOTO_TINY) {
    ParShaTot = ParDiff * (LSha / Lai);
    *ParSha = ParShaTot * PHOTO_PAR_PER_WATT / LSha;
  }
  else {
    *ParSha = 0.0f;
  }

  /* --- Photosynthetic capacity per fraction ------------------------------ */
  VcSunInt = Vcmax25 * (1.0f - (float)exp(-(Kn + Kb) * Lai)) / (Kn + Kb);
  VcShaInt = VcTot - VcSunInt;
  if (VcShaInt < 0.0f)
    VcShaInt = 0.0f;

  *VcmaxSun = (LSun > PHOTO_TINY) ? (VcSunInt / LSun) : 0.0f;
  *VcmaxSha = (LSha > PHOTO_TINY) ? (VcShaInt / LSha) : 0.0f;

  *LaiSun = LSun;
  *LaiSha = LSha;
}

/*****************************************************************************
 PhotoCanopyConductance()
 Run the leaf kernel twice (sunlit, shaded) and aggregate.  The return value
 is a canopy conductance in m/s; the caller inverts it to the resistance
 that the Penman-Monteith term in EvapoTranspiration() expects.
 *****************************************************************************/
float PhotoCanopyConductance(float Vcmax25, float G1, float Lai, float SinAlt,
  float ParBeam, float ParDiff, float Tair, float Vpd, float Ca, float Press,
  float Gb, float Beta, float Dormancy, float *AnCanopy)
{
  float LaiSun, LaiSha, ParSun, ParSha, VcmaxSun, VcmaxSha;
  float AnSun, GsSun, CiSun;
  float AnSha, GsSha, CiSha;
  float GsCanopyMol, AnTot;

  AnSun = 0.0f;  GsSun = 0.0f;  CiSun = 0.0f;
  AnSha = 0.0f;  GsSha = 0.0f;  CiSha = 0.0f;

  if (Lai < PHOTO_TINY) {
    if (AnCanopy != NULL)
      *AnCanopy = 0.0f;
    return 0.0f;
  }

  PhotoTwoLeafPartition(Lai, SinAlt, ParBeam, ParDiff, Vcmax25,
    &LaiSun, &LaiSha, &ParSun, &ParSha, &VcmaxSun, &VcmaxSha);

  /* Each big leaf keeps the class water-use strategy but carries its own
     photosynthetic capacity from the canopy nitrogen profile. */
  
  if (LaiSun > PHOTO_TINY)
    PhotoLeafFlux(VcmaxSun, G1, ParSun, Tair, Vpd, Ca, Press, Gb, Beta, Dormancy,
      &AnSun, &GsSun, &CiSun);

  if (LaiSha > PHOTO_TINY)
    PhotoLeafFlux(VcmaxSha, G1, ParSha, Tair, Vpd, Ca, Press, Gb, Beta, Dormancy,
      &AnSha, &GsSha, &CiSha);

  /* Conductances add in parallel across leaf area */
  GsCanopyMol = GsSun * LaiSun + GsSha * LaiSha;
  AnTot = AnSun * LaiSun + AnSha * LaiSha;

  if (AnCanopy != NULL)
    *AnCanopy = AnTot;

  return GsCanopyMol * PhotoMolarToVelocity(Tair, Press);
}
