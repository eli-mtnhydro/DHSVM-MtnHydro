/*****************************************************************************
  RootHydraulics.c

  Macroscopic root water uptake.  Every function maps to one numbered equation
  in roothydraulics.h; the harness suite `phase3` tests them individually
  before testing them together.

  No DHSVM dependencies: math.h, stdlib.h, string.h and its own header.
*****************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "roothydraulics.h"

/* Conversion between m of head and MPa for water at ~20 degC.
   rho g = 9810 Pa/m, so 1 m head = 9.81e-3 MPa. */
#define ROOT_M_TO_MPA   (9810.0 / 1.0e6)

/* M_PI is POSIX, not C99; DHSVM must build with -std=c99 on MSVC too. */
#ifndef ROOT_PI
#define ROOT_PI 3.14159265358979323846
#endif
#define ROOT_MPA_TO_M   (1.0e6 / 9810.0)

/* ========================================================================= */
/* Brooks-Corey matric flux potential                                        */
/* ========================================================================= */

/* Phi(h) = int_{-inf}^{h} K(h') dh'.

   Brooks-Corey unsaturated conductivity, in the exponent DHSVM already uses
   in UnsaturatedFlow.c:
       K(h) = Ks (Hae/|h|)^m,    m = 3 lambda + 2
   Substituting u = |h| and integrating from infinity down to u,
       Phi = Ks Hae^m u^(1-m) / (m-1)
   which converges for m > 1, true for all lambda > 0.

   Above air entry the soil is saturated and K = Ks, so Phi continues
   linearly.  Analytic throughout: no quadrature, no table. */
float RootMatricFluxPotential(float h, float Ks, float Hae, float Lambda)
{
  double m, u, PhiAe;

  if (Ks <= 0.0f) return 0.0f;
  if (Hae <= 0.0f) Hae = 0.01f;
  if (Lambda <= 0.0f) Lambda = 0.1f;

  m = 3.0 * (double)Lambda + 2.0;
  if (m <= 1.0 + ROOT_TINY) m = 1.0 + 1.0e-6;

  u = fabs((double)h);
  if (u < (double)Hae) {
    /* Wetter than air entry: saturated branch. */
    PhiAe = (double)Ks * (double)Hae / (m - 1.0);
    return (float)(PhiAe + (double)Ks * ((double)Hae - u));
  }

  return (float)((double)Ks * pow((double)Hae, m) * pow(u, 1.0 - m) / (m - 1.0));
}

/* ========================================================================= */
/* Perirhizal geometry                                                       */
/* ========================================================================= */

/* [LEI25] Eqn 24:
       B = 2(rho^2 - 1) / [ (1 - 0.53 rho)^2 + 2 rho^2 ln(0.53 rho) ]
   The 0.53 is [VLI06].  B diverges as rho -> 1 (no perirhizal zone at all),
   so rho is floored. */
float RootGeometryFactor(float Rho)
{
  double r, c, Num, Den;

  r = (double)Rho;
  if (r < ROOT_RHO_MIN) r = ROOT_RHO_MIN;

  c = ROOT_VANLIER_FACTOR * r;
  if (c < ROOT_TINY) return 0.0f;

  Num = 2.0 * (r * r - 1.0);
  Den = (1.0 - c) * (1.0 - c) + 2.0 * r * r * log(c);

  if (fabs(Den) < ROOT_TINY)
    return 0.0f;

  /* B must be positive; the expression can go negative for rho very close
     to 1/0.53 where the log changes sign. */
  if (Num / Den <= 0.0)
    return 0.0f;

  return (float)(Num / Den);
}

/* [LEI25] Eqns 29 + 30.  Perirhizal volume per unit root length is the
   reciprocal of root length density, so
       a_prhiz = sqrt( 1/(pi RLD) + a_root^2 ).  */
float RootPerirhizalRadius(float Rld, float Aroot)
{
  double v;

  if (Aroot <= 0.0f) Aroot = 1.0e-4f;
  if (Rld <= 0.0f)   return 1.0e3f * Aroot;   /* no roots: unbounded zone   */

  v = 1.0 / (ROOT_PI * (double)Rld) + (double)Aroot * (double)Aroot;
  return (float)sqrt(v);
}

/* ========================================================================= */
/* Perirhizal conductance and interface potential                            */
/* ========================================================================= */

/* [LEI25] Eqn 23:  K_prhiz = [Phi(h_s) - Phi(h_sr)] / (H_s - H_sr).
   As Hsr -> Hs this is 0/0 and the limit is dPhi/dh = K(h), so the
   unsaturated conductivity at Hs is returned instead. */
float RootPerirhizalK(float Hs, float Hsr, float Ks, float Hae, float Lambda)
{
  float dH, PhiS, PhiSr;
  double m, u;

  dH = Hs - Hsr;

  if ((float)fabs((double)dH) < 1.0e-6f) {
    /* Limit: K(h) itself. */
    if (Ks <= 0.0f) return 0.0f;
    if (Hae <= 0.0f) Hae = 0.01f;
    if (Lambda <= 0.0f) Lambda = 0.1f;
    m = 3.0 * (double)Lambda + 2.0;
    u = fabs((double)Hs);
    if (u < (double)Hae) return Ks;
    return (float)((double)Ks * pow((double)Hae / u, m));
  }

  PhiS  = RootMatricFluxPotential(Hs,  Ks, Hae, Lambda);
  PhiSr = RootMatricFluxPotential(Hsr, Ks, Hae, Lambda);

  return (PhiS - PhiSr) / dH;
}

/* [LEI25] Eqn 27:
       Hsr = (a_root k_r Hx + B K_prhiz Hs) / (a_root k_r + B K_prhiz)
   with K_prhiz itself a function of Hsr.  Solved by damped fixed-point
   iteration, which is what [LEI25] Algorithm 1 does at the system level.

   All potentials here are in m of head.  Hsr is bracketed between Hx and Hs
   by construction of the weighted mean, so the iteration cannot run away. */
float RootInterfacePotential(float Hs, float Hx, float ArootKr, float Bgeom,
                             float Ks, float Hae, float Lambda)
{
  float Hsr, HsrNew, Kp, Wr, Wp, Denom;
  int it;

  if (ArootKr <= 0.0f || Bgeom <= 0.0f)
    return Hs;                       /* no resolved interface: Hsr = Hs     */

  if (Hx >= Hs)
    return Hs;                       /* no uptake, or redistribution        */

  Hsr = Hs;                          /* start at bulk soil                  */

  for (it = 0; it < ROOT_FIXPOINT_MAX; it++) {
    Kp = RootPerirhizalK(Hs, Hsr, Ks, Hae, Lambda);
    if (Kp < 0.0f) Kp = 0.0f;

    Wr = ArootKr;
    Wp = Bgeom * Kp;
    Denom = Wr + Wp;
    /* Both weights are SI conductivities (m/s), of order 1e-13 for real
       root systems in drying soil, so the guard has to be a true zero test:
       the former ROOT_TINY (1e-12) silently returned Hs = Hsr and switched
       the whole stage off exactly where it matters. */
    if (Denom <= 0.0f)
      return Hs;

    HsrNew = (Wr * Hx + Wp * Hs) / Denom;

    /* Keep it in the physical bracket [Hx, Hs]. */
    if (HsrNew > Hs) HsrNew = Hs;
    if (HsrNew < Hx) HsrNew = Hx;

    /* Damping: K_prhiz falls very steeply with Hsr in dry soil, and an
       undamped step overshoots and oscillates. */
    HsrNew = 0.5f * (HsrNew + Hsr);

    if ((float)fabs((double)(HsrNew - Hsr)) <
        (float)(ROOT_FIXPOINT_TOL * ROOT_MPA_TO_M)) {
      Hsr = HsrNew;
      break;
    }
    Hsr = HsrNew;
  }

  return Hsr;
}

/* ========================================================================= */
/* Setup                                                                     */
/* ========================================================================= */

void RootZoneInit(ROOTZONE *R, int NLayers)
{
  int i;

  if (R == NULL) return;
  memset(R, 0, sizeof(*R));

  if (NLayers < 1) NLayers = 1;
  if (NLayers > ROOT_MAXLAYERS) NLayers = ROOT_MAXLAYERS;
  R->N = NLayers;

  for (i = 0; i < NLayers; i++) {
    R->SUF[i]     = 1.0f / (float)NLayers;
    R->Hs[i]      = 0.0f;
    R->Hsr[i]     = 0.0f;
    R->Bgeom[i]   = 0.0f;
    R->ArootKr[i] = 0.0f;
    R->KsVert[i]  = 1.0e-6f;
    R->Lambda[i]  = 0.2f;
    R->HaeM[i]    = 0.3f;
  }

  R->Krs   = 0.0f;
  R->Kcomp = 0.0f;
  R->UsePerirhizal = 0;
  R->HcollarMin = -8.0f;      /* MPa; override per vegetation class         */
}

void RootSetSUFFromFractions(ROOTZONE *R, const float *RootFract)
{
  float Sum = 0.0f;
  int i;

  if (R == NULL || RootFract == NULL) return;

  for (i = 0; i < R->N; i++)
    Sum += (RootFract[i] > 0.0f) ? RootFract[i] : 0.0f;

  if (Sum < ROOT_TINY) {
    for (i = 0; i < R->N; i++) R->SUF[i] = 1.0f / (float)R->N;
    return;
  }

  for (i = 0; i < R->N; i++)
    R->SUF[i] = ((RootFract[i] > 0.0f) ? RootFract[i] : 0.0f) / Sum;
}

void RootSetPerirhizal(ROOTZONE *R, const float *Rld, float Aroot,
                       float Kradial)
{
  float Aprhiz, Rho;
  int i;

  if (R == NULL) return;
  if (Aroot <= 0.0f) Aroot = 1.0e-4f;

  for (i = 0; i < R->N; i++) {
    float rld = (Rld != NULL && Rld[i] > 0.0f) ? Rld[i] : 0.0f;

    Aprhiz = RootPerirhizalRadius(rld, Aroot);
    Rho    = Aprhiz / Aroot;

    R->Bgeom[i]   = RootGeometryFactor(Rho);
    R->ArootKr[i] = Aroot * Kradial;
  }
}

float RootRadialFromKrs(float KrsGround, float Aroot, float RootLength)
{
  double KrsSI, Area;

  if (KrsGround <= 0.0f || Aroot <= 0.0f || RootLength <= 0.0f) return 0.0f;

  /* mmol/m2/s/MPa -> (m3/m2/s) per m of head:
     x 18.015e-9 m3 per mmol, / (m of head per MPa). */
  KrsSI = (double)KrsGround * 18.015e-9 * ROOT_M_TO_MPA;
  Area  = 2.0 * ROOT_PI * (double)Aroot * (double)RootLength;   /* m2/m2 */
  return (float)(KrsSI / Area);
}

float RootDisconnectPotential(float ArootKr, float Bgeom, float Ks,
                              float Hae, float Lambda)
{
  double m, ratio, h;

  if (ArootKr <= 0.0f || Bgeom <= 0.0f || Ks <= 0.0f) return 0.0f;
  if (Hae <= 0.0f) Hae = 0.01f;
  if (Lambda <= 0.0f) Lambda = 0.1f;

  /* B Ks (Hae/|h|)^m = a k_r  ->  |h| = Hae (B Ks / (a k_r))^(1/m) */
  m = 3.0 * (double)Lambda + 2.0;
  ratio = (double)Bgeom * (double)Ks / (double)ArootKr;
  if (ratio <= 1.0) return -(float)Hae * (float)ROOT_M_TO_MPA;
  h = (double)Hae * pow(ratio, 1.0 / m);          /* m of head            */
  return -(float)(h * ROOT_M_TO_MPA);
}

/* Departure R3.  The root system is given ROOT_KRS_FRAC of the total
   soil-to-canopy resistance at saturation, so
       1/Krs = f/(1/Kx) ... i.e. Krs = Kx (1-f)/f.
   With f = 0.5 (Sperry et al. 2017 root:stem:leaf = 2:1:1) this is
   Krs = KxMax.  Fallback only; prefer a configured ROOT CONDUCTANCE. */
float RootKrsFromXylem(float KxMax)
{
  if (KxMax <= 0.0f) return 0.0f;
  return KxMax * (1.0f - (float)ROOT_KRS_FRAC) / (float)ROOT_KRS_FRAC;
}

/* ========================================================================= */
/* The macroscopic model                                                     */
/* ========================================================================= */

float RootEffectivePotential(const ROOTZONE *R)
{
  double Heff = 0.0;
  int i;

  if (R == NULL) return 0.0f;

  for (i = 0; i < R->N; i++)
    Heff += (double)R->SUF[i] *
            (double)(R->UsePerirhizal ? R->Hsr[i] : R->Hs[i]);

  /* Brooks-Corey runs to -infinity as saturation -> 0.  The plant is long
     past hydraulic failure well before that, so bound it: an unbounded Heff
     only propagates inf/NaN into the stomatal schemes. */
  if (Heff < (double)ROOT_PSI_FLOOR) Heff = (double)ROOT_PSI_FLOOR;

  return (float)Heff;
}

/* [LEI25] Eqn 20.  One subtraction and one divide -- this is the whole
   soil-to-collar stage, and it is why the stomatal schemes can call it
   inside their inner loop without noticing. */
float RootCollarPotential(const ROOTZONE *R, float Heff, float Tup)
{
  if (R == NULL || R->Krs <= ROOT_TINY)
    return Heff;
  return Heff - Tup / R->Krs;
}

/* [COU12] compensation, the reduced form of [VAN21] Eqn 17 (departure R2):
       Q_i = SUF_i Tup + Kcomp SUF_i (H_i - Heff)

   sum(Q_i) = Tup for any Kcomp, because sum(SUF_i (H_i - Heff)) = 0 by the
   definition of Heff.  Mass conservation is therefore structural. */
void RootUptakeDistribution(const ROOTZONE *R, float Heff, float Tup,
                            float *Q, float *Efrac)
{
  float Hi, Qi, Comp, MaxComp = 0.0f, Sum = 0.0f, Pos = 0.0f;
  float Ql[ROOT_MAXLAYERS];
  int i;

  if (R == NULL) return;

  for (i = 0; i < R->N; i++) {
    Hi = R->UsePerirhizal ? R->Hsr[i] : R->Hs[i];
    Comp = R->Kcomp * R->SUF[i] * (Hi - Heff);
    Qi = R->SUF[i] * Tup + Comp;
    Ql[i] = Qi;
    if (Q != NULL) Q[i] = Qi;
    Sum += Qi;
    if ((float)fabs((double)Comp) > MaxComp)
      MaxComp = (float)fabs((double)Comp);
  }

  if (Efrac == NULL) return;

  /* Efrac is a WEIGHTING, and DHSVM multiplies potential transpiration by it.
     Q_i/Tup is the mathematically right fraction but it is unusable as a
     weight whenever the compensation flux is large next to net uptake: as
     Tup -> 0 the compensatory terms stay finite and opposite-signed, their sum
     is Tup, and Q_i/Tup runs to +/- 1e7.  That is hydraulic redistribution
     expressing itself, not an error -- but it must not reach the water
     balance as a multiplier.

     So: use Q_i/Tup only where net uptake dominates the compensation.  Where
     it does not, fall back to SUF, which is bounded, sums to 1, and is the
     uptake distribution at uniform soil potential -- exactly the right
     answer when there is no net flux to distribute.  The true fluxes are
     still in Q[] for diagnostics and for hydraulic redistribution. */
  if ((float)fabs((double)Sum) < ROOT_TINY ||
      MaxComp > (float)fabs((double)Sum)) {
    for (i = 0; i < R->N; i++) Efrac[i] = R->SUF[i];
    return;
  }

  for (i = 0; i < R->N; i++) {
    Efrac[i] = Ql[i] / Sum;
    if (Efrac[i] < 0.0f) Efrac[i] = 0.0f;   /* HR is reported via Q[], not  */
    Pos += Efrac[i];                        /* fed back as a negative weight */
  }
  if (Pos > ROOT_TINY) {
    for (i = 0; i < R->N; i++) Efrac[i] /= Pos;
  }
  else {
    for (i = 0; i < R->N; i++) Efrac[i] = R->SUF[i];
  }
}

void RootSolve(ROOTZONE *R, float Tup, const float *HsrPrev, ROOTSOLUTION *S)
{
  float Heff, Hcollar, HxM, HsM, HsrM;
  float HsrOld[ROOT_MAXLAYERS];
  float MaxChange;
  int i, it;

  if (R == NULL || S == NULL) return;

  memset(S, 0, sizeof(*S));
  S->NLayers = R->N;
  S->Converged = 1;
  S->Iterations = 0;

  /* --- perirhizal stage off: closed form ------------------------------- */
  if (!R->UsePerirhizal) {
    for (i = 0; i < R->N; i++) R->Hsr[i] = R->Hs[i];
    Heff = RootEffectivePotential(R);
    Hcollar = RootCollarPotential(R, Heff, Tup);

    /* Dirichlet switch: the collar cannot go below the critical potential. */
    if (Hcollar < R->HcollarMin) {
      Hcollar = R->HcollarMin;
      Tup = R->Krs * (Heff - Hcollar);
      if (Tup < 0.0f) Tup = 0.0f;
      S->Limited = 1;
    }

    S->Heff = Heff;
    S->Hcollar = Hcollar;
    S->Tactual = Tup;
    RootUptakeDistribution(R, Heff, Tup, S->Q, S->Efrac);
    return;
  }

  /* --- perirhizal stage on: [LEI25] Algorithm 1 ------------------------ */
  for (i = 0; i < R->N; i++)
    R->Hsr[i] = (HsrPrev != NULL) ? HsrPrev[i] : R->Hs[i];

  Heff = RootEffectivePotential(R);
  Hcollar = RootCollarPotential(R, Heff, Tup);

  for (it = 0; it < ROOT_FIXPOINT_MAX; it++) {
    memcpy(HsrOld, R->Hsr, sizeof(float) * (size_t)R->N);

    /* Step 1: Hsr from the current xylem potential.
       In the parallel-root reduction every layer's root xylem sits at the
       collar potential, so Hx_i = Hcollar.  [LEI25] Sect. 2.6, model Cxx --
       the reduction they show reproduces the full architecture well. */
    HxM = Hcollar * (float)ROOT_MPA_TO_M;

    for (i = 0; i < R->N; i++) {
      HsM = R->Hs[i] * (float)ROOT_MPA_TO_M;
      HsrM = RootInterfacePotential(HsM, HxM, R->ArootKr[i], R->Bgeom[i],
                                    R->KsVert[i], R->HaeM[i], R->Lambda[i]);
      R->Hsr[i] = HsrM * (float)ROOT_M_TO_MPA;
      if (R->Hsr[i] < (float)ROOT_PSI_FLOOR) R->Hsr[i] = (float)ROOT_PSI_FLOOR;
      if (R->Hsr[i] > R->Hs[i]) R->Hsr[i] = R->Hs[i];
    }

    /* Step 2: Heff and collar from the updated interface potentials, with
       the Neumann -> Dirichlet switch.  This is what bounds the iteration:
       without it a lower Hsr lowers Heff, which lowers Hcollar, which lowers
       Hsr further, and the fixed point diverges in dry soil. */
    Heff = RootEffectivePotential(R);
    Hcollar = RootCollarPotential(R, Heff, Tup);
    if (Hcollar < R->HcollarMin) {
      Hcollar = R->HcollarMin;
      S->Limited = 1;
    }

    MaxChange = 0.0f;
    for (i = 0; i < R->N; i++) {
      float d = (float)fabs((double)(R->Hsr[i] - HsrOld[i]));
      if (d > MaxChange) MaxChange = d;
    }

    S->Iterations = it + 1;
    if (MaxChange < (float)ROOT_FIXPOINT_TOL)
      break;
  }

  if (S->Iterations >= ROOT_FIXPOINT_MAX &&
      ROOT_FIXPOINT_MAX > 1)
    S->Converged = 0;

  if (S->Limited) {
    Tup = R->Krs * (Heff - Hcollar);
    if (Tup < 0.0f) Tup = 0.0f;
  }

  S->Heff = Heff;
  S->Hcollar = Hcollar;
  S->Tactual = Tup;
  RootUptakeDistribution(R, Heff, Tup, S->Q, S->Efrac);
}
