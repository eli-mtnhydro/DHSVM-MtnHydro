/*****************************************************************************
  roothydraulics.h

  Macroscopic root water uptake for DHSVM-MtnHydro.

  Replaces the layered-rhizosphere scheme that previously lived inside
  PlantHydraulics.c.  That scheme was invented for this model; this one is the
  published macroscopic formalism, and the three departures it carried
  (emergent layered uptake, root radial conductance in series, non-monotone
  soil-plant conductance) are retired by construction rather than patched.

  SOURCES
  -------------------------------------------------------------------------
  [COU12] Couvreur V, Vanderborght J, Javaux M (2012) HESS 16:2957-2971.
          doi:10.5194/hess-16-2957-2012
          Krs / SUF / compensatory uptake.

  [VAN21] Vanderborght J et al. (2021) HESS 25:4835-4860.
          doi:10.5194/hess-25-4835-2021
          Exact upscaling of 3-D root architecture to Krs, SUF and a
          compensatory matrix; proof that the parallel-root reduction
          reproduces the full architecture.

  [LEI25] Leitner D, Schnepf A, Vanderborght J (2025) HESS 29:1759-1782.
          doi:10.5194/hess-29-1759-2025    CC-BY.
          Equation numbers cited below are from this paper, which restates
          [VAN21] and [VDB23] in one place and quantifies the accuracy and
          cost of each upscaling step.

  [VDB23] Vanderborght J et al. (2023/2024) Vadose Zone J 23:e20273.
          doi:10.1002/vzj2.20273
          Perirhizal resistance coupled to the macroscopic root model.

  [SCH08] Schröder Tup et al. (2008) — steady-rate perirhizal model, the basis
          of the K_prhiz used here.

  [VLI06] Van Lier Q de J et al. (2006) — the 0.53 factor in the geometry
          term B: the ratio between the radial distance from the root surface
          at which water content equals the average perirhizal water content,
          and the perirhizal radius.

  [BC64]  Brooks RH, Corey AT (1964) — retention and conductivity relations,
          which DHSVM already uses, and which give the matric flux potential
          in closed form (see RootMatricFluxPotential).

  GOVERNING EQUATIONS  ([LEI25] numbering)
  -------------------------------------------------------------------------
  Eqn 19   Tup       = Krs (Heff - Hcollar)
  Eqn 20   Hcollar = (Krs Heff - Tup) / Krs
           Heff    = SUF^Tup . Hsr          (SUF-weighted mean SOIL-ROOT
                                           INTERFACE potential, not bulk soil)

  Compensation, [COU12] reduced form of [VAN21] Eqn 17 with C7 ~ I:
           Q_i     = SUF_i Tup + Kcomp SUF_i (Hsr_i - Heff)
           sum(Q_i) = Tup exactly, for any Kcomp -- mass conservation is
           structural, not enforced.

  Perirhizal stage, [LEI25] Eqns 21-27:
  Eqn 21   q_r    = (Hsr - Hx) / r1,   r1 = (2 a_root pi l_root k_r)^-1
  Eqn 22   q_sr   = (Hs - Hsr) / r2,   r2 = (2 a_root pi l_root K_prhiz B)^-1
  Eqn 23   K_prhiz = [Phi(h_s) - Phi(h_sr)] / (Hs - Hsr)
  Eqn 24   B      = 2(rho^2 - 1) / [(1 - 0.53 rho)^2 + 2 rho^2 ln(0.53 rho)]
  Eqn 25   rho    = a_prhiz / a_root
  Eqn 27   Hsr    = (a_root k_r Hx + B K_prhiz Hs) / (a_root k_r + B K_prhiz)
  Eqn 29   a_prhiz = sqrt( vol / (pi l_root) + a_root^2 )

  K_prhiz depends on Hsr, so Eqn 27 is implicit and is solved by the
  fixed-point iteration of [LEI25] Algorithm 1, warm-started from the
  previous timestep.

  WHY THIS RETIRES THREE DEPARTURES
  -------------------------------------------------------------------------
  1. Emergent layered uptake is now [COU12]/[VAN21], not an invention.
  2. Winner-take-all is structurally impossible: Heff is a SUF-weighted mean,
     so a single wet layer can only dominate to the extent its SUF is large.
     The ad hoc root radial conductance added to suppress it is unnecessary.
  3. Total uptake and its distribution are SEPARABLE ([VAN21]): the
     redistribution does not depend on total uptake or on collar potential.
     Heff and SUF are computed once per pixel per timestep, before the
     stomatal solve, and never enter its iteration.  That is what makes the
     scheme affordable, and it is why the soil-plant conductance seen by the
     stomatal scheme is monotone again.

  ZERO-DIFF FALLBACK
  -------------------------------------------------------------------------
  With SUF = normalized RootFract and Kcomp = 0, Q_i = SUF_i Tup reduces
  exactly to DHSVM's prescribed-root-fraction uptake.  Phase III is therefore
  a strict generalization with an exact fallback, which the harness checks.

  DEPARTURES (see docs/provenance.md sec. 7)
  -------------------------------------------------------------------------
  R1  The papers work in cm of head and cm3/d.  DHSVM works in m of head and
      MPa.  Conversion only; the perirhizal stage is computed in head units
      and handed over in MPa.
  R2  [COU12]'s scalar Kcomp replaces [VAN21]'s full compensatory matrix.
      Couvreur showed Kcomp ~ Krs for many root systems; the default here is
      Kcomp = Krs.  At three soil layers a full matrix is not identifiable.
  R3  Krs is a config parameter.  Where it is not supplied it is derived from
      the xylem conductance by a resistance-partition fraction, which is the
      pre-Phase-III assumption and remains the least-constrained number in
      the scheme.  It is exposed, not buried.
*****************************************************************************/

/* NAMING HAZARD.  DHSVM's brent.h does  #define T 1e-5  and is included
   BEFORE massenergy.h in RootBrent.c, SnowInterception.c and SnowMelt.c.
   Any identifier named T in a header reachable from massenergy.h is silently
   replaced by a numeric constant, and the resulting error is reported at
   brent.h rather than here.  The transpiration argument is Tup for that
   reason.  Do not "tidy" it back to T. */
#ifndef ROOTHYDRAULICS_H
#define ROOTHYDRAULICS_H

#define ROOT_MAXLAYERS      10

/* [VLI06] via [LEI25] Eqn 24 */
#define ROOT_VANLIER_FACTOR 0.53

/* [LEI25] Algorithm 1 */
#define ROOT_FIXPOINT_MAX   30
#define ROOT_FIXPOINT_TOL   1.0e-3    /* MPa; 1e-4 left the damped interface
                                          iteration at its cap near the
                                          disconnection point                */

#define ROOT_TINY           1.0e-12
#define ROOT_RHO_MIN        1.05      /* a_prhiz/a_root floor; B -> inf at 1 */
#define ROOT_PSI_FLOOR      -100.0    /* MPa, numerical bound                */

/* Fraction of total soil-to-canopy resistance carried by the root system
   when Krs is not supplied (departure R3).  Sperry et al. (2017) and
   Venturas et al. (2018) partition root:stem:leaf resistance 2:1:1, i.e. the
   root system carries about half; Tyree & Ewers (1991) likewise put roughly
   half of whole-plant resistance below ground.  With 0.5, Krs = KxMax.
   (The pre-refactor value 0.05 gave the roots 5% and made the soil-collar
   stage invisible.)  Prefer the configured ROOT CONDUCTANCE. */
#define ROOT_KRS_FRAC       0.5

/* Fine root radius default (m).  Jackson et al. (1997) give mean fine root
   diameters of ~0.3 mm across biomes; overridden per class by FINE ROOT
   RADIUS.  Enters only the perirhizal geometry (rho = a_prhiz/a_root). */
#define ROOT_RADIUS_DEFAULT 1.5e-4

typedef struct {
  int   N;                              /* soil layers in the root zone      */

  /* --- root system, macroscopic ([COU12], [VAN21]) --------------------- */
  float SUF[ROOT_MAXLAYERS];            /* standard uptake fraction, sums 1  */
  float Krs;                            /* root system conductance
                                           (mmol/m2 GROUND/s/MPa).  The
                                           vegetation table holds the
                                           leaf-specific value; the seam
                                           multiplies by LAI before the
                                           ground-basis RootSolve().        */
  float Kcomp;                          /* compensatory conductance, same
                                           units; 0 disables compensation    */

  /* --- per-layer soil state ------------------------------------------- */
  float Hs[ROOT_MAXLAYERS];             /* bulk soil water potential (MPa)   */
  float Hsr[ROOT_MAXLAYERS];            /* soil-root interface potl (MPa)    */

  /* --- perirhizal geometry and soil hydraulics ([LEI25] 21-29) --------- */
  float ArootKr[ROOT_MAXLAYERS];        /* a_root * k_r  (m/s)               */
  float Bgeom[ROOT_MAXLAYERS];          /* geometry factor B, Eqn 24         */
  float KsVert[ROOT_MAXLAYERS];         /* saturated conductivity (m/s)      */
  float Lambda[ROOT_MAXLAYERS];         /* Brooks-Corey pore size index      */
  float HaeM[ROOT_MAXLAYERS];           /* air entry, m of head (positive)   */

  int   UsePerirhizal;                  /* 0 = Hsr := Hs                     */

  /* Critical collar potential.  [LEI25] Sect. 2.3: "the boundary condition
     will automatically be switched between Neumann and Dirichlet, ensuring
     that the root collar potential cannot be below a critical potential
     where we assume the plant's wilting point."
     Without this the perirhizal fixed point runs away in dry soil: a lower
     Hsr lowers Heff, which lowers Hcollar, which lowers Hsr again.  With a
     prescribed Tup that the soil cannot supply there is no solution, and the
     switch is what makes the problem well posed. */
  float HcollarMin;                     /* MPa, negative                     */
} ROOTZONE;

/* Result of one solve. */
typedef struct {
  float Heff;                           /* SUF-weighted interface potl (MPa) */
  float Hcollar;                        /* root collar potential (MPa)       */
  float Q[ROOT_MAXLAYERS];              /* uptake per layer (mmol/m2/s)      */
  float Efrac[ROOT_MAXLAYERS];          /* Q_i / Tup, sums to 1                */
  float Tactual;                        /* uptake actually supplied
                                           (mmol/m2/s); < Tup when limited     */
  int   Limited;                        /* 1 if the Dirichlet switch fired   */
  int   NLayers;
  int   Iterations;                     /* perirhizal fixed-point count      */
  int   Converged;
} ROOTSOLUTION;

/* --------------------------------------------------------------------------
   Building blocks -- each maps to one numbered equation, and each is
   separately testable.
   -------------------------------------------------------------------------- */

/* Brooks-Corey matric flux potential Phi(h) = integral of K(h') dh' from
   -infinity to h, in closed form.  h in m of head (negative), Hae positive.
   Units m2/s.  [BC64] with K = Ks (Hae/|h|)^(3 lambda + 2). */
float RootMatricFluxPotential(float h, float Ks, float Hae, float Lambda);

/* Geometry factor B.  [LEI25] Eqn 24. */
float RootGeometryFactor(float Rho);

/* Outer perirhizal radius from root length density.  [LEI25] Eqn 29 with
   Eqn 30: the perirhizal volume per unit root length is 1/RLD. */
float RootPerirhizalRadius(float Rld, float Aroot);

/* Average perirhizal conductance.  [LEI25] Eqn 23.  Returns the limiting
   dPhi/dH as Hs -> Hsr rather than dividing by zero. */
float RootPerirhizalK(float Hs, float Hsr, float Ks, float Hae, float Lambda);

/* Soil-root interface potential for one layer at a known xylem potential.
   [LEI25] Eqn 27, solved for the implicit K_prhiz dependence. */
float RootInterfacePotential(float Hs, float Hx, float ArootKr, float Bgeom,
                             float Ks, float Hae, float Lambda);

/* --------------------------------------------------------------------------
   The macroscopic model
   -------------------------------------------------------------------------- */

/* Zero the struct and set safe defaults. */
void RootZoneInit(ROOTZONE *R, int NLayers);

/* SUF from prescribed root fractions, normalized to sum to 1.  This is the
   zero-diff default (see header note). */
void RootSetSUFFromFractions(ROOTZONE *R, const float *RootFract);

/* Perirhizal geometry for every layer, from root length density and fine
   root radius.  Rld in m root per m3 soil, Aroot in m. */
void RootSetPerirhizal(ROOTZONE *R, const float *Rld, float Aroot,
                       float Kradial);

/* Root radial conductivity k_r (m/s per m of head, i.e. 1/s) implied by a
   ground-basis root system conductance Krs (mmol/m2 ground/s/MPa) for a
   radially limited root system: Krs = 2 pi a_root L k_r, with L the root
   length per m2 ground (m/m2).  Using this instead of a separate k_r keeps
   the perirhizal stage ([LEI25] Eqn 27) consistent with the macroscopic
   Krs the rest of the model runs on. */
float RootRadialFromKrs(float KrsGround, float Aroot, float RootLength);

/* Bulk soil water potential (MPa, negative) at which the perirhizal
   conductance B K(h) has fallen to the root radial conductance a_root k_r,
   i.e. the interface sits midway between bulk soil and root xylem.  A
   startup diagnostic: where, for this soil and this root system, the
   rhizosphere starts to disconnect.  Returns 0 if undefined. */
float RootDisconnectPotential(float ArootKr, float Bgeom, float Ks,
                              float Hae, float Lambda);

/* Krs when it is not configured (departure R3): the root system is assigned
   ROOT_KRS_FRAC of the total soil-to-canopy resistance at saturation. */
float RootKrsFromXylem(float KxMax);

/* Effective (SUF-weighted) potential.  Uses Hsr when the perirhizal stage is
   on, Hs otherwise. */
float RootEffectivePotential(const ROOTZONE *R);

/* Collar potential for a demanded transpiration.  [LEI25] Eqn 20.
   ANALYTIC -- this is the call the stomatal schemes make inside their solve,
   and it costs one subtraction and one divide. */
float RootCollarPotential(const ROOTZONE *R, float Heff, float Tup);

/* Full solve: given demanded Tup (mmol/m2 ground/s), find consistent Hsr,
   Heff, Hcollar and the per-layer uptake.  When UsePerirhizal is 0 this is
   closed-form; when it is 1 it runs [LEI25] Algorithm 1.

   HsrPrev may be NULL; if given it warm-starts the iteration from the
   previous timestep, as [LEI25] recommend. */
void RootSolve(ROOTZONE *R, float Tup, const float *HsrPrev, ROOTSOLUTION *S);

/* Uptake distribution at a known Heff and Tup.  [COU12] compensation.
   Separated out because it does NOT depend on Hcollar or on the stomatal
   solve, which is the property that makes the scheme cheap. */
void RootUptakeDistribution(const ROOTZONE *R, float Heff, float Tup,
                            float *Q, float *Efrac);

#endif
