/*****************************************************************************
  stomatalscheme.h

  Stomatal conductance schemes for DHSVM-MtnHydro.  Each is translated from
  its own source; none shares equations with another beyond the FvCB kernel
  in Photosynthesis.c and the hydraulics in PlantHydraulics.c /
  RootHydraulics.c.

  SOURCES
  -------------------------------------------------------------------------
  [MED11] Medlyn BE et al. (2011) Glob Change Biol 17:2134-2144.
          Unified stomatal optimization.  Implemented in Photosynthesis.c
          (PhotoLeafFlux); this file only wraps it as a scheme.

  [SPE17] Sperry JS et al. (2017) Plant Cell Environ 40:816-830.
          Profit maximization: maximize (beta - theta), beta = A/Amax,
          theta = (Kmax - K)/(Kmax - Kcrit), with K = dE/dPleaf along the
          supply function.

  [WAN20] Wang Y et al. (2020) New Phytologist 227:311-325.
          ProfitMax2: maximize An (Ecrit - E)/Ecrit.  Restated as Eqn 17 in
          Sabot et al. (2022) JAMES 14:e2021MS002761.

  [ELL18] Eller CB et al. (2018) Phil Trans R Soc B 373:20170315.
          SOX solved by numerical iteration over ci.  This is the reference
          solution.

  [ELL20] Eller CB et al. (2020) New Phytologist 226:1622-1637.
          doi:10.1111/nph.16419.  The semi-analytical SOX.  Equation numbers
          below are from this paper and its Notes S1/S2.

  [JON22] Jones S, Eller CB, Cox PM (2022) Front Environ Sci 10:970266.
          Head-to-head comparison of the three SOX solvers.

  SOX EQUATIONS, transcribed
  -------------------------------------------------------------------------
  Eqn 1    maximize  A[ci(gs)] K[Psi_m(gs)]
  Eqn 2    K(Psi)  = 1 / [1 + (Psi/Psi_50)^a]
  Eqn 3    d(AK)/dgs = 0
  Eqn 4    gs = 0.5 (dA/dci) ( sqrt( 4 xi/(dA/dci) + 1 ) - 1 )
  Eqn 5    xi = 2 / [ (1/K)(dK/dPsi_m) r_p 1.6 D ]
  Eqn 6    r_p = r_p,min / K(Psi_pd)
  S1.1a    Psi_m  = (Psi_pd + Psi_c)/2
  S1.2     Psi_pd = Psi_r - h g rho 1e-6
  S1.8     E      = 1.6 gs D
  S1.9     Psi_m  = Psi_pd - r_p 1.6 gs D / 2
  S1.17b   gs^2 + gs (dA/dci) - xi (dA/dci) = 0     (Eqn 4 is its + root)
  S2.1     dA/dci = [A(ca) - A(ci,col)] / (ca - ci,col)
  S2.2     dA/dci = 0  =>  gs = A(ci,col) / (ca - ci,col)
  S2.9     (dK/dPsi_m)(1/K) = [K(Psi_pd) - K((Psi_pd+Psi_50)/2)]
                              / [Psi_pd - (Psi_pd+Psi_50)/2] / K(Psi_pd)

  UNITS.  gs here is conductance to CO2 in mol m-2 s-1; D is the leaf-to-air
  vapour MOLE FRACTION deficit (mol/mol), so E = 1.6 gs D comes out in
  mol H2O m-2 s-1.  dA/dci in umol m-2 s-1 per umol/mol is numerically
  mol m-2 s-1, the same units as gs and xi -- which is why Eqn 4 is
  dimensionally consistent.  r_p,min is m2 s MPa mol-1 H2O.

  DEPARTURES (see docs/provenance.md sec. 7)
  -------------------------------------------------------------------------
  S1  SOX is formulated on Collatz photosynthesis; DHSVM uses FvCB.  The
      translation is confined to PhotoCiColimit()/PhotoDAdCi() (departure D1
      in photosynthesis.h) -- this file uses their output unchanged.
  S2  [ELL20] treat the soil-plant-atmosphere pathway as ONE hydraulic
      compartment with a single r_p,min: "SOX only requires the hydraulic
      parameters rpmin, Psi50 and a".  DHSVM resolves root and xylem stages
      separately, so r_p,min is the series combination,
          1/r_p,min = 1/(1/Krs + 1/KxMax).
      Psi_pd maps to the Phase III effective potential Heff, which is a
      SUF-weighted root-zone mean -- a better-defined version of [ELL20]'s
      "mean soil water potential in the root zone".  The gravity term of
      Notes S1.2 is carried by the xylem object (HYDXYLEM.Pgrav), so SOX,
      ProfitMax and Medlyn-PHS see the same hydrostatic drop.
  S5  Water stress in the Medlyn-based schemes multiplies Vcmax
      (PHOTO_BETA_MODE 0), which is what CLM4.5 BTRAN and CLM5-PHS
      (Kennedy et al. 2019 Eqn 2) do; the slope form (De Kauwe et al. 2015)
      is retained behind PHOTO_BETA_MODE 1.  Before this refactor the slope
      form ran under the CLM5-PHS label.
  S3  [SPE17]'s theta uses K = dE/dPleaf along the supply function.  Before
      Phase III that quantity was non-monotone and the model substituted the
      xylem conductance (HYD_COST_XYLEM).  With Krs/SUF the series
      conductance is monotone again, so Sperry's original definition is
      restored and the substitute is deleted.
  S4  Leaf-basis to air-basis VPD conversion before Penman-Monteith.  A unit
      consistency requirement of coupling a leaf-level scheme to a canopy PM
      expression, not an invention.  Documented, not removed.

  ON THE LOW-LIGHT APERTURE PROBLEM
  -------------------------------------------------------------------------
  Buckley (2017) showed that Wolf et al. (2016) and [SPE17] predict
  unregulated stomatal aperture when the cost function is zero.  [ELL20]
  Notes S2.9 avoid it by evaluating the vulnerability gradient at Psi_50, so
  the cost stays positive even in the flat part of the curve; [WAN20]
  ProfitMax2 avoids it by weighting carbon gain by proximity to hydraulic
  failure.  Both are implemented here.  The pre-Phase-IV HYD_ABSOLUTE_GAIN
  flag was an unattributed patch for the same defect and is deleted.
*****************************************************************************/

/* NAMING HAZARD.  DHSVM's brent.h does  #define T 1e-5  and is pulled in
   through massenergy.h, so no identifier here may be named T.  The STOMTRAIT
   parameter is Tr for that reason -- not style. */
#ifndef STOMATALSCHEME_H
#define STOMATALSCHEME_H

#include "photosynthesis.h"
#include "planthydraulics.h"

enum StomScheme {
  STOM_JARVIS = 0,      /* handled in CanopyResistance.c, listed for parity */
  STOM_MEDLYN,          /* [MED11]                                          */
  STOM_PROFITMAX,       /* [SPE17]                                          */
  STOM_PROFITMAX2,      /* [WAN20]                                          */
  STOM_SOX              /* [ELL18] / [ELL20]                                */
};

enum StomSoxSolver {
  STOM_SOX_NUMERICAL = 0,   /* [ELL18], the reference solution              */
  STOM_SOX_ANALYTICAL       /* [ELL20] Eqns 4-5                             */
};

/* Points scanned when a scheme searches a curve. */
#define STOM_NSCAN        200   /* ci points in the [ELL18] reference scan */
#define STOM_MIN_GS       1.0e-6f      /* mol CO2/m2/s                      */
#define STOM_TINY         1.0e-12f

typedef struct {
  float Vcmax25;      /* umol/m2/s                                          */
  float G1;           /* Medlyn slope (kPa^0.5)                             */
  float G0;           /* Medlyn intercept (mol H2O/m2/s)                    */
  float Dormancy;     /* 0-1                                                */
  float Gmax;         /* cap on conductance to H2O (mol/m2/s); <=0 = none   */
  float Gb;           /* boundary layer conductance to H2O (mol/m2/s)       */
  float CanopyHeight; /* m; informational -- gravity comes from the xylem
                         object's Pgrav so all schemes share it            */
  float RpMin;        /* m2 s MPa mol-1 H2O; <=0 derives from Krs and KxMax */
  float JmaxRatio;    /* Jmax25/Vcmax25; <=0 uses PHOTO_JMAXRATIO           */
  float Rd25Ratio;    /* Rd25/Vcmax25;   <=0 uses PHOTO_RD25RATIO           */
  float LeafWidth;    /* m; <=0 uses PHOTO_LEAF_WIDTH                       */
} STOMTRAIT;

typedef struct {
  float Tair;         /* degC                                               */
  float VpdAir;       /* kPa, AIR basis                                     */
  float Press;        /* Pa                                                 */
  float Ca;           /* umol/mol                                           */
  float ParAbs;       /* umol/m2 leaf/s                                     */
  float Rabs;         /* W/m2 LEAF absorbed by THIS leaf, from
                         PhotoLeafRabs() -- never the canopy net radiation  */
  float Wind;         /* m/s                                                */
} STOMMET;

typedef struct {
  float Gs;           /* conductance to CO2 (mol/m2/s)                      */
  float Gw;           /* conductance to H2O (mol/m2/s)                      */
  float An;           /* umol/m2/s                                          */
  float Ci;           /* umol/mol                                           */
  float E;            /* mmol H2O/m2/s                                      */
  float Ecrit;        /* mmol H2O/m2/s                                      */
  float PsiLeaf;      /* MPa                                                */
  float PsiCollar;    /* MPa                                                */
  float PsiPd;        /* MPa, predawn (SOX)                                 */
  float Tleaf;        /* degC                                               */
  float VpdLeaf;      /* kPa                                                */
  float Objective;    /* scheme's objective at the optimum                  */
  int   Failed;
  int   Iterations;
} STOMSOLUTION;

void StomInit(STOMSOLUTION *S, const STOMMET *M);

/* [MED11].  No hydraulics; Beta is the external soil-moisture multiplier. */
void StomMedlyn(const STOMTRAIT *Tr, const STOMMET *M, float Beta,
                STOMSOLUTION *S);

/* [SPE17] and [WAN20].  Both scan the Phase III supply curve.  Which
   objective is used is the only difference between them. */
void StomProfitMax(const STOMTRAIT *Tr, const STOMMET *M,
                   const HYDSUPPLY *Sup, int Variant, STOMSOLUTION *S);

/* [ELL18] / [ELL20].  Solver selects numerical iteration or the
   semi-analytical Eqns 4-5. */
void StomSOX(const STOMTRAIT *Tr, const STOMMET *M, const HYDXYLEM *X,
             float Heff, float Krs, int Solver, STOMSOLUTION *S);

/* [MED11] stomatal model with plant hydraulic stress, i.e. CLM5-PHS
   (Kennedy et al. 2019, JAMES 11:485).  Beta comes from the LEAF WATER
   POTENTIAL through the vulnerability curve, not from soil moisture:
   beta = K(psi_leaf)/K(0), and attenuates Vcmax (Kennedy Eqn 2; see
   PHOTO_BETA_MODE).  gs and psi_leaf are mutually dependent, so the pair is
   iterated to a fixed point -- CLM5 Newton-solves the same coupling.

   This is the MEDLYN + KRSSUF cell of the option matrix. */
void StomMedlynPHS(const STOMTRAIT *Tr, const STOMMET *M, const HYDXYLEM *X,
                   float Heff, float Krs, STOMSOLUTION *S);

/* Derive KxMax from Vcmax25 via the [SPE17] Ci/Ca coordination.

   The coordination is DEFINED through the profit-maximization optimum, so
   this necessarily runs the optimizer -- which is why it lives here and not
   in PlantHydraulics.c.  Startup cost only: once per vegetation class, never
   per pixel per timestep.  The class's own photosynthetic traits, curve
   form, canopy height (gravity) and the run's CO2 all enter, so that the
   coordination point is the one the run will actually experience.

   CAVEAT.  kmax therefore depends on ProfitMax even when the run uses SOX or
   Medlyn-PHS.  That is Sperry's construction, not an artifact, but it is a
   reason to prefer a MEASURED kmax where one exists (XYLEM CONDUCTANCE).
   Eller et al. (2020) Table 2 give r_p,min by PFT: 1.5-8 mmol-1 m2 s MPa
   observed, i.e. total soil-plant conductance of order 0.1-0.7 mmol m-2 s-1
   MPa-1; the coordination gives 2-5x more for evergreen PFTs.
   Returns <= 0 if no coordination point exists in the bracket. */
float StomKmaxFromVcmax(const PHOTOTRAIT *P, float P50, float VulnShape,
                        int VulnForm, float CiCaTarget, float CanopyHeight,
                        float Ca, float Press);

/* The [SPE17] reference state -- PAR 2000 umol/m2/s, 25 degC, D = 1 kPa,
   wind 2 m/s, saturated soil -- solved with ProfitMax for a BUILT xylem
   object (with its gravity term) and a leaf-basis Krs.  Returns the
   ProfitMax leaf solution (gs, E, An, Ci, psi_leaf) and, through *Ecrit,
   the leaf-basis critical transpiration.  This is what the coordination
   bisects on; it is exposed so startup can report where a class sits
   (Ci/Ca, gs, E/Ecrit) whether KxMax was derived or configured. */
void  StomReferencePoint(const PHOTOTRAIT *P, const HYDXYLEM *X, float Krs,
                         float Ca, float Press, STOMSOLUTION *S,
                         float *Ecrit);

/* Series soil-plant resistance, departure S2. */
float StomRpMin(float Krs, float KxMax);

/* [ELL20] Eqn 6 and Notes S2.9, exposed for testing. */
float StomSoxRp(float RpMin, float PsiPd, const HYDCURVE *V);
float StomSoxKGradient(float PsiPd, const HYDCURVE *V);

#endif
