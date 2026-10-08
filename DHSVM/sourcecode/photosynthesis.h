/*****************************************************************************
  photosynthesis.h

  DHSVM-MtnHydro leaf photosynthesis, temperature kinetics and leaf energy
  balance.  Rebuilt in Phase II directly from published sources; see
  docs/provenance.md for the full dossier.

  SOURCES (one block, one source, no silent hybrids)
  -------------------------------------------------------------------------
  [FvCB]  Farquhar GD, von Caemmerer S, Berry JA (1980) Planta 149:78-90.
          doi:10.1007/BF00386231
          The Rubisco- and RuBP-limited carboxylation rates.

  [BER]   Bernacchi CJ, Singsaas EL, Pimentel C, Portis AR, Long SP (2001)
          Plant Cell Environ 24:253-260.
          Rubisco kinetics temperature responses (Kc, Ko, Gamma*, Rd).
          Values below are the CLM4.5/CLM5 tabulation of Bernacchi, i.e.
          Kc25 = 404.9 umol/mol, Ko25 = 278.4 mmol/mol, Gamma*25 = 42.75
          umol/mol, each multiplied by 101325 Pa to give Pa.  NOTE: rpmodel
          and some other implementations use the direct Bernacchi values
          (39.97 Pa / 27480 Pa), which differ from these by a factor 1.027.
          Both trace to the same paper.  Do not "fix" one to the other.

  [MED02] Medlyn BE et al. (2002) Plant Cell Environ 25:1167-1179.
          doi:10.1046/j.1365-3040.2002.00891.x
          The FUNCTIONAL FORM of the peaked Arrhenius response.

  [KUM19] Kumarathunge DP et al. (2019) New Phytologist 222:768-784.
          doi:10.1111/nph.15668
          The PARAMETER VALUES for the peaked Arrhenius response, as adopted
          as defaults in plantecophys >= v1.4 (EaV, delsC, EdVC, EaJ, delsJ,
          EdVJ).  plantecophys is GPL: its R source was used only as a
          numerical oracle during verification, never transcribed.

  [MED11] Medlyn BE et al. (2011) Glob Change Biol 17:2134-2144.
          doi:10.1111/j.1365-2486.2010.02375.x
          Unified stomatal optimization.  Analytic limit at g0 = 0 is
          Ci/Cs = g1/(g1 + sqrt(D)) -- note Cs, not Ca.

  [DPF97] de Pury DGG, Farquhar GD (1997) Plant Cell Environ 20:537-557.
          Two-leaf sunlit/shaded canopy scaling.

  [CN98]  Campbell GS, Norman JM (1998) An Introduction to Environmental
          Biophysics, 2nd ed., Springer.  Ch. 7 (boundary layer conductance,
          Table 7.6) and Ch. 14 (leaf energy budget, isothermal form).

  [KEN19] Kennedy D et al. (2019) J Adv Model Earth Syst 11:485-513.
          doi:10.1029/2018MS001500.  CLM5 plant hydraulic stress: the water
          stress factor attenuates Vcmax (their Eqn 2), not the Medlyn slope.

  [DEK15] De Kauwe MG et al. (2015) Biogeosciences 12:7503-7518.
          The alternative: stress applied to the Medlyn g1 slope (CABLE).

  [MAK04] Makela A, Hari P, Berninger F, Hanninen H, Nikinmaa E (2004)
          Tree Physiol 24:369-376; Kolari P et al. (2007) Tellus B 59:542-552.
          Delayed-temperature state of acclimation, dS/dt = (T - S)/tau,
          tau ~ 8-14 d (190-330 h).

  [ELL20] Eller CB et al. (2020) New Phytologist 226:1622-1637.
          doi:10.1111/nph.16419
          Notes S2: the co-limitation Ci (Eqn S2.6b) and the secant estimate
          of dA/dCi (Eqn S2.1) that SOX requires.

  [COL91] Collatz GJ, Ball JT, Grivet C, Berry JA (1991) Agric For Meteorol
          54:107-136.  Co-limitation smoothing (quadratic form).

  DEPARTURES (see docs/provenance.md sec. 7)
  -------------------------------------------------------------------------
  D1  SOX is formulated on the Collatz photosynthesis model; DHSVM uses FvCB.
      PhotoCiColimit() therefore translates Eqn S2.6b from Collatz symbols to
      FvCB symbols.  The mapping is documented at the function.
  D2  Absorbed PAR arrives in W/m2 (DHSVM's Rp) and is converted with a fixed
      4.57 umol/J.  This is a unit convention, not physics.
  D3  Leaf energy balance.  Each big leaf is given its OWN absorbed radiation
      (PhotoLeafRabs): absorbed PAR scaled to total shortwave with a leaf
      NIR:PAR absorptance ratio, plus the isothermal longwave term eps*sigma*
      Tair^4, so the residual is the [CN98] isothermal net radiation form.
      Sensible heat leaves both faces of the leaf (PHOTO_LEAF_SIDES).  This
      replaces the pre-refactor practice of handing the leaf the CANOPY net
      radiation (already net of longwave emission) and subtracting emission a
      second time, which cooled leaves by 5-8 K whenever Rn was small.
  D4  Two-leaf diffuse split follows [DPF97]: sunlit leaves receive the
      share of absorbed diffuse PAR given by the diffuse extinction kd and
      the beam extinction kb, not their LAI fraction.

  WHAT CHANGED IN PHASE II
  -------------------------------------------------------------------------
  1. There is now ONE FvCB kernel.  The duplicate Arrhenius/peaked/FvCB set
     in PlantHydraulics.c is deleted and its callers use this file.
  2. The diffusion balance closes on NET assimilation everywhere.  The old
     HydAssimilationAtGc() closed on gross, drawing Ci down by an extra
     Rd/Gc.  See PhotoAssimilationAtGc().
  3. The leaf energy balance moved here from PlantHydraulics.c and is now
     solved against an explicit residual rather than a linearization.
  4. Every constant carries its source.

  WHAT CHANGED IN THE PHASE-VI REFACTOR
  -------------------------------------------------------------------------
  5. Per-class photosynthetic traits (PHOTOTRAIT): Jmax25:Vcmax25, Rd25:
     Vcmax25 and leaf width are no longer compile-time constants.  The old
     entry points remain and use the defaults below.
  6. Water stress (beta) multiplies Vcmax, as in CLM4.5/CLM5 [KEN19]; the
     De Kauwe slope form is kept behind PHOTO_BETA_MODE.
  7. Leaf energy balance rebuilt per departure D3; two-leaf diffuse split
     per D4.  Atmospheric CO2 is a run constant read from the config
     ([CONSTANTS] ATMOSPHERIC CO2), PHOTO_CA is only its fallback.
*****************************************************************************/

#ifndef PHOTOSYNTHESIS_H
#define PHOTOSYNTHESIS_H

/* ========================================================================= */
/* Physical constants                                                        */
/* ========================================================================= */

#define PHOTO_RGAS          8.3145   /* universal gas constant (J/mol/K)      */
#define PHOTO_TFRZ          273.15   /* 0 degC in K                           */
#define PHOTO_P0            101325.0 /* standard sea-level pressure (Pa)      */
#define PHOTO_VONKARMAN     0.4

/* Photon flux per unit visible energy flux.  4.57 umol/J is the standard
   conversion for the 400-700 nm waveband under daylight. */
#define PHOTO_PAR_PER_WATT  4.57

#define PHOTO_PA_TO_KPA     0.001
#define PHOTO_KPA_TO_PA     1000.0

/* Diffusivity ratios of H2O to CO2 */
#define PHOTO_H2O_CO2_STOM  1.6      /* stomatal pore, molecular diffusion    */
#define PHOTO_H2O_CO2_BL    1.37     /* leaf boundary layer, turbulent        */

/* ========================================================================= */
/* Rubisco kinetics at 25 degC -- [BER], CLM tabulation                      */
/* ========================================================================= */

#define PHOTO_KC25          41.03    /* M-M constant for CO2 (Pa)             */
#define PHOTO_KO25          28210.0  /* M-M constant for O2 (Pa)              */
#define PHOTO_GAMMASTAR25   4.332    /* CO2 compensation pt w/o Rd (Pa)       */
#define PHOTO_O2_MOLFRAC    0.209    /* O2 mole fraction of dry air (-)       */

#define PHOTO_EA_KC         79430.0  /* activation energy, Kc (J/mol)   [BER] */
#define PHOTO_EA_KO         36380.0  /* activation energy, Ko (J/mol)   [BER] */
#define PHOTO_EA_GAMMASTAR  37830.0  /* activation energy, Gamma* (J/mol)[BER]*/
#define PHOTO_EA_RD         46390.0  /* activation energy, Rd (J/mol)   [BER] */

/* ========================================================================= */
/* Vcmax / Jmax temperature response                                         */
/*   form   [MED02]                                                          */
/*   values [KUM19] via plantecophys >= 1.4 defaults                         */
/* ========================================================================= */

#define PHOTO_EA_VCMAX      58550.0  /* EaV   (J/mol)                         */
#define PHOTO_HD_VCMAX      200000.0 /* EdVC  (J/mol)                         */
#define PHOTO_DS_VCMAX      629.26   /* delsC (J/mol/K)                       */
#define PHOTO_EA_JMAX       29680.0  /* EaJ   (J/mol)                         */
#define PHOTO_HD_JMAX       200000.0 /* EdVJ  (J/mol)                         */
#define PHOTO_DS_JMAX       631.88   /* delsJ (J/mol/K)                       */

/* Jmax25:Vcmax25 DEFAULT.  Temperate value; [MED02] review; Kattge & Knorr
   (2007) give 2.59 - 0.035*Tgrowth.  Overridden per class by the config key
   JMAX VCMAX RATIO (PHOTOTRAIT.JmaxRatio). */
#define PHOTO_JMAXRATIO     1.70

/* Rd25 as a fraction of Vcmax25, DEFAULT.  [COL91] give 0.015 for C3.
   Overridden per class by RD VCMAX RATIO (PHOTOTRAIT.Rd25Ratio). */
#define PHOTO_RD25RATIO     0.015

/* How the water-stress factor beta enters the leaf model.
     0  beta multiplies Vcmax (and therefore An); Jmax and Rd unstressed.
        This is CLM4.5 BTRAN and CLM5-PHS [KEN19] Eqn 2.  DEFAULT.
     1  beta multiplies the Medlyn slope 1.6(1+g1/sqrt(D)) [DEK15].
   Both are cited forms; 0 is what the option matrix in docs/provenance.md
   sec. 5 advertises (CLM4.5 BTRAN, CLM5-PHS). */
#define PHOTO_BETA_MODE     0

/* ========================================================================= */
/* Electron transport                                                        */
/*                                                                           */
/* Non-rectangular hyperbola.  ALPHA is the quantum yield of electron        */
/* transport per absorbed photon and THETA the curvature.                    */
/*                                                                           */
/* PROVENANCE (docs/provenance.md item 3, resolved): 0.24 and 0.85 are the    */
/* plantecophys::Photosyn() defaults (alpha = 0.24, theta = 0.85), i.e. the   */
/* same source family as the temperature response above.  plantecophys does   */
/* not cite a primary paper for either; 0.24 is the Medlyn et al. (2002)      */
/* quantum yield 0.3 reduced for ~0.8 absorptance, 0.85 a common curvature.   */
/* ========================================================================= */

#define PHOTO_ALPHA         0.24     /* mol e- / mol absorbed photon [plantecophys default] */
#define PHOTO_THETA         0.85     /* curvature of J vs PAR (-)    [plantecophys default] */

/* Product/export-limited rate for C3:  Ae = k Vcmax,  k = 0.5.
   [COL91]; reproduced verbatim in [ELL20] Notes S2 Eqn S2.3c ("k is a
   constant equal to 0.5 for C3 plants").  The pre-Phase-II file carried this
   with the comment "Noah-MP's WE = 0.5*Vcmax for C3" -- Noah-MP does use it,
   but the source is Collatz, and citing Noah-MP hid that.
   Ae is INDEPENDENT of Ci, which is what makes the SOX co-limitation Ci
   (Eqn S2.6b) well posed; dropping it changes ci,col materially. */
#define PHOTO_KEXPORT       0.5

/* Co-limitation smoothing, [COL91] quadratic form, applied in two stages:
     stage 1   Ac with Aj
     stage 2   that result with Ae
   Values are the CLM4.5/CLM5 defaults.  JULES uses 0.83 and 0.93
   (Jones et al. 2022); the choice shifts A by a few percent near the
   transitions and is a legitimate tuning point, so both are recorded.

   BEHAVIOUR CHANGE, PHASE II.  The pre-Phase-II Photosynthesis.c took a hard
   three-way minimum (the theta -> 1 limit) while PlantHydraulics.c smoothed
   at 0.98.  The two paths therefore disagreed with each other, and unifying
   them necessarily moves one of them.  Smoothing is the citable form. */
/* DEFAULT IS 1.0 = HARD MINIMUM.  Measured with the Phase 0 harness: at
   theta = 1 the rebuilt kernel reproduces the pre-Phase-II Medlyn path to
   0.001% median / 0.06% max in Rc, i.e. Phase II is behaviour-neutral.  At
   the CLM5 values (0.98, 0.95) canopy resistance shifts by 11% median /
   47% max and transpiration by 4.7% median / 29% max, which would silently
   invalidate existing calibration and would contaminate every Phase III-IV
   diff.  Smoothing is the better-supported form and should be adopted -- but
   as its own measured experiment in Phase VI, not as a side effect here. */
#define PHOTO_COLIM1        1.00     /* set 0.98 for CLM5, 0.83 for JULES     */
#define PHOTO_COLIM2        1.00     /* set 0.95 for CLM5, 0.93 for JULES     */

/* ========================================================================= */
/* Ambient conditions and numerical control                                  */
/* ========================================================================= */

#define PHOTO_CA            420.0    /* atmospheric CO2 (umol/mol): FALLBACK
                                        only.  The model reads [CONSTANTS]
                                        ATMOSPHERIC CO2 into ATMOS_CO2 and
                                        passes it in as Ca everywhere.      */
#define PHOTO_GB            2.0      /* leaf boundary layer conductance to
                                        H2O (mol/m2/s), when not computed     */

#define PHOTO_CI_TOL        1.0e-6   /* Ci convergence tolerance, relative     */
#define PHOTO_MAX_ITER      60       /* bisection cap on the Ci solve         */
#define PHOTO_AGC_NEWTON    4        /* Newton steps refining A(Gc)           */
#define PHOTO_CICOL_CAP     3.0      /* ci,col bound, as a multiple of Ca      */
#define PHOTO_TINY          1.0e-10
#define PHOTO_MIN_GS        1.0e-6   /* mol H2O/m2/s                          */
#define PHOTO_MAX_GS        2.0      /* mol H2O/m2/s, sanity cap              */

/* Cold acclimation.  [MAK04] delayed-temperature state of acclimation,
   dS/dt = (Tair - S)/tau, with tau = 200 h inside the 190-330 h range fitted
   by Makela et al. (2004) and Kolari et al. (2007) for Scots pine.  The
   linear ramp T0..T1 is this model's reading of their capacity-vs-S curve
   (item 4 closed as a documented departure: the end points are ours).  The
   state variable is advanced in MassEnergyBalance.c; only the thresholds
   live here. */
#define PHOTO_T0            -4.0     /* lower temperature (degC)              */
#define PHOTO_T1            6.0      /* upper temperature (degC)              */
#define PHOTO_TAU           720000.0 /* time constant (s), ~200 h             */

/* ========================================================================= */
/* Two-leaf canopy -- [DPF97]                                                */
/* ========================================================================= */

#define PHOTO_G_SPHERICAL   0.5      /* leaf projection, spherical angle dist */
#define PHOTO_KD_DIFFUSE    0.78     /* diffuse PAR extinction, [DPF97]       */
#define PHOTO_KN_NITROGEN   0.30     /* canopy N / Vcmax extinction coeff     */
#define PHOTO_MIN_SINALT    0.02     /* floor on sin(solar altitude)          */

/* ========================================================================= */
/* Leaf energy balance -- [CN98]                                             */
/* ========================================================================= */

#define PHOTO_LEAF_WIDTH    0.01     /* characteristic leaf width (m), DEFAULT;
                                        per class via LEAF WIDTH             */
#define PHOTO_LEAF_SIDES    2.0      /* leaf faces exchanging sensible heat.
                                        [CN98] Table 7.6 is per side.        */
#define PHOTO_NIR_PER_PAR   0.30     /* absorbed NIR : absorbed PAR (energy).
                                        Leaf absorptance ~0.85 PAR, ~0.25 NIR
                                        for equal incident energy -> 0.25/0.85.
                                        Departure D3.                         */
#define PHOTO_D_FACTOR      0.72     /* d = 0.72 * leaf width, [CN98] Tab 7.6 */
#define PHOTO_BL_COEF       0.135    /* gHa = 0.135*sqrt(u/d), [CN98] Tab 7.6 */
#define PHOTO_BL_TURB       1.4      /* outdoor turbulence enhancement [CN98] */
#define PHOTO_EMISSIVITY    0.97
#define PHOTO_STEFAN        5.67e-8  /* W/m2/K4                               */
#define PHOTO_CP_MOLAR      29.3     /* heat capacity of dry air (J/mol/K)    */
#define PHOTO_WIND_MIN      0.1      /* m/s, floor for the boundary layer     */
#define PHOTO_TLEAF_TOL     1.0e-4   /* degC                                  */
#define PHOTO_TLEAF_SPAN    25.0     /* degC bracket half-width about Tair    */
#define PHOTO_MIN_VPD       0.05     /* kPa                                   */

/* ========================================================================= */
/* Temperature-dependent kinetics, assembled once per (Tleaf, PAR, Vcmax25)  */
/* ========================================================================= */

typedef struct {
  float Vcmax;      /* umol/m2/s, at Tleaf, after dormancy scaling           */
  float Jmax;       /* umol/m2/s                                             */
  float J;          /* umol e-/m2/s, after the light response                */
  float Ae;         /* umol/m2/s, export-limited rate, k*Vcmax  [COL91]      */
  float Rd;         /* umol/m2/s                                             */
  float Kc;         /* Pa                                                    */
  float Ko;         /* Pa                                                    */
  float Km;         /* Pa,  Kc*(1 + O/Ko)                                    */
  float GammaStar;  /* Pa                                                    */
  float KmMol;      /* umol/mol, Km converted at the ambient pressure        */
  float GStarMol;   /* umol/mol                                              */
  float Press;      /* Pa, the pressure the mole-fraction forms assume       */
} PHOTOKIN;

/* Per-class photosynthetic traits.  Fill with PhotoTraitDefaults() and then
   override what the vegetation table supplies. */
typedef struct {
  float Vcmax25;    /* umol/m2 leaf/s, top of canopy                         */
  float JmaxRatio;  /* Jmax25 / Vcmax25 (-)                                  */
  float Rd25Ratio;  /* Rd25 / Vcmax25 (-)                                    */
  float LeafWidth;  /* characteristic leaf dimension (m)                     */
} PHOTOTRAIT;

/* ========================================================================= */
/* API                                                                       */
/* ========================================================================= */

void PhotoTraitDefaults(PHOTOTRAIT *P, float Vcmax25);

/* Conductance unit conversion, molar (mol/m2/s) -> velocity (m/s). */
float PhotoMolarToVelocity(float Tair, float Press);

/* Saturation vapour pressure (kPa) and its slope (kPa/K).
   Argument is Tc, not T: brent.h does  #define T 1e-5  and is included ahead
   of massenergy.h in three DHSVM source files. */
float PhotoSatVaporPressure(float Tc);
float PhotoSatVaporSlope(float Tc);

/* Assemble the temperature- and light-dependent kinetics.
   Dormancy in [0,1] scales Vcmax and Jmax but NOT Rd.  Beta in [0,1] is the
   water-stress factor and, with PHOTO_BETA_MODE 0, scales Vcmax only
   ([KEN19] Eqn 2); with mode 1 it is ignored here and applied to the slope
   in the Medlyn solve. */
void PhotoKineticsTrait(const PHOTOTRAIT *P, float Tleaf, float ParAbs,
                        float Dormancy, float Beta, float Press, PHOTOKIN *K);

/* Default-trait, unstressed wrapper (kept for the harness and old callers). */
void PhotoKinetics(float Tleaf, float ParAbs, float Vcmax25, float Dormancy,
                   float Press, PHOTOKIN *K);

/* Gross assimilation at a given intercellular CO2 mole fraction.
   Ci in umol/mol.  Returns Ag (umol/m2/s); Ac and Aj optional. */
float PhotoAssimilationAtCi(const PHOTOKIN *K, float Ci,
                            float *AcOut, float *AjOut);

/* Gross assimilation and its analytic slope dAg/dCi at the same Ci. */
void PhotoAgAndSlope(const PHOTOKIN *K, float Ci, float *AgOut,
                     float *SlopeOut);

/* NET assimilation at a known TOTAL leaf conductance to CO2.
   Gc in mol CO2/m2/s (stomatal and boundary layer already in series),
   Ca in umol/mol.  Analytic; no iteration.

   The diffusion balance closes on NET assimilation:  An = Gc (Ca - Ci).
   This is the correction to the pre-Phase-II HydAssimilationAtGc(), which
   solved its quadratics without Rd and then set Ci = Ca - Agross/Gc, drawing
   Ci down by an extra Rd/Gc.  The error was negligible in bright light and
   large at low conductance -- about 18 umol/mol of Ci at PAR = 100,
   T = 25 degC.  Used by ProfitMax and by the SOX A(gs) evaluation. */
float PhotoAssimilationAtGc(const PHOTOKIN *K, float Gc, float Ca,
                            float *CiOut);

/* Ci at co-limitation: the Ci above which further increases no longer raise
   assimilation because it has become light-limited.  [ELL20] Eqn S2.6b,
   translated from Collatz to FvCB symbols (departure D1). */
float PhotoCiColimit(const PHOTOKIN *K, float Ca);

/* dA/dCi as the secant between Ca and the co-limitation Ci.
   [ELL20] Eqn S2.1.  Units umol CO2 /m2/s per umol/mol.
   Returns 0 when the leaf is already light-limited at Ca, which is the
   signal SOX uses to stop opening stomata (Eqn S2.2). */
float PhotoDAdCi(const PHOTOKIN *K, float Ca, float *CiColOut, float *AColOut);

/* Coupled FvCB + Medlyn solve at fixed leaf temperature.
   Vpd is LEAF-TO-AIR (kPa).  Beta in [0,1] is the water-stress factor (see
   PHOTO_BETA_MODE).  Any output may be NULL. */
void PhotoLeafFluxTrait(const PHOTOTRAIT *P, float G1, float G0, float ParAbs,
                        float Tleaf, float Vpd, float Ca, float Press,
                        float Gb, float Beta, float Dormancy,
                        float *An, float *Gs, float *Ci);

/* Default-trait wrapper. */
void PhotoLeafFlux(float Vcmax25, float G1, float G0, float ParAbs,
                   float Tleaf, float Vpd, float Ca, float Press, float Gb,
                   float Beta, float Dormancy,
                   float *An, float *Gs, float *Ci);

/* Leaf boundary layer conductance to heat (mol/m2/s).  [CN98] Table 7.6. */
float PhotoBoundaryLayerH(float Wind, float LeafWidth);

/* Absorbed radiation of ONE leaf (W/m2 LEAF) for the energy balance:
   absorbed PAR (umol/m2 leaf/s) scaled to total shortwave with
   PHOTO_NIR_PER_PAR, plus the isothermal longwave term eps*sigma*Tair^4.
   Departure D3.  This -- not the canopy net radiation -- is what
   PhotoLeafEnergyBalance() expects as Rabs. */
float PhotoLeafRabs(float ParAbsLeaf, float Tair);

/* Leaf temperature and leaf-to-air VPD from the energy balance.
   E is transpiration (mol H2O/m2 leaf/s), Vpd is AIR vapour pressure deficit
   (kPa), Rabs is the leaf's absorbed radiation from PhotoLeafRabs() (W/m2
   leaf), LeafWidth in m.  Solved against an explicit residual, so the result
   can be checked by evaluating the balance. */
void PhotoLeafEnergyBalanceW(float E, float Tair, float Vpd, float Rabs,
                             float Wind, float Press, float LeafWidth,
                             float *Tleaf, float *VpdLeaf);

/* Default-leaf-width wrapper. */
void PhotoLeafEnergyBalance(float E, float Tair, float Vpd, float Rabs,
                            float Wind, float Press,
                            float *Tleaf, float *VpdLeaf);

/* Residual of the leaf energy balance (W/m2 leaf), exposed for verification:
   Rabs - eps sigma Tl^4 - PHOTO_LEAF_SIDES cp gHa (Tl - Ta) - lambda E. */
float PhotoEnergyResidualW(float Tleaf, float E, float Tair, float Rabs,
                           float Wind, float Press, float LeafWidth);
float PhotoEnergyResidual(float Tleaf, float E, float Tair, float Vpd,
                          float Rabs, float Wind, float Press);

/* Two-leaf sunlit/shaded partitioning.  [DPF97]
   ParBeam and ParDiff are absorbed visible flux (W/m2 GROUND);
   ParSun and ParSha come back as absorbed PAR (umol/m2 LEAF/s).  Beam goes
   to sunlit leaves; diffuse is shared using the [DPF97] sunlit fraction of
   the diffuse absorption profile, kd/(kd+kb) (1-e^-(kd+kb)L)/(1-e^-kd L). */
void PhotoTwoLeafPartition(float Lai, float SinAlt, float ParBeam,
                           float ParDiff, float Vcmax25,
                           float *LaiSun, float *LaiSha,
                           float *ParSun, float *ParSha,
                           float *VcmaxSun, float *VcmaxSha);

/* Canopy conductance (m/s) from the Medlyn scheme, ready to be inverted into
   the resistance EvapoTranspiration() expects.  AnCanopy optional. */
float PhotoCanopyConductance(float Vcmax25, float G1, float G0, float Lai,
                             float SinAlt, float ParBeam, float ParDiff,
                             float Tair, float Vpd, float Ca, float Press,
                             float Gb, float Beta, float Dormancy,
                             float *AnCanopy);

#endif
