
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "massenergy.h"
#include "constants.h"
#include "photosynthesis.h"
#include "planthydraulics.h"
#include "roothydraulics.h"
#include "stomatalscheme.h"
#include "data.h"
#include "functions.h"

/*****************************************************************************
  CanopyResistance()
*****************************************************************************/
float CanopyResistance(float LAI, float RsMin, float RsMax, float Rpc,
  float VpdThres, float MoistThres, float WP,
  float TSoil, float SoilMoisture, float Vpd, float Rp)
{
  float MoistFactor;	/* multiplier for resistance due to soil moisture feed-back */
  float Resistance;		/* Canopy resistance (s/m) */
  float RpFactor;		/* multiplier for resistance due to light level feed-back */
  float TFactor;		/* multiplier for resistance due to soil temperaure feed-back */
  float VpdFactor;		/* multiplier for resistance due to vapor pressure deficit feed-back */

  if (TSoil <= 0) {
    Resistance = DHSVM_HUGE;
    return Resistance;
  }

  /* for OBS */
  TFactor = 1.0 / (0.176 + 0.0770 * TSoil - 0.0018 * TSoil * TSoil);

  if (TFactor <= 0) {
    Resistance = DHSVM_HUGE;
    return Resistance;
  }

  /* equation 14, Wigmosta et al [1994] */

  if (Vpd >= VpdThres) {
    Resistance = DHSVM_HUGE;
    return Resistance;
  }
  else
    VpdFactor = 1.0 / (1 - Vpd / VpdThres);

  /* equation 15, Wigmosta et al [1994 */

  RpFactor = 1.0 / ((RsMin / RsMax + Rp / Rpc) / (1 + Rp / Rpc));

  /* equation 16, Wigmosta et al [1994] */

  if (SoilMoisture <= WP) {
    Resistance = DHSVM_HUGE;
    return Resistance;
  }
  else if (SoilMoisture < MoistThres)
    MoistFactor = (MoistThres - WP) / (SoilMoisture - WP);
  else
    MoistFactor = 1.0;

  Resistance = TFactor * VpdFactor * RpFactor * MoistFactor * RsMin / LAI;

  return Resistance;

}

/*****************************************************************************
  Phase V.  CanopyResistancePhoto() and CanopyResistanceSperry() are replaced
  by the single CanopyResistanceScheme() below, dispatching on the two
  orthogonal options.  Jarvis (Wigmosta et al. 1994) above is unchanged.

  This is the seam: DHSVM units above (m/s, m of head, volumetric moisture),
  physiology units below (MPa, mol/m2/s, mmol/m2/s).  Conversion lives here.

  Phase-VI refactor, what changed at the seam:
  * Each big leaf gets its own absorbed radiation (PhotoLeafRabs) for the
    leaf energy balance; the canopy net radiation is no longer handed to the
    leaf (photosynthesis.h departure D3).
  * Leaf-specific hydraulic conductances (KxMax, Krs, Kcomp, per m2 leaf) are
    multiplied by LAI for the ground-basis root solve, so RootSolve() sees the
    same plant the supply curves describe.
  * Atmospheric CO2 comes from [CONSTANTS] ATMOSPHERIC CO2.
  * Rc is floored at RsMin/LAI, the DHSVM convention, and the stomatal cap
    Gmax (MAXIMUM STOMATAL CONDUCTANCE, default 1/RsMin) is passed to the
    schemes.
  * Out->EConsistent tells EvapoTranspiration() when Penman-Monteith must
    reproduce the scheme's own transpiration (CanopyResistanceFromFlux).
*****************************************************************************/

/* ========================================================================= */
/* Unit conversions, DHSVM <-> physiology                                    */
/* ========================================================================= */

#define CR_M_TO_MPA      (9810.0f * 1.0e-6f)   /* m of head -> MPa          */
#define CR_MMOL_TO_MS(T,P) \
  (PhotoMolarToVelocity((T), (P)) * 0.001f)    /* mmol/m2/s -> m/s of vapour */

/* Transpiration mmol H2O /m2/s  ->  m/s of liquid water.
   1 mol H2O = 18.015 g = 18.015e-6 m3. */
#define CR_MMOL_TO_MPS   (18.015e-9f)

/* Perirhizal fixed point between the stomatal and root solves: maximum
   number of stomatal re-evaluations per timestep, and the change in Heff
   (MPa) below which it is considered settled. */
#define CR_PERIRHIZAL_PASSES 3
#define CR_HEFF_TOL          0.02f

/* DHSVM stores bubbling pressure as a positive head in m; the physiology
   wants a negative air-entry potential in MPa. */
static float CrHeadToMPa(float HeadM)
{
  float h = (HeadM > 0.0f) ? HeadM : 0.01f;
  return -h * CR_M_TO_MPA;
}

/* Brooks-Corey, matching the exponent DHSVM already uses in
   UnsaturatedFlow.c:  psi = psi_ae (theta/theta_sat)^(-1/lambda). */
static float CrSoilPotential(float Moist, float Porosity, float PsatMPa,
                             float Lambda)
{
  float Sat;

  if (Porosity < 1.0e-9f || Lambda < 1.0e-9f) return PsatMPa;

  Sat = Moist / Porosity;
  if (Sat > 1.0f) Sat = 1.0f;
  if (Sat < 0.001f) Sat = 0.001f;

  if (PsatMPa > 0.0f) PsatMPa = -PsatMPa;

  return PsatMPa * (float)pow((double)Sat, -1.0 / (double)Lambda);
}

/* ========================================================================= */
/* Root-zone water potential -- SCHEME-INDEPENDENT                           */
/* ========================================================================= */

/*****************************************************************************
  CanopyRootZonePotential()

  Root-fraction-weighted soil water potential over the root zone (MPa).

  This is a property of the SOIL STATE and the root distribution, not of the
  stomatal scheme: Brooks-Corey applied to the same moisture DHSVM already
  tracks.  It is therefore defined for Jarvis and for the empirical-beta
  paths just as much as for the hydraulic ones, and it is what makes those
  runs comparable.

  Volumetric moisture is not comparable across soil types or depths, and
  Jarvis and Medlyn-with-beta both limit transpiration on volumetric moisture
  relative to MoistThres/WP.  Expressing the same state as a potential puts
  every scheme on one axis, so "how much does this scheme transpire at a
  given soil water potential" becomes a question that can actually be asked.

  With hydraulics on, the seam overwrites this with the SUF-weighted
  effective potential (Vanderborght et al. 2021 Eqn 7).  The two are
  identical whenever SUF defaults to RootFract, which is the shipped
  configuration; they diverge only if SUF is configured separately.
*****************************************************************************/
float CanopyRootZonePotential(int Layer, VEGTABLE *VType, SOILTABLE *SType,
                              SOILPIX *LocalSoil, float *Moist)
{
  float Psi = 0.0f, Wt = 0.0f, Psat, f;
  int i, n;

  n = VType->NSoilLayers;
  if (n > ROOT_MAXLAYERS) n = ROOT_MAXLAYERS;

  for (i = 0; i < n; i++) {
    f = VType->RootFract[Layer][i];
    if (f <= 0.0f) continue;
    Psat = CrHeadToMPa(SType->Press[i]);
    Psi += f * CrSoilPotential(Moist[i], LocalSoil->Porosity[i], Psat,
                               SType->PoreDist[i]);
    Wt += f;
  }

  if (Wt <= 0.0f) return 0.0f;
  Psi /= Wt;

  /* Brooks-Corey runs to -infinity as saturation -> 0; bound it the same way
     RootEffectivePotential() does so the two agree. */
  if (Psi < (float)ROOT_PSI_FLOOR) Psi = (float)ROOT_PSI_FLOOR;

  return Psi;
}

/* ========================================================================= */
/* Root zone assembly                                                        */
/* ========================================================================= */

static void CrBuildRootZone(ROOTZONE *R, int Layer, VEGTABLE *VType,
                            SOILTABLE *SType, SOILPIX *LocalSoil,
                            float *Moist)
{
  float RootFract[ROOT_MAXLAYERS];
  float Psat;
  int i, n;

  n = VType->NSoilLayers;
  if (n > ROOT_MAXLAYERS) n = ROOT_MAXLAYERS;

  RootZoneInit(R, n);

  for (i = 0; i < n; i++) {
    Psat = CrHeadToMPa(SType->Press[i]);
    R->Hs[i]     = CrSoilPotential(Moist[i], LocalSoil->Porosity[i], Psat,
                                   SType->PoreDist[i]);
    /* Perirhizal stage: matrix-scale Ks from RHIZOSPHERE CONDUCTIVITY when
       given, else the model's (drainage-scale) KsVert. */
    R->KsVert[i] = (SType->KsRhizo != NULL && SType->KsRhizo[i] > 0.0f)
                   ? SType->KsRhizo[i] : LocalSoil->KsVert[i];
    R->Lambda[i] = SType->PoreDist[i];
    R->HaeM[i]   = SType->Press[i];
    RootFract[i] = VType->RootFract[Layer][i];
  }

  /* SUF defaults to the prescribed root fractions.  With Kcomp = 0 this
     reduces EXACTLY to DHSVM's pre-Phase-III uptake, which is the fallback
     the harness checks (suite phase3, R2). */
  RootSetSUFFromFractions(R, RootFract);

  R->Krs   = VType->Krs[Layer];
  R->Kcomp = VType->Kcomp[Layer];
  R->UsePerirhizal = 0;                  /* Phase V: simplest version       */
  R->HcollarMin    = VType->PsiCollarMin[Layer];
}

/* ========================================================================= */
/* The scheme dispatcher                                                     */
/* ========================================================================= */

/*****************************************************************************
  CanopyResistanceScheme()

  Returns canopy resistance (s/m) for the Penman-Monteith term in
  EvapoTranspiration().  Fills Out with the diagnostics DHSVM reports and,
  critically, with Efrac (per-layer uptake shares) and Tsupply (the
  transpiration the soil can actually deliver).

  StomScheme:  MEDLYN | PROFITMAX | PROFITMAX2 | SOX
  Hydraulics:  HYD_NONE | HYD_KRSSUF   (MEDLYN is the only scheme valid with
                                        HYD_NONE; validated at startup)
*****************************************************************************/
float CanopyResistanceScheme(int StomScheme, int Hydraulics,
  float Lai, int Layer, float Rp, float NetRad,
  PIXMET *Met, VEGTABLE *VType, SOILTABLE *SType, SOILPIX *LocalSoil,
  float *Moist, float Dormancy, CANOPYHYD *Out)
{
  ROOTZONE Root;
  ROOTSOLUTION RSol;
  STOMTRAIT Trait;
  STOMMET M;
  STOMSOLUTION Sun, Sha;
  HYDSUPPLY Supply;
  float FBeam, ParBeam, ParDiff;
  float LaiSun, LaiSha, ParSun, ParSha, VcSun, VcSha;
  float Heff, GwTot, AnTot, ETot, GsMs, Rc, RcFloor;
  float Beta, BetaLayer, Sum, KrsLeaf;
  int i, n, Solver, Pass;

  n = VType->NSoilLayers;
  if (n > ROOT_MAXLAYERS) n = ROOT_MAXLAYERS;

  /* Inert defaults: any early return leaves the caller with a safe state. */
  memset(Out, 0, sizeof(*Out));
  Out->NLayers = n;
  Out->Rc      = VType->RsMax[Layer];
  Out->Tleaf   = Met->Tair;
  Out->VpdLeaf = Met->Vpd * PHOTO_PA_TO_KPA;
  Out->Tsupply = -1.0f;                  /* < 0 means "no supply limit"     */
  for (i = 0; i < n; i++) Out->Efrac[i] = VType->RootFract[Layer][i];

  if (Lai <= 0.0f || VType->Vcmax25[Layer] <= 0.0f)
    return VType->RsMax[Layer];

  /* --- meteorology, in physiology units -------------------------------- */
  memset(&M, 0, sizeof(M));
  M.Tair   = Met->Tair;
  M.VpdAir = Met->Vpd * PHOTO_PA_TO_KPA;
  M.Press  = Met->Press;
  M.Ca     = (ATMOS_CO2 > 0.0f) ? ATMOS_CO2 : (float)PHOTO_CA;
  M.Wind   = Met->Wind;
  /* M.Rabs is set per big leaf below (PhotoLeafRabs).  NetRad -- the canopy
     net radiation per m2 ground, already net of longwave emission -- stays
     with Penman-Monteith in EvapoTranspiration(); handing it to a single
     leaf that then subtracts its own emission was the pre-refactor bug. */
  (void)NetRad;

  /* --- two-leaf partition, shared by every scheme ----------------------- */
  /* Layer 0 is the topmost canopy and sees the direct beam; layer 1 sits
     under an overstory where transmitted light is nearly all diffuse. */
  if (Layer == 0 && (Met->SinBeam + Met->SinDiffuse) > 0.0f)
    FBeam = Met->SinBeam / (Met->SinBeam + Met->SinDiffuse);
  else
    FBeam = 0.0f;

  ParBeam = Rp * FBeam;
  ParDiff = Rp * (1.0f - FBeam);

  PhotoTwoLeafPartition(Lai, Met->SineSolarAltitude, ParBeam, ParDiff,
    VType->Vcmax25[Layer], &LaiSun, &LaiSha, &ParSun, &ParSha,
    &VcSun, &VcSha);

  memset(&Trait, 0, sizeof(Trait));
  Trait.G1           = VType->G1[Layer];
  Trait.G0           = VType->G0[Layer];
  Trait.Dormancy     = Dormancy;
  /* Leaf boundary-layer conductance to water vapour from the same wind and
     leaf width the energy balance uses ([CN98] Table 7.6, x1.08 for vapour
     vs heat), rather than the fixed PHOTO_GB. */
  Trait.Gb           = 1.08f * PhotoBoundaryLayerH(Met->Wind,
                                                   VType->LeafWidth[Layer]);
  if (Trait.Gb < 0.1f) Trait.Gb = 0.1f;
  Trait.CanopyHeight = VType->Height[Layer];   /* informational; gravity is
                                                  in VType->Xylem.Pgrav     */
  Trait.RpMin        = 0.0f;                   /* derived from Krs, KxMax   */
  Trait.Gmax         = VType->Gmax[Layer];     /* mol/m2 leaf/s; 1/RsMin by
                                                  default                   */
  Trait.JmaxRatio    = VType->JmaxRatio[Layer];
  Trait.Rd25Ratio    = VType->Rd25Ratio[Layer];
  Trait.LeafWidth    = VType->LeafWidth[Layer];

  /* Leaf-specific root conductance (mmol/m2 leaf/s/MPa) feeds the per-leaf
     supply curves; the ground-basis root solve gets Krs x LAI below. */
  KrsLeaf = VType->Krs[Layer];

  GwTot = 0.0f; AnTot = 0.0f; ETot = 0.0f;

  /* Initialize BOTH leaf solutions before any branch runs.  At night
     PhotoTwoLeafPartition() returns LaiSun = 0, so the sunlit call is
     skipped -- and anything read out of Sun afterwards is uninitialized
     stack.  That produced leaf potentials of ~5e22 MPa in the first output
     with these diagnostics, on exactly the hours with zero radiation. */
  StomInit(&Sun, &M);
  StomInit(&Sha, &M);

  /* ===================================================================== */
  /* MEDLYN without hydraulics: the empirical beta path                     */
  /* ===================================================================== */
  if (Hydraulics == HYD_NONE) {

    /* Report the root-zone soil water potential even though this path does
       not use it.  EvapoTranspiration() computes it before the branch and
       then copies Out->Heff over the top, so without this the empirical-beta
       runs wrote a flat zero into PsiSoil while Jarvis, which never enters
       this function, kept the good value.  PsiSoil is the one axis on which
       the schemes can be compared, and 0 MPa reads as saturated. */
    Out->Heff = CanopyRootZonePotential(Layer, VType, SType, LocalSoil,
                                        Moist);

    Beta = 0.0f;
    for (i = 0; i < n; i++) {
      BetaLayer = (Moist[i] - SType->WP[i]) /
                  (VType->MoistThres[Layer] - SType->WP[i]);
      if (BetaLayer > 1.0f) BetaLayer = 1.0f;
      if (BetaLayer < 0.0f) BetaLayer = 0.0f;
      Out->Efrac[i] = VType->RootFract[Layer][i] * BetaLayer;
      Beta += Out->Efrac[i];
    }
    Sum = Beta;
    for (i = 0; i < n; i++)
      Out->Efrac[i] = (Sum > 0.0f) ? Out->Efrac[i] / Sum
                                   : VType->RootFract[Layer][i];
    Out->Beta = Beta;

    /* CLM4.5-style empirical stress: Beta attenuates Vcmax inside the leaf
       model (PHOTO_BETA_MODE 0). */
    if (LaiSun > 1.0e-6f) {
      Trait.Vcmax25 = VcSun; M.ParAbs = ParSun;
      M.Rabs = PhotoLeafRabs(ParSun, Met->Tair);
      StomMedlyn(&Trait, &M, Beta, &Sun);
      GwTot += Sun.Gw * LaiSun; AnTot += Sun.An * LaiSun;
      ETot  += Sun.E  * LaiSun;
    }
    if (LaiSha > 1.0e-6f) {
      Trait.Vcmax25 = VcSha; M.ParAbs = ParSha;
      M.Rabs = PhotoLeafRabs(ParSha, Met->Tair);
      StomMedlyn(&Trait, &M, Beta, &Sha);
      GwTot += Sha.Gw * LaiSha; AnTot += Sha.An * LaiSha;
      ETot  += Sha.E  * LaiSha;
    }
  }

  /* ===================================================================== */
  /* Hydraulic schemes                                                      */
  /* ===================================================================== */
  else {

    CrBuildRootZone(&Root, Layer, VType, SType, LocalSoil, Moist);

    /* The root model works per m2 GROUND (it receives the canopy total
       transpiration), the supply curves per m2 LEAF: a canopy of LAI leaves
       in parallel, each on its own leaf-specific pathway, has a ground-basis
       root system conductance of Krs x LAI.  Pre-refactor the leaf value was
       used on the ground basis, which understated the supply limit by a
       factor LAI. */
    Root.Krs   = KrsLeaf * Lai;
    Root.Kcomp = VType->Kcomp[Layer] * Lai;

    /* --- perirhizal (rhizosphere) stage, Leitner et al. 2025 Eqns 21-27 --
       On when the class has a root length index.  Geometry per soil layer
       from the root length density (root length index x root fraction /
       layer thickness) and the fine root radius; the root radial
       conductivity is the one implied by Krs itself, so the interface
       balance and the macroscopic root model describe the same roots.  Bulk
       soil potential Hs stays what it was; what the plant now sees is the
       soil-root INTERFACE potential Hsr, which falls below Hs as the
       perirhizal soil dries and its conductivity collapses. */
    if (VType->RootLengthIndex[Layer] > 0.0f) {
      float Rld[ROOT_MAXLAYERS], Aroot, Kr, Ltot, d;
      Ltot  = VType->RootLengthIndex[Layer];
      Aroot = VType->RootRadius[Layer];
      for (i = 0; i < Root.N; i++) {
        d = VType->RootDepth[i];
        if (d < 0.01f) d = 0.01f;
        Rld[i] = Ltot * VType->RootFract[Layer][i] / d;
      }
      Kr = RootRadialFromKrs(Root.Krs, Aroot, Ltot);
      RootSetPerirhizal(&Root, Rld, Aroot, Kr);
      Root.UsePerirhizal = 1;
    }
    for (i = 0; i < Root.N; i++) Root.Hsr[i] = Root.Hs[i];

    Heff = RootEffectivePotential(&Root);
    Out->Heff = Heff;
    Out->Beta = 1.0f;                     /* no empirical stress factor     */

    /* With the perirhizal stage on, the interface potential depends on the
       uptake and the uptake on the interface potential, so the stomatal
       solve and the root solve are iterated to a fixed point on Heff (at
       most CR_PERIRHIZAL_PASSES passes; one pass, i.e. the old behaviour,
       when the stage is off). */
    for (Pass = 0; Pass < CR_PERIRHIZAL_PASSES; Pass++) {
    GwTot = 0.0f; AnTot = 0.0f; ETot = 0.0f;

    if (StomScheme == MEDLYN) {

      /* CLM5-PHS (Kennedy et al. 2019): Medlyn stomata, beta from leaf water
         potential through the vulnerability curve, attenuating Vcmax
         (Kennedy Eqn 2; PHOTO_BETA_MODE 0).  Conductance is the state
         variable here, so Penman-Monteith runs from Rc (EConsistent = 0). */
      if (LaiSun > 1.0e-6f) {
        Trait.Vcmax25 = VcSun; M.ParAbs = ParSun;
        M.Rabs = PhotoLeafRabs(ParSun, Met->Tair);
        StomMedlynPHS(&Trait, &M, &VType->Xylem[Layer], Heff, KrsLeaf, &Sun);
        if (!Sun.Failed) {
          GwTot += Sun.Gw * LaiSun; AnTot += Sun.An * LaiSun;
          ETot  += Sun.E  * LaiSun;
        }
      }
      if (LaiSha > 1.0e-6f) {
        Trait.Vcmax25 = VcSha; M.ParAbs = ParSha;
        M.Rabs = PhotoLeafRabs(ParSha, Met->Tair);
        StomMedlynPHS(&Trait, &M, &VType->Xylem[Layer], Heff, KrsLeaf, &Sha);
        if (!Sha.Failed) {
          GwTot += Sha.Gw * LaiSha; AnTot += Sha.An * LaiSha;
          ETot  += Sha.E  * LaiSha;
        }
      }
      /* Ecrit depends only on the hydraulics, not on either leaf, so take
         it from the supply function directly (per leaf) and scale to the
         ground basis. */
      Out->Ecrit = HydEcrit(&VType->Xylem[Layer], Heff, KrsLeaf)
                   * Lai * CR_MMOL_TO_MPS;
    }
    else if (StomScheme == PROFITMAX || StomScheme == PROFITMAX2) {

      HydBuildSupply(&VType->Xylem[Layer], Heff, KrsLeaf, &Supply);
      if (Supply.N < 3) return VType->RsMax[Layer];
      Out->Ecrit = Supply.Ecrit * Lai * CR_MMOL_TO_MPS;
      Out->EConsistent = 1;               /* E solved with its own leaf
                                             energy balance                 */

      if (LaiSun > 1.0e-6f) {
        Trait.Vcmax25 = VcSun; M.ParAbs = ParSun;
        M.Rabs = PhotoLeafRabs(ParSun, Met->Tair);
        StomProfitMax(&Trait, &M, &Supply,
                      (StomScheme == PROFITMAX2) ? STOM_PROFITMAX2
                                                 : STOM_PROFITMAX, &Sun);
        if (!Sun.Failed) {
          GwTot += Sun.Gw * LaiSun; AnTot += Sun.An * LaiSun;
          ETot  += Sun.E  * LaiSun;
        }
      }
      if (LaiSha > 1.0e-6f) {
        Trait.Vcmax25 = VcSha; M.ParAbs = ParSha;
        M.Rabs = PhotoLeafRabs(ParSha, Met->Tair);
        StomProfitMax(&Trait, &M, &Supply,
                      (StomScheme == PROFITMAX2) ? STOM_PROFITMAX2
                                                 : STOM_PROFITMAX, &Sha);
        if (!Sha.Failed) {
          GwTot += Sha.Gw * LaiSha; AnTot += Sha.An * LaiSha;
          ETot  += Sha.E  * LaiSha;
        }
      }
    }
    else {   /* SOX */

      /* Numerical iteration is the default.  The semi-analytical form
         (Eller 2020 Eqns 4-5) leaves stomata ~32% open in deep shade where
         the numerical solver closes; see docs/provenance.md sec. 11. */
      Solver = STOM_SOX_NUMERICAL;
      Out->EConsistent = 1;

      if (LaiSun > 1.0e-6f) {
        Trait.Vcmax25 = VcSun; M.ParAbs = ParSun;
        M.Rabs = PhotoLeafRabs(ParSun, Met->Tair);
        StomSOX(&Trait, &M, &VType->Xylem[Layer], Heff, KrsLeaf, Solver,
                &Sun);
        if (!Sun.Failed) {
          GwTot += Sun.Gw * LaiSun; AnTot += Sun.An * LaiSun;
          ETot  += Sun.E  * LaiSun;
        }
      }
      if (LaiSha > 1.0e-6f) {
        Trait.Vcmax25 = VcSha; M.ParAbs = ParSha;
        M.Rabs = PhotoLeafRabs(ParSha, Met->Tair);
        StomSOX(&Trait, &M, &VType->Xylem[Layer], Heff, KrsLeaf, Solver,
                &Sha);
        if (!Sha.Failed) {
          GwTot += Sha.Gw * LaiSha; AnTot += Sha.An * LaiSha;
          ETot  += Sha.E  * LaiSha;
        }
      }
      Out->Ecrit = HydEcrit(&VType->Xylem[Layer], Heff, KrsLeaf)
                   * Lai * CR_MMOL_TO_MPS;
    }

    /* --- uptake distribution and the supply limit ---------------------- */
    /* Without the perirhizal stage RootSolve is called once, at the
       transpiration the canopy scheme predicts: Heff and SUF do not then
       depend on the stomatal solve (Vanderborght et al. 2021).  With it,
       RootSolve returns the interface-weighted Heff for this uptake, and the
       pass loop above re-runs the stomata at that Heff until it settles. */
    RootSolve(&Root, ETot, Root.UsePerirhizal ? Root.Hsr : NULL, &RSol);

    if (!Root.UsePerirhizal) break;
    if ((float)fabs((double)(RSol.Heff - Heff)) < CR_HEFF_TOL) break;
    Heff = RSol.Heff;
    }   /* end of perirhizal pass loop */

    for (i = 0; i < n; i++) Out->Efrac[i] = RSol.Efrac[i];
    Out->PsiRoot = RSol.Hcollar;
    Out->Heff    = RSol.Heff;

    /* If the Dirichlet switch fired, the soil cannot deliver what the
       canopy asked for.  Report the limit in DHSVM units so
       EvapoTranspiration() can cap the Penman-Monteith flux with it. */
    if (RSol.Limited)
      Out->Tsupply = RSol.Tactual * CR_MMOL_TO_MPS;

    /* Canopy-level diagnostics are the LAI-weighted mean of the two big
       leaves, not the sunlit leaf alone.  Same weighting already used for
       conductance and assimilation, and it stays defined when either leaf
       has no area -- which is the whole of the night for the sunlit one. */
    {
      float wSun = (LaiSun > 1.0e-6f && !Sun.Failed) ? LaiSun : 0.0f;
      float wSha = (LaiSha > 1.0e-6f && !Sha.Failed) ? LaiSha : 0.0f;
      float wTot = wSun + wSha;

      if (wTot > 1.0e-6f) {
        Out->PsiLeaf = (Sun.PsiLeaf * wSun + Sha.PsiLeaf * wSha) / wTot;
        Out->Tleaf   = (Sun.Tleaf   * wSun + Sha.Tleaf   * wSha) / wTot;
        Out->VpdLeaf = (Sun.VpdLeaf * wSun + Sha.VpdLeaf * wSha) / wTot;
        if (StomScheme == MEDLYN)
          Out->Beta  = (Sun.Objective * wSun + Sha.Objective * wSha) / wTot;
      }
      else {
        /* Neither leaf produced a solution.  Report the state the plant is
           actually in rather than leaving zeros that read as "unstressed":
           the collar potential, and no hydraulic stress relief. */
        Out->PsiLeaf = RSol.Hcollar;
        Out->Tleaf   = Met->Tair;
        Out->VpdLeaf = Met->Vpd * PHOTO_PA_TO_KPA;
        if (StomScheme == MEDLYN)
          Out->Beta = HydVulnerability(RSol.Hcollar,
                                       &VType->Xylem[Layer].Curve);
      }
    }

    /* Two derived diagnostics, computed here because this is the only place
       that has both the operating potential and the class's own curve.

       PLC is the fraction of xylem conductivity lost at the operating leaf
       potential, the quantity vulnerability curves are measured in.
       SafetyMargin is PsiLeaf - P50: positive means the plant is operating
       on the safe side of the 50% loss point.  Choat et al. (2012) frame the
       global drought-vulnerability picture in exactly this variable, and it
       is what a thinning experiment should move if reduced stand density
       genuinely relieves hydraulic stress. */
    {
      float Kfrac = HydVulnerability(Out->PsiLeaf, &VType->Xylem[Layer].Curve);
      if (Kfrac < 0.0f) Kfrac = 0.0f;
      if (Kfrac > 1.0f) Kfrac = 1.0f;
      Out->PLC = 100.0f * (1.0f - Kfrac);
      Out->SafetyMargin = Out->PsiLeaf - VType->P50[Layer];
    }
  }

  /* --- back to DHSVM units --------------------------------------------- */
  Out->An      = AnTot;
  Out->Escheme = ETot * CR_MMOL_TO_MPS;   /* m/s, per m2 ground              */

  /* Conductance-basis resistance, floored at RsMin/LAI as DHSVM's Jarvis
     path is (Wigmosta et al. 1994 eq. 15 with all factors at 1).  For the
     schemes with EConsistent set, EvapoTranspiration() replaces this with
     the resistance that makes Penman-Monteith return Escheme. */
  RcFloor = (Lai > 0.0f) ? VType->RsMin[Layer] / Lai : VType->RsMax[Layer];
  GsMs = GwTot * PhotoMolarToVelocity(Met->Tair, Met->Press);
  if (GsMs <= 1.0f / VType->RsMax[Layer])
    Rc = VType->RsMax[Layer];
  else
    Rc = 1.0f / GsMs;
  if (Rc < RcFloor) Rc = RcFloor;
  if (Rc > VType->RsMax[Layer]) Rc = VType->RsMax[Layer];

  Out->Rc = Rc;
  return Rc;
}

/*****************************************************************************
  CanopyResistanceFromFlux()

  DHSVM's transpiration is
      E = (Slope + Gamma) / (Slope + Gamma (1 + Rc/Ra)) * EPot
  (EvapoTranspiration.c).  ProfitMax, ProfitMax2 and SOX have already solved
  their transpiration with a leaf energy balance and the leaf-to-air VPD;
  re-evaluating Penman-Monteith with the air VPD at Rc = 1/gc gave a different
  number, and at low VPD a much larger one.  Inverting the expression gives
  the Rc for which the water balance removes exactly what the optimizer
  decided:
      Rc = Ra (Slope + Gamma) (EPot/E - 1) / Gamma.
  E >= EPot cannot be honoured (the atmosphere cannot take it) and returns
  RcMin; E <= 0 returns RcMax.  docs/provenance.md departure #1.
*****************************************************************************/
float CanopyResistanceFromFlux(float Eflux, float EPot, float Slope,
                               float Gamma, float Ra, float RcMin, float RcMax)
{
  float Rc;

  if (RcMin < 0.0f) RcMin = 0.0f;
  if (RcMax < RcMin) RcMax = RcMin;
  if (Eflux <= 0.0f || EPot <= 0.0f || Gamma <= 0.0f || Ra <= 0.0f)
    return RcMax;
  if (Eflux >= EPot)
    return RcMin;

  Rc = Ra * (Slope + Gamma) * (EPot / Eflux - 1.0f) / Gamma;

  if (Rc < RcMin) Rc = RcMin;
  if (Rc > RcMax) Rc = RcMax;
  return Rc;
}
