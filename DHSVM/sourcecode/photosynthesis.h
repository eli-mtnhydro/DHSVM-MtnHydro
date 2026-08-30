/*
 * SUMMARY:      Photosynthesis.h - Coupled photosynthesis / stomatal conductance
 * USAGE:        Part of DHSVM
 *
 * DESCRIPTION:  Farquhar-von Caemmerer-Berry (FvCB) photosynthesis model
 *               coupled to the Medlyn et al. (2011) unified stomatal
 *               optimization (USO) conductance model.
 *               Provides a replacement for the multiplicative Jarvis
 *               formulation in CanopyResistance().
 *
 *               Structure follows Noah-MP's STOMATA subroutine
 *               (module_sf_noahmplsm.F90); temperature responses use Arrhenius /
 *               peaked-Arrhenius forms rather than Noah-MP's Q10 forms (see note
 *               below).  Stomatal conductance uses Medlyn rather than Ball-Berry.
 *
 * REFERENCES:
 *   Farquhar, G.D., von Caemmerer, S., Berry, J.A. (1980) Planta 149:78-90.
 *   Medlyn, B.E. et al. (2011) Global Change Biology 17:2134-2144.
 *   Lin, Y.-S. et al. (2015) Nature Climate Change 5:459-464.
 *   Bernacchi, C.J. et al. (2001) Plant Cell Environ 24:253-259.
 *   Medlyn, B.E. et al. (2002) Plant Cell Environ 25:1167-1179.
 *   de Pury, D.G.G., Farquhar, G.D. (1997) Plant Cell Environ 20:537-557.
 *   De Kauwe, M.G. et al. (2015) Geosci. Model Dev. 8:431-452.
 *   Niu, G.-Y. et al. (2011) J. Geophys. Res. 116:D12109.  [Noah-MP]
 *
 * NOTE ON TEMPERATURE RESPONSES:
 *   Noah-MP uses Q10 forms (AKC=2.1, AKO=1.2, AVCMX=2.4) for Kc, Ko and Vcmax.
 *   We use Arrhenius (Kc, Ko, Gamma*, Rd) and peaked Arrhenius (Vcmax, Jmax)
 *   instead, for two reasons: (1) they are what plantecophys, CLM and Bonan's
 *   reference code use, so the kernel can be validated by direct diff; and
 *   (2) Q10 forms have no high-temperature deactivation and drift badly at the
 *   cold end, which matters for boreal and subalpine sites.  This is a
 *   deliberate deviation from Noah-MP -- everything else follows its structure.
 *
 */

#ifndef PHOTOSYNTHESIS_H
#define PHOTOSYNTHESIS_H

/*****************************************************************************
 Physical constants and unit conversions
 *****************************************************************************/

#define PHOTO_RGAS          8.3145   /* universal gas constant (J/mol/K)       */
#define PHOTO_TFRZ          273.15   /* 0 degC in K                            */
#define PHOTO_P0            101325.0 /* standard sea-level pressure (Pa)       */

/* Visible shortwave (W/m2) -> photosynthetic photon flux (umol/m2/s).
   McCree (1972) solar-spectrum quantum ratio.  Noah-MP uses 4.6. */
#define PHOTO_PAR_PER_WATT  4.57

/* Vapor pressure deficit: DHSVM carries Pa, Medlyn g1 is fitted in kPa^0.5 */
#define PHOTO_PA_TO_KPA     0.001

/* Diffusivity ratios relative to CO2 (Farquhar & Sharkey 1982) */
#define PHOTO_H2O_CO2_STOM  1.6      /* stomatal pathway (molecular diffusion) */
#define PHOTO_H2O_CO2_BL    1.37     /* leaf boundary layer (turbulent)        */

/* Additional parameters */
#define PHOTO_JMAXRATIO  1.7       /* 1.7 is common temperate value (Medlyn et al. 2002) */
#define PHOTO_RD25RATIO  0.015     /* 0.015 of Vcmax25 (Collatz et al. 1991) */
#define PHOTO_G0         0.01      /* 0.01 mol/m2/s is the SiB2/CABLE C3 default */
#define PHOTO_T0         -4.0      /* lower temperature (C) for dormancy */
#define PHOTO_T1         6.0       /* upper temperature (C) for dormancy */
#define PHOTO_TAU        720000.0  /* time (seconds) for dormancy--default 200 hours */
#define PHOTO_CA         420.0     /* (umol/mol) */
#define PHOTO_GB         2.0       /* (mol/m^2/s) leaf boundary layer conductance */

/*****************************************************************************
 FvCB kinetic constants
 Kc, Ko and Gamma* are partial-pressure quantities.  They are stored here in
 Pa at 25 degC and standard pressure, and the kernel works in Pa internally,
 so elevation is handled exactly.
 Values are Bernacchi et al. (2001) converted from umol/mol at 101.325 kPa.
 Noah-MP's equivalents (KC25=30 Pa, KO25=3.0e4 Pa) are in the same range.
 *****************************************************************************/

#define PHOTO_KC25          41.03    /* Rubisco M-M constant for CO2 (Pa)      */
#define PHOTO_KO25          28210.0  /* Rubisco M-M constant for O2 (Pa)       */
#define PHOTO_GAMMASTAR25   4.332    /* CO2 compensation pt w/o Rd (Pa)        */
#define PHOTO_O2_MOLFRAC    0.209    /* O2 mole fraction of dry air (-)        */

#define PHOTO_EA_KC         79430.0  /* activation energy, Kc (J/mol)          */
#define PHOTO_EA_KO         36380.0  /* activation energy, Ko (J/mol)          */
#define PHOTO_EA_GAMMASTAR  37830.0  /* activation energy, Gamma* (J/mol)      */
#define PHOTO_EA_RD         46390.0  /* activation energy, Rd (J/mol)          */

/* Peaked-Arrhenius parameters for Vcmax and Jmax (Medlyn et al. 2002).
   These match the plantecophys Photosyn() defaults. */
#define PHOTO_EA_VCMAX      58550.0  /* activation energy, Vcmax (J/mol)       */
#define PHOTO_HD_VCMAX      200000.0 /* deactivation energy, Vcmax (J/mol)     */
#define PHOTO_DS_VCMAX      629.26   /* entropy term, Vcmax (J/mol/K)          */
#define PHOTO_EA_JMAX       29680.0  /* activation energy, Jmax (J/mol)        */
#define PHOTO_HD_JMAX       200000.0 /* deactivation energy, Jmax (J/mol)      */
#define PHOTO_DS_JMAX       631.88   /* entropy term, Jmax (J/mol/K)           */

/* Electron transport.  ALPHA is the quantum yield of electron transport per
   absorbed photon; ALPHA/4 is the maximum quantum yield of CO2 assimilation
   (0.24/4 = 0.06), which is exactly Noah-MP's QE25 and the plantecophys
   default.  THETA is the curvature of the non-rectangular hyperbola. */
#define PHOTO_ALPHA         0.24     /* mol electrons / mol absorbed photons   */
#define PHOTO_THETA         0.85     /* curvature of J vs. absorbed PAR (-)    */

/* Solver control */
#define PHOTO_CI_TOL        1.0e-6   /* Ci convergence tolerance (mol/mol)     */
#define PHOTO_MAX_ITER      10       /* Brent iteration cap (typically <15)    */
#define PHOTO_TINY          1.0e-10  /* guard against division by zero         */

/* Leaf-level bounds.  MIN_GS is a floor for numerical safety only; the
   physical floor is G0 from the vegetation table. */
#define PHOTO_MIN_GS        1.0e-6   /* mol H2O/m2/s                           */
#define PHOTO_MAX_GS        2.0      /* mol H2O/m2/s (sanity cap)              */

/* Two-leaf canopy scaling (de Pury & Farquhar 1997) */
#define PHOTO_G_SPHERICAL   0.5      /* leaf projection, spherical angle dist. */
#define PHOTO_KN_NITROGEN   0.30     /* canopy N/Vcmax extinction coeff (-)    */
#define PHOTO_MIN_SINALT    0.02     /* floor on sin(solar altitude), limits   */
                                     /* k_b as the sun approaches the horizon  */

/* Convert conductance from molar (mol/m2/s) to velocity (m/s) units.
   Requires air temperature and pressure */
float PhotoMolarToVelocity(float Tair, float Press);

/* Core leaf kernel (self-contained)
   Solves the coupled FvCB / diffusion / Medlyn system at fixed leaf
   temperature by bracketed root-finding on Ci.

   Inputs:
     Vcmax25  Maximum rate of carboxylation by the Rubisco enzyme
     G1       Slope of stomatal conductance with respect to VPD
     ParAbs   Absorbed PAR (umol/m2/s, per unit leaf area)
     Tleaf    Leaf temperature (degC)
     Vpd      Leaf-to-air vapor pressure deficit (kPa)
     Ca       Atmospheric CO2 mole fraction (umol/mol)
     Press    Air pressure (Pa)
     Gb       Leaf boundary layer conductance to H2O (mol/m2/s)
     Beta     Soil moisture stress factor, 0-1 (applies to the slope term)
     Dormancy Temperature hysteresis factor, 0-1 (scales Vcmax)
   Outputs (any may be NULL):
     An       Net assimilation (umol/m2/s)
     Gs       Stomatal conductance to H2O (mol/m2/s)
     Ci       Intercellular CO2 mole fraction (umol/mol)                       */

void PhotoLeafFlux(float Vcmax25, float G1, float ParAbs, float Tleaf, float Vpd,
  float Ca, float Press, float Gb, float Beta, float Dormancy,
  float *An, float *Gs, float *Ci);

/* Two-leaf (sunlit/shaded) canopy partitioning after de Pury & Farquhar (1997)

   Inputs:
     Lai        canopy leaf area index (m2/m2)
     SinAlt     sine of solar altitude = cos(solar zenith); <=0 means night
     ParBeam    absorbed visible flux attributable to the direct beam (W/m2)
     ParDiff    absorbed visible flux attributable to diffuse radiation (W/m2)
     Vcmax25    top-of-canopy Vcmax25 (umol/m2/s)
   Outputs:
     LaiSun/LaiSha        sunlit and shaded leaf area index (m2/m2)
     ParSun/ParSha        mean absorbed PAR per unit leaf area (umol/m2/s)
     VcmaxSun/VcmaxSha    mean Vcmax25 per unit leaf area (umol/m2/s)          */
void PhotoTwoLeafPartition(float Lai, float SinAlt, float ParBeam, float ParDiff,
  float Vcmax25, float *LaiSun, float *LaiSha, float *ParSun, float *ParSha,
  float *VcmaxSun, float *VcmaxSha);

/* Full canopy conductance: runs the leaf kernel for the sunlit and shaded
   fractions and returns the LAI-weighted canopy conductance in m/s, ready to
   be inverted into the resistance that EvapoTranspiration() expects.
   Returns canopy conductance (m/s); optionally reports canopy An (umol/m2/s). */
float PhotoCanopyConductance(float Vcmax25, float G1, float Lai, float SinAlt,
  float ParBeam, float ParDiff, float Tair, float Vpd, float Ca, float Press,
  float Gb, float Beta, float Dormancy, float *AnCanopy);

#endif
