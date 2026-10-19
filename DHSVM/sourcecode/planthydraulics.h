/*****************************************************************************
  planthydraulics.h

  Xylem hydraulics: vulnerability curves, the Kirchhoff supply function from
  root collar to canopy, and the hydraulic failure threshold.

  Phase III scope note.  This header is much smaller than its predecessor
  because two things left: the duplicate FvCB kernel and leaf energy balance
  went to Photosynthesis.c (Phase II), and the whole layered rhizosphere went
  to RootHydraulics.c (Phase III).  What remains is one stage of the SPAC,
  with one source per equation.

  See PlantHydraulics.c for the source list and the structure of the supply
  function.
*****************************************************************************/

#ifndef PLANTHYDRAULICS_H
#define PLANTHYDRAULICS_H

#define HYD_TINY          1.0e-12f

/* "Physiological zero" conductance, as a fraction of maximum, defining
   Pcrit.  [SPE17] use a small fixed fraction; 5e-4 reproduces their quoted
   Pcrit of ca. -4 MPa for the default b=2, c=3 curve. */
#define HYD_KCRIT_FRAC    0.0005

/* The sigmoidal curve of [ELL20] Eqn 2 / CLM5-PHS has a long algebraic tail:
   at K = 5e-4 its Pcrit would be |psi50| * 1999^(1/a), i.e. -38 MPa for
   a = 3.  Neither Eller nor Kennedy define a Pcrit for it, so a separate,
   larger fraction bounds the supply curve for that form.  Numerical bound,
   not physiology; the Weibull/Sperry benchmarks are unaffected. */
#define HYD_KCRIT_FRAC_SIGMOIDAL 0.05

/* Default Weibull shape.  [SPE17] run c = 3 with b = 2 MPa; that pair is
   what their published Pcrit of ca. -4 MPa refers to and what the harness
   benchmarks reproduce. */
#define HYD_WEIBULL_C_DEFAULT  3.0

/* Ci/Ca under favourable conditions, used to derive KxMax from Vcmax25
   ([SPE17] coordination hypothesis).  A broad C3 generalization: sensitivity
   testing gives -30%/+52% in transpiration for 0.60/0.80, so this is the
   dominant control on magnitude and should be treated as a calibration
   target, not a constant. */
#define HYD_CICA_TARGET   0.70

#define HYD_NCUMULANT     512    /* Kirchhoff grid; built once per class     */
#define HYD_NSUPPLY       120    /* points on the supply curve               */
#define HYD_MAXSUPPLY     256

enum HydCurveForm { HYD_WEIBULL = 0, HYD_SIGMOIDAL = 1 };

/* Vulnerability curve.  Weibull: P1 = b (MPa), P2 = c.
   Sigmoidal: P1 = psi50 (MPa, negative), P2 = a. */
typedef struct {
  int   Form;
  float P1;
  float P2;
} HYDCURVE;

/* Xylem: the curve plus its Kirchhoff cumulant.  Built once per vegetation
   class at startup; it depends only on traits, not on soil state.

   UNITS.  KxMax, Krs and every E on the supply curve are per m2 of LEAF
   area (the coordination that derives KxMax is a leaf-scale calculation, and
   [ELL20] r_p,min is on a leaf-area basis).  The seam multiplies by LAI where
   a ground-area quantity is needed (RootSolve, Ecrit diagnostics). */
typedef struct {
  HYDCURVE Curve;
  float KxMax;                  /* mmol/m2 LEAF/s/MPa                        */
  float Pcrit;                  /* MPa                                       */
  float Pgrav;                  /* MPa, rho g h for the canopy height: the
                                   hydrostatic drop between collar and
                                   canopy.  [SPE17] canopy pressure and
                                   [ELL20] Notes S1.2 both carry it.  Set
                                   with HydSetCanopyHeight(); 0 if unset.  */
  float dP;                     /* grid step (negative)                      */
  float Fcrit;                  /* cumulant at Pcrit                         */
  int   N;
  float P[HYD_NCUMULANT];
  float F[HYD_NCUMULANT];
} HYDXYLEM;

/* Discretized supply curve for one pixel and timestep. */
typedef struct {
  int   N;
  float Heff;                   /* effective soil potential (MPa)            */
  float Krs;                    /* root system conductance                   */
  float Pcrit;
  float Ecrit;                  /* mmol/m2/s                                 */
  float Kcmax, Kcrit;
  float E[HYD_MAXSUPPLY];
  float Pleaf[HYD_MAXSUPPLY];
  float Pcollar[HYD_MAXSUPPLY];
  float Kc[HYD_MAXSUPPLY];
} HYDSUPPLY;

float HydWeibull(float Psi, float B, float C);
float HydSigmoidal(float Psi, float Psi50, float A);
float HydVulnerability(float Psi, const HYDCURVE *V);
float HydP50ToWeibullB(float P50, float C);
void  HydCurveWeibullFromP50(HYDCURVE *V, float P50, float C);
void  HydCurveSigmoidal(HYDCURVE *V, float P50, float A);
float HydCriticalPressure(const HYDCURVE *V);

void  HydBuildCumulant(HYDXYLEM *X, const HYDCURVE *V, float KxMax);
void  HydSetCanopyHeight(HYDXYLEM *X, float HeightM);
float HydCumulant(const HYDXYLEM *X, float Psi);
float HydCumulantInverse(const HYDXYLEM *X, float Ftarget);

float HydLeafPressure(const HYDXYLEM *X, float Heff, float Krs, float E,
                      float *CollarOut, int *Failed);
float HydEcrit(const HYDXYLEM *X, float Heff, float Krs);
void  HydBuildSupply(const HYDXYLEM *X, float Heff, float Krs, HYDSUPPLY *S);


#endif
