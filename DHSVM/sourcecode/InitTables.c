
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "Calendar.h"
#include "constants.h"
#include "fileio.h"
#include "getinit.h"
#include "planthydraulics.h"
#include "roothydraulics.h"
#include "stomatalscheme.h"
#include "photosynthesis.h"

/*******************************************************************************
  PrintRhizosphereSummary()

  One line per vegetation class/layer with a root length index, per soil
  type: the perirhizal geometry and the bulk soil potential at which the
  rhizosphere starts to disconnect (B K(h) = a_root k_r).  Uses the class's
  peak monthly LAI and the layer holding most of its roots.  Startup
  diagnostic only.
*******************************************************************************/
static void PrintRhizosphereSummary(OPTIONSTRUCT *Options, SOILTABLE *SType,
                                    int NSoil, VEGTABLE *VType, int NVeg)
{
  int i, j, k, m, imax;
  float LaiMax, Ltot, Rld, Aroot, Aprhiz, Rho, B, KrsGround, Kr, Hae, Ks,
        Lambda, Psi, d;

  for (i = 0; i < NVeg; i++) {
    for (j = 0; j < VType[i].NVegLayers; j++) {
      if (VType[i].RootLengthIndex[j] <= 0.0f) {
        printf("Veg %d layer %d: perirhizal stage OFF (no ROOT LENGTH INDEX)\n",
               i + 1, j);
        continue;
      }
      LaiMax = 0.0f;
      for (m = 0; m < 12; m++)
        if (VType[i].LAIMonthly[j][m] > LaiMax) LaiMax = VType[i].LAIMonthly[j][m];
      if (j == 0 && VType[i].OverStory) LaiMax *= VEG_LAI_ADJ;   /* as InitTerrainMaps applies it */
      if (LaiMax <= 0.0f) LaiMax = 1.0f;

      imax = 0;
      for (k = 1; k < VType[i].NSoilLayers; k++)
        if (VType[i].RootFract[j][k] > VType[i].RootFract[j][imax]) imax = k;

      Ltot  = VType[i].RootLengthIndex[j];
      Aroot = VType[i].RootRadius[j];
      d = VType[i].RootDepth[imax]; if (d < 0.01f) d = 0.01f;
      Rld    = Ltot * VType[i].RootFract[j][imax] / d;
      Aprhiz = RootPerirhizalRadius(Rld, Aroot);
      Rho    = Aprhiz / Aroot;
      B      = RootGeometryFactor(Rho);
      KrsGround = VType[i].Krs[j] * LaiMax;
      Kr     = RootRadialFromKrs(KrsGround, Aroot, Ltot);

      printf("Veg %d layer %d perirhizal: root length %.1f km/m2, radius %.2f mm, "
             "root surface index %.2f m2/m2, k_r=%.2e 1/s (from Krs x peak LAI %.2f), "
             "layer %d RLD=%.0f m/m3, rho=%.1f, B=%.2f\n",
             i + 1, j, Ltot / 1000.0f, Aroot * 1000.0f,
             (float)(2.0 * 3.14159265 * Aroot * Ltot), Kr, LaiMax, imax,
             Rld, Rho, B);

      for (k = 0; k < NSoil; k++) {
        if (imax >= SType[k].NLayers) continue;
        Ks     = SType[k].Ks[imax];
        /* Mirror InitTerrainMaps: with VERTICAL KSAT SOURCE = ANISOTROPY the
           layer's vertical Ks is the depth-averaged lateral conductivity
           over the layer divided by the anisotropy, not the table value. */
        if (Options->UseKsatAnisotropy) {
          float Top = 0.0f, Bot, Tr;
          for (m = 0; m < imax; m++) Top += VType[i].RootDepth[m];
          Bot = Top + d;
          Tr = CalcTransmissivity(Bot, Top, SType[k].KsLat, SType[k].KsLatExp,
                                  SType[k].DepthThresh);
          Ks = Tr / d / SType[k].KsAnisotropy;
        }
        if (SType[k].KsRhizo != NULL && SType[k].KsRhizo[imax] > 0.0f)
          Ks = SType[k].KsRhizo[imax];
        Lambda = SType[k].PoreDist[imax];
        Hae    = SType[k].Press[imax];          /* m of head, positive     */
        Psi    = RootDisconnectPotential(Aroot * Kr, B, Ks, Hae, Lambda);
        printf("   soil type %d (%s): %s=%.2e m/s lambda=%.2f air entry=%.2f m -> "
               "rhizosphere carries half the soil-to-root resistance at psi_soil ~ %.2f MPa\n",
               k + 1, SType[k].Desc,
               (SType[k].KsRhizo != NULL && SType[k].KsRhizo[imax] > 0.0f)
                 ? "Ks(rhizosphere key)" : "KsVert(model)",
               Ks, Lambda, Hae, Psi);
      }
    }
  }
}

/*******************************************************************************/
/*				  InitTables()                                 */
/*******************************************************************************/
void InitTables(int StepsPerDay, LISTPTR Input, OPTIONSTRUCT *Options, 
  MAPSIZE *Map, SOILTABLE **SType, LAYER *Soil, VEGTABLE **VType, LAYER *Veg,
  LAKETABLE **LType, TIMESTRUCT *Time)
{
  printf("Initializing tables\n");

  if ((Soil->NTypes = InitSoilTable(Options, SType, Input, Soil,
    Options->Infiltration, Time)) == 0)
    ReportError("Input Options File", 8);

  if ((Veg->NTypes = InitVegTable(VType, Input, Options, Veg)) == 0)
    ReportError("Input Options File", 8);

  if (Options->Hydraulics == HYD_KRSSUF)
    PrintRhizosphereSummary(Options, *SType, Soil->NTypes, *VType, Veg->NTypes);
  
  if (Options->LakeDynamics) {
    if ((Map->NumLakes = InitLakeTable(LType, Input, Options)) == 0)
      ReportError("Input Options File", 8);
  } else {
    Map->NumLakes = 0;
  }

  InitSatVaporTable();
}

/********************************************************************************
Function Name: InitSoilTable()

Purpose      : Initialize the soil lookup table
Processes most of the following section in InFileName:
[SOILS]

Required     :
SOILTABLE **SType - Pointer to lookup table
LISTPTR Input     - Pointer to linked list with input info
LAYER *Soil       - Pointer to structure with soil layer information

Returns      : Number of soil layers

Modifies     : SoilTable and Soil

Comments     :
********************************************************************************/
int InitSoilTable(OPTIONSTRUCT *Options, SOILTABLE ** SType,
  LISTPTR Input, LAYER * Soil, int InfiltOption, TIMESTRUCT *Time)
{
  const char *Routine = "InitSoilTable";
  int i;			/* counter */
  int j;			/* counter */
  int NSoils;			/* Number of soil types */
  char KeyName[soil_last_key + 1][BUFSIZE + 1];
  char *KeyStr[] = {
    "SOIL DESCRIPTION",
    "LATERAL CONDUCTIVITY",
    "EXPONENTIAL DECREASE",
    "DEPTH THRESHOLD",
    "VERTICAL ANISOTROPY",
    "MAXIMUM INFILTRATION",
    "CAPILLARY DRIVE",
    "DEEP FLUX",
    "SURFACE ALBEDO",
    "MANNINGS N",
    "NUMBER OF SOIL LAYERS",
    "POROSITY",
    "PORE SIZE DISTRIBUTION",
    "BUBBLING PRESSURE",
    "FIELD CAPACITY",
    "WILTING POINT",
    "BULK DENSITY",
    "VERTICAL CONDUCTIVITY",
    "THERMAL CONDUCTIVITY",
    "THERMAL CAPACITY",
    "RHIZOSPHERE CONDUCTIVITY"   /* optional, m/s per layer; matrix Ks for
                                    the perirhizal stage.  Absent or <=0 =
                                    the model KsVert                       */
  };
  char SectionName[] = "SOILS";
  char VarStr[soil_last_key + 1][BUFSIZE + 1];


  /* Get the number of different soil types */
  GetInitString(SectionName, "NUMBER OF SOIL TYPES", "", VarStr[0],
    (unsigned long)BUFSIZE, Input);
  if (!CopyInt(&NSoils, VarStr[0], 1))
    ReportError("NUMBER OF SOIL TYPES", 51);

  if (NSoils == 0)
    return NSoils;

  if (!(Soil->NLayers = (int *)calloc(NSoils, sizeof(int))))
    ReportError((char *)Routine, 1);

  if (!(*SType = (SOILTABLE *)calloc(NSoils, sizeof(SOILTABLE))))
    ReportError((char *)Routine, 1);

  /********** Read information and allocate memory for each soil type *********/

  Soil->MaxLayers = 0;

  for (i = 0; i < NSoils; i++) {

    /* Read the key-entry pairs from the input file */
    for (j = 0; j <= soil_last_key; j++) {
      sprintf(KeyName[j], "%s %d", KeyStr[j], i + 1);
      GetInitString(SectionName, KeyName[j], "", VarStr[j],
        (unsigned long)BUFSIZE, Input);
    }

    /* Assign the entries to the appropriate variables */
    if (IsEmptyStr(VarStr[soil_description]))
      ReportError(KeyName[soil_description], 51);

    strcpy((*SType)[i].Desc, VarStr[soil_description]);
    (*SType)[i].Index = i;

    if (!CopyFloat(&((*SType)[i].KsLat), VarStr[lateral_ks], 1))
      ReportError(KeyName[lateral_ks], 51);

    if (!CopyFloat(&((*SType)[i].KsLatExp), VarStr[exponent], 1))
      ReportError(KeyName[exponent], 51);

    if (!CopyFloat(&((*SType)[i].DepthThresh), VarStr[depth_thresh], 1))
      ReportError(KeyName[depth_thresh], 51);
    
    if (!CopyFloat(&((*SType)[i].KsAnisotropy), VarStr[anisotropy], 1))
      ReportError(KeyName[exponent], 51);
    
    if (!CopyFloat(&((*SType)[i].MaxInfiltrationRate), VarStr[max_infiltration], 1))
      ReportError(KeyName[max_infiltration], 51);

    if (InfiltOption == DYNAMIC) {
      if (!CopyFloat(&((*SType)[i].G_Infilt), VarStr[capillary_drive], 1))
        ReportError(KeyName[capillary_drive], 51);
    }
    else (*SType)[i].G_Infilt = NOT_APPLICABLE;

    if (!CopyFloat(&((*SType)[i].DeepFlux), VarStr[deepflux], 1))
      ReportError(KeyName[deepflux], 51);
    /* Convert m/yr (config) to m/timestep (used in DistributeSatflow) */
    (*SType)[i].DeepFlux /= (DAYPYEAR * Time->NDaySteps);

    if (!CopyFloat(&((*SType)[i].Albedo), VarStr[soil_albedo], 1))
      ReportError(KeyName[soil_albedo], 51);
    
    if (!CopyInt(&(*SType)[i].NLayers, VarStr[number_of_layers], 1))
      ReportError(KeyName[number_of_layers], 51);
    Soil->NLayers[i] = (*SType)[i].NLayers;

    if (Soil->NLayers[i] > Soil->MaxLayers)
      Soil->MaxLayers = Soil->NLayers[i];
    
    if (!CopyFloat(&((*SType)[i].Manning), VarStr[manning], 1))
      ReportError(KeyName[manning], 51);
    
    /* allocate memory for the soil layers */
    if (!((*SType)[i].Porosity = (float *)calloc((*SType)[i].NLayers,
      sizeof(float))))
      ReportError((char *)Routine, 1);
    if (!((*SType)[i].PoreDist = (float *)calloc((*SType)[i].NLayers,
      sizeof(float))))
      ReportError((char *)Routine, 1);
    if (!((*SType)[i].Press = (float *)calloc((*SType)[i].NLayers,
      sizeof(float))))
      ReportError((char *)Routine, 1);
    if (!((*SType)[i].FCap = (float *)calloc((*SType)[i].NLayers,
      sizeof(float))))
      ReportError((char *)Routine, 1);
    if (!((*SType)[i].WP = (float *)calloc((*SType)[i].NLayers,
      sizeof(float))))
      ReportError((char *)Routine, 1);
    if (!((*SType)[i].Dens = (float *)calloc((*SType)[i].NLayers,
      sizeof(float))))
      ReportError((char *)Routine, 1);
    if (!((*SType)[i].Ks = (float *)calloc((*SType)[i].NLayers,
      sizeof(float))))
      ReportError((char *)Routine, 1);
    if (!((*SType)[i].KhDry = (float *)calloc((*SType)[i].NLayers,
      sizeof(float))))
      ReportError((char *)Routine, 1);
    if (!((*SType)[i].KhSol = (float *)calloc((*SType)[i].NLayers,
      sizeof(float))))
      ReportError((char *)Routine, 1);
    if (!((*SType)[i].Ch = (float *)calloc((*SType)[i].NLayers,
      sizeof(float))))
      ReportError((char *)Routine, 1);

    if (!CopyFloat((*SType)[i].Porosity, VarStr[porosity], (*SType)[i].NLayers))
      ReportError(KeyName[porosity], 51);

    if (!CopyFloat((*SType)[i].PoreDist, VarStr[pore_size],
      (*SType)[i].NLayers))
      ReportError(KeyName[pore_size], 51);

    if (!CopyFloat((*SType)[i].Press, VarStr[bubbling_pressure],
      (*SType)[i].NLayers))
      ReportError(KeyName[bubbling_pressure], 51);

    if (!CopyFloat((*SType)[i].FCap, VarStr[field_capacity],
      (*SType)[i].NLayers))
      ReportError(KeyName[field_capacity], 51);

    if (!CopyFloat((*SType)[i].WP, VarStr[wilting_point], (*SType)[i].NLayers))
      ReportError(KeyName[wilting_point], 51);

    if (!CopyFloat((*SType)[i].Dens, VarStr[bulk_density], (*SType)[i].NLayers))
      ReportError(KeyName[bulk_density], 51);

    if (!CopyFloat((*SType)[i].Ks, VarStr[vertical_ks], (*SType)[i].NLayers))
      ReportError(KeyName[vertical_ks], 51);

    if (!CopyFloat((*SType)[i].KhSol, VarStr[solids_thermal], (*SType)[i].NLayers))
      ReportError(KeyName[solids_thermal], 51);

    if (!CopyFloat((*SType)[i].Ch, VarStr[thermal_capacity], (*SType)[i].NLayers))
      ReportError(KeyName[thermal_capacity], 51);

    /* Rhizosphere (matrix) conductivity, optional.  DHSVM's lateral and
       vertical Ks are effective hillslope-drainage values (macropores,
       pipes); flow through the last millimetres of soil to a root at
       -0.5 MPa is matrix flow, two to three orders slower.  The perirhizal
       stage needs the latter. */
    if (!((*SType)[i].KsRhizo = (float *)calloc((*SType)[i].NLayers, sizeof(float))))
      ReportError((char *)Routine, 1);
    if (IsEmptyStr(VarStr[rhizosphere_ks])) {
      for (j = 0; j < (*SType)[i].NLayers; j++) (*SType)[i].KsRhizo[j] = -1.0f;
    }
    else if (!CopyFloat((*SType)[i].KsRhizo, VarStr[rhizosphere_ks], (*SType)[i].NLayers))
      ReportError(KeyName[rhizosphere_ks], 51);
  }

  for (i = 0; i < NSoils; i++)
    for (j = 0; j < (*SType)[i].NLayers; j++) {
      (*SType)[i].KhDry[j] = CalcKhDry((*SType)[i].Dens[j]);
      if (((*SType)[i].Porosity[j] < (*SType)[i].FCap[j])
        || ((*SType)[i].Porosity[j] < (*SType)[i].WP[j])
        || ((*SType)[i].FCap[j] < (*SType)[i].WP[j]))
        ReportError((*SType)[i].Desc, 11);
    }

  return NSoils;
}

/********************************************************************************
Function Name: InitVegTable()

Purpose      : Initialize the vegetation lookup table
Processes most of the following section in the input file:
[VEGETATION]

Required     :
VEGTABLE **VType - Pointer to lookup table
LISTPTR Input    - Pointer to linked list with input info
LAYER *Veg       - Pointer to structure with veg layer information

Returns      : Number of vegetation types

Modifies     : VegTable and Veg

Comments     :
********************************************************************************/
/*****************************************************************************
  Optional per-class keys.  An absent key (empty string from GetInitString)
  fills the layers with Default; a present key must supply one value per
  vegetation layer, exactly like the required keys.
*****************************************************************************/
static void ReadOptionalFloats(float *Dest, char *Str, char *Key, int N,
                               float Default)
{
  int k;

  if (IsEmptyStr(Str)) {
    for (k = 0; k < N; k++) Dest[k] = Default;
    return;
  }
  if (!CopyFloat(Dest, Str, N))
    ReportError(Key, 51);
}

/* XYLEM CURVE: one token per layer, WEIBULL (default) or SIGMOIDAL. */
static void ReadOptionalCurveForm(int *Dest, char *Str, char *Key, int N)
{
  char Buf[BUFSIZE + 1];
  char *Tok;
  int k = 0;

  for (k = 0; k < N; k++) Dest[k] = HYD_WEIBULL;
  if (IsEmptyStr(Str)) return;

  strncpy(Buf, Str, BUFSIZE);
  Buf[BUFSIZE] = '\0';
  Tok = strtok(Buf, " \t,");
  for (k = 0; k < N; k++) {
    if (Tok == NULL) ReportError(Key, 51);
    if      (strncmp(Tok, "SIG", 3) == 0 || strncmp(Tok, "sig", 3) == 0)
      Dest[k] = HYD_SIGMOIDAL;
    else if (strncmp(Tok, "WEI", 3) == 0 || strncmp(Tok, "wei", 3) == 0)
      Dest[k] = HYD_WEIBULL;
    else
      ReportError(Key, 51);
    Tok = strtok(NULL, " \t,");
  }
}

int InitVegTable(VEGTABLE **VType, LISTPTR Input, OPTIONSTRUCT *Options, LAYER *Veg)
{
  const char *Routine = "InitVegTable";
  int i;			/* Counter */
  int j;			/* Counter */
  int k;      /* counter */
  int y;      /* counter */
  float impervious;	/* flag to check whether impervious layers are specified */

  int NVegs;		/* Number of vegetation types */

  char KeyName[veg_last_key + 1][BUFSIZE + 1];
  char *KeyStr[] = {
    "VEGETATION DESCRIPTION",
    "OVERSTORY PRESENT",
    "UNDERSTORY PRESENT",
    "FRACTIONAL COVERAGE",
    "HEMI FRACT COVERAGE",
    "TRUNK SPACE",
    "AERODYNAMIC ATTENUATION",
    "RADIATION ATTENUATION",
    "DIFFUSE RADIATION ATTENUATION",
    "CLUMPING FACTOR",
    "LEAF ANGLE A",
    "LEAF ANGLE B",
    "SCATTERING PARAMETER",
    "MAX SNOW INT CAPACITY",
    "MASS RELEASE DRIP RATIO",
    "SNOW INTERCEPTION EFF",
    "IMPERVIOUS FRACTION",
    "DETENTION FRACTION",
    "DETENTION DECAY",
    "HEIGHT",
    "MAXIMUM RESISTANCE",
    "MINIMUM RESISTANCE",
    "MAXIMUM CARBOXYLATION",
    "STOMATAL SLOPE",
    "STOMATAL INTERCEPT",
    "XYLEM PRESSURE",
    "MOISTURE THRESHOLD",
    "VAPOR PRESSURE DEFICIT",
    "RPC",
    "NUMBER OF ROOT ZONES",
    "ROOT ZONE DEPTHS",
    "OVERSTORY ROOT FRACTION",
    "UNDERSTORY ROOT FRACTION",
    "MONTHLY LIGHT EXTINCTION",
    "CANOPY VIEW ADJ FACTOR",
    "OVERSTORY MONTHLY LAI",
    "UNDERSTORY MONTHLY LAI",
    "OVERSTORY MONTHLY ALB",
    "UNDERSTORY MONTHLY ALB",
    /* optional physiology / hydraulics keys; absent = derived default */
    "XYLEM SHAPE",                 /* Weibull c or sigmoidal a per layer   */
    "XYLEM CURVE",                 /* WEIBULL | SIGMOIDAL per layer        */
    "XYLEM CONDUCTANCE",           /* KxMax, mmol/m2 leaf/s/MPa; <=0 derive*/
    "ROOT CONDUCTANCE",            /* Krs, mmol/m2 leaf/s/MPa; <=0 derive  */
    "ROOT COMPENSATION CONDUCTANCE",/* Kcomp; <0 = Krs, 0 = off            */
    "COLLAR PRESSURE MINIMUM",     /* MPa; >=0 -> 1.5 x P50                */
    "MAXIMUM STOMATAL CONDUCTANCE",/* mol/m2 leaf/s; <=0 -> 1/RsMin        */
    "LEAF WIDTH",                  /* m                                    */
    "JMAX VCMAX RATIO",            /* Jmax25/Vcmax25                       */
    "CICA TARGET",                 /* Ci/Ca for the kmax coordination      */
    "RD VCMAX RATIO",              /* Rd25/Vcmax25                         */
    "ROOT LENGTH INDEX",           /* km fine root / m2 ground; 0 = no
                                      perirhizal stage                     */
    "FINE ROOT RADIUS"             /* m                                    */
  };
  char SectionName[] = "VEGETATION";
  char VarStr[veg_last_key + 1][BUFSIZE + 1];
  float maxLAI;

  /* Get the number of different vegetation types */
  GetInitString(SectionName, "NUMBER OF VEGETATION TYPES", "", VarStr[0],
    (unsigned long)BUFSIZE, Input);
  if (!CopyInt(&NVegs, VarStr[0], 1))
    ReportError("NUMBER OF VEGETATION TYPES", 51);

  if (NVegs == 0)
    return NVegs;

  if (!(Veg->NLayers = (int *)calloc(NVegs, sizeof(int))))
    ReportError((char *)Routine, 1);

  if (!(*VType = (VEGTABLE *)calloc(NVegs, sizeof(VEGTABLE))))
    ReportError((char *)Routine, 1);

  /******* Read information and allocate memory for each vegetation type ******/

  Veg->MaxLayers = 0;
  impervious = 0.0;
  for (i = 0; i < NVegs; i++) {

    /* Read the key-entry pairs from the input file */
    for (j = 0; j <= veg_last_key; j++) {
      sprintf(KeyName[j], "%s %d", KeyStr[j], i + 1);
      GetInitString(SectionName, KeyName[j], "", VarStr[j],
        (unsigned long)BUFSIZE, Input);
    }

    /* Assign the entries to the appropriate variables */
    if (IsEmptyStr(VarStr[veg_description]))
      ReportError(KeyName[veg_description], 51);
    strcpy((*VType)[i].Desc, VarStr[veg_description]);
    MakeKeyString(VarStr[veg_description]);	/* basically makes the string all
                                            uppercase and removed spaces so
                                            it is easier to compare */
    (*VType)[i].Index = i;

    (*VType)[i].NVegLayers = 0;

    if (strncmp(VarStr[overstory], "TRUE", 4) == 0) {
      (*VType)[i].OverStory = TRUE;
      ((*VType)[i].NVegLayers)++;
    }
    else if (strncmp(VarStr[overstory], "FALSE", 5) == 0)
      (*VType)[i].OverStory = FALSE;
    else
      ReportError(KeyName[overstory], 51);

    if (strncmp(VarStr[understory], "TRUE", 4) == 0) {
      (*VType)[i].UnderStory = TRUE;
      ((*VType)[i].NVegLayers)++;
    }
    else if (strncmp(VarStr[understory], "FALSE", 5) == 0)
      (*VType)[i].UnderStory = FALSE;
    else
      ReportError(KeyName[understory], 51);

    Veg->NLayers[i] = (*VType)[i].NVegLayers;
    if ((*VType)[i].NVegLayers > Veg->MaxLayers)
      Veg->MaxLayers = (*VType)[i].NVegLayers;

    if (!CopyInt(&(*VType)[i].NSoilLayers, VarStr[number_of_root_zones], 1))
      ReportError(KeyName[number_of_root_zones], 51);

    if (!CopyFloat(&((*VType)[i].ImpervFrac), VarStr[imperv_frac], 1))
      ReportError(KeyName[imperv_frac], 51);
    impervious += (*VType)[i].ImpervFrac;

    if ((*VType)[i].ImpervFrac > 0) {
      if (!CopyFloat(&((*VType)[i].DetentionFrac), VarStr[detention_frac], 1))
        ReportError(KeyName[detention_frac], 51);
      if (!CopyFloat(&((*VType)[i].DetentionDecay), VarStr[detention_decay], 1))
        ReportError(KeyName[detention_decay], 51);
    }
    else {
      (*VType)[i].DetentionFrac = 0.;
      (*VType)[i].DetentionDecay = 0.;
    }

    /* allocate memory for the vegetation layers */
    if (!((*VType)[i].Fract = (float *)calloc((*VType)[i].NVegLayers,
      sizeof(float))))
      ReportError((char *)Routine, 1);

    if (Options->CanopyRadAtt == VARIABLE) {
      if (!((*VType)[i].HemiFract = (float *)calloc((*VType)[i].NVegLayers,
        sizeof(float))))
        ReportError((char *)Routine, 1);
    }
    else {
      (*VType)[i].HemiFract = NULL;
    }

    if (!((*VType)[i].Height = (float *)calloc((*VType)[i].NVegLayers, sizeof(float))))
      ReportError((char *)Routine, 1);

    if (!((*VType)[i].RsMax = (float *)calloc((*VType)[i].NVegLayers, sizeof(float))))
      ReportError((char *)Routine, 1);

    if (!((*VType)[i].RsMin = (float *)calloc((*VType)[i].NVegLayers, sizeof(float))))
      ReportError((char *)Routine, 1);

    if (!((*VType)[i].Vcmax25 = (float *)calloc((*VType)[i].NVegLayers, sizeof(float))))
      ReportError((char *)Routine, 1);

    if (!((*VType)[i].G1 = (float *)calloc((*VType)[i].NVegLayers, sizeof(float))))
      ReportError((char *)Routine, 1);

    if (!((*VType)[i].G0 = (float *)calloc((*VType)[i].NVegLayers, sizeof(float))))
      ReportError((char *)Routine, 1);

    if (!((*VType)[i].P50 = (float *)calloc((*VType)[i].NVegLayers, sizeof(float))))
      ReportError((char *)Routine, 1);

    if (!((*VType)[i].KxMax = (float *)calloc((*VType)[i].NVegLayers, sizeof(float))))
      ReportError((char *)Routine, 1);
    
    if (!((*VType)[i].Gmax = (float *)calloc((*VType)[i].NVegLayers, sizeof(float))))
      ReportError((char *)Routine, 1);

    if (!((*VType)[i].VulnShape = (float *)calloc((*VType)[i].NVegLayers, sizeof(float))))
      ReportError("InitVegTable", 1);
    if (!((*VType)[i].VulnForm = (int *)calloc((*VType)[i].NVegLayers, sizeof(int))))
      ReportError("InitVegTable", 1);
    if (!((*VType)[i].Krs = (float *)calloc((*VType)[i].NVegLayers, sizeof(float))))
      ReportError("InitVegTable", 1);
    if (!((*VType)[i].Kcomp = (float *)calloc((*VType)[i].NVegLayers, sizeof(float))))
      ReportError("InitVegTable", 1);
    if (!((*VType)[i].PsiCollarMin = (float *)calloc((*VType)[i].NVegLayers, sizeof(float))))
      ReportError("InitVegTable", 1);
    if (!((*VType)[i].Xylem = (HYDXYLEM *)calloc((*VType)[i].NVegLayers, sizeof(HYDXYLEM))))
      ReportError("InitVegTable", 1);
    if (!((*VType)[i].LeafWidth = (float *)calloc((*VType)[i].NVegLayers, sizeof(float))))
      ReportError("InitVegTable", 1);
    if (!((*VType)[i].JmaxRatio = (float *)calloc((*VType)[i].NVegLayers, sizeof(float))))
      ReportError("InitVegTable", 1);
    if (!((*VType)[i].Rd25Ratio = (float *)calloc((*VType)[i].NVegLayers, sizeof(float))))
      ReportError("InitVegTable", 1);
    if (!((*VType)[i].CiCaTarget = (float *)calloc((*VType)[i].NVegLayers, sizeof(float))))
      ReportError("InitVegTable", 1);
    if (!((*VType)[i].RootLengthIndex = (float *)calloc((*VType)[i].NVegLayers, sizeof(float))))
      ReportError("InitVegTable", 1);
    if (!((*VType)[i].RootRadius = (float *)calloc((*VType)[i].NVegLayers, sizeof(float))))
      ReportError("InitVegTable", 1);
    if (!((*VType)[i].WeibullB = (float *)calloc((*VType)[i].NVegLayers, sizeof(float))))
      ReportError((char *)Routine, 1);

    if (!((*VType)[i].MoistThres = (float *)calloc((*VType)[i].NVegLayers, sizeof(float))))
      ReportError((char *)Routine, 1);

    if (!((*VType)[i].VpdThres = (float *)calloc((*VType)[i].NVegLayers, sizeof(float))))
      ReportError((char *)Routine, 1);

    if (!((*VType)[i].Rpc = (float *)calloc((*VType)[i].NVegLayers, sizeof(float))))
      ReportError((char *)Routine, 1);

    if (!((*VType)[i].Albedo = (float *)calloc(((*VType)[i].NVegLayers + 1), sizeof(float))))
      ReportError((char *)Routine, 1);

    if (!((*VType)[i].MaxInt = (float *)calloc((*VType)[i].NVegLayers, sizeof(float))))

      ReportError((char *)Routine, 1);
    if (!((*VType)[i].LAI = (float *)calloc((*VType)[i].NVegLayers, sizeof(float))))
      ReportError((char *)Routine, 1);

    if (!((*VType)[i].RootFract = (float **)calloc((*VType)[i].NVegLayers, sizeof(float *))))
      ReportError((char *)Routine, 1);

    for (j = 0; j < (*VType)[i].NVegLayers; j++) {
      if (!((*VType)[i].RootFract[j] =
        (float *)calloc((*VType)[i].NSoilLayers, sizeof(float))))
        ReportError((char *)Routine, 1);
    }
    if (!((*VType)[i].RootDepth = (float *)calloc((*VType)[i].NSoilLayers, sizeof(float))))
      ReportError((char *)Routine, 1);

    if (!((*VType)[i].LAIMonthly = (float **)calloc((*VType)[i].NVegLayers, sizeof(float *))))
      ReportError((char *)Routine, 1);
    for (j = 0; j < (*VType)[i].NVegLayers; j++) {
      if (!((*VType)[i].LAIMonthly[j] = (float *)calloc(12, sizeof(float))))
        ReportError((char *)Routine, 1);
    }

    if (!((*VType)[i].AlbedoMonthly = (float **)calloc((*VType)[i].NVegLayers, sizeof(float *))))
      ReportError((char *)Routine, 1);
    for (j = 0; j < (*VType)[i].NVegLayers; j++) {
      if (!((*VType)[i].AlbedoMonthly[j] = (float *)calloc(12, sizeof(float))))
        ReportError((char *)Routine, 1);
    }

    /* assign the entries to the appropriate variables */
    /* allocation of zero memory is not supported on some
    compilers */
    if ((*VType)[i].OverStory == TRUE) {
      if (!CopyFloat(&((*VType)[i].Fract[0]), VarStr[fraction], 1))
        ReportError(KeyName[fraction], 51);
      if (Options->CanopyRadAtt == VARIABLE) {
        if (!CopyFloat(&((*VType)[i].HemiFract[0]), VarStr[hemifraction], 1))
          ReportError(KeyName[hemifraction], 51);
        if (!CopyFloat(&((*VType)[i].ClumpingFactor), VarStr[clumping_factor], 1))
          ReportError(KeyName[clumping_factor], 51);
        if (!CopyFloat(&((*VType)[i].LeafAngleA), VarStr[leaf_angle_a], 1))
          ReportError(KeyName[leaf_angle_a], 51);
        if (!CopyFloat(&((*VType)[i].LeafAngleB), VarStr[leaf_angle_b], 1))
          ReportError(KeyName[leaf_angle_b], 51);
        if (!CopyFloat(&((*VType)[i].Scat), VarStr[scat], 1))
          ReportError(KeyName[scat], 51);
        (*VType)[i].Atten = NOT_APPLICABLE;
      }
      else if (Options->CanopyRadAtt == FIXED && Options->ImprovRadiation == FALSE) {
        if (!CopyFloat(&((*VType)[i].Atten), VarStr[beam_attn], 1))
          ReportError(KeyName[beam_attn], 51);
        (*VType)[i].ClumpingFactor = NOT_APPLICABLE;
        (*VType)[i].Scat = NOT_APPLICABLE;
        (*VType)[i].LeafAngleA = NOT_APPLICABLE;
        (*VType)[i].LeafAngleB = NOT_APPLICABLE;
        (*VType)[i].Taud = NOT_APPLICABLE;
      }
      else if (Options->ImprovRadiation == TRUE) {
        if (!CopyFloat(&((*VType)[i].Taud), VarStr[diff_attn], 1))
          ReportError(KeyName[diff_attn], 51);
        (*VType)[i].Atten = NOT_APPLICABLE;
        (*VType)[i].ClumpingFactor = NOT_APPLICABLE;
        (*VType)[i].Scat = NOT_APPLICABLE;
        (*VType)[i].LeafAngleA = NOT_APPLICABLE;
        (*VType)[i].LeafAngleB = NOT_APPLICABLE;
      }

      if (!CopyFloat(&((*VType)[i].Trunk), VarStr[trunk_space], 1))
        ReportError(KeyName[trunk_space], 51);

      if (!CopyFloat(&((*VType)[i].Cn), VarStr[aerodynamic_att], 1))
        ReportError(KeyName[aerodynamic_att], 51);

      if (!CopyFloat(&((*VType)[i].MaxSnowInt), VarStr[snow_int_cap], 1))
        ReportError(KeyName[snow_int_cap], 51);

      if (!CopyFloat(&((*VType)[i].MDRatio), VarStr[mass_drip_ratio], 1))
        ReportError(KeyName[mass_drip_ratio], 51);

      if (!CopyFloat(&((*VType)[i].SnowIntEff), VarStr[snow_int_eff], 1))
        ReportError(KeyName[snow_int_eff], 51);

      if (!CopyFloat((*VType)[i].RootFract[0], VarStr[overstory_fraction],
        (*VType)[i].NSoilLayers))
        ReportError(KeyName[overstory_fraction], 51);

      if (!CopyFloat((*VType)[i].LAIMonthly[0], VarStr[overstory_monlai], 12))
        ReportError(KeyName[overstory_monlai], 51);

      maxLAI = -9999;
      for (k = 0; k < 12; k++) {
        if ((*VType)[i].LAIMonthly[0][k] > maxLAI)
          maxLAI = (*VType)[i].LAIMonthly[0][k];
      }

      if (!CopyFloat((*VType)[i].AlbedoMonthly[0], VarStr[overstory_monalb], 12))
        ReportError(KeyName[overstory_monalb], 51);

      if ((*VType)[i].UnderStory == TRUE) {
        (*VType)[i].Fract[1] = 1.0;
        if (!CopyFloat((*VType)[i].RootFract[1], VarStr[understory_fraction],
          (*VType)[i].NSoilLayers))
          ReportError(KeyName[understory_fraction], 51);

        if (!CopyFloat((*VType)[i].LAIMonthly[1], VarStr[understory_monlai], 12))
          ReportError(KeyName[understory_monlai], 51);

        if (!CopyFloat((*VType)[i].AlbedoMonthly[1], VarStr[understory_monalb], 12))
          ReportError(KeyName[understory_monalb], 51);
      }
    }
    else {
      if ((*VType)[i].UnderStory == TRUE) {
        (*VType)[i].Fract[0] = 1.0;
        if (!CopyFloat((*VType)[i].RootFract[0], VarStr[understory_fraction],
          (*VType)[i].NSoilLayers))
          ReportError(KeyName[understory_fraction], 51);

        if (!CopyFloat((*VType)[i].LAIMonthly[0], VarStr[understory_monlai], 12))
          ReportError(KeyName[understory_monlai], 51);

        if (!CopyFloat((*VType)[i].AlbedoMonthly[0], VarStr[understory_monalb], 12))
          ReportError(KeyName[understory_monalb], 51);
      }
      (*VType)[i].Trunk = NOT_APPLICABLE;
      (*VType)[i].Cn = NOT_APPLICABLE;
      (*VType)[i].Atten = NOT_APPLICABLE;
      (*VType)[i].ClumpingFactor = NOT_APPLICABLE;
    }

    if (!CopyFloat((*VType)[i].Height, VarStr[height], (*VType)[i].NVegLayers))
      ReportError(KeyName[height], 51);

    if (!CopyFloat((*VType)[i].RsMax, VarStr[max_resistance],
      (*VType)[i].NVegLayers))
      ReportError(KeyName[max_resistance], 51);

    if (!CopyFloat((*VType)[i].RsMin, VarStr[min_resistance],
      (*VType)[i].NVegLayers))
      ReportError(KeyName[min_resistance], 51);

    if (Options->StomScheme != JARVIS) {
      if (!CopyFloat((*VType)[i].G1, VarStr[stomatal_slope],
                     (*VType)[i].NVegLayers))
        ReportError(KeyName[stomatal_slope], 51);
      
      if (!CopyFloat((*VType)[i].G0, VarStr[stomatal_intercept],
                     (*VType)[i].NVegLayers))
        ReportError(KeyName[stomatal_intercept], 51);
    }
    
    if (Options->Hydraulics == HYD_KRSSUF) {
      if (!CopyFloat((*VType)[i].P50, VarStr[xylem_p50],
                     (*VType)[i].NVegLayers))
        ReportError(KeyName[xylem_p50], 51);
    }
    
    if (Options->StomScheme != JARVIS) {
      if (!CopyFloat((*VType)[i].Vcmax25, VarStr[max_carboxylation],
                     (*VType)[i].NVegLayers))
        ReportError(KeyName[max_carboxylation], 51);
    }

    if (Options->StomScheme != JARVIS) {
      for (j = 0; j < (*VType)[i].NVegLayers; j++)
        if ((*VType)[i].Vcmax25[j] <= 0.0)
          ReportError("InitVegTable: Vcmax25 must be > 0", 51);
    }
    if (Options->Hydraulics == HYD_KRSSUF) {
      for (j = 0; j < (*VType)[i].NVegLayers; j++)
        if ((*VType)[i].P50[j] >= 0.0)
          ReportError("InitVegTable: P50 must be negative (MPa)", 51);
    }

    /* ---- optional per-class physiology keys (Phase-VI refactor) --------
       Every key below may be omitted; the value then falls back to the
       documented default or derivation.  Photosynthetic traits are read for
       every scheme but Jarvis, hydraulic traits only with KRSSUF. */
    if (Options->StomScheme != JARVIS) {
      int n = (*VType)[i].NVegLayers;

      ReadOptionalFloats((*VType)[i].LeafWidth, VarStr[leaf_width],
                         KeyName[leaf_width], n, (float)PHOTO_LEAF_WIDTH);
      ReadOptionalFloats((*VType)[i].JmaxRatio, VarStr[jmax_vcmax_ratio],
                         KeyName[jmax_vcmax_ratio], n, (float)PHOTO_JMAXRATIO);
      ReadOptionalFloats((*VType)[i].Rd25Ratio, VarStr[rd_vcmax_ratio],
                         KeyName[rd_vcmax_ratio], n, (float)PHOTO_RD25RATIO);
      ReadOptionalFloats((*VType)[i].Gmax, VarStr[max_stomatal_conductance],
                         KeyName[max_stomatal_conductance], n, -1.0f);

      for (j = 0; j < n; j++) {
        if ((*VType)[i].LeafWidth[j] <= 0.0f)
          ReportError("InitVegTable: LEAF WIDTH must be > 0 (m)", 51);
        if ((*VType)[i].JmaxRatio[j] <= 0.0f)
          (*VType)[i].JmaxRatio[j] = (float)PHOTO_JMAXRATIO;
        if ((*VType)[i].Rd25Ratio[j] <= 0.0f)
          (*VType)[i].Rd25Ratio[j] = (float)PHOTO_RD25RATIO;

        /* Cap on stomatal conductance.  Default: the class's minimum stomatal
           resistance, which DHSVM already reads for Jarvis, converted from
           s/m to mol H2O/m2 leaf/s at 20 degC and standard pressure.  It is
           the same physiological quantity (gsmax = 1/RsMin). */
        if ((*VType)[i].Gmax[j] <= 0.0f) {
          float RsMin = (*VType)[i].RsMin[j];
          (*VType)[i].Gmax[j] = (RsMin > 0.0f)
            ? (1.0f / RsMin) / PhotoMolarToVelocity(20.0f, (float)PHOTO_P0)
            : 0.0f;
        }
      }
    }

    if (Options->Hydraulics == HYD_KRSSUF) {
      int n = (*VType)[i].NVegLayers;

      ReadOptionalFloats((*VType)[i].VulnShape, VarStr[xylem_shape],
                         KeyName[xylem_shape], n, -1.0f);
      ReadOptionalCurveForm((*VType)[i].VulnForm, VarStr[xylem_curve],
                            KeyName[xylem_curve], n);
      ReadOptionalFloats((*VType)[i].KxMax, VarStr[xylem_conductance],
                         KeyName[xylem_conductance], n, -1.0f);
      ReadOptionalFloats((*VType)[i].Krs, VarStr[root_conductance],
                         KeyName[root_conductance], n, -1.0f);
      ReadOptionalFloats((*VType)[i].Kcomp, VarStr[root_compensation],
                         KeyName[root_compensation], n, -1.0f);
      ReadOptionalFloats((*VType)[i].PsiCollarMin, VarStr[collar_pressure_min],
                         KeyName[collar_pressure_min], n, 0.0f);
      ReadOptionalFloats((*VType)[i].CiCaTarget, VarStr[cica_target],
                         KeyName[cica_target], n, (float)HYD_CICA_TARGET);

      /* Perirhizal (rhizosphere) stage, Leitner et al. 2025 / Vanderborght
         et al. 2021: on for a layer when its root length index is > 0.
         The config gives km of fine root per m2 ground (Jackson et al. 1997
         report values of order 3-8 km/m2 for temperate and boreal forests);
         stored in m/m2. */
      ReadOptionalFloats((*VType)[i].RootLengthIndex, VarStr[root_length_index],
                         KeyName[root_length_index], n, 0.0f);
      ReadOptionalFloats((*VType)[i].RootRadius, VarStr[fine_root_radius],
                         KeyName[fine_root_radius], n, (float)ROOT_RADIUS_DEFAULT);
      for (j = 0; j < n; j++) {
        if ((*VType)[i].RootLengthIndex[j] < 0.0f)
          (*VType)[i].RootLengthIndex[j] = 0.0f;
        (*VType)[i].RootLengthIndex[j] *= 1000.0f;          /* km -> m   */
        if ((*VType)[i].RootRadius[j] <= 0.0f)
          (*VType)[i].RootRadius[j] = (float)ROOT_RADIUS_DEFAULT;
      }
    }
    
    /* Derive the hydraulic traits once per vegetation class.

       Two things are frozen here rather than recomputed per pixel per
       timestep.  HydKmaxFromVcmax() runs a bisection that builds a supply
       curve and solves the optimization at every step; and the Kirchhoff
       cumulant (Sperry & Love 2015) is a 512-point integral of the
       vulnerability curve.  Both depend only on traits, not on soil state,
       so both belong here.  That is a large part of why the hydraulic
       schemes are affordable at 3 m and hourly.

       Computed at standard pressure: kmax varies under 5% from sea level to
       3000 m, well inside the uncertainty in P50 itself.  This freezes only
       the COORDINATION; the elevation dependence of CO2 diffusion, the
       mol-to-m/s conversion and the mole-fraction kinetics all still use the
       live LocalMet->Press. */
    if (Options->Hydraulics == HYD_KRSSUF) {
      HYDCURVE Curve;

      for (j = 0; j < (*VType)[i].NVegLayers; j++) {

        PHOTOTRAIT Ptrait;
        int KxConfigured;

        /* Vulnerability curve.  Weibull is the default so the Sperry et al.
           (2017) benchmarks keep reproducing; the sigmoidal form of Eller
           (2020) Eqn 2 / CLM5-PHS is selected per class by XYLEM CURVE.
           Shape defaults to c = 3 (Weibull) or a = 3 (sigmoidal). */
        if ((*VType)[i].VulnShape[j] <= 0.0)
          (*VType)[i].VulnShape[j] = (float)HYD_WEIBULL_C_DEFAULT;

        if ((*VType)[i].VulnForm[j] == HYD_SIGMOIDAL)
          HydCurveSigmoidal(&Curve, (*VType)[i].P50[j],
                            (*VType)[i].VulnShape[j]);
        else
          HydCurveWeibullFromP50(&Curve, (*VType)[i].P50[j],
                                 (*VType)[i].VulnShape[j]);

        (*VType)[i].WeibullB[j] = Curve.P1;

        PhotoTraitDefaults(&Ptrait, (*VType)[i].Vcmax25[j]);
        Ptrait.JmaxRatio = (*VType)[i].JmaxRatio[j];
        Ptrait.Rd25Ratio = (*VType)[i].Rd25Ratio[j];
        Ptrait.LeafWidth = (*VType)[i].LeafWidth[j];

        /* KxMax (per m2 leaf): configured (XYLEM CONDUCTANCE), or derived
           from the Ci/Ca coordination at the class's CICA TARGET (default
           0.70, Sperry et al. 2017).  Sensitivity testing shows this is the
           dominant control on transpiration magnitude (-30%/+52% for
           0.60/0.80) -- see docs/provenance.md sec. 3.  Treat it as a
           calibration target, and prefer a measured value where one exists. */
        KxConfigured = ((*VType)[i].KxMax[j] > 0.0);
        if ((*VType)[i].KxMax[j] <= 0.0) {
          (*VType)[i].KxMax[j] =
            StomKmaxFromVcmax(&Ptrait, (*VType)[i].P50[j],
                              (*VType)[i].VulnShape[j], (*VType)[i].VulnForm[j],
                              (*VType)[i].CiCaTarget[j], (*VType)[i].Height[j],
                              ATMOS_CO2, (float)PHOTO_P0);

          /* Refuse rather than run.  The previous derivation could silently
             return its bracket bound, which produced a plant roughly four
             orders of magnitude too resistant and a canopy that never opened
             -- and looked like a plausible drought response. */
          if ((*VType)[i].KxMax[j] <= 0.0)
            ReportError("InitVegTable: no Ci/Ca coordination point for this "
                        "Vcmax25/P50 pair; specify KxMax explicitly", 51);
        }

        /* Krs (per m2 leaf): configured (ROOT CONDUCTANCE), or the
           ROOT_KRS_FRAC fallback -- roots carry half the soil-to-canopy
           resistance (Sperry et al. 2017), so Krs = KxMax.  Still the
           least-constrained number in the scheme (open item 5). */
        if ((*VType)[i].Krs[j] <= 0.0)
          (*VType)[i].Krs[j] = RootKrsFromXylem((*VType)[i].KxMax[j]);

        /* Couvreur et al. (2012): Kcomp ~ Krs for many root systems, which is
           the default when ROOT COMPENSATION CONDUCTANCE is absent or
           negative.  An explicit 0 disables compensation and reduces the
           scheme EXACTLY to prescribed-root-fraction uptake.  (Before the
           refactor the unset value was 0, so compensation was silently
           off.) */
        if ((*VType)[i].Kcomp[j] < 0.0)
          (*VType)[i].Kcomp[j] = (*VType)[i].Krs[j];

        /* Dirichlet switch (Leitner et al. 2025 sec. 2.3).  1.5 x P50 is a
           placeholder, not a measurement. */
        if ((*VType)[i].PsiCollarMin[j] >= 0.0)
          (*VType)[i].PsiCollarMin[j] = 1.5f * (*VType)[i].P50[j];

        HydBuildCumulant(&(*VType)[i].Xylem[j], &Curve,
                         (*VType)[i].KxMax[j]);
        /* Hydrostatic drop over the canopy height, shared by every scheme. */
        HydSetCanopyHeight(&(*VType)[i].Xylem[j], (*VType)[i].Height[j]);

        /* Plausibility check against Eller et al. (2020) Fig. 1, which uses
           r_p,min of 1-2 mmol-1 m2 s MPa.  A total soil-plant conductance
           outside 0.02-50 mmol m-2 s-1 MPa-1 means something is wrong with
           the traits, and the symptom (a canopy that never opens, or one
           that never closes) is easy to mistake for physics. */
        {
          float KrsL  = (*VType)[i].Krs[j];
          float Ktot  = 1.0f / (1.0f / KrsL + 1.0f / (*VType)[i].KxMax[j]);
          if (Ktot < 0.02f || Ktot > 50.0f)
            printf("WARNING: veg class %d layer %d has soil-plant "
                   "conductance %.4g mmol/m2/s/MPa, outside the plausible "
                   "0.02-50 range (Eller et al. 2020 Fig. 1).  "
                   "KxMax=%.4g Krs=%.4g\n",
                   i, j, Ktot, (*VType)[i].KxMax[j], KrsL);

          /* One line per class and layer so a run's derived traits are on
             record with its output. */
          printf("Veg %d layer %d hydraulics: %s P50=%.2f shape=%.2f Pcrit=%.2f "
                 "MPa | KxMax=%.3f Krs=%.3f Kcomp=%.3f mmol/m2leaf/s/MPa | "
                 "Kplant=%.3f | PsiCollarMin=%.2f | gsmax=%.3f mol/m2/s | "
                 "Jmax:Vcmax=%.2f Rd:Vcmax=%.3f leaf width=%.3f m | "
                 "rho g h=%.3f MPa\n",
                 i + 1, j,
                 ((*VType)[i].VulnForm[j] == HYD_SIGMOIDAL) ? "sigmoidal" : "Weibull",
                 (*VType)[i].P50[j], (*VType)[i].VulnShape[j],
                 (*VType)[i].Xylem[j].Pcrit, (*VType)[i].KxMax[j], KrsL,
                 (*VType)[i].Kcomp[j], Ktot, (*VType)[i].PsiCollarMin[j],
                 (*VType)[i].Gmax[j], (*VType)[i].JmaxRatio[j],
                 (*VType)[i].Rd25Ratio[j], (*VType)[i].LeafWidth[j],
                 (*VType)[i].Xylem[j].Pgrav);

          /* Where the class sits at the Sperry reference state (PAR 2000,
             25 degC, D 1 kPa, wet soil), whether KxMax was derived (then
             Ci/Ca equals the target) or configured (then this is the Ci/Ca
             the configured kmax implies). */
          {
            STOMSOLUTION Ref;
            float EcritRef;
            StomReferencePoint(&Ptrait, &(*VType)[i].Xylem[j], KrsL,
                               ATMOS_CO2, (float)PHOTO_P0, &Ref, &EcritRef);
            if (!Ref.Failed)
              printf("Veg %d layer %d reference state (%s KxMax): "
                     "gs(H2O)=%.3f gs(CO2)=%.3f mol/m2/s  E=%.2f mmol/m2 leaf/s  Ecrit=%.2f  "
                     "E/Ecrit=%.2f  Ci/Ca=%.3f  psi_leaf=%.2f MPa  An=%.1f umol/m2/s\n",
                     i + 1, j, KxConfigured ? "configured" : "derived",
                     Ref.Gw, Ref.Gs, Ref.E, EcritRef,
                     (EcritRef > 0.0f) ? Ref.E / EcritRef : 0.0f,
                     Ref.Ci / ATMOS_CO2, Ref.PsiLeaf, Ref.An);
            else
              printf("Veg %d layer %d reference state: ProfitMax found no "
                     "solution at the reference conditions\n", i + 1, j);
          }
        }
      }
    }

    if (!CopyFloat((*VType)[i].MoistThres, VarStr[moisture_threshold],
      (*VType)[i].NVegLayers))
      ReportError(KeyName[moisture_threshold], 51);

    if (!CopyFloat((*VType)[i].VpdThres, VarStr[vpd], (*VType)[i].NVegLayers))
      ReportError(KeyName[vpd], 51);

    if (!CopyFloat((*VType)[i].Rpc, VarStr[rpc], (*VType)[i].NVegLayers))
      ReportError(KeyName[rpc], 51);

    if (!CopyFloat((*VType)[i].RootDepth, VarStr[root_zone_depth],
      (*VType)[i].NSoilLayers))
      ReportError(KeyName[root_zone_depth], 51);

    /* Calculate the wind speed profiles and the aerodynamical resistances
    for each layer.  The values are normalized for a reference height wind
    speed of 1 m/s, and are adjusted each timestep using actual reference
    height wind speeds */
    CalcAerodynamic((*VType)[i].NVegLayers, (*VType)[i].OverStory,
      (*VType)[i].Cn, (*VType)[i].Height, (*VType)[i].Trunk,
      (*VType)[i].U, &((*VType)[i].USnow), (*VType)[i].Ra,
      &((*VType)[i].RaSnow));

    /* Run the improved radiation scheme in which the tree height, solar altitude and fractional coverage
    are all taken into account into the radiation calculation */
    if (Options->ImprovRadiation == TRUE) {
      if ((*VType)[i].OverStory == TRUE) {
        if (!CopyFloat((*VType)[i].MonthlyExtnCoeff, VarStr[monextn], 12))
          ReportError(KeyName[monextn], 51);
        if (!CopyFloat(&((*VType)[i].VfAdjust), VarStr[vf_adj], 1))
          ReportError(KeyName[vf_adj], 51);
        (*VType)[i].Vf = (*VType)[i].Fract[0] * (*VType)[i].VfAdjust;
      }
      else {
        if ((*VType)[i].UnderStory == TRUE) {
          for (k = 0; k < 12; k++)
            (*VType)[i].MonthlyExtnCoeff[k] = 0;
          (*VType)[i].VfAdjust = 1.0;
          /* assuming 100% coverage if understory=TRUE & overstory=FALSE*/
          (*VType)[i].Vf = (*VType)[i].Fract[0] * (*VType)[i].VfAdjust;
        }
      }
    }
    
    (*VType)[i].TotalDepth = 0.0;
    for (y = 0; y < (*VType)[i].NSoilLayers; y++)
      (*VType)[i].TotalDepth += (*VType)[i].RootDepth[y];
    
  } /* end of the VEG TYPE loop */

  if (impervious) {
    GetInitString(SectionName, "IMPERVIOUS SURFACE ROUTING FILE", "", VarStr[0],
      (unsigned long)BUFSIZE, Input);
    if (IsEmptyStr(VarStr[0]))
      ReportError("IMPERVIOUS SURFACE ROUTING FILE", 51);
    strcpy(Options->ImperviousFilePath, VarStr[veg_description]);
  }

  return NVegs;
}

/********************************************************************************
InitLakeTable()
********************************************************************************/
int InitLakeTable(LAKETABLE **LType, 
  LISTPTR Input, OPTIONSTRUCT *Options)
{
  const char *Routine = "InitLakeTable";
  int i;			/* counter */
  int j;			/* counter */
  int NumLakes;			/* Number of unique lakes */
  int OutId;
  char KeyName[lake_exponent + 1][BUFSIZE + 1];
  char *KeyStr[] = {
    "LAKE NAME",
    "OUTFLOW CHANNEL",
    "POWER LAW SCALE",
    "POWER LAW EXPONENT"
  };
  char SectionName[] = "TERRAIN";
  char VarStr[lake_exponent + 1][BUFSIZE + 1];

  /* Get the number of unique lakes */
  GetInitString(SectionName, "NUMBER OF LAKES", "", VarStr[0],
    (unsigned long)BUFSIZE, Input);
  if (!CopyInt(&NumLakes, VarStr[0], 1))
    ReportError("NUMBER OF LAKES", 51);

  if (NumLakes == 0)
    return NumLakes;

  if (!(*LType = (LAKETABLE *)calloc(NumLakes, sizeof(LAKETABLE))))
    ReportError((char *)Routine, 1);

  /********** Read information and allocate memory for each lake *********/

  for (i = 0; i < NumLakes; i++) {

    /* Read the key-entry pairs from the input file */
    for (j = 0; j <= lake_exponent; j++) {
      sprintf(KeyName[j], "%s %d", KeyStr[j], i + 1);
      GetInitString(SectionName, KeyName[j], "", VarStr[j],
        (unsigned long)BUFSIZE, Input);
    }

    /* Assign the entries to the appropriate variables */
    if (IsEmptyStr(VarStr[lake_name]))
      ReportError(KeyName[lake_name], 51);

    strcpy((*LType)[i].Name, VarStr[lake_name]);
    (*LType)[i].Index = i+1;
    
    if (!CopyInt(&(OutId), VarStr[lake_outlet], 1))
      ReportError(KeyName[lake_outlet], 51);
    (*LType)[i].OutletID = OutId;
    
    if (!CopyFloat(&((*LType)[i].PowLawScale), VarStr[lake_scale], 1))
      ReportError(KeyName[lake_scale], 51);

    if (!CopyFloat(&((*LType)[i].PowLawExponent), VarStr[lake_exponent], 1))
      ReportError(KeyName[lake_exponent], 51);
  }
  
  return NumLakes;
}

