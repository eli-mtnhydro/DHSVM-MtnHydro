
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
#include "mapio.h"
#include "tables.h"

static float *TableFloats(TABLE *T, int Row, int *Columns, int N, float *Values);

/*******************************************************************************/
/*				  InitTables()                                 */
/*******************************************************************************/
void InitTables(int StepsPerDay, LISTPTR Input, OPTIONSTRUCT *Options,
  MAPSIZE *Map, SOILTABLE **SType, LAYER *Soil, VEGTABLE **VType, LAYER *Veg,
  LAKETABLE **LType, TIMESTRUCT *Time)
{
  char Path[BUFSIZE + 1];

  Message("Initializing tables\n");

  RequireConfigValue(Input, "INPUTS", "Soil Table", Path);
  Soil->NTypes = InitSoilTable(Options, SType, Path, Soil, Map->SoilLayers, Time);

  RequireConfigValue(Input, "INPUTS", "Vegetation Table", Path);
  Veg->NTypes = InitVegTable(VType, Path, Options, Veg, Map->SoilLayers);

  if (Options->LakeDynamics) {
    RequireConfigValue(Input, "INPUTS", "Lake Table", Path);
    Map->NumLakes = InitLakeTable(LType, Path, Options);
  }
  else {
    GetConfigValue(Input, "INPUTS", "Lake Table", Path);
    Map->NumLakes = 0;
  }

  InitSatVaporTable();
}

/* Values of the columns Columns[0..N-1] of a row */
static float *TableFloats(TABLE *T, int Row, int *Columns, int N, float *Values)
{
  int k;

  for (k = 0; k < N; k++)
    Values[k] = TableFloat(T, Row, Columns[k]);
  return Values;
}

/********************************************************************************
Function Name: InitSoilTable()

Purpose      : Initialize the soil lookup table from the soil class table
               (one row per class, Class = row number). Every class has the
               L root zone layers and the deep layer.

Returns      : Number of soil classes
********************************************************************************/
int InitSoilTable(OPTIONSTRUCT *Options, SOILTABLE **SType, const char *Path,
  LAYER *Soil, int SoilLayers, TIMESTRUCT *Time)
{
  const char *Routine = "InitSoilTable";
  TABLE T;
  int i;
  int NSoils;
  int NLayers = SoilLayers + 1;	/* root zone layers and the deep layer */
  int CName, CDepthThresh, CAnisotropy, CMaxInfilt, CCapillary, CDeepFlux;
  int CAlbedo, CManning;
  int *CPoreDist, *CPress, *CWP, *CKs;

  ReadTable(Path, &T);
  TableCheckRows(&T, "Class", TABLE_MAXROWS);
  NSoils = T.NRows;

  CPoreDist = calloc(NLayers, sizeof(int));
  CPress = calloc(NLayers, sizeof(int));
  CWP = calloc(NLayers, sizeof(int));
  CKs = calloc(NLayers, sizeof(int));
  if (!CPoreDist || !CPress || !CWP || !CKs)
    ReportError((char *)Routine, 1);

  CName = TableColumn(&T, "Name", "-", TRUE);
  CDepthThresh = TableColumn(&T, "DepthThreshold", "m", TRUE);
  CAnisotropy = TableColumn(&T, "VerticalAnisotropy", "-", TRUE);
  CDeepFlux = TableColumn(&T, "DeepFlux", "m/yr", TRUE);
  CAlbedo = TableColumn(&T, "Albedo", "fraction", TRUE);
  TableLayerColumns(&T, "PoreSizeDistribution", MAP_SOIL, SoilLayers, "-", TRUE, CPoreDist);
  TableLayerColumns(&T, "BubblingPressure", MAP_SOIL, SoilLayers, "m", TRUE, CPress);
  TableLayerColumns(&T, "WiltingPoint", MAP_SOIL, SoilLayers, "m3/m3", TRUE, CWP);
  CMaxInfilt = CCapillary = CManning = -1;
  if (!Options->UseKsatAnisotropy) {
    CMaxInfilt = TableColumn(&T, "MaxInfiltration", "m/s", TRUE);
    TableLayerColumns(&T, "VerticalConductivity", MAP_SOIL, SoilLayers, "m/s", TRUE, CKs);
  }
  else {
    TableKnown(&T, "MaxInfiltration");
    TableKnownLayers(&T, "VerticalConductivity", MAP_SOIL, SoilLayers);
  }
  if (Options->Infiltration == DYNAMIC)
    CCapillary = TableColumn(&T, "CapillaryDrive", "m", TRUE);
  else
    TableKnown(&T, "CapillaryDrive");
  if (Options->Routing)
    CManning = TableColumn(&T, "ManningsN", "-", TRUE);
  else
    TableKnown(&T, "ManningsN");
  TableWarnUnknown(&T);

  if (!(Soil->NLayers = (int *)calloc(NSoils, sizeof(int))))
    ReportError((char *)Routine, 1);

  if (!(*SType = (SOILTABLE *)calloc(NSoils, sizeof(SOILTABLE))))
    ReportError((char *)Routine, 1);

  Soil->MaxLayers = SoilLayers;

  for (i = 0; i < NSoils; i++) {
    if (strlen(TableText(&T, i, CName)) > BUFSIZE)
      ReportErrorCode(106, "%s: line %d: column Name is longer than %d characters", Path,
                      i + 3, BUFSIZE);
    strcpy((*SType)[i].Desc, TableText(&T, i, CName));
    (*SType)[i].Index = i;
    (*SType)[i].NLayers = SoilLayers;
    Soil->NLayers[i] = SoilLayers;

    (*SType)[i].DepthThresh = TableFloat(&T, i, CDepthThresh);
    (*SType)[i].KsAnisotropy = TableFloat(&T, i, CAnisotropy);

    if (!Options->UseKsatAnisotropy)
      (*SType)[i].MaxInfiltrationRate = TableFloat(&T, i, CMaxInfilt);

    if (Options->Infiltration == DYNAMIC)
      (*SType)[i].G_Infilt = TableFloat(&T, i, CCapillary);
    else
      (*SType)[i].G_Infilt = NOT_APPLICABLE;

    /* m/yr as read; ApplyTableParameters() adds the offset and converts it
       to m per time step */
    (*SType)[i].DeepFlux = TableFloat(&T, i, CDeepFlux);

    (*SType)[i].Albedo = TableFloat(&T, i, CAlbedo);

    if (Options->Routing)
      (*SType)[i].Manning = TableFloat(&T, i, CManning);

    /* Root zone layers and the deep layer */
    if (!((*SType)[i].PoreDist = (float *)calloc(NLayers, sizeof(float))) ||
        !((*SType)[i].Press = (float *)calloc(NLayers, sizeof(float))) ||
        !((*SType)[i].WP = (float *)calloc(NLayers, sizeof(float))) ||
        !((*SType)[i].Ks = (float *)calloc(NLayers, sizeof(float))))
      ReportError((char *)Routine, 1);

    TableFloats(&T, i, CPoreDist, NLayers, (*SType)[i].PoreDist);
    TableFloats(&T, i, CPress, NLayers, (*SType)[i].Press);
    TableFloats(&T, i, CWP, NLayers, (*SType)[i].WP);
    if (!Options->UseKsatAnisotropy)
      TableFloats(&T, i, CKs, NLayers, (*SType)[i].Ks);
  }

  free(CPoreDist);
  free(CPress);
  free(CWP);
  free(CKs);
  FreeTable(&T);
  return NSoils;
}

/********************************************************************************
Function Name: InitVegTable()

Purpose      : Initialize the vegetation lookup table from the vegetation
               class table (one row per class, Class = row number). The
               values of a class with an overstory and an understory are
               held in layer order: overstory first.

Returns      : Number of vegetation classes
********************************************************************************/
int InitVegTable(VEGTABLE **VType, const char *Path, OPTIONSTRUCT *Options, LAYER *Veg,
                 int SoilLayers)
{
  const char *Routine = "InitVegTable";
  TABLE T;
  int i;			/* Counter */
  int j;			/* Counter */
  int k;      /* counter */
  int y;      /* counter */
  int NVegs;		/* Number of vegetation types */
  int Layer;
  int CName, COver, CUnder, CImperv, CDetFrac, CDetDecay;
  int CTrunk, CAeroAtt, CMaxSnowInt, CDrip, CSnowEff;
  int CHemi, CClump, CLeafA, CLeafB, CScat, CAtten, CDiffAtten, CViewAdj;
  int *CRootDepth, *CRootFract[2];
  int CExtn[12], CLAI[12], CLAIFrac[12], CAlbedo[2][12];
  int CRsMax[2], CRsMin[2], CMoist[2], CVpd[2], CRpc[2];
  const char *Story[2] = { "Overstory", "Understory" };
  char Name[BUFSIZE + 1];
  int s;
  int Impervious = FALSE;

  CRootDepth = calloc(SoilLayers, sizeof(int));
  CRootFract[0] = calloc(SoilLayers, sizeof(int));
  CRootFract[1] = calloc(SoilLayers, sizeof(int));
  if (!CRootDepth || !CRootFract[0] || !CRootFract[1])
    ReportError((char *)Routine, 1);

  ReadTable(Path, &T);
  TableCheckRows(&T, "Class", TABLE_MAXROWS);
  NVegs = T.NRows;

  CName = TableColumn(&T, "Name", "-", TRUE);
  COver = TableColumn(&T, "Overstory", "-", TRUE);
  CUnder = TableColumn(&T, "Understory", "-", TRUE);
  CImperv = TableColumn(&T, "ImperviousFraction", "fraction", TRUE);
  TableLayerColumns(&T, "RootDepth", MAP_ONE, SoilLayers, "m", TRUE, CRootDepth);

  for (i = 0; i < NVegs; i++) {
    if (TableFloat(&T, i, CImperv) > 0)
      Impervious = TRUE;
  }
  CDetFrac = CDetDecay = -1;
  if (Impervious) {
    CDetFrac = TableColumn(&T, "DetentionFraction", "fraction", TRUE);
    CDetDecay = TableColumn(&T, "DetentionDecay", "-", TRUE);
  }
  else {
    TableKnown(&T, "DetentionFraction");
    TableKnown(&T, "DetentionDecay");
  }

  /* Overstory columns */
  CTrunk = TableColumn(&T, "TrunkSpace", "fraction", TRUE);
  CAeroAtt = TableColumn(&T, "AerodynamicAttenuation", "-", TRUE);
  CMaxSnowInt = TableColumn(&T, "MaxSnowInterception", "m", TRUE);
  CDrip = TableColumn(&T, "MassReleaseDripRatio", "-", TRUE);
  CSnowEff = TableColumn(&T, "SnowInterceptionEff", "fraction", TRUE);
  TableLayerColumns(&T, "OverstoryLAIFraction", MAP_MONTH, SoilLayers, "fraction", TRUE, CLAIFrac);
  CHemi = CClump = CLeafA = CLeafB = CScat = CAtten = CDiffAtten = CViewAdj = -1;
  if (Options->CanopyRadAtt == VARIABLE) {
    CHemi = TableColumn(&T, "HemiFractCoverage", "-", TRUE);
    CClump = TableColumn(&T, "ClumpingFactor", "-", TRUE);
    CLeafA = TableColumn(&T, "LeafAngleA", "-", TRUE);
    CLeafB = TableColumn(&T, "LeafAngleB", "-", TRUE);
    CScat = TableColumn(&T, "ScatteringParameter", "-", TRUE);
  }
  else {
    TableKnown(&T, "HemiFractCoverage");
    TableKnown(&T, "ClumpingFactor");
    TableKnown(&T, "LeafAngleA");
    TableKnown(&T, "LeafAngleB");
    TableKnown(&T, "ScatteringParameter");
  }
  if (Options->CanopyRadAtt == FIXED && Options->ImprovRadiation == FALSE)
    CAtten = TableColumn(&T, "RadiationAttenuation", "-", TRUE);
  else
    TableKnown(&T, "RadiationAttenuation");
  if (Options->CanopyRadAtt == FIXED && Options->ImprovRadiation == TRUE)
    CDiffAtten = TableColumn(&T, "DiffuseAttenuation", "-", TRUE);
  else
    TableKnown(&T, "DiffuseAttenuation");
  if (Options->ImprovRadiation == TRUE) {
    TableLayerColumns(&T, "LightExtinction", MAP_MONTH, SoilLayers, "1/m", TRUE, CExtn);
    CViewAdj = TableColumn(&T, "CanopyViewAdj", "-", TRUE);
  }
  else {
    TableKnownLayers(&T, "LightExtinction", MAP_MONTH, SoilLayers);
    TableKnown(&T, "CanopyViewAdj");
  }

  /* Columns of each story: overstory (s = 0), understory (s = 1) */
  for (s = 0; s < 2; s++) {
    snprintf(Name, sizeof(Name), "%sRootFraction", Story[s]);
    TableLayerColumns(&T, Name, MAP_ONE, SoilLayers, "fraction", TRUE, CRootFract[s]);
    snprintf(Name, sizeof(Name), "%sAlbedo", Story[s]);
    TableLayerColumns(&T, Name, MAP_MONTH, SoilLayers, "fraction", TRUE, CAlbedo[s]);
    snprintf(Name, sizeof(Name), "%sMaxResistance", Story[s]);
    CRsMax[s] = TableColumn(&T, Name, "s/m", TRUE);
    snprintf(Name, sizeof(Name), "%sMinResistance", Story[s]);
    CRsMin[s] = TableColumn(&T, Name, "s/m", TRUE);
    snprintf(Name, sizeof(Name), "%sMoistureThreshold", Story[s]);
    CMoist[s] = TableColumn(&T, Name, "m3/m3", TRUE);
    snprintf(Name, sizeof(Name), "%sVpdThreshold", Story[s]);
    CVpd[s] = TableColumn(&T, Name, "Pa", TRUE);
    snprintf(Name, sizeof(Name), "%sRpc", Story[s]);
    CRpc[s] = TableColumn(&T, Name, "W/m2", TRUE);
  }
  TableLayerColumns(&T, "UnderstoryLAI", MAP_MONTH, SoilLayers, "m2/m2", TRUE, CLAI);
  TableWarnUnknown(&T);

  if (!(Veg->NLayers = (int *)calloc(NVegs, sizeof(int))))
    ReportError((char *)Routine, 1);

  if (!(*VType = (VEGTABLE *)calloc(NVegs, sizeof(VEGTABLE))))
    ReportError((char *)Routine, 1);

  /******* Read information and allocate memory for each vegetation type ******/

  Veg->MaxLayers = 0;
  Options->Impervious = FALSE;
  for (i = 0; i < NVegs; i++) {

    if (strlen(TableText(&T, i, CName)) > BUFSIZE)
      ReportErrorCode(106, "%s: line %d: column Name is longer than %d characters", Path,
                      i + 3, BUFSIZE);
    strcpy((*VType)[i].Desc, TableText(&T, i, CName));
    (*VType)[i].Index = i;

    (*VType)[i].NVegLayers = 0;

    if (TableInt(&T, i, COver) != 0 && TableInt(&T, i, COver) != 1)
      ReportErrorCode(106, "%s: line %d: column Overstory is %s, expected 0 or 1", Path, i + 3,
                      TableText(&T, i, COver));
    if (TableInt(&T, i, CUnder) != 0 && TableInt(&T, i, CUnder) != 1)
      ReportErrorCode(106, "%s: line %d: column Understory is %s, expected 0 or 1", Path, i + 3,
                      TableText(&T, i, CUnder));
    (*VType)[i].OverStory = (unsigned char) TableInt(&T, i, COver);
    (*VType)[i].UnderStory = (unsigned char) TableInt(&T, i, CUnder);
    if ((*VType)[i].OverStory == TRUE)
      ((*VType)[i].NVegLayers)++;
    if ((*VType)[i].UnderStory == TRUE)
      ((*VType)[i].NVegLayers)++;

    Veg->NLayers[i] = (*VType)[i].NVegLayers;
    if ((*VType)[i].NVegLayers > Veg->MaxLayers)
      Veg->MaxLayers = (*VType)[i].NVegLayers;

    (*VType)[i].NSoilLayers = SoilLayers;

    (*VType)[i].ImpervFrac = TableFloat(&T, i, CImperv);
    if ((*VType)[i].ImpervFrac > 0.0)
      Options->Impervious = TRUE;

    if ((*VType)[i].ImpervFrac > 0) {
      (*VType)[i].DetentionFrac = TableFloat(&T, i, CDetFrac);
      (*VType)[i].DetentionDecay = TableFloat(&T, i, CDetDecay);
    }
    else {
      (*VType)[i].DetentionFrac = 0.;
      (*VType)[i].DetentionDecay = 0.;
    }

    /* allocate memory for the vegetation layers */
    if (Options->CanopyRadAtt == VARIABLE) {
      if (!((*VType)[i].HemiFract = (float *)calloc((*VType)[i].NVegLayers,
        sizeof(float))))
        ReportError((char *)Routine, 1);
    }
    else {
      (*VType)[i].HemiFract = NULL;
    }

    if (!((*VType)[i].RsMax = (float *)calloc((*VType)[i].NVegLayers, sizeof(float))))
      ReportError((char *)Routine, 1);

    if (!((*VType)[i].RsMin = (float *)calloc((*VType)[i].NVegLayers, sizeof(float))))
      ReportError((char *)Routine, 1);

    if (!((*VType)[i].MoistThres = (float *)calloc((*VType)[i].NVegLayers, sizeof(float))))
      ReportError((char *)Routine, 1);

    if (!((*VType)[i].VpdThres = (float *)calloc((*VType)[i].NVegLayers, sizeof(float))))
      ReportError((char *)Routine, 1);

    if (!((*VType)[i].Rpc = (float *)calloc((*VType)[i].NVegLayers, sizeof(float))))
      ReportError((char *)Routine, 1);

    if (!((*VType)[i].Albedo = (float *)calloc(((*VType)[i].NVegLayers + 1), sizeof(float))))
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

    /* assign the entries to the appropriate variables; the overstory LAI
       is a map times the monthly class fractions */
    if ((*VType)[i].OverStory == TRUE) {
      if (Options->CanopyRadAtt == VARIABLE) {
        (*VType)[i].HemiFract[0] = TableFloat(&T, i, CHemi);
        (*VType)[i].ClumpingFactor = TableFloat(&T, i, CClump);
        (*VType)[i].LeafAngleA = TableFloat(&T, i, CLeafA);
        (*VType)[i].LeafAngleB = TableFloat(&T, i, CLeafB);
        (*VType)[i].Scat = TableFloat(&T, i, CScat);
        (*VType)[i].Atten = NOT_APPLICABLE;
      }
      else if (Options->CanopyRadAtt == FIXED && Options->ImprovRadiation == FALSE) {
        (*VType)[i].Atten = TableFloat(&T, i, CAtten);
        (*VType)[i].ClumpingFactor = NOT_APPLICABLE;
        (*VType)[i].Scat = NOT_APPLICABLE;
        (*VType)[i].LeafAngleA = NOT_APPLICABLE;
        (*VType)[i].LeafAngleB = NOT_APPLICABLE;
        (*VType)[i].Taud = NOT_APPLICABLE;
      }
      else if (Options->ImprovRadiation == TRUE) {
        (*VType)[i].Taud = TableFloat(&T, i, CDiffAtten);
        (*VType)[i].Atten = NOT_APPLICABLE;
        (*VType)[i].ClumpingFactor = NOT_APPLICABLE;
        (*VType)[i].Scat = NOT_APPLICABLE;
        (*VType)[i].LeafAngleA = NOT_APPLICABLE;
        (*VType)[i].LeafAngleB = NOT_APPLICABLE;
      }

      (*VType)[i].Trunk = TableFloat(&T, i, CTrunk);
      (*VType)[i].Cn = TableFloat(&T, i, CAeroAtt);
      (*VType)[i].MaxSnowInt = TableFloat(&T, i, CMaxSnowInt);
      (*VType)[i].MDRatio = TableFloat(&T, i, CDrip);
      (*VType)[i].SnowIntEff = TableFloat(&T, i, CSnowEff);
      TableFloats(&T, i, CRootFract[0], (*VType)[i].NSoilLayers, (*VType)[i].RootFract[0]);
      TableFloats(&T, i, CAlbedo[0], 12, (*VType)[i].AlbedoMonthly[0]);
      TableFloats(&T, i, CLAIFrac, 12, (*VType)[i].LAIFraction);

      if ((*VType)[i].UnderStory == TRUE) {
        TableFloats(&T, i, CRootFract[1], (*VType)[i].NSoilLayers, (*VType)[i].RootFract[1]);
        TableFloats(&T, i, CLAI, 12, (*VType)[i].LAIMonthly[1]);
        TableFloats(&T, i, CAlbedo[1], 12, (*VType)[i].AlbedoMonthly[1]);
      }
    }
    else {
      if ((*VType)[i].UnderStory == TRUE) {
        TableFloats(&T, i, CRootFract[1], (*VType)[i].NSoilLayers, (*VType)[i].RootFract[0]);
        TableFloats(&T, i, CLAI, 12, (*VType)[i].LAIMonthly[0]);
        TableFloats(&T, i, CAlbedo[1], 12, (*VType)[i].AlbedoMonthly[0]);
      }
      (*VType)[i].Trunk = NOT_APPLICABLE;
      (*VType)[i].Cn = NOT_APPLICABLE;
      (*VType)[i].Atten = NOT_APPLICABLE;
      (*VType)[i].ClumpingFactor = NOT_APPLICABLE;
    }

    /* Per layer values in layer order, overstory first */
    Layer = 0;
    for (s = 0; s < 2; s++) {
      if ((s == 0 && (*VType)[i].OverStory == TRUE) || (s == 1 && (*VType)[i].UnderStory == TRUE)) {
        (*VType)[i].RsMax[Layer] = TableFloat(&T, i, CRsMax[s]);
        (*VType)[i].RsMin[Layer] = TableFloat(&T, i, CRsMin[s]);
        (*VType)[i].MoistThres[Layer] = TableFloat(&T, i, CMoist[s]);
        (*VType)[i].VpdThres[Layer] = TableFloat(&T, i, CVpd[s]);
        (*VType)[i].Rpc[Layer] = TableFloat(&T, i, CRpc[s]);
        Layer++;
      }
    }

    TableFloats(&T, i, CRootDepth, (*VType)[i].NSoilLayers, (*VType)[i].RootDepth);

    /* Run the improved radiation scheme in which the tree height, solar altitude and fractional coverage
    are all taken into account into the radiation calculation */
    if (Options->ImprovRadiation == TRUE) {
      if ((*VType)[i].OverStory == TRUE) {
        TableFloats(&T, i, CExtn, 12, (*VType)[i].MonthlyExtnCoeff);
        (*VType)[i].VfAdjust = TableFloat(&T, i, CViewAdj);
      }
      else {
        if ((*VType)[i].UnderStory == TRUE) {
          for (k = 0; k < 12; k++)
            (*VType)[i].MonthlyExtnCoeff[k] = 0;
          (*VType)[i].VfAdjust = 1.0;
        }
      }
    }

    (*VType)[i].TotalDepth = 0.0;
    for (y = 0; y < (*VType)[i].NSoilLayers; y++)
      (*VType)[i].TotalDepth += (*VType)[i].RootDepth[y];

  } /* end of the VEG TYPE loop */

  free(CRootDepth);
  free(CRootFract[0]);
  free(CRootFract[1]);
  FreeTable(&T);
  return NVegs;
}

/********************************************************************************
InitLakeTable()

The lake table: one row per lake (Lake = row number, the Lake.ID value),
its name, the segment it pours into and the scale and exponent of its
outflow power law.
********************************************************************************/
int InitLakeTable(LAKETABLE **LType, const char *Path, OPTIONSTRUCT *Options)
{
  const char *Routine = "InitLakeTable";
  TABLE T;
  int i;
  int NumLakes;
  int CName, COutlet, CScale, CExponent;

  ReadTable(Path, &T);
  TableCheckRows(&T, "Lake", TABLE_MAXROWS);
  NumLakes = T.NRows;
  CName = TableColumn(&T, "Name", "-", TRUE);
  COutlet = TableColumn(&T, "OutletSegment", "-", TRUE);
  CScale = TableColumn(&T, "OutflowScale", "-", TRUE);
  CExponent = TableColumn(&T, "OutflowExponent", "-", TRUE);
  TableWarnUnknown(&T);

  if (!(*LType = (LAKETABLE *)calloc(NumLakes, sizeof(LAKETABLE))))
    ReportError((char *)Routine, 1);

  for (i = 0; i < NumLakes; i++) {
    if (strlen(TableText(&T, i, CName)) > BUFSIZE)
      ReportErrorCode(106, "%s: line %d: column Name is longer than %d characters", Path,
                      i + 3, BUFSIZE);
    strcpy((*LType)[i].Name, TableText(&T, i, CName));
    (*LType)[i].Index = i+1;
    (*LType)[i].Outlet = TableInt(&T, i, COutlet);
    (*LType)[i].PowLawScale = TableFloat(&T, i, CScale);
    (*LType)[i].PowLawExponent = TableFloat(&T, i, CExponent);
  }
  strcpy(Options->LakeTablePath, Path);

  FreeTable(&T);
  return NumLakes;
}
