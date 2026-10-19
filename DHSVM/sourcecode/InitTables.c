
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
static int InitOverstoryTable(OVERTABLE **OType, const char *Path, OPTIONSTRUCT *Options,
                              int SoilLayers);
static int InitUnderstoryTable(UNDERTABLE **UType, const char *Path, OPTIONSTRUCT *Options,
                               int SoilLayers);

/*******************************************************************************/
/*				  InitTables()                                 */
/*******************************************************************************/
void InitTables(int StepsPerDay, LISTPTR Input, OPTIONSTRUCT *Options,
  MAPSIZE *Map, SOILTABLE **SType, LAYER *Soil, VEGTABLES *Veg,
  LAKETABLE **LType, TIMESTRUCT *Time)
{
  char Path[BUFSIZE + 1];

  Message("Initializing tables\n");

  RequireConfigValue(Input, "INPUTS", "Soil Table", Path);
  Soil->NTypes = InitSoilTable(Options, SType, Path, Soil, Map->SoilLayers, Time);

  RequireConfigValue(Input, "INPUTS", "Overstory Table", Path);
  Veg->NOver = InitOverstoryTable(&(Veg->Over), Path, Options, Map->SoilLayers);

  RequireConfigValue(Input, "INPUTS", "Understory Table", Path);
  Veg->NUnder = InitUnderstoryTable(&(Veg->Under), Path, Options, Map->SoilLayers);

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
  TableCheckRows(&T, "Class", 1, TABLE_MAXROWS);
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

/* Columns of the story parameters of a vegetation table */
typedef struct {
  int LAIFraction[12];
  int *RootFract;
  int Albedo[12];
  int RsMax, RsMin, Moist, Vpd, Rpc;
} STORYCOLUMNS;

static void StoryColumns(TABLE *T, int SoilLayers, STORYCOLUMNS *C)
{
  if (!(C->RootFract = calloc(SoilLayers, sizeof(int))))
    ReportError("StoryColumns", 1);
  TableLayerColumns(T, "LAIFraction", MAP_MONTH, SoilLayers, "fraction", TRUE, C->LAIFraction);
  TableLayerColumns(T, "RootFraction", MAP_ONE, SoilLayers, "fraction", TRUE, C->RootFract);
  TableLayerColumns(T, "Albedo", MAP_MONTH, SoilLayers, "fraction", TRUE, C->Albedo);
  C->RsMax = TableColumn(T, "MaxResistance", "s/m", TRUE);
  C->RsMin = TableColumn(T, "MinResistance", "s/m", TRUE);
  C->Moist = TableColumn(T, "MoistureThreshold", "m3/m3", TRUE);
  C->Vpd = TableColumn(T, "VpdThreshold", "Pa", TRUE);
  C->Rpc = TableColumn(T, "Rpc", "W/m2", TRUE);
}

/* The story parameters of row i */
static void ReadStory(TABLE *T, int i, STORYCOLUMNS *C, int SoilLayers, STORYTABLE *Story)
{
  if (!(Story->RootFract = (float *) calloc(SoilLayers, sizeof(float))))
    ReportError("ReadStory", 1);
  TableFloats(T, i, C->LAIFraction, 12, Story->LAIFraction);
  TableFloats(T, i, C->RootFract, SoilLayers, Story->RootFract);
  TableFloats(T, i, C->Albedo, 12, Story->AlbedoMonthly);
  Story->RsMax = TableFloat(T, i, C->RsMax);
  Story->RsMin = TableFloat(T, i, C->RsMin);
  Story->MoistThres = TableFloat(T, i, C->Moist);
  Story->VpdThres = TableFloat(T, i, C->Vpd);
  Story->Rpc = TableFloat(T, i, C->Rpc);
}

/* The name of row i */
static void ReadName(TABLE *T, int i, int CName, char *Desc)
{
  if (strlen(TableText(T, i, CName)) > BUFSIZE)
    ReportErrorCode(106, "%s: line %d: column Name is longer than %d characters", T->Path,
                    i + 3, BUFSIZE);
  strcpy(Desc, TableText(T, i, CName));
}

/********************************************************************************
Function Name: InitOverstoryTable()

Purpose      : Initialize the overstory lookup table from the overstory class
               table (one row per class, Class = row number, Veg.OverstoryClass
               0 = no overstory)

Returns      : Number of overstory classes
********************************************************************************/
static int InitOverstoryTable(OVERTABLE **OType, const char *Path, OPTIONSTRUCT *Options,
                              int SoilLayers)
{
  const char *Routine = "InitOverstoryTable";
  TABLE T;
  int i;
  int NRows;
  int CName, CTrunk, CAeroAtt, CMaxSnowInt, CDrip, CSnowEff, CExtn;
  STORYCOLUMNS CStory;
  OVERTABLE *Over;

  ReadTable(Path, &T);
  TableCheckRows(&T, "Class", 0, TABLE_MAXROWS);
  NRows = T.NRows;

  CName = TableColumn(&T, "Name", "-", TRUE);
  CTrunk = TableColumn(&T, "TrunkSpace", "fraction", TRUE);
  CAeroAtt = TableColumn(&T, "AerodynamicAttenuation", "-", TRUE);
  CMaxSnowInt = TableColumn(&T, "MaxSnowInterception", "m", TRUE);
  CDrip = TableColumn(&T, "MassReleaseDripRatio", "-", TRUE);
  CSnowEff = TableColumn(&T, "SnowInterceptionEff", "fraction", TRUE);
  StoryColumns(&T, SoilLayers, &CStory);
  CExtn = TableColumn(&T, "LightExtinction", "-", TRUE);
  TableWarnUnknown(&T);

  if (!(*OType = (OVERTABLE *) calloc((NRows > 0) ? NRows : 1, sizeof(OVERTABLE))))
    ReportError((char *)Routine, 1);

  for (i = 0; i < NRows; i++) {
    Over = &((*OType)[i]);
    ReadName(&T, i, CName, Over->Desc);
    ReadStory(&T, i, &CStory, SoilLayers, &(Over->Story));

    Over->Trunk = TableFloat(&T, i, CTrunk);
    Over->Cn = TableFloat(&T, i, CAeroAtt);
    Over->MaxSnowInt = TableFloat(&T, i, CMaxSnowInt);
    Over->MDRatio = TableFloat(&T, i, CDrip);
    Over->SnowIntEff = TableFloat(&T, i, CSnowEff);
    Over->ExtnCoeff = TableFloat(&T, i, CExtn);
  }

  free(CStory.RootFract);
  FreeTable(&T);
  return NRows;
}

/********************************************************************************
Function Name: InitUnderstoryTable()

Purpose      : Initialize the understory lookup table from the understory
               class table (one row per class, Class = row number,
               Veg.UnderstoryClass 0 = no understory). A class holds the
               impervious fraction of its cells whatever their understory
               cover.

Returns      : Number of understory classes
********************************************************************************/
static int InitUnderstoryTable(UNDERTABLE **UType, const char *Path, OPTIONSTRUCT *Options,
                               int SoilLayers)
{
  const char *Routine = "InitUnderstoryTable";
  TABLE T;
  int i;
  int NRows;
  int CName, CAtten, CImperv, CDetFrac, CDetDecay;
  STORYCOLUMNS CStory;
  UNDERTABLE *Under;

  ReadTable(Path, &T);
  TableCheckRows(&T, "Class", 0, TABLE_MAXROWS);
  NRows = T.NRows;

  CName = TableColumn(&T, "Name", "-", TRUE);
  StoryColumns(&T, SoilLayers, &CStory);
  CAtten = TableColumn(&T, "RadiationAttenuation", "-", TRUE);
  CImperv = TableColumn(&T, "ImperviousFraction", "fraction", TRUE);

  /* Some class with an impervious fraction */
  Options->Impervious = FALSE;
  for (i = 0; i < NRows; i++) {
    if (TableFloat(&T, i, CImperv) > 0.0)
      Options->Impervious = TRUE;
  }
  CDetFrac = CDetDecay = -1;
  if (Options->Impervious) {
    CDetFrac = TableColumn(&T, "DetentionFraction", "fraction", TRUE);
    CDetDecay = TableColumn(&T, "DetentionDecay", "-", TRUE);
  }
  else {
    TableKnown(&T, "DetentionFraction");
    TableKnown(&T, "DetentionDecay");
  }
  TableWarnUnknown(&T);

  if (!(*UType = (UNDERTABLE *) calloc((NRows > 0) ? NRows : 1, sizeof(UNDERTABLE))))
    ReportError((char *)Routine, 1);

  for (i = 0; i < NRows; i++) {
    Under = &((*UType)[i]);
    ReadName(&T, i, CName, Under->Desc);
    ReadStory(&T, i, &CStory, SoilLayers, &(Under->Story));

    Under->RadAtten = TableFloat(&T, i, CAtten);
    Under->ImpervFrac = TableFloat(&T, i, CImperv);
    if (Under->ImpervFrac > 0) {
      Under->DetentionFrac = TableFloat(&T, i, CDetFrac);
      Under->DetentionDecay = TableFloat(&T, i, CDetDecay);
    }
    else {
      Under->DetentionFrac = 0.;
      Under->DetentionDecay = 0.;
    }
  }

  free(CStory.RootFract);
  FreeTable(&T);
  return NRows;
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
  TableCheckRows(&T, "Lake", 1, TABLE_MAXROWS);
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
