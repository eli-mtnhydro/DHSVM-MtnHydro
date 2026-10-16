
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "fileio.h"
#include "functions.h"
#include "constants.h"
#include "getinit.h"
#include "sizeofnt.h"
#include "slopeaspect.h"
#include "mapio.h"

static void *ReadMap(LISTPTR Input, const char *Name, int NumberType, const char *Units,
                     int LayerSet, MAPSIZE *Map);

/*****************************************************************************
  ReadMap()

  Allocates and reads every layer of the input map Name ([INPUTS] key);
  the header is checked first
*****************************************************************************/
static void *ReadMap(LISTPTR Input, const char *Name, int NumberType, const char *Units,
                     int LayerSet, MAPSIZE *Map)
{
  void *Values;
  size_t NValues = (size_t) Map->NY * Map->NX * MapLayerCount(LayerSet, Map->SoilLayers);

  if (!(Values = calloc(NValues, SizeOfNumberType(NumberType))))
    ReportError("ReadMap", 1);
  ReadInputMap(Input, Name, NumberType, Units, LayerSet, Map, Values);
  return Values;
}

/*****************************************************************************
  ReadInputMaps()

  Reads every input map the run uses, as given (the values are checked by
  ValidateInputs(), then ApplyParameters() applies the parameters). The
  monthly precipitation pattern stays open: its layers are read by
  ValidateInputs() and, a month at a time, by InitNewMonth().
*****************************************************************************/
void ReadInputMaps(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map, INPUTMAPS *Raw)
{
  char Path[BUFSIZE + 1];

  Message("\nReading input maps\n");

  memset(Raw, 0, sizeof(INPUTMAPS));
  Raw->Elevation = ReadMap(Input, "Terrain.Elevation", NC_FLOAT, "m", MAP_ONE, Map);
  Raw->Mask = ReadMap(Input, "Terrain.Mask", NC_BYTE, "-", MAP_ONE, Map);

  /* Horizon angles: the map is open until InitShading() has read it */
  Options->HorizonAngle.Ptr = NULL;
  if (Options->Shading) {
    RequireConfigValue(Input, "INPUTS", "Terrain.HorizonAngle", Path);
    OpenInputMap(&(Options->HorizonAngle), Path, "Terrain.HorizonAngle", NC_FLOAT,
                 "degrees", MAP_AZIMUTH, Map);
    Raw->SkyView = ReadMap(Input, "Terrain.SkyView", NC_FLOAT, "fraction", MAP_ONE, Map);
  }

  Raw->SoilClass = ReadMap(Input, "Soil.Class", NC_BYTE, "-", MAP_ONE, Map);
  Raw->Depth = ReadMap(Input, "Soil.Depth", NC_FLOAT, "m", MAP_ONE, Map);
  Raw->LateralConductivity =
    ReadMap(Input, "Soil.LateralConductivity", NC_FLOAT, "m/s", MAP_ONE, Map);
  Raw->ExponentialDecrease =
    ReadMap(Input, "Soil.ExponentialDecrease", NC_FLOAT, "1/m", MAP_ONE, Map);
  Raw->Porosity = ReadMap(Input, "Soil.Porosity", NC_FLOAT, "m3/m3", MAP_SOIL, Map);
  Raw->FieldCapacity = ReadMap(Input, "Soil.FieldCapacity", NC_FLOAT, "m3/m3", MAP_SOIL, Map);
  Raw->VegClass = ReadMap(Input, "Veg.Class", NC_BYTE, "-", MAP_ONE, Map);
  Raw->OverstoryCover = ReadMap(Input, "Veg.OverstoryCover", NC_FLOAT, "fraction", MAP_ONE, Map);
  Raw->OverstoryLAI = ReadMap(Input, "Veg.OverstoryLAI", NC_FLOAT, "m2/m2", MAP_ONE, Map);
  Raw->OverstoryHeight = ReadMap(Input, "Veg.OverstoryHeight", NC_FLOAT, "m", MAP_ONE, Map);
  Raw->UnderstoryHeight = ReadMap(Input, "Veg.UnderstoryHeight", NC_FLOAT, "m", MAP_ONE, Map);

  Options->PrecipPattern.Ptr = NULL;
  if (Options->Prism) {
    RequireConfigValue(Input, "INPUTS", "Met.PrecipPattern", Path);
    OpenInputMap(&(Options->PrecipPattern), Path, "Met.PrecipPattern", NC_FLOAT,
                 "-", MAP_MONTH, Map);
  }
  if (Options->SnowPattern)
    Raw->SnowPattern = ReadMap(Input, "Met.SnowPattern", NC_FLOAT, "-", MAP_ONE, Map);
  if (Options->MeltScale)
    Raw->MeltScale = ReadMap(Input, "Snow.MeltScale", NC_FLOAT, "-", MAP_ONE, Map);
  if (Options->LakeDynamics)
    Raw->LakeID = ReadMap(Input, "Lake.ID", NC_BYTE, "-", MAP_ONE, Map);
  if (!IsEmptyStr(Options->InitialSWEPath))
    Raw->InitialSWE = ReadMap(Input, "Snow.InitialSWE", NC_FLOAT, "m", MAP_ONE, Map);
}

/*****************************************************************************
  InitTopoMap()

  Elevation and basin mask of every cell, lake ID of the cells inside the
  mask; lake areas; slope, aspect and flow directions
*****************************************************************************/
void InitTopoMap(OPTIONSTRUCT *Options, MAPSIZE *Map, INPUTMAPS *Raw,
  TOPOPIX ***TopoMap, LAKETABLE *LType)
{
  const char *Routine = "InitTopoMap";
  int i;			/* Counter */
  int x;			/* Counter */
  int y;			/* Counter */
  float Area; /* Area of a lake */

  Message("\nInitializing terrain maps\n");

  if (!(*TopoMap = (TOPOPIX **)calloc(Map->NY, sizeof(TOPOPIX *))))
    ReportError((char *)Routine, 1);
  for (y = 0; y < Map->NY; y++) {
    if (!((*TopoMap)[y] = (TOPOPIX *)calloc(Map->NX, sizeof(TOPOPIX))))
      ReportError((char *)Routine, 1);
  }

  for (y = 0, i = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++, i++) {
      (*TopoMap)[y][x].Dem = Raw->Elevation[i];
      (*TopoMap)[y][x].Mask = Raw->Mask[i];
      (*TopoMap)[y][x].LakeID = 0;
      if (Options->LakeDynamics && INBASIN((*TopoMap)[y][x].Mask))
        (*TopoMap)[y][x].LakeID = Raw->LakeID[i];
    }
  }
  free(Raw->Elevation);
  free(Raw->Mask);
  Raw->Elevation = NULL;
  Raw->Mask = NULL;
  if (Raw->LakeID != NULL) {
    free(Raw->LakeID);
    Raw->LakeID = NULL;
  }

  /* Calculate the area of each lake */
  if (Options->LakeDynamics) {
    Message("\nSetting up lakes\n");
    for (i = 0; i < Map->NumLakes; i++) {
      Area = 0.0;

      for (y = 0; y < Map->NY; y++) {
        for (x = 0; x < Map->NX; x++) {
          if (INBASIN((*TopoMap)[y][x].Mask) && (*TopoMap)[y][x].LakeID == i+1)
            Area += 1.0;
        }
      }

      LType[i].Area = Area * (Map->DX * Map->DY);
      LType[i].Storage = 0.0;
      LType[i].Inflow = 0.0;
      LType[i].Outflow = 0.0;

      Message("%20s: %4.2f km^2\t\tQ ~ %3.2f*L^%3.2f\n",
             LType[i].Name, LType[i].Area/(1000*1000),
             LType[i].PowLawScale, LType[i].PowLawExponent);
    }
    Message("Total: %d lakes\n\n",Map->NumLakes);
  }

  /* Calculate slope, aspect, magnitude of subsurface flow gradient, and
     fraction of flow flowing in each direction based on the land surface
     slope. */
  ElevationSlopeAspect(Map, *TopoMap, Options->MultiFlowDir);
}
