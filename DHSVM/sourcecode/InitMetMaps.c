
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "constants.h"
#include "data.h"
#include "DHSVMerror.h"
#include "fileio.h"
#include "functions.h"
#include "getinit.h"
#include "rad.h"
#include "mapio.h"

 /*****************************************************************************
   InitMetMaps()
 *****************************************************************************/
void InitMetMaps(MAPSIZE *Map, OPTIONSTRUCT *Options, INPUTMAPS *Raw,
  float ***PrismMap, float ***SnowPatternMap, float ***SnowPatternMapBase,
  EVAPPIX ***EvapMap, PRECIPPIX ***PrecipMap, PIXRAD ***RadMap,
  SOILPIX **SoilMap, LAYER *Soil, VEGPIX **VegMap,
  LAYER *Veg, TOPOPIX **TopoMap)
{
  Message("Initializing meteorological maps\n");

  InitEvapMap(Map, EvapMap, SoilMap, Soil, VegMap, Veg, TopoMap);
  InitPrecipMap(Map, PrecipMap, VegMap, Veg, TopoMap);

  /* Monthly precipitation pattern: the map is open (ReadInputMaps()), the
     month's layer is read by InitNewMonth() */
  if (Options->Prism == TRUE)
    InitPrismMap(Map->NY, Map->NX, PrismMap);

  if (Options->SnowPattern == TRUE)
    InitSnowPatternMap(Raw, SnowPatternMap, SnowPatternMapBase, Map);

  InitRadMap(Map, RadMap);
}


/*****************************************************************************
  InitEvapMap()
*****************************************************************************/
void InitEvapMap(MAPSIZE *Map, EVAPPIX ***EvapMap, SOILPIX **SoilMap,
  LAYER *Soil, VEGPIX **VegMap, LAYER *Veg,
  TOPOPIX **TopoMap)
{
  const char *Routine = "InitEvapMap";
  int i;			/* counter */
  int x;			/* counter */
  int y;			/* counter */
  int NSoil;			/* Number of soil layers for current pixel */
  int NVeg;			/* Number of veg layers for current pixel */

  if (DEBUG)
    Message("Initializing evaporation map\n");

  if (!(*EvapMap = (EVAPPIX **)calloc(Map->NY, sizeof(EVAPPIX *))))
    ReportError((char *)Routine, 1);

  for (y = 0; y < Map->NY; y++) {
    if (!((*EvapMap)[y] = (EVAPPIX *)calloc(Map->NX, sizeof(EVAPPIX))))
      ReportError((char *)Routine, 1);
  }

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
        NVeg = Veg->NLayers[(VegMap[y][x].Veg - 1)];
        NSoil = Soil->NLayers[(SoilMap[y][x].Soil - 1)];
        
        if (!((*EvapMap)[y][x].EPot =
          (float *)calloc(NVeg + 1, sizeof(float))))
          ReportError((char *)Routine, 1);

        if (!((*EvapMap)[y][x].EAct =
          (float *)calloc(NVeg + 1, sizeof(float))))
          ReportError((char *)Routine, 1);

        if (!((*EvapMap)[y][x].EInt = (float *)calloc(NVeg, sizeof(float))))
          ReportError((char *)Routine, 1);

        if (!((*EvapMap)[y][x].ESoil =
          (float **)calloc(NVeg, sizeof(float *))))
          ReportError((char *)Routine, 1);

        for (i = 0; i < NVeg; i++) {
          if (!((*EvapMap)[y][x].ESoil[i] =
            (float *)calloc(NSoil, sizeof(float))))
            ReportError((char *)Routine, 1);
        }
      }
    }
  }
}

/*****************************************************************************
  InitPrecipMap()
*****************************************************************************/
void InitPrecipMap(MAPSIZE * Map, PRECIPPIX *** PrecipMap, VEGPIX ** VegMap,
  LAYER * Veg, TOPOPIX ** TopoMap)
{
  const char *Routine = "InitPrecipMap";
  int x;			/* counter */
  int y;			/* counter */
  int NVeg;			/* Number of veg layers at current pixel */

  if (DEBUG)
    Message("Initializing precipitation map\n");

  if (!(*PrecipMap = (PRECIPPIX **)calloc(Map->NY, sizeof(PRECIPPIX *))))
    ReportError((char *)Routine, 1);

  for (y = 0; y < Map->NY; y++) {
    if (!((*PrecipMap)[y] = (PRECIPPIX *)calloc(Map->NX, sizeof(PRECIPPIX))))
      ReportError((char *)Routine, 1);
  }

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
        NVeg = Veg->NLayers[(VegMap[y][x].Veg - 1)];
        if (!((*PrecipMap)[y][x].IntRain =
          (float *)calloc(NVeg, sizeof(float))))
          ReportError((char *)Routine, 1);
      }
    }
  }

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
        NVeg = Veg->NLayers[(VegMap[y][x].Veg - 1)];
        if (!((*PrecipMap)[y][x].IntSnow =
          (float *)calloc(NVeg, sizeof(float))))
          ReportError((char *)Routine, 1);
      }
    }
  }

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      (*PrecipMap)[y][x].SumPrecip = 0.0;
      (*PrecipMap)[y][x].SnowAccum = 0.0;
      (*PrecipMap)[y][x].SnowMelt = 0.0;
      if (INBASIN(TopoMap[y][x].Mask))
        (*PrecipMap)[y][x].PrecipStart = TRUE;
    }
  }
}

/******************************************************************************
  InitRadMap()
******************************************************************************/
void InitRadMap(MAPSIZE *Map, PIXRAD ***RadMap)
{
  const char *Routine = "InitRadMap";
  int y;			/* counter */

  if (DEBUG)
    Message("Initializing radiation map\n");

  if (!(*RadMap = (PIXRAD **)calloc(Map->NY, sizeof(PIXRAD *))))
    ReportError((char *)Routine, 1);

  for (y = 0; y < Map->NY; y++) {
    if (!((*RadMap)[y] = (PIXRAD *)calloc(Map->NX, sizeof(PIXRAD))))
      ReportError((char *)Routine, 1);
  }
}

/******************************************************************************/
/*			       InitPRISMMap                             */
/******************************************************************************/

void InitPrismMap(int NY, int NX, float ***PrismMap)
{
  const char *Routine = "InitPRISMMap";
  int x;			/* counter */
  int y;			/* counter */

  if (!((*PrismMap) = (float **)calloc(NY, sizeof(float *))))
    ReportError((char *)Routine, 1);

  for (y = 0; y < NY; y++) {
    if (!((*PrismMap)[y] = (float *)calloc(NX, sizeof(float))))
      ReportError((char *)Routine, 1);
  }

  for (y = 0; y < NY; y++) {
    for (x = 0; x < NX; x++) {
      (*PrismMap)[y][x] = 1.0;
    }
  }

}

/******************************************************************************/
/*			       InitSnowPatternMap                             */
/******************************************************************************/

void InitSnowPatternMap(INPUTMAPS *Raw, float ***SnowPatternMap, float ***SnowPatternMapBase,
                        MAPSIZE *Map)
{
  const char *Routine = "InitSnowPatternMap";
  int x, y, i;

  if (!((*SnowPatternMap) = (float **)calloc(Map->NY, sizeof(float *))))
    ReportError((char *)Routine, 1);
  if (!((*SnowPatternMapBase) = (float **)calloc(Map->NY, sizeof(float *))))
    ReportError((char *)Routine, 1);
  
  for (y = 0; y < Map->NY; y++) {
    if (!((*SnowPatternMap)[y] = (float *)calloc(Map->NX, sizeof(float))))
      ReportError((char *)Routine, 1);
  }
  for (y = 0; y < Map->NY; y++) {
    if (!((*SnowPatternMapBase)[y] = (float *)calloc(Map->NX, sizeof(float))))
      ReportError((char *)Routine, 1);
  }
  
  /* The snow pattern map (read by ReadInputMaps()) */
  for (y = 0, i = 0; y < Map->NY; y++)
    for (x = 0; x < Map->NX; x++, i++)
      (*SnowPatternMapBase)[y][x] = Raw->SnowPattern[i];

  free(Raw->SnowPattern);
  Raw->SnowPattern = NULL;
}
