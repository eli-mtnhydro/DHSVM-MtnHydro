
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "constants.h"
#include "data.h"
#include "DHSVMerror.h"
#include "fileio.h"
#include "functions.h"
#include "rad.h"
#include "sizeofnt.h"
#include "varid.h"

 /*****************************************************************************
   InitMetMaps()
 *****************************************************************************/
void InitMetMaps(LISTPTR Input, int NDaySteps, MAPSIZE *Map, 
  OPTIONSTRUCT *Options,
  float ***PrismMap, float ***SnowPatternMap, float ***SnowPatternMapBase,
  unsigned char ****ShadowMap, float ***SkyViewMap,
  EVAPPIX ***EvapMap, PRECIPPIX ***PrecipMap, float ***PptMultiplierMap,
  float ***MeltMultiplierMap, PIXRAD ***RadMap,
  SOILPIX **SoilMap, LAYER *Soil, VEGPIX **VegMap,
  LAYER *Veg, TOPOPIX **TopoMap)
{
  int y, x;

  printf("Initializing meteorological maps\n");

  InitEvapMap(Map, EvapMap, SoilMap, Soil, VegMap, Veg, TopoMap);
  InitPrecipMap(Map, PrecipMap, VegMap, Veg, TopoMap);
  InitMultiplierMaps(Options, Map, PptMultiplierMap, MeltMultiplierMap);                                                            

  if (Options->Prism == TRUE)
    InitPrismMap(Map->NY, Map->NX, PrismMap);
  if (Options->SnowPattern == TRUE)
    InitSnowPatternMap(SnowPatternMap, SnowPatternMapBase, Map, Options);
  if (Options->Shading == TRUE)
    InitShadeMap(Options, NDaySteps, Map, ShadowMap, SkyViewMap);
  
  if (!((*SkyViewMap) = (float **)calloc(Map->NY, sizeof(float *))))
    ReportError("InitMetMaps()", 1);
  for (y = 0; y < Map->NY; y++) {
    if (!((*SkyViewMap)[y] = (float *)calloc(Map->NX, sizeof(float))))
      ReportError("InitMetMaps()", 1);
  }
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      (*SkyViewMap)[y][x] = 1.0;
    }
  }
  
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
    printf("Initializing evaporation map\n");

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
    printf("Initializing precipitation map\n");

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
    printf("Initializing radiation map\n");

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

void InitSnowPatternMap(float ***SnowPatternMap, float ***SnowPatternMapBase,
                        MAPSIZE *Map, OPTIONSTRUCT *Options)
{
  const char *Routine = "InitSnowPatternMap";
  int x, y, i;
  char FileName[BUFSIZE + 1];
  char VarName[BUFSIZE + 1];
  int NumberType;
  float *Array = NULL;
  
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
  
  printf("\nReading in snow pattern map\n");
  sprintf(FileName, "%s", Options->SnowPatternDataPath);
  GetVarName(207, 0, VarName);
  GetVarNumberType(207, &NumberType);
  if (!(Array = (float *)calloc(Map->NY * Map->NX, sizeof(float))))
    ReportError((char *)Routine, 1);
  Read2DMatrix(FileName, Array, NumberType, Map, 0, VarName, 0);
  
  for (y = 0, i = 0; y < Map->NY; y++)
    for (x = 0; x < Map->NX; x++, i++)
      (*SnowPatternMapBase)[y][x] = Array[i];
  
  free(Array);
}

/******************************************************************************/
/*				  InitShadeMap                                */
/******************************************************************************/
void InitShadeMap(OPTIONSTRUCT * Options, int NDaySteps, MAPSIZE *Map,
  unsigned char ****ShadowMap, float ***SkyViewMap)
{
  const char *Routine = "InitShadeMap";
  char VarName[BUFSIZE + 1];	/* Variable name */
  int x;			/* counter */
  int y;			/* counter */
  int n;
  int NumberType;
  float *Array = NULL;

  if (!((*ShadowMap) =
    (unsigned char ***)calloc(NDaySteps, sizeof(unsigned char **))))
    ReportError((char *)Routine, 1);
  for (n = 0; n < NDaySteps; n++) {
    if (!((*ShadowMap)[n] =
      (unsigned char **)calloc(Map->NY, sizeof(unsigned char *))))
      ReportError((char *)Routine, 1);
    for (y = 0; y < Map->NY; y++) {
      if (!((*ShadowMap)[n][y] =
        (unsigned char *)calloc(Map->NX, sizeof(unsigned char))))
        ReportError((char *)Routine, 1);
    }
  }

  if (!((*SkyViewMap) = (float **)calloc(Map->NY, sizeof(float *))))
    ReportError((char *)Routine, 1);
  for (y = 0; y < Map->NY; y++) {
    if (!((*SkyViewMap)[y] = (float *)calloc(Map->NX, sizeof(float))))
      ReportError((char *)Routine, 1);
  }

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      (*SkyViewMap)[y][x] = 1.0;
    }
  }

  GetVarName(305, 0, VarName);
  GetVarNumberType(305, &NumberType);
  if (!(Array = (float *)calloc(Map->NY * Map->NX, sizeof(float))))
    ReportError((char *)Routine, 1);
  Read2DMatrix(Options->SkyViewDataPath, Array, NumberType, Map, 0,
    VarName, 0);
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      (*SkyViewMap)[y][x] = Array[y * Map->NX + x];
    }
  }

  free(Array);
}

/*****************************************************************************
 InitMultiplierMaps()
*****************************************************************************/
void InitMultiplierMaps(OPTIONSTRUCT * Options, MAPSIZE *Map,
                        float ***PptMultiplierMap, float ***MeltMultiplierMap)
{
  const char *Routine = "InitMultiplierMaps";
  char VarName[BUFSIZE + 1];
  int i;			/* counter */
  int x;			/* counter */
  int y;			/* counter */
  int NumberType;		/* number type */
  float *Array;
  
  /* Precip multiplier */
  
  if (!((*PptMultiplierMap) = (float **)calloc(Map->NY, sizeof(float *))))
    ReportError((char *)Routine, 1);
  for (y = 0; y < Map->NY; y++) {
    if (!((*PptMultiplierMap)[y] = (float *)calloc(Map->NX, sizeof(float))))
      ReportError((char *)Routine, 1);
  }
  if (PRECIP_MULTIPLIER > NA) {
    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*PptMultiplierMap)[y][x] = PRECIP_MULTIPLIER;
      }
    }
  }
  else if (!IsEmptyStr(Options->PrecipMultiplierMapPath)) {
    /* Read the map path */
    GetVarName(100, 0, VarName);
    GetVarNumberType(100, &NumberType);
    if (!(Array = (float *)calloc(Map->NX * Map->NY, SizeOfNumberType(NumberType))))
      ReportError((char *)Routine, 1);
    Read2DMatrix(Options->PrecipMultiplierMapPath, Array, NumberType, Map, 0, VarName, 0);

    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*PptMultiplierMap)[y][x] = Array[i];
      }
    }
    free(Array);
  }
  else {
    printf("No valid input of precipitation multiplier ...\n");
    exit(88);
  }

  /* Snow melt multiplier */
  
  if (!((*MeltMultiplierMap) = (float **)calloc(Map->NY, sizeof(float *))))
    ReportError((char *)Routine, 1);
  for (y = 0; y < Map->NY; y++) {
    if (!((*MeltMultiplierMap)[y] = (float *)calloc(Map->NX, sizeof(float))))
      ReportError((char *)Routine, 1);
  }
  if (SNOWMELT_MULTIPLIER > NA) {
    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*MeltMultiplierMap)[y][x] = SNOWMELT_MULTIPLIER;
      }
    }
  }
  else if (!IsEmptyStr(Options->SnowMeltMultiplierMapPath)) {
    /* Read the map path */
    GetVarName(101, 0, VarName);
    GetVarNumberType(101, &NumberType);
    if (!(Array = (float *)calloc(Map->NX * Map->NY, SizeOfNumberType(NumberType))))
      ReportError((char *)Routine, 1);
    Read2DMatrix(Options->SnowMeltMultiplierMapPath, Array, NumberType, Map, 0, VarName, 0);

    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*MeltMultiplierMap)[y][x] = Array[i];
      }
    }
    free(Array);
  }
  else {
    printf("No valid input of snow melt multiplier ...\n");
    exit(88);
  }

}                                                               