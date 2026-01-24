
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"
#include "fifobin.h"
#include "fileio.h"
#include "rad.h"
#include "slopeaspect.h"
#include "sizeofnt.h"
#include "varid.h"

 /*****************************************************************************
   InitNewMonth()
   At the start of a new month, read the new radiation files
   (diffuse and direct beam), and potentially a new LAI value.
 *****************************************************************************/
void InitNewMonth(TIMESTRUCT *Time, OPTIONSTRUCT *Options, MAPSIZE *Map,
  TOPOPIX **TopoMap, float **PrismMap, float **SnowPatternMap, float **SnowPatternMapBase, unsigned char ***ShadowMap, 
  INPUTFILES *InFiles, int NVegs, VEGTABLE *VType, int NStats,
  METLOCATION *Stat, char *Path, VEGPIX ***VegMap, SNOWPIX **SnowMap)
{
  const char *Routine = "InitNewMonth";
  char FileName[BUFSIZE * 2 + 5];
  char VarName[BUFSIZE + 1];	/* Variable name */
  int i;
  int j, jj;
  int y, x;
  float a, b, l;
  int NumberType;
  float *Array = NULL;
  unsigned char *Array1 = NULL;
  
  if (DEBUG)
    printf("Initializing new month\n");

  /* If PRISM precipitation fields are being used to interpolate the
     observed precipitation fields, then read in the new months field */

  if (Options->Prism == TRUE) {
    printf("reading in new PRISM field for month %d \n", Time->Current.Month);
    sprintf(FileName, "%s.%02d.%s", Options->PrismDataPath,
      Time->Current.Month, Options->PrismDataExt);
    GetVarName(205, 0, VarName);
    GetVarNumberType(205, &NumberType);
    if (!(Array = (float *)calloc(Map->NY * Map->NX, sizeof(float))))
      ReportError((char *)Routine, 1);
    Read2DMatrix(FileName, Array, NumberType, Map, 0, VarName, 0);

    for (y = 0, i = 0; y < Map->NY; y++)
      for (x = 0; x < Map->NX; x++, i++)
        PrismMap[y][x] = Array[i];

    free(Array);
    
    /* Re-weight snow pattern with fractional amount of current precip pattern */
    if (Options->SnowPattern == TRUE) {
      for (y = 0; y < Map->NY; y++) {
        for (x = 0; x < Map->NX; x++) {
          SnowPatternMap[y][x] = (SnowPatternMapBase[y][x] * SNOWPAT_WEIGHT) +
                                 (PrismMap[y][x] * (1.0 - SNOWPAT_WEIGHT));
        }
      }
      for (i = 0; i < NStats; i++) {
        Stat[i].SnowPattern = (Stat[i].SnowPatternBase * SNOWPAT_WEIGHT) +
                              (Stat[i].PrismPrecip[Time->Current.Month - 1] * (1.0 - SNOWPAT_WEIGHT));
      }
    } /* End snow pattern re-weighting */
    
  }

  if (Options->Shading == TRUE) {
    printf("reading in new shadow map for month %d \n", Time->Current.Month);
    sprintf(FileName, "%s.%02d.%s", Options->ShadingDataPath,
      Time->Current.Month, Options->ShadingDataExt);
    GetVarName(304, 0, VarName);
    GetVarNumberType(304, &NumberType);
    if (!(Array1 = (unsigned char *)calloc(Map->NY * Map->NX, sizeof(unsigned char))))
      ReportError((char *)Routine, 1);
    for (i = 0; i < Time->NDaySteps; i++) {
	  /* if computational time step is finer than hourly, make the shade factor equal within
	  the hourly interval */
	  if (Time->NDaySteps > 24) {
		jj = round(i / (Time->NDaySteps / 24));
		Read2DMatrix(FileName, Array1, NumberType, Map, jj, VarName, jj);
	  }
	  else   
      Read2DMatrix(FileName, Array1, NumberType, Map, i, VarName, i);
      for (y = 0; y < Map->NY; y++) {
        for (x = 0; x < Map->NX; x++) {
          ShadowMap[i][y][x] = Array1[y * Map->NX + x];
        }
      }
    }
    free(Array1);
  }

  printf("changing LAI, albedo and diffuse transmission parameters\n");

  for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
        for (j = 0; j < VType[(*VegMap)[y][x].Veg - 1].NVegLayers; j++) {
          (*VegMap)[y][x].LAI[j] = (*VegMap)[y][x].LAIMonthly[j][Time->Current.Month - 1];
          /*Due to LAI and FC change, have to change MaxInt to spatial as well*/
          (*VegMap)[y][x].MaxInt[j] = (*VegMap)[y][x].LAI[j] * (*VegMap)[y][x].Fract[j] * LAI_WATER_MULTIPLIER;
        }
      }
		}
	}
  for (i = 0; i < NVegs; i++) {
    if (Options->ImprovRadiation) {
      if (VType[i].OverStory == TRUE) {
        VType[i].ExtnCoeff = VType[i].MonthlyExtnCoeff[Time->Current.Month - 1];
      }
      else
        VType[i].ExtnCoeff = 0.;
    }
    
    /* Update vegetation albedo and corresponding understory albedo for snowpack */
    for (j = 0; j < VType[i].NVegLayers; j++) {
      VType[i].Albedo[j] = VType[i].AlbedoMonthly[j][Time->Current.Month - 1];
    }
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
        if (INBASIN(TopoMap[y][x].Mask)) {
          if (VType[(*VegMap)[y][x].Veg - 1].UnderStory == TRUE) {
            if (VType[(*VegMap)[y][x].Veg - 1].OverStory == TRUE) {
              SnowMap[y][x].AlbedoGround = VType[(*VegMap)[y][x].Veg - 1].Albedo[1];
            } else {
              SnowMap[y][x].AlbedoGround = VType[(*VegMap)[y][x].Veg - 1].Albedo[0];
            }
          }
        }
      }
    }
    
	if (Options->CanopyRadAtt == VARIABLE) {
      if (VType[i].OverStory) {
        a = VType[i].LeafAngleA;
        b = VType[i].LeafAngleB;
        l = VType[i].LAI[0] / VType[i].ClumpingFactor;
        if (l == 0)
          VType[i].Taud = 1.0;
        else
          VType[i].Taud = exp(-b * l) * ((1 - a * l) * exp(-a * l) +
            (a * l) * (a * l) * evalexpint(1, a * l));
      }
      else {
        VType[i].Taud = 0.0;
      }
    }
  }
}


/*****************************************************************************
  Function name: InitNewDay()

  Purpose      : Initialize the Earth-Sun geometry variables at the beginning
                 of each day

  Required     :
    int DayOfYear           - day of year (January 1 = 1)
    SOLARGEOMETRY *SolarGeo - structure with information about Earth-Sun
                              geometry

  Comments     : To be excuted at the beginning of each new day
*****************************************************************************/
void InitNewDay(int DayOfYear, SOLARGEOMETRY * SolarGeo)
{
  SolarDay(DayOfYear, SolarGeo->Longitude, SolarGeo->Latitude,
    SolarGeo->StandardMeridian, &(SolarGeo->NoonHour),
    &(SolarGeo->Declination), &(SolarGeo->HalfDayLength),
    &(SolarGeo->Sunrise), &(SolarGeo->Sunset),
    &(SolarGeo->TimeAdjustment), &(SolarGeo->SunEarthDistance));
}


/*****************************************************************************
  Function name: InitNewStep()

  Purpose      : Initialize Earth-Sun geometry and meteorological data at the
                 beginning of each timestep

  Required     :
    MAPSIZE Map              - Structure with information about location
    TIMESTRUCT Time          - Structure with time information
    int FlowGradient         - Type of FlowGradient calculation
    int NStats               - Number of meteorological stations
    METLOCATION *Stat        - Structure with information about the
                               meteorological stations in or near the study
                               area
    SOLARGEOMETRY *SolarGeo  - structure with information about Earth-Sun
                               geometry
    SOILPIX **SoilMap        - structure with soil information

  Comments     : To be executed at the beginning of each time step
*****************************************************************************/
void InitNewStep(INPUTFILES *InFiles, MAPSIZE *Map, TIMESTRUCT *Time,
                 int NSoilLayers, OPTIONSTRUCT *Options, int NStats,
                 METLOCATION *Stat, SOLARGEOMETRY *SolarGeo,
                 TOPOPIX **TopoMap, SOILPIX **SoilMap)
{
  
  /* Calculate variables related to the position of the sun above the
     horizon, this is only necessary if shading is TRUE */

  SolarHour(SolarGeo->Latitude,
            (Time->DayStep + 1) * ((float)Time->Dt) / SECPHOUR,
            ((float)Time->Dt) / SECPHOUR, SolarGeo->NoonHour,
            SolarGeo->Declination, SolarGeo->Sunrise, SolarGeo->Sunset,
            SolarGeo->TimeAdjustment, SolarGeo->SunEarthDistance,
            &(SolarGeo->SineSolarAltitude), &(SolarGeo->DayLight),
            &(SolarGeo->SolarTimeStep), &(SolarGeo->SunMax),
            &(SolarGeo->SolarAzimuth));

  GetMetData(Options, Time, NSoilLayers, NStats, SolarGeo->SunMax, Stat);
}

/*****************************************************************************
   InitNewWaterYear()
   At the start of a new water year, re-initiate the SWE stats maps 
 *****************************************************************************/
void InitNewWaterYear(TIMESTRUCT *Time, OPTIONSTRUCT *Options, MAPSIZE *Map,
                TOPOPIX **TopoMap, SNOWPIX **SnowMap, PRECIPPIX **PrecipMap)
{
  int y, x;
  if (DEBUG)
    printf("Initializing new water year \n");
  
  
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      PrecipMap[y][x].SumPrecip = 0;
      PrecipMap[y][x].SnowAccum = 0;
    }
  }
  
  /* If PRISM precipitation fields are being used to interpolate the
     observed precipitation fields, then read in the new months field */

  if (Options->SnowStats == TRUE) {
    printf("resetting SWE stats map %d \n", Time->Current.Year);
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
        if (INBASIN(TopoMap[y][x].Mask)) {
          SnowMap[y][x].MaxSwe = 0.0;
          SnowMap[y][x].MaxSweDate = 0;
          SnowMap[y][x].MeltOutDate = 0;
        }
      }
    }
  }
}
