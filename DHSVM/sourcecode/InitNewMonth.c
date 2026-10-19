
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"
#include "fileio.h"
#include "rad.h"
#include "slopeaspect.h"
#include "sizeofnt.h"
#include "mapio.h"

 /*****************************************************************************
   InitNewMonth()
   At the start of a new month, read the new radiation files
   (diffuse and direct beam), and potentially a new LAI value.
 *****************************************************************************/
void InitNewMonth(TIMESTRUCT *Time, OPTIONSTRUCT *Options, MAPSIZE *Map,
  TOPOPIX **TopoMap, float **PrismMap, float **SnowPatternMap, float **SnowPatternMapBase,
  INPUTFILES *InFiles, VEGTABLES *Veg, int NStats,
  METLOCATION *Stat, char *Path, VEGPIX ***VegMap, SNOWPIX **SnowMap,
  SOILTABLE *SType, SOILPIX **SoilMap)
{
  const char *Routine = "InitNewMonth";
  int i;
  int y, x;
  float *Array = NULL;
  OVERTABLE *Over;
  UNDERTABLE *Under;

  if (DEBUG)
    Message("Initializing new month\n");

  /* If PRISM precipitation fields are being used to interpolate the
     observed precipitation fields, then read in the new months field */

  if (Options->Prism == TRUE) {
    Message("reading in new PRISM field for month %d \n", Time->Current.Month);
    if (!(Array = (float *)calloc(Map->NY * Map->NX, sizeof(float))))
      ReportError((char *)Routine, 1);
    ReadMapLayer(&(Options->PrecipPattern), Time->Current.Month - 1, Array, Map);

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

  Message("changing LAI, albedo and diffuse transmission parameters\n");

  for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask))
        SetCellMonth(Options, Veg, &((*VegMap)[y][x]), Time->Current.Month);
		}
	}

  /* Albedo of the month of every overstory and understory class */
  for (i = 0; i < Veg->NOver; i++) {
    Over = &(Veg->Over[i]);
    Over->Story.Albedo = Over->Story.AlbedoMonthly[Time->Current.Month - 1];
  }
  for (i = 0; i < Veg->NUnder; i++) {
    Under = &(Veg->Under[i]);
    Under->Story.Albedo = Under->Story.AlbedoMonthly[Time->Current.Month - 1];
  }

  /* The ground albedo under the snowpack */
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask))
        SnowMap[y][x].AlbedoGround = GroundAlbedo(Veg, &((*VegMap)[y][x]),
          SType[SoilMap[y][x].Soil - 1].Albedo, Time->Current.Month);
    }
  }
}


/*****************************************************************************
  SetCellMonth()

  The LAI of the month of each vegetation layer of a mask cell, its
  interception capacity, the shortwave transmittance of the understory from
  its LAI, and the diffuse transmittance of the overstory from its LAI
  (hemispherical integration of the beam transmittance, 2 E3(k LAI))
*****************************************************************************/
void SetCellMonth(OPTIONSTRUCT *Options, VEGTABLES *Veg, VEGPIX *Cell, int Month)
{
  OVERTABLE *Over = OverstoryRow(Veg, Cell);
  int j;
  float l;

  for (j = 0; j < Cell->NVegLayers; j++) {
    Cell->LAI[j] = Cell->LAIMonthly[j][Month - 1];
    /*Due to LAI and FC change, have to change MaxInt to spatial as well*/
    Cell->MaxInt[j] = Cell->LAI[j] * Cell->Fract[j] * LAI_WATER_MULTIPLIER;
  }

  if (Cell->UnderStory == TRUE)
    Cell->TauU = exp(-UnderstoryRow(Veg, Cell)->RadAtten * Cell->LAI[CanopyIndex(Cell, 1)]);
  else
    Cell->TauU = 1.0;

  if (Cell->OverStory) {
    l = Over->ExtnCoeff * Cell->LAI[0];
    if (l == 0)
      Cell->Taud = 1.0;
    else
      Cell->Taud = (1 - l) * exp(-l) + l * l * evalexpint(1, l);
  }
  else {
    Cell->Taud = 0.0;
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
                 OPTIONSTRUCT *Options, int NStats,
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

  GetMetData(Options, Time, NStats, SolarGeo->SunMax, Stat);
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
    Message("Initializing new water year \n");
  
  
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      PrecipMap[y][x].SumPrecip = 0;
      PrecipMap[y][x].SnowAccum = 0;
      PrecipMap[y][x].SnowMelt = 0;
    }
  }
  
  /* Reset the snow statistics of the water year */
  Message("resetting SWE stats map %d \n", Time->Current.Year);
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
