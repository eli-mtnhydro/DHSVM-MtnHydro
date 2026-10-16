
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "snow.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"
#include "rad.h"

/*****************************************************************************
Function name: MakeLocalMetData()

Purpose      : Generates meteorological for each individual cell

Required     :
int y 
int x
MAPSIZE Map
int DayStep
int NStats
METLOCATION *Stat
uchar *MetWeights
float LocalElev
RADCLASSPIX *RadMap 
PRECIPPIX *PrecipMap

Comments     :
Reference: Shuttleworth, W.J., Evaporation,  In: Maidment, D. R. (ed.),
Handbook of hydrology,  1993, McGraw-Hill, New York, etc..
*****************************************************************************/
PIXMET MakeLocalMetData(int y, int x, MAPSIZE *Map, int DayStep, int NDaySteps,
                        OPTIONSTRUCT *Options, int NStats,
                        METLOCATION *Stat, uchar *MetWeights,
                        float LocalElev, PIXRAD *RadMap,
                        PRECIPPIX *PrecipMap, float **PrismMap, float **SnowPatternMap,
                        SNOWPIX *LocalSnow,
                        float precipMultiplier, int Month, double ShadeFactor,
                        float SkyView)
{
  float CurrentWeight;		/* weight for current station */
  float Temp;			/* Temporary variable */
  float WeightSum = 0.0;		/* sum of the weights */
  int i;			/* counter */
  float TempLapseRate;
  PIXMET LocalMet;		/* local met data */
  float ContribPrecip, ContribSnow, ContribRain;
  float SnowfallRate;		/* snowfall of the step as a daily rate (m/d) */

  LocalMet.Tair = 0.0;
  LocalMet.Rh = 0.0;
  LocalMet.Wind = 0.0;
  LocalMet.Sin = 0.0;
  LocalMet.SinBeam = 0.0;
  LocalMet.SinDiffuse = 0.0;
  LocalMet.Lin = 0.0;
  TempLapseRate = 0.0;

  WeightSum = 0.0;
  for (i = 0; i < NStats; i++) {
    WeightSum += (float) MetWeights[i];
  }
  for (i = 0; i < NStats; i++) {
    CurrentWeight = ((float) MetWeights[i]) / WeightSum;
    LocalMet.Tair += CurrentWeight *
      LapseT(Stat[i].Data.Tair, Stat[i].Elev, LocalElev,
             Stat[i].Data.TempLapse);
    LocalMet.Rh += CurrentWeight * Stat[i].Data.Rh;
    LocalMet.Wind += CurrentWeight * Stat[i].Data.Wind;
    LocalMet.Lin += CurrentWeight * Stat[i].Data.Lin;
    LocalMet.Sin += CurrentWeight * Stat[i].Data.Sin;
    
    LocalMet.SinBeam += CurrentWeight * Stat[i].Data.SinBeamObs;
    LocalMet.SinDiffuse += CurrentWeight * Stat[i].Data.SinDiffuseObs;
    
    TempLapseRate += CurrentWeight * Stat[i].Data.TempLapse;
  }
  LocalMet.Tair += TEMPERATURE_OFFSET;
  
  /* Optional additional elevation-dependent bias */
  LocalMet.Tair += LAPSE_RATE_BIAS * (LocalElev - LAPSE_BIAS_ELEV);
  
  /* WORK IN PROGRESS, taken from old DHSVM version */
  /* Air pressure */
  /* In rare cases - i.e. when the lapse rate has a different sign for 
   different met stations - you can end up with a TemplapseRate of 0.0
   This will result in a crash, so a check was put in (Jul 28, 1997 - Bart
   Nijssen).  It is somewhat awkward to interpolate lapse rates anyway, so
   a better way of doing this would be welcome */
  if (TempLapseRate != 0.0) {
    Temp = 9.8067 / (TempLapseRate * 287.0);
    LocalMet.Press = 101300. * pow(((288.0 - TempLapseRate * LocalElev) / 288.0), Temp);
  }
  else
    LocalMet.Press = 101300.;

  /* Terrain shading: the beam on the cell (its slope, aspect and horizon)
     relative to the beam on a flat, unshaded surface, the diffuse scaled by
     the sky view factor; where beam and diffuse together would exceed the
     solar constant the beam is reduced */
  if (Options->Shading == TRUE) {
    LocalMet.SinBeam *= ShadeFactor;
    LocalMet.SinDiffuse *= SkyView;

    if (LocalMet.SinBeam + LocalMet.SinDiffuse > SOLARCON)
      LocalMet.SinBeam = SOLARCON - LocalMet.SinDiffuse;
  }

  RadMap->BeamIn = LocalMet.SinBeam;
  RadMap->DiffuseIn = LocalMet.SinDiffuse;
  RadMap->Tair = LocalMet.Tair;

  /* Store the VIC incoming shortwave radiatio without topo or canopy shading */
  LocalMet.VICSin = LocalMet.Sin;

  /* the incoming shortwave radiation adjusted for shading */
  LocalMet.Sin = RadMap->BeamIn + RadMap->DiffuseIn;

  if (Options->Prism == FALSE) {
    PrecipMap->Precip = 0.0;
    PrecipMap->SnowFall = 0.0;
    PrecipMap->RainFall = 0.0;
    for (i = 0; i < NStats; i++) {
      CurrentWeight = ((float) MetWeights[i]) / WeightSum;
      PrecipMap->Precip += CurrentWeight * Stat[i].Data.Precip * precipMultiplier;
    }
  }
  else if (Options->Prism == TRUE) {
    PrecipMap->Precip = 0.0;
    PrecipMap->SnowFall = 0.0;
    PrecipMap->RainFall = 0.0;
    for (i = 0; i < NStats; i++) {
      CurrentWeight = ((float) MetWeights[i]) / WeightSum;
      /* note that X = position from left  boundary, ie # of columns */
      /* note that Y = position from upper boundary, ie # of rows   */
      if (Options->SnowPattern == TRUE) {
        /* Separate weighting of snow and rain */
        ContribPrecip = CurrentWeight * Stat[i].Data.Precip;
        /* First need to partition rain vs. snow */
        if (ContribPrecip > 0.0 && LocalMet.Tair < LocalSnow->Ts) {
          if (LocalMet.Tair > LocalSnow->Tr)
            ContribSnow = ContribPrecip * (LocalSnow->Ts - LocalMet.Tair) / (LocalSnow->Ts - LocalSnow->Tr);
          else
            ContribSnow = ContribPrecip;
        } else {
          ContribSnow = 0.0;
        }
        ContribRain = ContribPrecip - ContribSnow;
        
        /* Apply snow pattern weighting */
        ContribSnow *= (SnowPatternMap[y][x] / Stat[i].SnowPattern);
        PrecipMap->SnowFall += ContribSnow;
        
        /* Apply rain pattern weighting (here PRISM == rain pattern) */
        ContribRain *= (PrismMap[y][x] / Stat[i].PrismPrecip[Month - 1]);
        PrecipMap->RainFall += ContribRain;
        
        PrecipMap->Precip += (ContribSnow + ContribRain);
      } /* End of snow pattern weighting */
        else
          PrecipMap->Precip += (CurrentWeight * Stat[i].Data.Precip) * (PrismMap[y][x] / Stat[i].PrismPrecip[Month - 1]);
    }
    PrecipMap->Precip *= precipMultiplier;
    PrecipMap->SnowFall *= precipMultiplier;
    PrecipMap->RainFall *= precipMultiplier;
  }

  /* Separate precipitation into rainfall and snowfall by air temperature
  (with a snow pattern the split was made per station above) */
  if (Options->SnowPattern == FALSE) {
    if (PrecipMap->Precip > 0.0 && LocalMet.Tair < LocalSnow->Ts) {
      if (LocalMet.Tair > LocalSnow->Tr)
        PrecipMap->SnowFall = PrecipMap->Precip *
        (LocalSnow->Ts - LocalMet.Tair) / (LocalSnow->Ts - LocalSnow->Tr);
      else
        PrecipMap->SnowFall = PrecipMap->Precip;
    }
    else
      PrecipMap->SnowFall = 0.0;
    PrecipMap->RainFall = PrecipMap->Precip - PrecipMap->SnowFall;
  }

  /* Local heat of vaporization, Eq. 4.2.1, Shuttleworth (1993) */
  LocalMet.Lv = 2501000 - 2361 * LocalMet.Tair;

  /* Psychrometric constant */
  LocalMet.Gamma = CP * LocalMet.Press / (EPS * LocalMet.Lv);

  /* Saturated vapor pressure, Eq. 4.2.2, Shuttleworth (1993) */
  LocalMet.Es = SatVaporPressure(LocalMet.Tair);

  /* Slope of vapor pressure curve, Eq. 4.2.3, Shuttleworth (1993) */
  LocalMet.Slope = 4098.0 * LocalMet.Es /
    ((237.3 + LocalMet.Tair) * (237.3 + LocalMet.Tair));

  /* Actual vapor pressure */
  LocalMet.Eact = LocalMet.Es * (LocalMet.Rh / 100.);

  /* Vapor pressure deficit */
  LocalMet.Vpd = LocalMet.Es - LocalMet.Eact;

  /* Air density, Eq. 4.2.4 Shuttleworth (1993) */
  LocalMet.AirDens = 0.003486 * LocalMet.Press / (275 + LocalMet.Tair);

  /* Snow albedo as a function of days since last snow */
  if (LocalSnow->HasSnow) {
    if (PrecipMap->SnowFall > 0.0 && LocalSnow->TSurf < 0.0)
      LocalSnow->AccumSeason = TRUE;
    else if (fequal(LocalSnow->TSurf, 0.0))
      LocalSnow->AccumSeason = FALSE;
    
    /* The step's snowfall as a daily rate r: above the reset rate R it resets
       the days since snowfall; 0 < r <= R scales them by 1 - r/R per day,
       (1 - r/R)^(1/NDaySteps) per step */
    SnowfallRate = PrecipMap->SnowFall * (float) NDaySteps;
    if (SnowfallRate > MIN_SNOW_RESET_ALBEDO)
      LocalSnow->LastSnow = 0.0;
    else if (SnowfallRate > 0.0)
      LocalSnow->LastSnow *= pow(1.0 - SnowfallRate / MIN_SNOW_RESET_ALBEDO, 1.0 / NDaySteps);
    else
      LocalSnow->LastSnow += (1.0 / (float) NDaySteps);
    
    LocalSnow->Albedo = CalcSnowAlbedo(LocalSnow, PrecipMap->SnowFall);
  }
  /* Snowfall on a snow-free cell starts a new pack at age 0; the ground below
     it is blended in by depth (CalcSnowAlbedo) */
  else if (PrecipMap->SnowFall > 0.0)
    LocalSnow->LastSnow = 0.0;
  else
    LocalSnow->LastSnow += (1.0 / (float) NDaySteps);

  return LocalMet;
}
