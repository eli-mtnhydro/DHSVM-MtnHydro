/*
 * SUMMARY:      SurfaceEvaporation.c - Calculate evaporation from various surfaces
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Eli Boardman (adapted from SoilEvaporation.c by Bart Nijssen)
 * ORG:          Mountain Hydrology LLC / University of Washington
 * E-MAIL:       eli.boardman@mountainhydrology.com / nijssen@u.washington.edu
 * ORIG-DATE:    Jun-2024 (replaces SoilEvaporation.c orig. Apr-1996)
 * DESCRIPTION:  Evaporation from surface ponding, stream channel, or bare soil
 * DESCRIP-END.
 * FUNCTIONS:    PondEvaporation()
 *               ChannelEvaporation()
 *               SoilEvaporation()
 * COMMENTS:
 * $Id:
 */

#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "DHSVMerror.h"
#include "massenergy.h"
#include "constants.h"
#include "settings.h"
#include "DHSVMChannel.h"

/*****************************************************************************
 PondEvaporation()
*****************************************************************************/
float PondEvaporation(int Dt, float Temp, float Slope, float Gamma,
                      float Lv, float AirDens, float Vpd, float NetRad, float LowerRa,
                      float Evapotranspiration, float *IExcess)
{
  float EPot;			/* Potential evaporation from open water during timestep (m) */
  float PondEvap;	/* Amount of evaporation from surface ponding (m) */
  
  EPot = (Slope * NetRad + AirDens * CP * Vpd / LowerRa) /
    (WATER_DENSITY * Lv * (Slope + Gamma)) * Dt;
  
  /* The potential evaporation rate accounts for the amount of moisture that
     the atmosphere can absorb.  If we do not account for the amount of
     evaporation from overlying evaporation, we can end up with the situation
     that all vegetation layers and the soil layer transpire/evaporate at the
     potential rate, resulting in an overprediction of the actual evaporation
     rate.  Thus we subtract the amount of evaporation that has already
     been calculated for overlying layers from the potential evaporation.
     Another mechanism that could be used to account for this would be to 
     decrease the vapor pressure deficit while going down through the canopy
     (not implemented here) */
  
  /* Prior snow vapor deposition should not increase PET */
  if (Evapotranspiration < 0.0)
    Evapotranspiration = 0.0;
  
  EPot -= Evapotranspiration;
  if (EPot < 0.0)
    EPot = 0.0;
  
  /* Water evaporates from surface ponding up to IExcess or EPot */
  if (EPot > *IExcess)
    PondEvap = *IExcess;
  else
    PondEvap = EPot;
  
  if (PondEvap < 0.0)
    PondEvap = 0.0;
  
  *IExcess -= PondEvap;
  
  return PondEvap;
}

/*****************************************************************************
 ChannelEvaporation()
 *****************************************************************************/
float ChannelEvaporation(int Dt, float DXDY, float Temp, float Slope, float Gamma,
                         float Lv, float AirDens, float Vpd, float NetRad, float LowerRa,
                         float Evapotranspiration, int x, int y, CHANNEL *ChannelData)
{
  float EPot;			/* Potential evaporation rate during timestep (m) */
  float EPotCell;  /* Unsatisfied atmospheric demand during timestep (m) */
  float ChannelEvap;	/* Amount of evaporation directly from the channel (m) */
  
  /* Maximum atmospheric demand if the whole grid cell is covered by water */
  EPot = (Slope * NetRad + AirDens * CP * Vpd / LowerRa) /
    (WATER_DENSITY * Lv * (Slope + Gamma)) * Dt;
  
  /* The potential evaporation rate accounts for the amount of moisture that
   the atmosphere can absorb.  If we do not account for the amount of
   evaporation from overlying evaporation, we can end up with the situation
   that all vegetation layers and the soil layer transpire/evaporate at the
   potential rate, resulting in an overprediction of the actual evaporation
   rate.  Thus we subtract the amount of evaporation that has already
   been calculated for overlying layers from the potential evaporation.
   Another mechanism that could be used to account for this would be to 
   decrease the vapor pressure deficit while going down through the canopy
   (not implemented here) */
  
  /* Prior snow vapor deposition should not increase PET */
  if (Evapotranspiration < 0.0)
    Evapotranspiration = 0.0;
  
  EPotCell = EPot;
  EPotCell -= Evapotranspiration;
  if (EPotCell < 0.0)
    EPotCell = 0.0;
  
  /* Water evaporates from each channel segment at the potential rate
   up to the limit imposed by the total grid cell EPot */
  
  if (EPotCell > 0.0)
    ChannelEvap = channel_grid_evaporation(ChannelData->stream_map, x, y,
                                           EPot, (EPotCell * DXDY)) / DXDY;
  else
    ChannelEvap = 0.0;
  
  return ChannelEvap;
}

/*****************************************************************************
 SoilEvaporation()
 *****************************************************************************/
float SoilEvaporation(int Dt, float Temp, float Slope, float Gamma, float Lv,
                      float AirDens, float Vpd, float NetRad, float RaSoil,
                      float Evapotranspiration, float Porosity, float FCap, float Ks,
                      float Press, float m, float RootDepth,
                      float *MoistContent, float Adjust)
{
  float DesorptionVolume;	/* Amount of water the soil can deliver to the atmosphere during a timestep (m) */
  float EPot;			/* Potential evaporation from soil during timestep (m) */
  float SoilEvap;		/* Amount of evaporation directly from the soil (m) */
  float SoilMoisture;   /* Amount of water in surface soil layer (m) */
  float MoistThrhld;    /* threshold that limits evap to maintain soil at a moisture level */
  float tmp;
  
  DesorptionVolume = Desorption(Dt, *MoistContent, Porosity, Ks, Press, m);
  
  /* Eq.4 Wigmosta et al [1994] */
  
  /* Calculate the density of pure water as a function of temperature.
   Thiesen, Scheel-Diesselhorst Equation (in Handbook of hydrology, fig
   11.1.1) */
  
  EPot = (Slope * NetRad + AirDens * CP * Vpd / RaSoil) /
    (WATER_DENSITY * Lv * (Slope + Gamma)) * Dt;
  
  /* The potential evaporation rate accounts for the amount of moisture that
   the atmosphere can absorb.  If we do not account for the amount of
   evaporation from overlying evaporation, we can end up with the situation
   that all vegetation layers and the soil layer transpire/evaporate at the
   potential rate, resulting in an overprediction of the actual evaporation
   rate.  Thus we subtract the amount of evaporation that has already
   been calculated for overlying layers from the potential evaporation.
   Another mechanism that could be used to account for this would be to 
   decrease the vapor pressure deficit while going down through the canopy
   (not implemented here) */
  
  /* Prior snow vapor deposition should not increase PET */
  if (Evapotranspiration < 0.0)
    Evapotranspiration = 0.0;
  
  EPot -= Evapotranspiration;
  if (EPot < 0.0)
    EPot = 0.0;
  
  /* Eq.8 Wigmosta et al [1994] */
  
  SoilEvap = MIN(EPot, DesorptionVolume);
  SoilEvap *= Adjust;
  SoilMoisture = *MoistContent * RootDepth * Adjust;
  
  MoistThrhld = FCap;
  tmp = MoistThrhld *RootDepth * Adjust;
  if (SoilEvap > SoilMoisture - tmp) {
    SoilEvap = SoilMoisture - tmp;
    *MoistContent = MoistThrhld;
  }
  else {
    SoilMoisture -= SoilEvap;
    *MoistContent = SoilMoisture / (RootDepth * Adjust);
  }
  return SoilEvap;
}
