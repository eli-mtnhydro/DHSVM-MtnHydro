
#ifdef SNOW_ONLY
  #define NO_SOIL
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "massenergy.h"
#include "snow.h"
#include "constants.h"
#include "soilmoisture.h"
#include "photosynthesis.h"
#include "Calendar.h"

 /*****************************************************************************
   Function name: MassEnergyBalance()

   Reference    :
     Epema, G.F. and H.T. Riezbos, 1983, Fall Velocity of waterdrops at different
 heights as a factor influencing erosivity of simulated rain. Rainfall simulation,
 Runoff and Soil Erosion. Catena suppl. 4, Braunschweig. Jan de Ploey (Ed), 1-17.

   Laws, J.o., and D.A. Parsons, 1943, the relation of raindrop size to intensity.
 Trans. Am. Geophys. Union, 24: 452-460.

 *****************************************************************************/
void MassEnergyBalance(OPTIONSTRUCT *Options, int y, int x,
  float DX, float DY, int Dt,
  int HeatFluxOption, int CanopyRadAttOption,
  int InfiltOption, int MaxSoilLayers, int MaxVegLayers, PIXMET *LocalMet,
  NETSTRUCT *LocalNetwork, PRECIPPIX *LocalPrecip, float SnowMeltMultiplier,
  VEGTABLE *LocalVType, VEGPIX *LocalVeg, SOILTABLE *SType,
  SOILPIX *LocalSoil, SNOWPIX *LocalSnow, PIXRAD *LocalRad,
  EVAPPIX *LocalEvap, PIXRAD *TotalRad, CHANNEL *ChannelData,
  float **skyview,
  VEGTABLE *VType, VEGPIX **VegMap, NETSTRUCT **Network, SOILPIX **SoilMap,
  TOPOPIX *LocalTopo, TOPOPIX **TopoMap, MAPSIZE *Map)
{
  float SurfaceWater;		/* Pixel average depth of water before infiltration is calculated (m) */
  float ChannelWater;       /* Precip that hits the channel */
  float Infiltration;		/* Infiltration into the top soil layer (m) */
  float Infiltrability;     /* Dynamic infiltration capacity (m/s)*/
  float B;                  /* Capillary drive and soil saturation deficit used
                               in dynamic infiltration calculation*/
  float LowerRa;		    /* Aerodynamic resistance for lower layer (s/m) */
  float MaxInfiltration;	/* Maximum infiltration into the top soil layer (m) */
  float NetRadiation;		/* Total Net long- and shortwave radiation for each veg layer (W/m2) */
  float PercArea;           /* Surface area of percolation corrected for
                               channel area, divided by the grid cell area (0-1)  */
  float Reference;			/* Reference height for sensible heat calculation (m) */
  float Roughness;			/* Roughness length (m) */
  float Rp;					/* radiation flux in visible part of the spectrum (W/m^2) */
  float UpperRa;		    /* Aerodynamic resistance for upper layer (s/m) */
  float UpperWind;			/* Wind for upper layer (m/s) */
  float SnowLongIn;			/* Incoming longwave radiation at snow surface (W/m2) */
  float SnowNetShort;		/* Net amount of short wave radiation at the snow surface (W/m2) */
  float SnowRa;				/* Aerodynamic resistance for snow */
  float SnowWind;		    /* Wind 2 m above snow */
  float Tsurf;				/* Surface temperature used in LongwaveBalance (C) */
  int   NVegLActual;		/* Number of vegetation layers above snow */
  int i, j;
  float R;                  /* Radius of canopy gap */
  double weight = 0.0;             /* ratio of canopy gap to the grid cell area*/
  float DryChannelEvap;

  /* Used in snow surface energy balance */
  float OldSnowTSurf;       /* Effective surface temperature at the end of the last timestep (C) */
  float SnowTmean;          /* Average snow surface temperature*/
  double Tmp;			    /* Temporary value */
  float Ls;			        /* Latent heat of sublimation (J/kg) */

  /* Edited by Zhuoran Duan zhuoran.duan@pnnl.gov 06/21/2006*/
  /*Add a function to modify soil moisture by add/extract SatFlow from previous time step*/
  DistributeSatflow(Dt, DX, DY, LocalSoil->SatFlow,
    SType->NLayers, LocalSoil->Depth, LocalVType->RootDepth,
    LocalSoil->Porosity, LocalSoil->FCap,
    LocalNetwork->Adjust, &(LocalSoil->TableDepth),
    &(LocalSoil->IExcess), LocalSoil->Moist);
  
  /* Calculate the number of vegetation layers above the snow.
  Note that veg cells with gap must have both over- and under-story as stipulated
  in InitTerrainMap.c, in which gapping is set to FALSE if no overstory regardless
  of canopy gap map value */
  NVegLActual = LocalVType->NVegLayers;
  if (LocalSnow->HasSnow == TRUE && LocalVType->UnderStory == TRUE)
    --NVegLActual;
  if (LocalVeg->Gapping > 0.0) {
    LocalVeg->Type[Opening].NVegLActual = LocalVType->NVegLayers - 1;
    if (LocalVeg->Type[Opening].HasSnow == TRUE && LocalVType->UnderStory == TRUE)
      --LocalVeg->Type[Opening].NVegLActual;

    /* initialize soil moisture */
    for (i = 0; i < CELL_PARTITION; i++) {
      for (j = 0; j <= MaxSoilLayers; j++)
        LocalVeg->Type[i].Moist[j] = LocalSoil->Moist[j];
      LocalVeg->Type[i].MeltEnergy = 0.0;
      LocalVeg->Type[i].MoistureFlux = 0.0;
      LocalVeg->Type[i].ETot = 0.;
    }

    R = LocalVeg->Gapping / 2.0;
    weight = (PI*R*R) / (DX*DY);
    if (weight > 1)
      ReportError("MassEnergyBalance()", 72);
  }

  /* initialize the total amount of evapotranspiration, and MeltEnergy */
  LocalEvap->ETot = 0.0;
  LocalVeg->MeltEnergy = 0.0;
  LocalVeg->MoistureFlux = 0.0;

  /* calculate the radiation balance for pixels */
  RadiationBalance(Options, HeatFluxOption, CanopyRadAttOption,
    LocalVType->OverStory, LocalVType->UnderStory, LocalMet->SineSolarAltitude,
    LocalMet->VICSin, LocalMet->Sin, LocalMet->SinBeam,
    LocalMet->SinDiffuse, LocalMet->Lin, LocalMet->Tair, LocalVeg->Tcanopy,
    LocalSoil->TSurf, SType->Albedo, LocalVType, LocalSnow, LocalRad, LocalVeg);

  /* if a gap is present, calculate radiation balance */
  if (Options->CanopyGapping && (LocalVeg->Gapping > 0.0)) {
    CanopyGapRadiation(&(LocalVeg->Type), LocalMet->SineSolarAltitude, LocalMet->Sin,
      LocalMet->SinBeam, LocalMet->SinDiffuse, LocalMet->Lin, LocalSoil->TSurf,
      LocalVeg->Tcanopy, SType->Albedo, LocalVType, LocalSnow, LocalRad, LocalVeg->Gapping, LocalVeg);

    if (LocalVeg->Type[Forest].HasSnow == TRUE)
      Tsurf = LocalSnow->TSurf;
    else if (HeatFluxOption == TRUE)
      Tsurf = LocalSoil->TSurf;
    else
      Tsurf = LocalMet->Tair;
    GapSurroundingLongRadiation(&(LocalVeg->Type[Forest]), LocalMet->Lin, LocalVeg->Vf,
      LocalVeg->Fract[0], LocalVeg->Type[Forest].Tcanopy, Tsurf);

    GapSurroundingShortRadiation(&(LocalVeg->Type[Forest]), LocalVType, LocalSnow,
      SType->Albedo, LocalMet->SineSolarAltitude, LocalMet->Sin, LocalVeg);
  }

  /* calculate the actual aerodynamic resistances and wind speeds */
  UpperWind = LocalVType->U[0] * LocalMet->Wind;
  UpperRa = LocalVType->Ra[0] / LocalMet->Wind;
  if (LocalVType->OverStory == TRUE)
    LowerRa = LocalVType->Ra[1] / LocalMet->Wind;
  else
    LowerRa = UpperRa;
  if (LocalVeg->Gapping > 0.0) {
    /* calculate the aerodynamic resistance */
    CalcCanopyGapAerodynamic(&(LocalVeg->Type), LocalVType->NVegLayers, LocalVType->Height);
  }

  /* calculate the amount of interception storage, and the amount of
     throughfall. Of course this only needs to be done if there is
     vegetation present. */
#ifndef NO_SNOW

  if (LocalVType->OverStory == TRUE &&
    (LocalPrecip->IntSnow[0] || LocalPrecip->SnowFall > 0.0)) {
    SnowInterception(Options, y, x, Dt, LocalVeg->Fract[0], LocalVeg->Vf,
      LocalVeg->LAI[0], LocalVeg->MaxInt[0], LocalVType->MaxSnowInt, LocalVType->MDRatio,
      LocalVType->SnowIntEff, UpperRa, LocalMet->AirDens,
      LocalMet->Eact, LocalMet->Lv, LocalRad, LocalMet->Press,
      LocalMet->Tair, LocalMet->Vpd, UpperWind,
      &(LocalPrecip->RainFall), &(LocalPrecip->SnowFall),
      &(LocalPrecip->IntRain[0]), &(LocalPrecip->IntSnow[0]),
      &(LocalPrecip->TempIntStorage),
      &(LocalSnow->CanopyVaporMassFlux), &(LocalVeg->Tcanopy),
      &(LocalVeg->MeltEnergy), LocalVeg->Height, LocalVType->UnderStory);
    LocalVeg->MoistureFlux -= LocalSnow->CanopyVaporMassFlux;
    /* Because we now have a new estimate of the canopy temperature we can
       recalculate the longwave balance */
    if (LocalSnow->HasSnow == TRUE)
      Tsurf = LocalSnow->TSurf;
    else if (HeatFluxOption == TRUE)
      Tsurf = LocalSoil->TSurf;
    else
      Tsurf = LocalMet->Tair;
    LongwaveBalance(Options, LocalVType->OverStory, LocalVeg->Fract[0], LocalVeg->Vf,
      LocalMet->Lin, LocalVeg->Tcanopy, Tsurf, LocalRad);
  }
  /* if no snow */
  else if (LocalVType->NVegLayers > 0) {
    LocalVeg->Tcanopy = LocalMet->Tair;
    LocalSnow->CanopyVaporMassFlux = 0.0;
    LocalPrecip->TempIntStorage = 0.0;
    InterceptionStorage(NVegLActual, LocalVeg->MaxInt, LocalVeg->Fract, LocalPrecip->IntRain,
      &(LocalPrecip->RainFall));
  }

  /* if snow is present, simulate the snow pack dynamics */
  if (LocalSnow->HasSnow || LocalPrecip->SnowFall > 0.0) {
    if (LocalVType->OverStory == TRUE) {
      SnowLongIn = LocalRad->LongIn[1];
      SnowNetShort = LocalRad->NetShort[1];
    }
    else {
      SnowLongIn = LocalRad->LongIn[0];
      SnowNetShort = LocalRad->NetShort[0];
    }

    SnowWind = LocalVType->USnow * LocalMet->Wind;
    SnowRa = LocalVType->RaSnow / LocalMet->Wind;

    OldSnowTSurf = LocalSnow->TSurf;
    LocalSnow->Outflow =
      SnowMelt(y, x, Dt, 2. + Z0_SNOW, 0.f, Z0_SNOW, SnowRa, LocalMet->AirDens,
        LocalMet->Eact, LocalMet->Lv, SnowNetShort, SnowLongIn,
        LocalMet->Press, LocalPrecip->RainFall, LocalPrecip->SnowFall,
        LocalMet->Tair, LocalMet->Vpd, SnowWind,
        &(LocalSnow->PackWater), &(LocalSnow->SurfWater),
        &(LocalSnow->Swq), &(LocalSnow->VaporMassFlux),
        &(LocalSnow->TPack), &(LocalSnow->TSurf), &(LocalVeg->MeltEnergy),
        SnowMeltMultiplier);

    /* Calculate the terms of the snow energy balance.  This is similar to the
    code in SnowPackEnergyBalance.c */
    SnowTmean = 0.5 * (OldSnowTSurf + LocalSnow->TSurf);
    /* Apply the stability correction to the aerodynamic resistance */
    if (SnowWind > 0.0)
      SnowRa /= StabilityCorrection(2.0f, 0.f, SnowTmean, LocalMet->Tair,
        SnowWind, Z0_SNOW);
    else
      SnowRa = DHSVM_HUGE;
    /* convert snow surface temperature from C to K */
    Tmp = SnowTmean + 273.15;
    /* net shortwave radiation */
    LocalSnow->Qsw = SnowNetShort;
    /* net longwave radiation */
    LocalSnow->Qlw = SnowLongIn - STEFAN * (Tmp * Tmp * Tmp * Tmp);
    /* sensible heat */
    LocalSnow->Qs = LocalMet->AirDens * CP * (LocalMet->Tair - SnowTmean) / SnowRa;

    /* Calculate latent heat flux */
    if (SnowTmean >= 0.0) {
      /* Melt conditions: use latent heat of vaporization */
      LocalSnow->Qe = LocalMet->Lv * LocalSnow->VaporMassFlux * WATER_DENSITY;
    }
    else {
      /* Accumulation: use latent heat of sublimation (Eq. 3.19, Bras 1990 */
      Ls = (677. - 0.07 * SnowTmean) * JOULESPCAL * GRAMSPKG;
      LocalSnow->Qe = Ls * LocalSnow->VaporMassFlux * WATER_DENSITY;
    }
    LocalSnow->Qe /= Dt;

    /* Calculate advected heat flux from rain */
    LocalSnow->Qp = (CH_WATER * LocalMet->Tair * LocalPrecip->RainFall) / Dt;

    /* Calculate change in cold content */
    LocalSnow->MeltEnergy = LocalVeg->MeltEnergy;
    /* end of snow energy tems */
    /* Rainfall was added to SurfWater of the snow pack and has to be set to zero */
    LocalPrecip->RainFall = 0.0;
    LocalVeg->MoistureFlux -= LocalSnow->VaporMassFlux;

    /* Because we now have a new estimate of the snow surface temperature we
       can recalculate the longwave balance */
    Tsurf = LocalSnow->TSurf;
    LongwaveBalance(Options, LocalVType->OverStory, LocalVeg->Fract[0],
      LocalVeg->Vf, LocalMet->Lin, LocalVeg->Tcanopy, Tsurf, LocalRad);
  }
  else {
    LocalSnow->Outflow = 0.0;
    LocalSnow->VaporMassFlux = 0.0;
    LocalSnow->Qe = 0.;
    LocalSnow->Qs = 0.;
    LocalSnow->Qp = 0.;
    LocalSnow->Qsw = 0.;
    LocalSnow->Qlw = 0.;
    LocalSnow->MeltEnergy = 0.;							   
  }

  /* Determine whether a snow pack is still present, or whether everything
  has melted */
  if (LocalSnow->Swq > 0.0)
    LocalSnow->HasSnow = TRUE;
  else
    LocalSnow->HasSnow = FALSE;

  /************ if canopy gap is present *************/
  if (LocalVeg->Gapping > 0.0) {

    /* calculate intercept rain/snow */
    CanopyGapInterception(Options, &(LocalVeg->Type), HeatFluxOption, y, x,
      Dt, NVegLActual, DX, DY, UpperRa, UpperWind, LocalVType, LocalSoil, LocalVeg,
      LocalSnow, LocalPrecip, LocalRad, LocalMet);

    /* calcuate outflow from snowpack */
    CanopyGapSnowMelt(Options, y, x, Dt, &(LocalVeg->Type), DX, DY, LocalVType,
      LocalVeg, LocalSnow, LocalPrecip, LocalRad, LocalMet);

    /* calcuate snow interception and melt for gap surroudings */
    CalcGapSurroudingIntercept(Options, Options->HeatFlux, y, x, Dt, NVegLActual, 
      &(LocalVeg->Type), LocalVType, LocalRad, LocalMet, UpperRa, UpperWind, LocalVeg);
  }

#endif

#ifndef NO_ET
  /* calculate the amount of evapotranspiration from each vegetation layer
     above the ground/soil surface.  Also calculate the total amount of
     evapotranspiration from the vegetation */
  
  LocalVeg->PhotoAcclim += (LocalMet->Tair - LocalVeg->PhotoAcclim) * Dt / PHOTO_TAU;
  LocalVeg->PhotoDormancy = (LocalVeg->PhotoAcclim - PHOTO_T0) / (PHOTO_T1 - PHOTO_T0);
  LocalVeg->PhotoDormancy = MAX(0.0, MIN(1.0, LocalVeg->PhotoDormancy));
  
  if (LocalVType->OverStory == TRUE) {
    Rp = VISFRACT * LocalRad->NetShort[0];
    if (Options->ImprovRadiation)
      NetRadiation = LocalRad->NetShort[0] +
      LocalRad->LongIn[0] - 2 * LocalVeg->Vf * LocalRad->LongOut[0];
    else
      NetRadiation = LocalRad->NetShort[0] +
      LocalRad->LongIn[0] - 2 * LocalVeg->Fract[0] * LocalRad->LongOut[0];
    LocalRad->NetRadiation[0] = NetRadiation;
    
    EvapoTranspiration(0, Options->ImprovRadiation, Dt, LocalMet, NetRadiation,
      Rp, LocalVType, SType, LocalVeg->MoistureFlux, LocalSoil->Moist, LocalSoil->Temp,
      &(LocalPrecip->IntRain[0]), LocalEvap->EPot, LocalEvap->EInt, LocalEvap->ESoil,
      LocalEvap->EAct, &(LocalEvap->ETot), LocalNetwork->Adjust, UpperRa, LocalVeg,
      Options->PhotoET);
    
    LocalVeg->MoistureFlux += LocalEvap->EAct[0] + LocalEvap->EInt[0];

    if (LocalSnow->HasSnow != TRUE && LocalVType->UnderStory == TRUE) {
      Rp = VISFRACT * LocalRad->NetShort[1];
      NetRadiation =
        LocalRad->NetShort[1] +
        LocalRad->LongIn[1] - LocalVeg->Fract[1] * LocalRad->LongOut[1];
      LocalRad->NetRadiation[1] = NetRadiation;
      
      EvapoTranspiration(1, Options->ImprovRadiation, Dt, LocalMet, NetRadiation,
        Rp, LocalVType, SType, LocalVeg->MoistureFlux, LocalSoil->Moist, LocalSoil->Temp,
        &(LocalPrecip->IntRain[1]), LocalEvap->EPot, LocalEvap->EInt, LocalEvap->ESoil,
        LocalEvap->EAct, &(LocalEvap->ETot), LocalNetwork->Adjust, LowerRa, LocalVeg,
        Options->PhotoET);
      
      LocalVeg->MoistureFlux += LocalEvap->EAct[1] + LocalEvap->EInt[1];
    }
    else if (LocalVType->UnderStory == TRUE) {
      LocalEvap->EAct[1] = 0.;
      LocalEvap->EInt[1] = 0.;
      LocalRad->NetRadiation[1] = 0.;
    }
  }				/* end if(LocalVType->OverStory == TRUE) */
  else if (LocalSnow->HasSnow != TRUE && LocalVType->UnderStory == TRUE) {
    Rp = VISFRACT * LocalRad->NetShort[0];
    NetRadiation =
      LocalRad->NetShort[0] +
      LocalRad->LongIn[0] - LocalVeg->Fract[0] * LocalRad->LongOut[0];
    
    EvapoTranspiration(0, Options->ImprovRadiation, Dt, LocalMet, NetRadiation,
      Rp, LocalVType, SType, LocalVeg->MoistureFlux, LocalSoil->Moist, LocalSoil->Temp,
      &(LocalPrecip->IntRain[0]), LocalEvap->EPot, LocalEvap->EInt, LocalEvap->ESoil,
      LocalEvap->EAct, &(LocalEvap->ETot), LocalNetwork->Adjust, LowerRa, LocalVeg,
      Options->PhotoET);
    
    LocalVeg->MoistureFlux += LocalEvap->EAct[0] + LocalEvap->EInt[0];
    LocalRad->NetRadiation[0] = NetRadiation;
    LocalRad->NetRadiation[1] = 0.;
  }
  else if (LocalVType->UnderStory == TRUE) {
    LocalEvap->EAct[0] = 0.;
    LocalEvap->EInt[0] = 0.;
    LocalRad->NetRadiation[0] = 0;
    LocalRad->NetRadiation[1] = 0.;
  }
  
  /* Calculate open water evaporation from surface ponding if present */
  if (LocalSnow->HasSnow != TRUE && LocalSoil->IExcess > 0.) {
    if (LocalVType->OverStory == TRUE && LocalVType->UnderStory == TRUE)
      NetRadiation = LocalRad->NetShort[1] + LocalRad->LongIn[1] - LocalVeg->Fract[1] * LocalRad->LongOut[1];
    else if (LocalVType->OverStory == TRUE && LocalVType->UnderStory != TRUE)
      NetRadiation = LocalRad->NetShort[1] + LocalRad->LongIn[1] - LocalRad->LongOut[1];
    else if (LocalVType->UnderStory == TRUE)
      NetRadiation = LocalRad->NetShort[0] + LocalRad->LongIn[0] - LocalVeg->Fract[0] * LocalRad->LongOut[0];
    else
      NetRadiation = LocalRad->NetShort[0] + LocalRad->LongIn[0] - LocalRad->LongOut[0];
    
    LocalEvap->EvapSoil =
      PondEvaporation(Dt, (DX*DY), LocalNetwork->Area,
                      LocalMet->Tair, LocalMet->Slope, LocalMet->Gamma,
                      LocalMet->Lv, LocalMet->AirDens, LocalMet->Vpd, NetRadiation, LowerRa,
                      LocalVeg->MoistureFlux, &(LocalSoil->IExcess));
  }
  else
    LocalEvap->EvapSoil = 0.0;
  
  /* Calculate soil evaporation from the upper soil layer if no snow is
     present and there is no understory */
  if (LocalSnow->HasSnow != TRUE && LocalVType->UnderStory != TRUE) {
    if (LocalVType->OverStory == TRUE) {
      NetRadiation =
        LocalRad->NetShort[1] + LocalRad->LongIn[1] - LocalRad->LongOut[1];
      LocalRad->NetRadiation[1] = NetRadiation;
    }
    else {
      NetRadiation =
        LocalRad->NetShort[0] + LocalRad->LongIn[0] - LocalRad->LongOut[0];
      LocalRad->NetRadiation[0] = NetRadiation;
      LocalRad->NetRadiation[1] = 0.;
    }
    /* New: SoilEvaporation is added to the pre-existing value of EvapSoil,
       which is either equal to SurfaceEvaporation or 0. */
    LocalEvap->EvapSoil +=
      SoilEvaporation(Dt,
                      LocalMet->Tair, LocalMet->Slope, LocalMet->Gamma,
                      LocalMet->Lv, LocalMet->AirDens, LocalMet->Vpd, NetRadiation, LowerRa,
                      (LocalVeg->MoistureFlux + LocalEvap->EvapSoil),
                      LocalSoil->Porosity[0], LocalSoil->FCap[0], LocalSoil->KsVert[0],
                      SType->Press[0], SType->PoreDist[0], LocalVType->RootDepth[0],
                      &(LocalSoil->Moist[0]), LocalNetwork->Adjust[0]);
  }
  
  LocalVeg->MoistureFlux += LocalEvap->EvapSoil;
  LocalEvap->ETot += LocalEvap->EvapSoil;
  
  /* Calculate open water and/or dry soil evaporation from stream channels */
  if (channel_grid_has_channel(ChannelData->stream_map, x, y)) {
    if (LocalVType->OverStory == TRUE && LocalVType->UnderStory == TRUE)
      NetRadiation = LocalRad->NetShort[1] + LocalRad->LongIn[1] - LocalVeg->Fract[1] * LocalRad->LongOut[1];
    else if (LocalVType->OverStory == TRUE && LocalVType->UnderStory != TRUE)
      NetRadiation = LocalRad->NetShort[1] + LocalRad->LongIn[1] - LocalRad->LongOut[1];
    else if (LocalVType->UnderStory == TRUE)
      NetRadiation = LocalRad->NetShort[0] + LocalRad->LongIn[0] - LocalVeg->Fract[0] * LocalRad->LongOut[0];
    else
      NetRadiation = LocalRad->NetShort[0] + LocalRad->LongIn[0] - LocalRad->LongOut[0];
    
    /* Calculate soil evaporation from dry channels if not snow-covered */
    if (LocalSnow->HasSnow != TRUE) {
      DryChannelEvap =
        ChannelSoilEvaporation(Dt, (DX*DY),
                               LocalMet->Tair, LocalMet->Slope, LocalMet->Gamma,
                               LocalMet->Lv, LocalMet->AirDens, LocalMet->Vpd, NetRadiation, LowerRa,
                               LocalVeg->MoistureFlux,
                               LocalSoil->Porosity, LocalSoil->FCap, LocalSoil->KsVert,
                               SType->Press, SType->PoreDist,
                               ((SType->NLayers == LocalNetwork->CutBankZone) ?
                                  LocalSoil->Depth - LocalVType->TotalDepth :
                                  LocalVType->RootDepth[LocalNetwork->CutBankZone]),
                               LocalSoil->Moist, LocalNetwork->Adjust,
                               x, y, ChannelData, LocalNetwork->CutBankZone);
      LocalEvap->EvapSoil += DryChannelEvap;
      LocalVeg->MoistureFlux += DryChannelEvap;
      LocalEvap->ETot += DryChannelEvap;
    }
    
    /* Calculate amount of potential evaporation from flowing channels
       Note that actual subtraction and water balance is handled during routing
       Note that it is assumed that flowing channels are always snow-free */
    ChannelEvaporation(Dt, (DX*DY),
                       LocalMet->Tair, LocalMet->Slope, LocalMet->Gamma,
                       LocalMet->Lv, LocalMet->AirDens, LocalMet->Vpd, NetRadiation, LowerRa,
                       LocalVeg->MoistureFlux, x, y, ChannelData);
  }
  else
    LocalEvap->EvapChannel = 0.0;
  
  /* with canopy gaps */
  if (LocalVeg->Gapping > 0.0) {

    CalcGapSurroudingET(Dt, &(LocalVeg->Type), SType, LocalVType, LocalRad, LocalMet,
      LocalSoil, LocalNetwork, UpperRa, LowerRa, LocalVeg,
      DX, DY, x, y, ChannelData, Options->PhotoET);

    /* update wind and aero resistance for gap opening */
    LowerRa = LocalVeg->Type[Opening].Ra[1] / LocalMet->Wind;

    CalcCanopyGapET(&(LocalVeg->Type), MaxSoilLayers, LocalVType, LocalVeg, SType,
      LocalSoil, LocalMet, LocalEvap, LocalNetwork, Dt, UpperRa, LowerRa,
      DX, DY, x, y, ChannelData, Options->PhotoET);

  }
#endif
  
  /* aggregate the gap and non-gap variables based on area weight*/
  if (LocalVeg->Gapping > 0.0)
    AggregateCanopyGap(&(LocalVeg->Type), LocalVeg, LocalSoil, LocalSnow,
		LocalEvap, LocalPrecip, LocalRad, weight, MaxSoilLayers, MaxVegLayers, LocalVType->NVegLayers);

  /* add the water that was not intercepted to the upper soil layer */

#ifndef NO_SOIL

  /* This has been modified so that PercArea for infiltration is calculated
     locally to account for the fact that some cells have streams.
     I am not sure if the old PercArea (which does not account for the fact
     that some cells have streams) needs to remain the same.
     Currently, the old PercArea is passed to UnsaturatedFlow */

  MaxInfiltration = 0.;
  ChannelWater = 0.;
  SurfaceWater = 0.;
  PercArea = 1.;
  
  /* ChannelWater is precipitation falling on the channel */
  if (channel_grid_has_channel(ChannelData->stream_map, x, y)) {
    PercArea = 1. - LocalNetwork->Area / (DX*DY);
    ChannelWater = (LocalNetwork->Area / (DX*DY)) * LocalPrecip->RainFall;
  }

  /* SurfaceWater is rain falling on the hillslope +
     snowmelt on the hillslope (there is no snowmelt on the channel) +
     existing IExcess */
  SurfaceWater = (PercArea * LocalPrecip->RainFall) + LocalSnow->Outflow + LocalSoil->IExcess;
  
  if (InfiltOption == STATIC)
    MaxInfiltration = (1. - LocalVType->ImpervFrac) * PercArea * LocalSoil->MaxInfiltrationRate * Dt;
  else { /* InfiltOption == DYNAMIC
        Dynamic Infiltration Capacity after Parlange and Smith 1978,
        as used in KINEROS and THALES */
    Infiltration = 0.0;
    if (SurfaceWater > 0.) {
      /* Infiltration is a function of the amount of water infiltrated since
     the storm started */
      if (LocalPrecip->PrecipStart) {
        LocalSoil->MoistInit = LocalSoil->Moist[0];
        LocalSoil->InfiltAcc = 0.0;
      }
      /* Check that the B parameter > 0 */
      if ((LocalSoil->InfiltAcc > 0.) && (LocalSoil->Porosity[0] > LocalSoil->MoistInit)) {
        B = (LocalSoil->Porosity[0] - LocalSoil->MoistInit) * (SType->G_Infilt + SurfaceWater);
        Infiltrability = LocalSoil->KsVert[0] * exp((LocalSoil->InfiltAcc) / B) /
          (exp((LocalSoil->InfiltAcc) / B) - 1.);
      }
      else
        Infiltrability = SurfaceWater / Dt;

      MaxInfiltration = Infiltrability * PercArea * (1. - LocalVType->ImpervFrac) * Dt;
      LocalPrecip->PrecipStart = FALSE;
    }/* end  if (SurfaceWater > 0.) */
    else
      LocalPrecip->PrecipStart = TRUE;

  } /* end Dynamic MaxInfiltration calculation */

  Infiltration = (1. - LocalVType->ImpervFrac) * SurfaceWater;

  if (Infiltration > MaxInfiltration)
    Infiltration = MaxInfiltration;
  
  LocalSoil->IExcess = SurfaceWater - Infiltration;

  /*Add water that hits the channel network to the channel network */
  if (ChannelWater > 0.) {
    channel_grid_inc_inflow(ChannelData->stream_map, x, y, ChannelWater * DX * DY);
    LocalSoil->ChannelInt += ChannelWater;
  }

  /* Calculate unsaturated soil water movement, and adjust soil water table depth */
  UnsaturatedFlow(Options, Dt, DX, DY, Infiltration,
    SType->NLayers, LocalSoil->Depth,
    LocalNetwork->Area, LocalVType->RootDepth, LocalSoil->KsVert, SType->KsAnisotropy,
    SType->PoreDist, LocalSoil->Porosity, LocalSoil->FCap, LocalSoil->Perc,
    LocalNetwork->PercArea, LocalNetwork->Adjust, LocalNetwork->CutBankZone,
    LocalNetwork->BankHeight, &(LocalSoil->TableDepth), &(LocalSoil->IExcess),
    LocalSoil->Moist, InfiltOption,
    VType, VegMap, Network, SoilMap,
    LocalTopo->CosSlope, LocalTopo->SinSlope, TopoMap,
    Map, y, x);

  /* Infiltration is updated in UnsaturatedFlow and accumulated
     below */
  if ((InfiltOption == DYNAMIC) && (SurfaceWater > 0.))
    LocalSoil->InfiltAcc += Infiltration;
#endif

  if (HeatFluxOption == TRUE) {
    if (LocalSnow->HasSnow == TRUE) {
      Reference = 2. + Z0_SNOW;
      Roughness = Z0_SNOW;
    }
    else {
      Reference = 2. + Z0_GROUND;
      Roughness = Z0_GROUND;
    }

    SensibleHeatFlux(y, x, Dt, LowerRa, Reference, 0.0f, Roughness,
      LocalMet, LocalRad->PixelNetShort, LocalRad->PixelLongIn,
      LocalVeg->MoistureFlux, SType->NLayers, LocalVType->RootDepth,
      SType, LocalVeg->MeltEnergy, LocalSoil);
    Tsurf = LocalSoil->TSurf;
    LongwaveBalance(Options, LocalVType->OverStory, LocalVeg->Fract[0],
      LocalVeg->Vf, LocalMet->Lin, LocalVeg->Tcanopy, Tsurf, LocalRad);
  }
  else
    NoSensibleHeatFlux(Dt, LocalMet, LocalVeg->MoistureFlux, LocalSoil);


  /* add the components of the radiation balance for the current pixel to
     the total */
  AggregateRadiation(MaxVegLayers, LocalVType->NVegLayers, LocalRad, TotalRad);
  
}
