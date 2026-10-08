
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"

/*****************************************************************************
  Aggregate()
  
  Calculate the average values for the different fluxes and state variables
  over the basin.  
  In the current implementation the local radiation
  elements are not stored for the entire area.  Therefore these components
  are aggregated in AggregateRadiation() inside MassEnergyBalance().
  
  The aggregated values are set to zero in the function RestAggregate,
  which is executed at the beginning of each time step.
*****************************************************************************/
void Aggregate(MAPSIZE *Map, OPTIONSTRUCT *Options, TOPOPIX **TopoMap,
	       LAYER *Soil, LAYER *Veg, VEGPIX **VegMap, EVAPPIX **Evap,
	       PRECIPPIX **Precip, PIXRAD **RadMap, SNOWPIX **Snow,
	       SOILPIX **SoilMap, AGGREGATED *Total, VEGTABLE *VType,
	       NETSTRUCT **Network, CHANNEL *ChannelData, 
         int Dt, int NDaySteps)
{
  int NPixels;		/* Number of pixels in the basin */
  int NSoilL;			/* Number of soil layers for current pixel */
  int NVegL;			/* Number of vegetation layers for current pixel */
  int NSnow;      /* Number of pixels with snow cover currently */
  int i;				/* counter */
  int j;				/* counter */
  int x;
  int y;
  float DeepDepth;		/* depth to bottom of lowest rooting zone */

  NPixels = 0;
  NSnow = 0;
  
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
  NPixels++;
  NSoilL = Soil->NLayers[SoilMap[y][x].Soil - 1];
  NVegL = Veg->NLayers[VegMap[y][x].Veg - 1];
  
  /* aggregate the evaporation data */
  Total->Evap.ETot += Evap[y][x].ETot;
  for (i = 0; i < NVegL; i++) {
	  Total->Evap.EPot[i] += Evap[y][x].EPot[i];
	  Total->Evap.EAct[i] += Evap[y][x].EAct[i];
	  Total->Evap.EInt[i] += Evap[y][x].EInt[i];
  }
  Total->Evap.EPot[Veg->MaxLayers] += Evap[y][x].EPot[NVegL];
  Total->Evap.EAct[Veg->MaxLayers] += Evap[y][x].EAct[NVegL];
  
  for (i = 0; i < NVegL; i++) {
	  for (j = 0; j < NSoilL; j++) {
		  Total->Evap.ESoil[i][j] += Evap[y][x].ESoil[i][j];
	  }
  }
  Total->Evap.EvapSoil += Evap[y][x].EvapSoil;
  Total->Evap.EvapChannel += Evap[y][x].EvapChannel;
  
  /* Sperry/Medlyn canopy diagnostics.  Rc is accumulated as CONDUCTANCE:
   resistances combine in parallel across a basin, so an arithmetic mean of
   Rc is dominated by the RsMax pixels and overstates basin resistance. */
  for (i = 0; i < NVegL; i++)
    Total->Veg.Rc[i] += (VegMap[y][x].Rc[i] > 0.0) ? 1.0 / VegMap[y][x].Rc[i] : 0.0;
  Total->Veg.PhotoDormancy += VegMap[y][x].PhotoDormancy;

  /* Plant hydraulic diagnostics, per canopy layer.

     Water potentials are INTENSIVE, so unlike a flux they cannot simply be
     summed over the basin and divided by NPixels: a pixel with no vegetation
     in a given layer contributes a potential of 0 MPa, which is the wettest
     possible value and would bias the basin mean toward "unstressed".  So
     each layer carries its own count of contributing pixels, incremented
     only where that layer actually has leaf area, and the divide below uses
     it.  NVegPix is a diagnostic counter, not model state. */
  for (i = 0; i < NVegL; i++) {
    Total->Veg.AnCanopy[i] += VegMap[y][x].AnCanopy[i];
    Total->Veg.Escheme[i]  += VegMap[y][x].Escheme[i];
    Total->Veg.Ecrit[i]    += VegMap[y][x].Ecrit[i];
    Total->Veg.SupplyLimited[i] += VegMap[y][x].SupplyLimited[i];

    /* PsiSoil is computed for every scheme, so it is averaged over the
       vegetated fraction whether or not the hydraulic scheme ran. */
    if (VegMap[y][x].LAI != NULL && VegMap[y][x].LAI[i] > 0.0)
      Total->Veg.PsiSoil[i] += VegMap[y][x].PsiSoil[i];

    if (VegMap[y][x].LAI != NULL && VegMap[y][x].LAI[i] > 0.0) {
      Total->NVegPix[i]++;
      Total->Veg.Tleaf[i]        += VegMap[y][x].Tleaf[i];
      Total->Veg.VpdLeaf[i]      += VegMap[y][x].VpdLeaf[i];
      Total->Veg.PsiRoot[i]      += VegMap[y][x].PsiRoot[i];
      Total->Veg.PsiLeaf[i]      += VegMap[y][x].PsiLeaf[i];
      Total->Veg.PLC[i]          += VegMap[y][x].PLC[i];
      Total->Veg.SafetyMargin[i] += VegMap[y][x].SafetyMargin[i];
      Total->Veg.HydStress[i]    += VegMap[y][x].HydStress[i];
    }
  }
  
  /* aggregate precipitation data */
  Total->Precip.Precip += Precip[y][x].Precip;
      Total->Precip.SnowFall += Precip[y][x].SnowFall;
  for (i = 0; i < NVegL; i++) {
	  Total->Precip.IntRain[i] += Precip[y][x].IntRain[i];
	  Total->Precip.IntSnow[i] += Precip[y][x].IntSnow[i];
	  Total->CanopyWater += Precip[y][x].IntRain[i] +
	  Precip[y][x].IntSnow[i];
  }

	/* aggregate radiation data */
	Total->Rad.Tair += RadMap[y][x].Tair;
  Total->Rad.ObsShortIn += RadMap[y][x].ObsShortIn;
  Total->Rad.BeamIn += RadMap[y][x].BeamIn;
  Total->Rad.DiffuseIn += RadMap[y][x].DiffuseIn;
  Total->Rad.PixelNetShort += RadMap[y][x].PixelNetShort;
  Total->NetRad += RadMap[y][x].NetRadiation[0] + RadMap[y][x].NetRadiation[1];

	/* aggregate snow data */
	if (Snow[y][x].HasSnow) {
	  NSnow++;
	  Total->Snow.Albedo += Snow[y][x].Albedo;
	  Total->Snow.LastSnow += Snow[y][x].LastSnow;
	  Total->Snow.PackWater += Snow[y][x].PackWater;
	  Total->Snow.TPack += Snow[y][x].TPack;
	  Total->Snow.SurfWater += Snow[y][x].SurfWater;
	  Total->Snow.TSurf += Snow[y][x].TSurf;
	  Total->Snow.ColdContent += Snow[y][x].ColdContent;
	  Total->Snow.Depth += Snow[y][x].Depth;
	  Total->Snow.Qe += Snow[y][x].Qe;
	  Total->Snow.Qs += Snow[y][x].Qs;
	  Total->Snow.Qsw += Snow[y][x].Qsw;
	  Total->Snow.Qlw += Snow[y][x].Qlw;
	  Total->Snow.Qp += Snow[y][x].Qp;
	  Total->Snow.MeltEnergy += Snow[y][x].MeltEnergy;
	}
  Total->Snow.Swq += Snow[y][x].Swq;
  Total->Snow.Melt += Snow[y][x].Outflow;
  Total->Snow.VaporMassFlux += Snow[y][x].VaporMassFlux;
  Total->Snow.CanopyVaporMassFlux += Snow[y][x].CanopyVaporMassFlux;

	if (VegMap[y][x].Gapping > 0.0) {
	  Total->Veg.Type[Opening].Qsw += VegMap[y][x].Type[Opening].Qsw;
	  Total->Veg.Type[Opening].Qlin += VegMap[y][x].Type[Opening].Qlin;
	  Total->Veg.Type[Opening].Qlw += VegMap[y][x].Type[Opening].Qlw;
	  Total->Veg.Type[Opening].Qe += VegMap[y][x].Type[Opening].Qe;
	  Total->Veg.Type[Opening].Qs += VegMap[y][x].Type[Opening].Qs;
	  Total->Veg.Type[Opening].Qp += VegMap[y][x].Type[Opening].Qp;
	  Total->Veg.Type[Opening].Swq += VegMap[y][x].Type[Opening].Swq;
	  Total->Veg.Type[Opening].MeltEnergy += VegMap[y][x].Type[Opening].MeltEnergy;
	}
	/* aggregate soil moisture data */
	Total->Soil.Depth += SoilMap[y][x].Depth;
	DeepDepth = 0.0;

	for (i = 0; i < NSoilL; i++) {
		Total->Soil.Moist[i] += SoilMap[y][x].Moist[i];
		if (SoilMap[y][x].Moist[i] <= 0.0)
		  SoilMap[y][x].Moist[i] = 0.0;
		Total->Soil.InterFlow[i] += SoilMap[y][x].InterFlow[i];
		
		Total->Soil.Perc[i] += SoilMap[y][x].Perc[i];
		Total->Soil.Temp[i] += SoilMap[y][x].Temp[i];
		Total->SoilWater += SoilMap[y][x].Moist[i] * VType[VegMap[y][x].Veg - 1].RootDepth[i] * Network[y][x].Adjust[i]; 
		DeepDepth += VType[VegMap[y][x].Veg - 1].RootDepth[i];
	}

	Total->Soil.Moist[Soil->MaxLayers] += SoilMap[y][x].Moist[NSoilL];
	Total->SoilWater += SoilMap[y][x].Moist[NSoilL] * (SoilMap[y][x].Depth - DeepDepth) * Network[y][x].Adjust[NSoilL];
	Total->Soil.TableDepth += SoilMap[y][x].TableDepth;

	if (SoilMap[y][x].TableDepth <= 0)
		(Total->Saturated)++;
	
	Total->Soil.WaterLevel += SoilMap[y][x].WaterLevel;
	Total->Soil.SatFlow += SoilMap[y][x].SatFlow;
	Total->Soil.DeepFlux += SoilMap[y][x].DeepFlux;
	Total->Soil.TSurf += SoilMap[y][x].TSurf;
	Total->Soil.Qnet += SoilMap[y][x].Qnet;
	Total->Soil.Qs += SoilMap[y][x].Qs;
	Total->Soil.Qe += SoilMap[y][x].Qe;
	Total->Soil.Qg += SoilMap[y][x].Qg;
	Total->Soil.Qst += SoilMap[y][x].Qst;
	Total->Soil.IExcess += SoilMap[y][x].IExcess;
	Total->Soil.DetentionStorage += SoilMap[y][x].DetentionStorage;
	
	if (Options->Infiltration == DYNAMIC)
		Total->Soil.InfiltAcc += SoilMap[y][x].InfiltAcc;
	
	Total->Soil.Runoff += SoilMap[y][x].Runoff;
	Total->ChannelInt += SoilMap[y][x].ChannelInt;
	Total->ChannelInfiltration += SoilMap[y][x].ChannelInfiltration;
	SoilMap[y][x].ChannelInt = 0.0;
	SoilMap[y][x].ChannelInfiltration = 0.0;
  }
  }
  }
  
  /* calculate average values for all quantities except the surface flow */
  
  /* Snow variables calculated for snow area only, except for mass variables */
  if (NSnow < 1)
    NSnow = 1;
  Total->Snow.HasSnow = (uchar) (((float) NSnow / (float) NPixels) * 255.0);
  
  /* average evaporation data */
  Total->Evap.ETot /= NPixels;
  for (i = 0; i < Veg->MaxLayers + 1; i++) {
    /* convert EPot from m/s to m */
    Total->Evap.EPot[i] /= NPixels;
    Total->Evap.EPot[i] *= Dt;
    Total->Evap.EAct[i] /= NPixels;
  }
  for (i = 0; i < Veg->MaxLayers; i++)
    Total->Evap.EInt[i] /= NPixels;
  for (i = 0; i < Veg->MaxLayers; i++) {
    for (j = 0; j < Soil->MaxLayers; j++) {
      Total->Evap.ESoil[i][j] /= NPixels;
    }
  }
  Total->Evap.EvapSoil /= NPixels;
  Total->Evap.EvapChannel /= NPixels;
  
  for (i = 0; i < Veg->MaxLayers; i++)
    Total->Veg.Rc[i] = (Total->Veg.Rc[i] > 0.0) ? (float)NPixels / Total->Veg.Rc[i] : DHSVM_HUGE;
  Total->Veg.PhotoDormancy /= NPixels;

  /* Extensive quantities over the whole basin; intensive ones over the
     vegetated fraction only (see the note where these are accumulated). */
  for (i = 0; i < Veg->MaxLayers; i++) {
    Total->Veg.AnCanopy[i] /= NPixels;
    Total->Veg.Escheme[i]  /= NPixels;
    Total->Veg.Ecrit[i]    /= NPixels;

    if (Total->NVegPix[i] > 0) {
      float N = (float) Total->NVegPix[i];
      Total->Veg.Tleaf[i]        /= N;
      Total->Veg.VpdLeaf[i]      /= N;
      Total->Veg.PsiSoil[i]      /= N;
      Total->Veg.PsiRoot[i]      /= N;
      Total->Veg.PsiLeaf[i]      /= N;
      Total->Veg.PLC[i]          /= N;
      Total->Veg.SafetyMargin[i] /= N;
      Total->Veg.HydStress[i]    /= N;
    }
  }
  
  /* average precipitation data */
  Total->Precip.Precip /= NPixels;
  Total->Precip.SnowFall /= NPixels;
  for (i = 0; i < Veg->MaxLayers; i++) {
    Total->Precip.IntRain[i] /= NPixels;
    Total->Precip.IntSnow[i] /= NPixels;
  }
  Total->CanopyWater /= NPixels;

  /* average radiation data */
  Total->Rad.Tair /= NPixels;
  Total->Rad.ObsShortIn /= NPixels;
  Total->Rad.PixelNetShort /= NPixels;
  Total->NetRad /= NPixels;
  Total->Rad.BeamIn /= NPixels;
  Total->Rad.DiffuseIn /= NPixels;
  for (i = 0; i < 2; i++) {
    Total->Rad.NetShort[i] /= NPixels;
	  Total->Rad.LongIn[i] /= NPixels;
	  Total->Rad.LongOut[i] /= NPixels;
  }

  /* average snow data */
  Total->Snow.Swq /= NPixels;
  Total->Snow.Melt /= NPixels;
  Total->Snow.PackWater /= NSnow;
  Total->Snow.TPack /= NSnow;
  Total->Snow.SurfWater /= NSnow;
  Total->Snow.TSurf /= NSnow;
  Total->Snow.ColdContent /= NSnow;
  Total->Snow.Albedo /= NSnow;
  Total->Snow.LastSnow /= NSnow;
  Total->Snow.Depth /= NSnow;
  Total->Snow.Qe /= NSnow;
  Total->Snow.Qs /= NSnow;
  Total->Snow.Qsw /= NSnow;
  Total->Snow.Qlw /= NSnow;
  Total->Snow.Qp /= NSnow;
  Total->Snow.MeltEnergy /= NSnow;
  Total->Snow.VaporMassFlux /= NPixels;
  Total->Snow.CanopyVaporMassFlux /= NPixels;

  if (TotNumGap > 0) {
	Total->Veg.Type[Opening].Qsw /= TotNumGap;
	Total->Veg.Type[Opening].Qlin /= TotNumGap;
	Total->Veg.Type[Opening].Qlw /= TotNumGap;
	Total->Veg.Type[Opening].Qe /= TotNumGap;
	Total->Veg.Type[Opening].Qs /= TotNumGap;
	Total->Veg.Type[Opening].Qp /= TotNumGap;
	Total->Veg.Type[Opening].Swq /= TotNumGap;
	Total->Veg.Type[Opening].MeltEnergy /= TotNumGap;
  }
  /* average soil moisture data */
  Total->Soil.Depth /= NPixels;
  for (i = 0; i < Soil->MaxLayers; i++) {
    Total->Soil.Moist[i] /= NPixels;
    Total->Soil.InterFlow[i] /= NPixels;
    Total->Soil.Perc[i] /= NPixels;
    Total->Soil.Temp[i] /= NPixels;
  }
  Total->Soil.Moist[Soil->MaxLayers] /= NPixels;
  Total->Soil.InterFlow[Soil->MaxLayers] /= NPixels;
  Total->Soil.TableDepth /= NPixels;
  Total->Soil.WaterLevel /= NPixels;
  Total->Soil.SatFlow /= NPixels;
  Total->Soil.DeepFlux /= NPixels;
  Total->Soil.TSurf /= NPixels;
  Total->Soil.Qnet /= NPixels;
  Total->Soil.Qs /= NPixels;
  Total->Soil.Qe /= NPixels;
  Total->Soil.Qg /= NPixels;
  Total->Soil.Qst /= NPixels;
  Total->Soil.IExcess /= NPixels;
  Total->Soil.DetentionStorage /= NPixels;
  
  if (Options->Infiltration == DYNAMIC)
    Total->Soil.InfiltAcc /= NPixels;

  Total->SoilWater /= NPixels;
  Total->Soil.Runoff /= NPixels;
  Total->ChannelInt /= NPixels;
  Total->ChannelInfiltration /= NPixels;
}
