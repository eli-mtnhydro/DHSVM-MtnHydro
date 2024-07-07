/*
* SUMMARY:      DistrubteSatflow.c - 
* USAGE:        Part of DHSVM
*
* AUTHOR:       Bart Nijssen and Mark Wigmosta (*)
* ORG:          University of Washington, Department of Civil Engineering
* E-MAIL:       nijssen@u.washington.edu
* ORIG-DATE:    Apr-1996
* DESCRIPTION:  Distribute the SatFlow from Previous Timestep
* DESCRIP-END.
* FUNCTIONS:    DistrubteSatflow()
* COMMENTS: (*) Mark Wigmosta, Batelle Pacific Northwest Laboratories,
*               ms_wigmosta@pnl.gov
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "constants.h"
#include "settings.h"
#include "functions.h"
#include "soilmoisture.h"

/*****************************************************************************/
void DistributeSatflow(int Dt, float DX, float DY, float SatFlow,
  int NSoilLayers, float TotalDepth, float *RootDepth,
  float *Porosity, float *FCap, float *Adjust,
  float *TableDepth, float *Runoff, float *Moist)

{
  float DeepLayerDepth;		/* depth of the layer below the deepest root layer */
  int i;			        /* counter */
  float DeepFCap;		/* field capacity of the layer below the deepest root layer */
  float DeepPorosity;		/* porosity of the layer below the deepest root layer */
  float DeepAvaWater;
  float AvaWater;
  float DeepWaterGap;
  float WaterGap;
  float ExtracWater;
  float DeepExtracWater;
  float Depth;

  DeepPorosity = Porosity[NSoilLayers - 1];
  DeepFCap = FCap[NSoilLayers - 1];
  
  DeepLayerDepth = TotalDepth;
  for (i = 0; i < NSoilLayers; i++)
    DeepLayerDepth -= RootDepth[i];

  /* Added 06/09/2016 by Zhuoran Duan(zhuoran.duan@pnnl.gov) */
  /* When calculating lateral outflow, available water was calculated using all 3 root zone layers 
  and the deep layer underneath but outflow was extracted from the bottom deep soil layer. This 
  might lead to negative soil moisture in deep soil while the layer above it remains saturated. This tends
  to happen more often in dry climate. In order to avoid negative deep soil moisture, here I add a loop 
  to redistribution of water extraction. The outflow starts from the (top) water table layer, extract the 
  excess water larger than field capacity. While there's no enough water from the current layer, extract 
  water from one layer below it until it reaches bottom layer. */
  
  /*New algorithm for SatFlow, remove water from top layer to bottom layer*/
  if (SatFlow < 0.0) {
    
    Depth = 0.0;
    for (i = 0; i < NSoilLayers && Depth < TotalDepth; i++) {
      AvaWater = 0.0;
      ExtracWater = 0.0;
      if (RootDepth[i] < (TotalDepth - Depth))
        Depth += RootDepth[i];
      else
        Depth = TotalDepth;

      /* if water table is in the ith root zone layer */
      if (Depth > *TableDepth) {
        /* fully saturated in the current layer */
        if ((Depth - *TableDepth) > RootDepth[i])
          AvaWater = (Porosity[i]-FCap[i]) * RootDepth[i] * Adjust[i];
        else
          AvaWater = (Moist[i]-FCap[i]) * RootDepth[i] * Adjust[i];
      }
      ExtracWater = (-SatFlow > AvaWater) ? -AvaWater : SatFlow;

      Moist[i] += ExtracWater / (RootDepth[i] * Adjust[i]);

      SatFlow -= ExtracWater;
      if (fequal(SatFlow, 0.0))
        break;
    }

    if (SatFlow < 0.0) {
      DeepAvaWater = 0.0;
      DeepExtracWater = 0.0;
      if (Depth < TotalDepth) {
        Depth = TotalDepth;
        if ((Depth - *TableDepth) > DeepLayerDepth)
          DeepAvaWater = (DeepPorosity - DeepFCap) * DeepLayerDepth * Adjust[NSoilLayers];
        else
          DeepAvaWater = (Moist[NSoilLayers] - DeepFCap) * DeepLayerDepth * Adjust[NSoilLayers];;
      }

      DeepExtracWater = (-SatFlow > DeepAvaWater) ? -DeepAvaWater : SatFlow;

      Moist[NSoilLayers] += DeepExtracWater / (DeepLayerDepth * Adjust[NSoilLayers]);

      SatFlow -= DeepExtracWater;
    }
  } else if (SatFlow > 0.0) {

    Depth = DeepLayerDepth;
    DeepWaterGap = 0.0;
    DeepExtracWater = 0.0;
    if (Depth > (TotalDepth - *TableDepth)) {
      DeepWaterGap = (DeepPorosity - Moist[NSoilLayers]) * DeepLayerDepth * Adjust[NSoilLayers];
      DeepExtracWater = (SatFlow > DeepWaterGap) ? DeepWaterGap : SatFlow;
      SatFlow -= DeepExtracWater;
      Moist[NSoilLayers] += DeepExtracWater / (DeepLayerDepth * Adjust[NSoilLayers]);
    }

    if (SatFlow > 0.0) {
      for (i = NSoilLayers - 1; i >= 0; i--) {
        WaterGap = 0.0;
        ExtracWater = 0.0;
        Depth += RootDepth[i];
        if (Depth > (TotalDepth - *TableDepth)) {
          WaterGap = (Porosity[i] - Moist[i]) * RootDepth[i] * Adjust[i];
          ExtracWater = (SatFlow > WaterGap) ? WaterGap : SatFlow;
          SatFlow -= ExtracWater;
          Moist[i] += ExtracWater / (RootDepth[i] * Adjust[i]);
        }
        if (fequal(SatFlow, 0.0))
          break;
      }
    }
  }

  if (SatFlow > 0.0)
    *Runoff += SatFlow;
}
