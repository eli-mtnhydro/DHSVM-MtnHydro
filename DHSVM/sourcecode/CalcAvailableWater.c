
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "soilmoisture.h"

 /****************************************************************************
   Function name: CalcAvailableWater()

   Purpose      : This routine calculates the amount of soil moisture
                  available for saturated flow below the groundwater table.
                  No flow is allowed to occur if the moisture content falls
                  below the field capacity.

   Required     :
     int NRootLayers  - Number of soil layers
     float TotalDepth - Total depth of the soil profile (m)
     float *RootDepth - Depth of each of the soil layers (m)
     float *Porosity  - Porosity of each soil layer
     float *FCap      - Field capacity of each soil layer
     float TableDepth - Depth of the water table below the surface (m)
     float *Adjust    - Correction for each layer for loss of soil storage
                        due to channel cut.  Multiplied with RootDepth
                        to give the layer thickness for use in calculating
                        soil moisture

   Returns      : float AvailableWater - Total water available for saturated flow

 *****************************************************************************/
float CalcAvailableWater(int NRootLayers, float TotalDepth, float *RootDepth,
  float *Porosity, float *FCap, float *Moist, float TableDepth,  float *Adjust)

{
  float AvailableWater;		/* amount of water available for movement (m) */

  float DeepFCap;		    /* field capacity of the layer below the  deepest root layer */

  float DeepLayerDepth;		/* depth of layer below deepest root zone layer */

  float DeepPorosity;		/* porosity of the layer below the deepest root layer */

  float Depth;			    /* depth below the ground surface (m) */
  int i;			        /* counter */
 
  DeepPorosity = Porosity[NRootLayers - 1];
  DeepFCap = FCap[NRootLayers - 1];
  
  DeepLayerDepth = TotalDepth;
  for (i = 0; i < NRootLayers; i++)
    DeepLayerDepth -= RootDepth[i];
  
  AvailableWater = 0.0;

  Depth = 0.0;
  for (i = 0; i < NRootLayers && Depth < TotalDepth; i++) {
    if (RootDepth[i] < (TotalDepth - Depth))
      Depth += RootDepth[i];
    else
      Depth = TotalDepth;
    if (Depth > TableDepth) {
      if ((Depth - TableDepth) > RootDepth[i])
        AvailableWater += (Porosity[i] - FCap[i]) * RootDepth[i] * Adjust[i];
      else
        AvailableWater += (Moist[i] - FCap[i]) * RootDepth[i] * Adjust[i];
    }
  }

  if (Depth < TotalDepth) {
    Depth = TotalDepth;

    if ((Depth - TableDepth) > DeepLayerDepth)
      AvailableWater += (DeepPorosity - DeepFCap) * DeepLayerDepth * Adjust[NRootLayers];
    else
      AvailableWater += (Moist[NRootLayers] - DeepFCap) * DeepLayerDepth * Adjust[NRootLayers];
    
  }
  
  if (AvailableWater < 0.0)
    AvailableWater = 0.0;
  
  return AvailableWater;
}
