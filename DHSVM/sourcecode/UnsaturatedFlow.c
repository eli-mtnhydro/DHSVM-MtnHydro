
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "constants.h"
#include "settings.h"
#include "functions.h"
#include "soilmoisture.h"

/*****************************************************************************
Function name: UnsaturatedFlow()

Purpose      : Calculate the unsaturated flow in the soil column, and
adjust the moisture in each soil layer

Sources:
Bras, R. A., Hydrology, an introduction to hydrologic science, Addisson
Wesley, Inc., Reading, etc., 1990.
Wigmosta, M. S., L. W. Vail, and D. P. Lettenmaier, A distributed
hydrology-vegetation model for complex terrain, Water Resour. Res.,
30(6), 1665-1679, 1994.

This function is based on Wigmosta et al [1994], and assumes a unit
hydraulic gradient in the unsaturated zone.  This implies a
steady-state situation and uniform moisture distribution.  This is not
unreasonable for situations where the groundwater is fairly deep
[Bras, 1990], but may be unrealistic for the boreal forest system, or
other similar situations where the groundwater table is relatively close
to the surface. The current formulation does not allow for upward
movement of water in the unsaturated zone, by capillary and matrix
forces.  No unsaturated flux is assumed to occur if the water content
drops below the field capacity.

The unsaturated hydraulic conductivity is calculated using the
Brooks-Corey equation see for example Wigmosta et al. [1994], eq.41

The calculated amount of drainage is averaged with the amount calculated
for the previous timestep, see eq. 42 [Wigmosta et al., 1994].

The residual moisture content is assumed to be zero (easy to adapt if
needed).

If the amount of soil moisture in a layer after drainage still
exceeds the porosity for that layer, this additional amount of water is
added to the drainage to the next layer.

CHANGES:

Changes have been made to account for the potental loss of soil storage
in a grid cell due to a channel.  Correction coefficents are
calculated in AdjustStorage() and CutBankGeometry()

Eli Boardman added unsaturated lateral gravimetric interflow (2026)
*****************************************************************************/
void UnsaturatedFlow(int Dt, float DX, float DY, float Infiltration,
  int NSoilLayers,
  float TotalDepth, float Area, float *RootDepth, float *Ks,
  float *PoreDist, float *Porosity, float *FCap,
  float *Perc, float *PercArea, float *Adjust,
  int CutBankZone, float BankHeight, float *TableDepth,
  float *IExcess, float *Moist, int InfiltOption,
  float *MoistDownhill, float *PorosityDownhill, float *InterFlowDownhill,
  float *RootDepthDownhill, float *AdjustDownhill, float LateralFrac)
{
  float DeepDrainage;		/* amount of drainage from the lowest root
                               zone to the layer below it (m) */
  float DeepLayerDepth;		/* depth of the layer below the deepest root layer */
  float Drainage;		    /* amount of water drained from each soil
                               layer during the current timestep */
  float Exponent;		    /* Brooks-Corey exponent */
  float FieldCapacity;		/* amount of water in soil at field capacity (m) */
  float MaxSoilWater;		/* maximum allowable amount of soil moiture in each layer (m) */
  float SoilWater;		    /* amount of water in each soil layer (m) */
  float LateralFracLayer; /* Fraction of total percolation that is routed laterally */
  float DownhillCapacity; /* Amount of water that neighbor cell's layer can accept */
  int i;			        /* counter */

  DeepLayerDepth = TotalDepth;
  for (i = 0; i < NSoilLayers; i++)
    DeepLayerDepth -= RootDepth[i];
  
  if (*TableDepth <= 0) { /* watertable above surface */
    *IExcess += Infiltration;
    if (InfiltOption == DYNAMIC) 
      Infiltration = 0.;
  }
  else {
    Moist[0] += Infiltration / (RootDepth[0] * Adjust[0]);
  }

  /* From top to bottom soil layer */
  for (i = 0; i < NSoilLayers; i++) {

    /* No movement if soil moisture is below field capacity */
    if (Moist[i] > FCap[i]) {
      Exponent = 2.0 / PoreDist[i] + 3.0;

      if (Moist[i] > Porosity[i])
        /* This can happen because the moisture content can exceed the
        porosity the way the algorithm is implemented */
        Drainage = Ks[i];
      else
        Drainage = Ks[i] * pow((double)(Moist[i]/Porosity[i]), (double)Exponent);
      /* Convert to m */
      Drainage *= Dt;

      /* Average across timesteps */
      Perc[i] = 0.5 * (Perc[i] + Drainage) * PercArea[i];

      MaxSoilWater = RootDepth[i] * Porosity[i] * Adjust[i];
      SoilWater = RootDepth[i] * Moist[i] * Adjust[i];
      FieldCapacity = RootDepth[i] * FCap[i] * Adjust[i];

      /* No unsaturated flow if the moisture content drops below field
      capacity */

      if ((SoilWater - Perc[i]) < FieldCapacity)
        Perc[i] = SoilWater - FieldCapacity;

      /* If the moisture content is greater than the porosity add the
      additional soil moisture to the percolation */
      SoilWater -= Perc[i];
      if (SoilWater > MaxSoilWater)
        Perc[i] += SoilWater - MaxSoilWater;

      /* Adjust the moisture content in the current layer, and the layer
      immediately below it */
      Moist[i] -= Perc[i] / (RootDepth[i] * Adjust[i]);
      if (i < (NSoilLayers - 1)) {
        /* Only allow enough lateral flow to raise neighboring layer to saturation */
        DownhillCapacity = (PorosityDownhill[i] - MoistDownhill[i]) * RootDepthDownhill[i] * AdjustDownhill[i];
        if (DownhillCapacity > (LateralFrac * Perc[i])) {
          LateralFracLayer = LateralFrac;
        } else if (Perc[i] > 1e-10) {
          LateralFracLayer = DownhillCapacity / Perc[i];
        } else {
          LateralFracLayer = 0.0;
        }
        
        Moist[i + 1] += (1.0 - LateralFracLayer) * Perc[i] / (RootDepth[i + 1] * Adjust[i + 1]);
        InterFlowDownhill[i] += LateralFracLayer * Perc[i] / (RootDepthDownhill[i] * AdjustDownhill[i]);
      }
    }
    else
      Perc[i] = 0.0;

    /* Convert back to straight 1-d flux */
    Perc[i] /= PercArea[i];
  }
  
  DeepDrainage = (Perc[NSoilLayers - 1] * PercArea[NSoilLayers - 1]);
  
  /* Only allow enough lateral flow to raise neighboring layer to saturation */
  DownhillCapacity = (PorosityDownhill[NSoilLayers - 1] - MoistDownhill[NSoilLayers - 1]) *
                      RootDepthDownhill[NSoilLayers - 1] * AdjustDownhill[NSoilLayers - 1];
  if (DownhillCapacity > (LateralFrac * DeepDrainage)) {
    LateralFracLayer = LateralFrac;
  } else if (DeepDrainage > 1e-10) {
    LateralFracLayer = DownhillCapacity / DeepDrainage;
  } else {
    LateralFracLayer = 0.0;
  }
  
  Moist[NSoilLayers] += (1.0 - LateralFracLayer) * DeepDrainage / (DeepLayerDepth * Adjust[NSoilLayers]);
  InterFlowDownhill[NSoilLayers - 1] += LateralFracLayer * DeepDrainage / (RootDepthDownhill[NSoilLayers - 1] * AdjustDownhill[NSoilLayers - 1]);

  /* Calculate the depth of the water table based on the soil moisture
  profile and adjust the soil moisture profile, to assure that the soil
  moisture is never more than the maximum allowed soil moisture amount,
  i.e. the porosity.  A negative water table depth means that the water is
  ponding on the surface.  This amount of water becomes surface IExcess */

  *TableDepth = WaterTableDepth(NSoilLayers, TotalDepth, RootDepth, Porosity, FCap, Adjust, Moist);
  
  if (*TableDepth < 0.0) {
    *IExcess += -(*TableDepth);
    if (InfiltOption == DYNAMIC) {
      if (Infiltration > -(*TableDepth))
        Infiltration += *TableDepth;
      else
        Infiltration = 0.;
    }

    *TableDepth = 0.0;
  }
}
