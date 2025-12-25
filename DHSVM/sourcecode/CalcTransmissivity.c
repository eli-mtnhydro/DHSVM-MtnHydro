
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "functions.h"

/*****************************************************************************
  Function name: CalcTransmissivity()

  Purpose      : Calculates the transmissivity through the saturated part of
                 the soil profile
                 
  Required     : 
    float SoilDepth  - Total soil depth in m
    float WaterTable - Depth of the water table below the soil surface in m
    float LateralKs  - Lateral hydraulic conductivity in m/s
    float KsExponent - Exponent that describes exponential decay of LateralKs
                       with depth below the soil surface

  Returns      : Transmissivity in m2/s

  Comments     :
    Source:
    Wigmosta, M. S., L. W. Vail, and D. P. Lettenmaier, A distributed 
      hydrology-vegetation model for complex terrain, Water Resour. Res.,
      30(6), 1665-1679, 1994.

    Based on:
    Beven, K. J., On subsurface stormflow:  An analysis of response times,
      Hydrolog. Sci. J., 4, 505-521, 1982.

    The hydraulic conductivity is assumed exponentially with depth, based on
    material in Beven [1982].
*****************************************************************************/
float CalcTransmissivity(float SoilDepth, float WaterTable, float LateralKs,
			 float KsExponent, float DepthThresh)
{
  float Transmissivity;		/* Transmissivity (m^2/s) */
  float TransThresh;

  if (fequal(KsExponent, 0.0))
    Transmissivity = LateralKs * (SoilDepth - WaterTable);
  else {
	/* A smaller value of WaterTable variables indicates a higher actual water table depth */
	if (WaterTable < DepthThresh) {
	  Transmissivity = (LateralKs / KsExponent) * (exp(-KsExponent * WaterTable) - exp(-KsExponent * SoilDepth));
	} else if (SoilDepth < DepthThresh) {
	  /* Water table > depth thresh, but soil depth < depth thresh, so no transmissivity */
	  Transmissivity = 0.0;
	} else {
	  /* Water table and soil depth both deeper than depth thresh,
	     transmissivity decreases linearly with water table depth */
	  TransThresh = (LateralKs / KsExponent) * (exp(-KsExponent * DepthThresh) - exp(-KsExponent * SoilDepth));
	  Transmissivity = ((SoilDepth-WaterTable) / (SoilDepth-DepthThresh)) * TransThresh;
	}
  }
  return Transmissivity;
}
