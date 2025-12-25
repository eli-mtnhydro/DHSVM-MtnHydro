
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "massenergy.h"
#include "constants.h"

/*****************************************************************************
  Function name: StabilityCorrection()

  Purpose      : Calculate atmospheric stability correction for non-neutral
                 conditions

  Required     :
    float Z          - Reference height (m)
    float d          - Displacement height (m)
    float TSurf      - Surface temperature (C)
    float Tair       - Air temperature (C)
    float Wind       - Wind speed (m/s)
    float Z0         - Roughness length (m)

  Returns      :
    float Correction - Multiplier for aerodynamic resistance

*****************************************************************************/
float StabilityCorrection(float Z, float d, float TSurf, float Tair,
			  float Wind, float Z0)
{
  float Correction;		/* Correction to aerodynamic resistance */
  float Ri;			/* Richardson's Number */
  float RiCr = 0.2;		/* Critical Richardson's Number */
  float RiLimit;		/* Upper limit for Richardson's Number */

  Correction = 1.0;

  /* Calculate the effect of the atmospheric stability using a Richardson 
     Number approach */

  if (TSurf != Tair) {

    /* Non-neutral conditions */

    Ri = G * (Tair - TSurf) * (Z - d) /
      (((Tair + 273.15) + (TSurf + 273.15)) / 2.0 * Wind * Wind);

    RiLimit = (Tair + 273.15) /
      (((Tair + 273.15) + (TSurf + 273.15)) / 2.0 * (log((Z - d) / Z0) + 5));

    if (Ri > RiLimit)
      Ri = RiLimit;

    if (Ri > 0.0)
      Correction = (1 - Ri / RiCr) * (1 - Ri / RiCr);

    else {
      if (Ri < -0.5)
	Ri = -0.5;

      Correction = sqrt(1 - 16 * Ri);
    }
  }

  return Correction;
}
