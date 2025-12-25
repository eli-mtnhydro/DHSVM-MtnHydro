
/*
 Calculate the maximum amount of moisture the soil can
 deliver to the atmosphere in one time step
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "massenergy.h"
#include "constants.h"

/*****************************************************************************
  Desorption()
*****************************************************************************/
float Desorption(int Dt, float MoistContent, float Porosity, float Ks,
		 float Press, float m)
{
  float Sorptivity;		/* sorptivity */
  float DesorptionVolume;	/* total desorption volume during timestep */

  /* Eq. 46, Wigmosta et al [1994] */

  if (MoistContent > Porosity)
    MoistContent = Porosity;

  Sorptivity = sqrt((double)
		    ((8 * Porosity * Ks * Press) /
		     (3.0 * (1 + 3 * m) * (1 + 4 * m)))) *
    pow((double) (MoistContent / Porosity), (double) (1.0 / (2.0 * m) + 2));

  /* Eq. 45, Wigmosta et al [1994] */
  
  DesorptionVolume = Sorptivity * sqrt((double) Dt);

  return DesorptionVolume;
}
