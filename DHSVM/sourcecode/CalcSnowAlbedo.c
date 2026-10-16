
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "constants.h"
#include "data.h"
#include "Calendar.h"
#include "functions.h"

/*****************************************************************************
  CalcSnowAlbedo()

  Source:
  Laramie, R. L., and J. C. Schaake, Jr., Simulation of the continuous
  snowmelt process, Ralph M. Parsons Laboratory, Mass. Inst. of Technol.,
  1972

  Snow albedo is calculated as a function of the number of days since the
  last observed snow fall. There are separete albedo curves for the freeze
  and thaw conditions.
 
  Laramie and Schaake (1972)
  Updated based on Storck (2000)
  Updated by Eli Boardman 2024 to fix issues and introduce reset threshold
*****************************************************************************/
float CalcSnowAlbedo(SNOWPIX *LocalSnow, float NewSnow)
{
  
  float Last, Albedo, GroundFrac, Depth;
  
  Last = LocalSnow->LastSnow;
  
  if (Last > (float) DAYPYEAR)
    Last = (float) DAYPYEAR;
  
  if (LocalSnow->AccumSeason == TRUE) {
    /* Accumulation season */
    Albedo = LocalSnow->amax * pow(LocalSnow->LamdaAcc, pow(Last, 0.58));
  } else {
    /* Melt season */
    Albedo = LocalSnow->amax * pow(LocalSnow->LamdaMelt, pow(Last, 0.46));
  }
  
  /* Snow (the step's snowfall NewSnow included) shallower than the reset rate
     over one day, in m of SWE, blends in the ground albedo by depth, as the
     step's snowfall blends into an older pack's days since snowfall */
  Depth = LocalSnow->Swq + NewSnow;
  if (Depth < MIN_SNOW_RESET_ALBEDO) {
    GroundFrac = 1.0 - Depth / MIN_SNOW_RESET_ALBEDO;
    Albedo = (GroundFrac * LocalSnow->AlbedoGround) + ((1.0 - GroundFrac) * Albedo);
  }
  
  /* Enforce constraints */
  if (LocalSnow->AccumSeason == TRUE && Albedo < LocalSnow->AccMin)
    Albedo = LocalSnow->AccMin;
  else if (Albedo < LocalSnow->MeltMin)
    Albedo = LocalSnow->MeltMin;
  
  return Albedo;
}
