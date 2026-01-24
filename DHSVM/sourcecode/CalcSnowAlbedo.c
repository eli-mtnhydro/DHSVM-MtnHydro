
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
float CalcSnowAlbedo(SNOWPIX *LocalSnow, int StepsPerDay)
{
  
  float Last, Albedo, GroundFrac;
  
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
  
  /* Average snow albedo with ground albedo for shallow snowpack */
  if (LocalSnow->Swq < MIN_SNOW_RESET_ALBEDO) {
    GroundFrac = 1.0 - LocalSnow->Swq / MIN_SNOW_RESET_ALBEDO;
    Albedo = (GroundFrac * LocalSnow->AlbedoGround) + ((1.0 - GroundFrac) * Albedo);
  }
  
  /* Enforce constraints */
  if (LocalSnow->AccumSeason == TRUE && Albedo < LocalSnow->AccMin)
    Albedo = LocalSnow->AccMin;
  else if (Albedo < LocalSnow->MeltMin)
    Albedo = LocalSnow->MeltMin;
  
  return Albedo;
}
