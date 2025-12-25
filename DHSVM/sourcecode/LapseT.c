
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "functions.h"
#include "constants.h"

/*****************************************************************************
  Function name: LapseT()

  Purpose      : Lapse temperature with elevation

  Required     :
    float Temp       - Air temperature at FromElev elevation (C)
    float FromElev   - Elevation to lapse from (m)
    float ToElev     - Elevation to lapse to (m)
    float LapseRate  - Lapse rate in (m/m)

  Returns      :
    float LapsedTemp - Air temperature at ToElev elevation (C)

*****************************************************************************/
float LapseT(float Temp, float FromElev, float ToElev, float LapseRate)
{
  float LapsedTemp;

  LapsedTemp = Temp + (ToElev - FromElev) * LapseRate;

  return LapsedTemp;
}

