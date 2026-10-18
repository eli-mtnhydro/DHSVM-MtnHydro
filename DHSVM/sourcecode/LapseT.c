#include <math.h>

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

/*****************************************************************************
  AirPressure()

  Air pressure (Pa) at elevation Elev (m) of the standard atmosphere with
  the temperature lapse rate LapseRate (C/m), 101300 Pa where it is 0
*****************************************************************************/
float AirPressure(float Elev, float LapseRate)
{
  float Exponent;

  if (LapseRate != 0.0) {
    Exponent = 9.8067 / (LapseRate * 287.0);
    return 101300. * pow(((288.0 - LapseRate * Elev) / 288.0), Exponent);
  }
  return 101300.;
}

