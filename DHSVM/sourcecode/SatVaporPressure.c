
#include <stdlib.h>
#include <math.h>
#include "lookuptable.h"

float CalcVaporPressure(float T);
static FLOATTABLE svp;		/* Table that contains saturated vapor 
				   pressures as a function of temperature 
				   in degrees C */

/*****************************************************************************
  Function name: InitSatVaporTable()

  Purpose      : Initialize lookup table for saturated vapor pressure as a 
                 function of temperature in degrees C.

  Comments     :  Table runs from -100 C to 100 C with an interval of 0.02 C
*****************************************************************************/
void InitSatVaporTable(void)
{
  InitFloatTable(30000L, -300., .02, CalcVaporPressure, &svp);
}

/*****************************************************************************
  Function name: CalcVaporPressure() 

  Purpose      : Calculates the saturated vapor pressure in Pa for a certain
                 temperature in degrees C 
  Required     : 
    float T    - Temperature (C)
  Returns      :
    float      - Saturated vapor pressure (Pa) 

    References: Shuttleworth, W.J., Evaporation,  In: Maidment, D. R. (ed.),
                  Handbook of hydrology,  1993, McGraw-Hill, New York, etc..
                Bras, R. A., Hydrology, an introduction to hydrologic
                  science, Addisson Wesley, Inc., Reading, etc., 1990.

*****************************************************************************/
float CalcVaporPressure(float T)
{
  float Pressure;

  Pressure = 610.78 * exp((double) ((17.269 * T) / (237.3 + T)));

  /* Calculate the saturated vapor pressure in the snow pack, 
     (Equation 3.32, Bras 1990) */

  if (T < 0.0)
    Pressure *= 1.0 + .00972 * T + .000042 * T * T;

  return Pressure;
}

/*****************************************************************************
  Function name: SatVaporPressure() - new version, using a lookup table

  Purpose      : Looks up the saturated vapor pressure in Pa for a certain
                 temperature in a table
*****************************************************************************/
float SatVaporPressure(float T)
{
  return FloatLookup(T, &svp);
}