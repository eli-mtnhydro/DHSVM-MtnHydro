
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"
#include "rad.h"

 /*****************************************************************************
   GetMetData()
 *****************************************************************************/
void GetMetData(OPTIONSTRUCT *Options, TIMESTRUCT *Time,
  int NStats, float SunMax, METLOCATION *Stat)
{
  int i;			/* counter */
  float Ea;			/* vapor pressure (Pa) */
  float Press;			/* air pressure (Pa) */

  if (DEBUG)
    Message("Reading all met data for current timestep\n");

  for (i = 0; i < NStats; i++)
    ReadMetRecord(&(Time->Current), &(Stat[i]));

  /* Specific humidity of each station from its relative humidity, air
     temperature and air pressure (its elevation and lapse rate) */
  for (i = 0; i < NStats; i++) {
    Ea = SatVaporPressure(Stat[i].Data.Tair) * (Stat[i].Data.Rh / 100.);
    Press = AirPressure(Stat[i].Elev, Stat[i].Data.TempLapse);
    Stat[i].Data.Q = EPS * Ea / (Press - (1.0 - EPS) * Ea);
  }

  for (i = 0; i < NStats; i++) {
    if (SunMax > 0.0) {
      Stat[i].Data.ClearIndex = Stat[i].Data.Sin / SunMax;
      SeparateRadiation(Stat[i].Data.Sin, Stat[i].Data.ClearIndex,
        &(Stat[i].Data.SinBeamObs),
        &(Stat[i].Data.SinDiffuseObs));
    }
    else {
      /* if sun is below horizon, then force all shortwave to zero */
      Stat[i].Data.Sin = 0.0;
      Stat[i].Data.SinBeamObs = 0.0;
      Stat[i].Data.SinDiffuseObs = 0.0;
    }
  }
}
