
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
void GetMetData(OPTIONSTRUCT *Options, TIMESTRUCT *Time, int NSoilLayers,
  int NStats, float SunMax, METLOCATION *Stat)
{
  int i;			/* counter */

  if (DEBUG)
    printf("Reading all met data for current timestep\n");

  for (i = 0; i < NStats; i++)
    ReadMetRecord(Options, &(Time->Current), NSoilLayers, &(Stat[i].MetFile), &(Stat[i].Data));

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
