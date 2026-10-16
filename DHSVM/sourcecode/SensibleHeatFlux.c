
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "massenergy.h"
#include "constants.h"
#include "functions.h"

/*****************************************************************************
  Function name: NoSensibleHeatFlux()
 
  Comments     : This function sets all the energy fluxes at the pixel level 
                 to 0.0, but calculates the evapotranspiration in W/m2
*****************************************************************************/
void NoSensibleHeatFlux(int Dt, PIXMET *LocalMet, float ETot, SOILPIX *LocalSoil)
{
  LocalSoil->TSurf = 0.0;

  LocalSoil->Ra = 0.0;

  LocalSoil->Qnet = 0.0;

  LocalSoil->Qs = 0.0;

  LocalSoil->Qe = -(LocalMet->Lv * ETot) / Dt * WATER_DENSITY;

  LocalSoil->Qg = 0.0;

  LocalSoil->Qst = 0.0;

  LocalSoil->Qrest = 0.0;
}
