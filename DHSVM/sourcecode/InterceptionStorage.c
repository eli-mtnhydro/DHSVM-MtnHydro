
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "massenergy.h"
#include "constants.h"

/*****************************************************************************
InterceptionStorage()
*****************************************************************************/
void InterceptionStorage(int NAct, float *MaxInt, float *Fract,
  float *Int, float *Precip)
{
  float Available;		/* Available storage */
  float Intercepted;	/* Amount of water intercepted during this timestep */
  int i;			    /* counter */

                        /* The precipitation is multiplied by the fractional coverage, since if the
                        vegetation covers only 10% of the grid cell, only 10% can be intercepted as a
                        maximum */
  for (i = 0; i < NAct; i++) {
    Available = MaxInt[i] - Int[i];
    if (Available > *Precip * Fract[i])
      Intercepted = (*Precip) * Fract[i];
    else
      Intercepted = Available;
    *Precip -= Intercepted;
    Int[i] += Intercepted;
  }
}

/*****************************************************************************
CanopyGapInterception()
*****************************************************************************/
void CanopyGapInterceptionStorage(int NAct, float *MaxInt,
  float *Fract, float *Int, float *Precip)
{
  float Available;		/* Available storage */
  float Intercepted;	/* Amount of water intercepted during this timestep */
  int i;			    /* counter */

  if (NAct >= 1) {
    i = 1;
    Available = MaxInt[i] - Int[i];
    if (Available > *Precip * Fract[i])
      Intercepted = (*Precip) * Fract[i];
    else
      Intercepted = Available;
    *Precip -= Intercepted;
    Int[i] += Intercepted;
  }

}
