
/* Used to calculate distance from meteorology stations for interpolation */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "settings.h"
#include "data.h"
#include "functions.h"

double CalcDistance(COORD * LocA, COORD * LocB)
{
  double Distance;
  double dN = (double) (LocA->N - LocB->N);
  double dE = (double) (LocA->E - LocB->E);
  Distance = sqrt(dN * dN + dE * dE);

  return Distance;
}
