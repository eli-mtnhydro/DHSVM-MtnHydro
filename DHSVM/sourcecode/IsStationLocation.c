
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "functions.h"

#define NOTSAME -1

/*****************************************************************************
function  : IsStationLocation()
input     : Current Location, address of structure with station data, and
            pointer to int indicating which station is at current location 
	    (if any)

This routine determines whether there is a station at the current location
*****************************************************************************/
uchar IsStationLocation(COORD * Loc, int NStats, METLOCATION * Station,
			int *WhichStation)
{
  int i;			/* Station Counter */

  for (i = 0, *WhichStation = NOTSAME; i < NStats &&
       *WhichStation == NOTSAME; i++) {
    if (Loc->N == Station[i].Loc.N && Loc->E == Station[i].Loc.E) {
      *WhichStation = i;
      return TRUE;
    }
  }

  return FALSE;
}
