
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"

/*****************************************************************************
  Function name: InitSnowMap()

  Purpose      : Initialize the snow information for each pixel in the basin

  Required     :
    MAPSIZE Map        - Size and location of the model area
    SNOWPIX ***SnowMap - Address of array with snow information 

*****************************************************************************/
void InitSnowMap(MAPSIZE *Map, SNOWPIX ***SnowMap, TIMESTRUCT *Time)
{
  const char *Routine = "InitSnowMap";
  int y;			/* counter */
  
  printf("Initializing snow map\n");

  if (!(*SnowMap = (SNOWPIX **) calloc(Map->NY, sizeof(SNOWPIX *))))
    ReportError((char *) Routine, 1);

  for (y = 0; y < Map->NY; y++) {
    if (!((*SnowMap)[y] = (SNOWPIX *) calloc(Map->NX, sizeof(SNOWPIX))))
      ReportError((char *) Routine, 1);
  }
}
