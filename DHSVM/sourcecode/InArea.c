
#include <stdio.h>
#include <stdlib.h>
#include "constants.h"
#include "settings.h"
#include "data.h"

/*****************************************************************************
  InArea()
*****************************************************************************/
uchar InArea(MAPSIZE * Map, COORD * Loc)
{
  if (Loc->N < 0 || Loc->N > (Map->NY - 1))
    return FALSE;
  else if (Loc->E < 0 || Loc->E > (Map->NX - 1))
    return FALSE;
  else
    return TRUE;
}
