
#ifndef SLOPEASPECT_H
#define SLOPEASPECT_H

#include "settings.h"
#include "data.h"

/* -------------------------------------------------------------
   available variables
   ------------------------------------------------------------- */
extern int xneighbor[NNEIGHBORS];
extern int yneighbor[NNEIGHBORS];

/* -------------------------------------------------------------
   available functions
   ------------------------------------------------------------- */
void ElevationSlopeAspect(MAPSIZE * Map, TOPOPIX ** TopoMap, int MultiFlowDir);
void HeadSlopeAspect(MAPSIZE * Map, TOPOPIX ** TopoMap, SOILPIX ** SoilMap,
  float **FlowGrad, unsigned char ***Dir, unsigned int **TotalDir, int MultiFlowDir, int x, int y);
void SnowSlopeAspect(MAPSIZE * Map, TOPOPIX ** TopoMap, SNOWPIX ** Snow,
  float **FlowGrad, unsigned char ***Dir, unsigned int **TotalDir, int MultiFlowDir);
int valid_cell(MAPSIZE * Map, int x, int y);
void quick(ITEM *OrderedCells, int count);
#endif

