
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "settings.h"
#include "data.h"
#include "functions.h"
#include "slopeaspect.h"
#include "DHSVMerror.h"

float temp_aspect[NNEIGHBORS] = {
  225., 180., 135., 90., 45., 0., 315., 270.
};

/* NNEIGHBORS used for redistribution of subsurface flow in topoindex
 * and for slope/aspect calculations, must equal 8. */
int xneighbor[NNEIGHBORS] = {
-1, 0, 1, 1, 1, 0, -1, -1
};

int yneighbor[NNEIGHBORS] = {
 1, 1, 1, 0, -1, -1, -1, 0
};

/* -------------------------------------------------------------
   valid_cell
   Checks to see if grid indices, x and y, are within the grid 
   defined by the specified Map
   ------------------------------------------------------------- */
int valid_cell(MAPSIZE *Map, int x, int y)
{
  return (x >= 0 && y >= 0 && x < Map->NX && y < Map->NY);
}

/* -------------------------------------------------------------
   slope_aspect
   Calculation of slope and aspect given elevations of cell and neighbors
   ------------------------------------------------------------- */
/* The slope is used by the D8 routing gradient, the kinematic overland
   routing and the shading, the aspect by the shading */
static void slope_aspect(float dx, float dy, float celev, float
			 nelev[NNEIGHBORS], float *slope, float *aspect)
{
  int n;
  float dzdx, dzdy;
  float *dummyelev;
  /* This dummy variable is added for calculation of elev difference,
  in which the elev of OUTSIDEBASIN cells (which is ZERO) is
  replaced by the elev of the central cell */

  /* allocate memory */
  if (!(dummyelev = (float*) calloc(NNEIGHBORS, sizeof(float))))
	  ReportError("slope_aspect( )", 1);

  for (n = 0; n < NNEIGHBORS; n++) {
      if (nelev[n] == OUTSIDEBASIN) {
		  dummyelev[n] = celev;
      }
	  else
		  dummyelev[n] = nelev[n];
    }
  dzdx = ((dummyelev[0] + 2 * dummyelev[7] + dummyelev[6]) -
	    (dummyelev[2] + 2 * dummyelev[3] + dummyelev[4])) / (8 * dx);
  dzdy = ((dummyelev[0] + 2 * dummyelev[1] + dummyelev[2]) -
	    (dummyelev[4] + 2 * dummyelev[5] + dummyelev[6])) / (8 * dy);

  *slope = sqrt(dzdx * dzdx + dzdy * dzdy);
  if (fequal(dzdx, 0.0) && fequal(dzdy, 0.0)) {
    *aspect = 0.0;
  }
  else {
	  /* convert from radian to degree */
	  *aspect = atan2(dzdx, dzdy) ;
  }
  free(dummyelev);
  return;
}

/* -------------------------------------------------------------
   flow_fractions
   Computes subsurface flow fractions given the slope and aspect 

   Comment: this function is considerably modified to avoid any
   out flow to the cells outside of the basin mask (Ning, 2013)
 
 Update May 2024 by Eli Boardman
 The multi-direction 4-neighbor algorithm used in the original version
 can cause problems when the 8-neighbor average slope from slope_aspect()
 conflicts with the true flow path due to local terrain heterogeneity,
 leading to infinite loops and unphysical surface ponding.
 However, the D8 algorithm does not allow flow to spread out,
 which leads to over-concentrated and unrealistic hillslope flow paths.
 This update implements a multi-direction 8-neighbor flow algorithm
 that routes water to all downhill cells in proportion to the local slope.
 Note that the effective flow width now has a circular rather than
 rectangular interpretation: if only one cell is downhill,
 the effective flow width is dx * PI/8, and if all 8 cells are downhill,
 the effective flow width is dx * PI.
 
------------------------------------------------------------- */
static void flow_fractions(float dx, float dy, float slope, float aspect,
			   float celev, float nelev[NDIRS], float *grad,
			   unsigned char dir[NDIRS], unsigned int *total_dir, int MultiFlowDir)
{
  float total_width, effective_width;
  int n;
  float drop[NDIRS];

    /* Flow goes to all downhill neighbors in proportion to the slope, or
       just to the single steepest */
    float total_drop;

    if (MultiFlowDir) {
      /* For MFD8 flow directions, water discharges to ALL of its downhill neighbors: 
         As with D8, this requires the DEM to be pre-filled for D8 routing scheme
         as any flat areas will confuse the model. */
      
      total_drop = 0.0;
      *grad = 0.0;
      effective_width = dx * PI / 8; // Per direction
      
      for (n = 0; n < NDIRS; n++) {
        /* Make sure flow is inside boundary */
        if (nelev[n] == (float) OUTSIDEBASIN){
          dir[n] = 0;
          drop[n] = 0.0;
        } else {
          /* Find all downhill cells */
          if (n == 0 || n == 2 || n == 4 || n == 6)
            drop[n] = (celev - nelev[n]) / sqrt(dx * dx + dy * dy);
          else
            drop[n] = (celev - nelev[n]) / dx;
          
          /* Add up total drop in all directions */
          if (drop[n] > 0.0) {
            total_drop += drop[n];
            *grad += drop[n] * effective_width;
          }
          else
            drop[n] = 0.0;
        }
      }
      
      if (total_drop <= 0.0)
        total_drop = 1.0; /* To avoid dividing by zero if no cells are lower */
      
      *total_dir = 0;
      for (n = 0; n < NDIRS; n++) {
        dir[n] = (int) ((drop[n] / total_drop) * 255.0 + 0.5);
        *total_dir += dir[n];
      }
    } else {
      /* For D8 flow directions, water discharges to ONE of its eight neighbors: 
       to one located in the direction of steepest descent. This requires the DEM
       to be pre-filled for D8 routing scheme as any flat areas will confuse the model. */
      
      float maxdrop = -9999; 
      int steepest = -9999;
      
      /*Determine flow direction based on deepest drop */
      for (n = 0; n < NDIRS; n++) {
        dir[n] = 0; /* only the steepest direction is set below */
        /*Make sure flow is inside boundary*/
        if (nelev[n] == (float) OUTSIDEBASIN){
          dir[n] = 0;
          drop[n] = 0;
        }
        else {
          /*Find steepest descent*/
          if ( n == 0 || n == 2 || n == 4 ||n == 6)
            drop[n] = (celev - nelev[n]) / sqrt( dx * dx + dy * dy);
          else 
            drop[n] = (celev - nelev[n]) / dx;
          
          if ((drop[n] < 0.0)&&(drop[n] > -0.001)){
            // printf("Reset minor negative flow slope from %f to 0.0\n", drop[n]);
            drop[n]=0.0;
          }
          
          if (drop[n] >= 0 && drop[n] > maxdrop){
            steepest = n;
            maxdrop = drop[n];
          }
        }
      }
      
      *total_dir = 0;
      if (steepest >= 0){
        dir[steepest] = 1.0;
        *total_dir += dir[steepest];
        
        /* This requires dx = dy */
        if (steepest == 0 || steepest == 2 || steepest == 4 ||steepest == 6)
          total_width = sqrt( dx * dx + dy * dy);
        else 
          total_width = dx; 
      }
      else{
        total_width = dx; 
      }
      *grad = slope * total_width;
    }
  return;
}
/* -------------------------------------------------------------
   ElevationSlopeAspect
   ------------------------------------------------------------- */
void ElevationSlopeAspect(MAPSIZE * Map, TOPOPIX ** TopoMap, int MultiFlowDir)
{
  const char *Routine = "ElevationSlopeAspect";
  int x, y, n, k, xn, yn;
  float neighbor_elev[NNEIGHBORS];
  int steepestdirection;
  float min, dzdx, dzdx_max, slope;

  /* Fill neighbor array */
  
  for (x = 0; x < Map->NX; x++) {
    for (y = 0; y < Map->NY; y++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
	/* Count the number of cells in the basin.  
	   Need this to allocate memory for
	   the new, smaller Elev[] and Coords[][].  */
	Map->NumCells++;
	for (n = 0; n < NNEIGHBORS; n++) {
	  xn = x + xneighbor[n];
	  yn = y + yneighbor[n];	  
	  if (valid_cell(Map, xn, yn)) {
	    neighbor_elev[n] = ((TopoMap[yn][xn].Mask) ? TopoMap[yn][xn].Dem : (float) OUTSIDEBASIN);
	  }
	  else {
	    neighbor_elev[n] = (float) OUTSIDEBASIN;
	  }
	}
	slope_aspect(Map->DX, Map->DY, TopoMap[y][x].Dem, neighbor_elev,
		     &(TopoMap[y][x].Slope), &(TopoMap[y][x].Aspect));

	/* Fill Dirs in TopoMap too */
	flow_fractions(Map->DX, Map->DY, TopoMap[y][x].Slope,
		       TopoMap[y][x].Aspect, TopoMap[y][x].Dem,
		       neighbor_elev, &(TopoMap[y][x].FlowGrad),
		       TopoMap[y][x].Dir, &(TopoMap[y][x].TotalDir), MultiFlowDir);
	
	/* If there is a sink, check again to see if there 
	   is a direction of steepest descent. Does not account 
	   for ties.*/
	if(TopoMap[y][x].TotalDir == 0) {
	  steepestdirection = -99;
	  min = DHSVM_HUGE;	       
	  for (n = 0; n < NDIRS; n++) {
	    xn = x + xdirection[n];
	    yn = y + ydirection[n];	  
	    if (valid_cell(Map, xn, yn)) {
	      if (INBASIN(TopoMap[yn][xn].Mask)) {
			  if(TopoMap[yn][xn].Dem < min) { 
				  min = TopoMap[yn][xn].Dem;
				  steepestdirection = n;
		}
		}
	  }
	  }
	  if(min < TopoMap[y][x].Dem) {
	    TopoMap[y][x].Dir[steepestdirection] = (int)(255.0 + 0.5);
	    TopoMap[y][x].TotalDir = (int)(255.0 + 0.5);
	  }
	  else {
	    /*  Last resort: set the Dir of the cell to the cell that is
		closest in elevation. This should only happen for the 
		basin outlet, unless the Dem wasn't filled. With no mask cell
		among the directions, every Dir stays 0 (TotalDir stays
		nonzero, so the routing fractions are 0). */
	    if (steepestdirection >= 0)
	      TopoMap[y][x].Dir[steepestdirection] = (int)(255.0 + 0.5);
	    TopoMap[y][x].TotalDir = (int)(255.0 + 0.5);
	  }
	  }
  }
  }
  } /* End of loop over y,x */
  /* Create a structure to hold elevations of only those cells
     within the basin and the y,x of those cells */
  if (!(Map->OrderedCells = (ITEM *) calloc(Map->NumCells, sizeof(ITEM))))
    ReportError((char *) Routine, 1);
  k = 0;
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      /* Save the elevation, y, and x in the ITEM structure */
      if (INBASIN(TopoMap[y][x].Mask)) {
        Map->OrderedCells[k].Rank = TopoMap[y][x].Dem;
        Map->OrderedCells[k].y = y;
        Map->OrderedCells[k].x = x;
        k++;
      }
    }
  }
  /* Sort Elev in descending order-- Elev.x and Elev.y hold indices */ 
  quick(Map->OrderedCells, Map->NumCells);
  
  /* Determine unsaturated flow directions and slopes (steepest descent) */
  for (x = 0; x < Map->NX; x++) {
    for (y = 0; y < Map->NY; y++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
        steepestdirection = -99;
        dzdx_max = -1.0 * DHSVM_HUGE;	       
        for (n = 0; n < NDIRS; n++) {
          xn = x + xdirection[n];
          yn = y + ydirection[n];
          if (valid_cell(Map, xn, yn)) {
            if (INBASIN(TopoMap[yn][xn].Mask)) {
              if (n == 0 || n == 2 || n == 4 ||n == 6)
                dzdx = (TopoMap[y][x].Dem - TopoMap[yn][xn].Dem) / (sqrt(2.0) * Map->DX);
              else 
                dzdx = (TopoMap[y][x].Dem - TopoMap[yn][xn].Dem) / Map->DX;
              if(dzdx > dzdx_max) { 
                dzdx_max = dzdx;
                steepestdirection = n;
              }
            }
          }
        }
        TopoMap[y][x].LateralDir = (uchar) steepestdirection;
        
        if (dzdx_max > 0.0) {
          slope = atan(dzdx_max);
          TopoMap[y][x].CosSlope = cos(slope);
          TopoMap[y][x].SinSlope = sin(slope);
        } else {
          TopoMap[y][x].CosSlope = 1.0;
          TopoMap[y][x].SinSlope = 0.0;
        }
        
  }
  }
  } /* End of loop over y,x */
  
  
  /* End of modifications to create ordered cell coordinates */
  return;
}

/* -------------------------------------------------------------
   QuickSort
   ------------------------------------------------------------- */

/**********************************************************************
        this subroutine starts the quick sort
**********************************************************************/

void quick(ITEM *OrderedCells, int count)
{
  qs(OrderedCells,0,count-1);
}

void qs(ITEM *item, int left, int right)
/**********************************************************************
        this is the quick sort subroutine - it returns the values in
        an array from high to low.
**********************************************************************/
{
  register int i,j;
  ITEM x,y;

  i=left;
  j=right;
  x=item[(left+right)/2];

  do {
    while(item[i].Rank<x.Rank && i<right) i++;
    while(x.Rank<item[j].Rank && j>left) j--;

    if (i<=j) {
      y=item[i];
      item[i]=item[j];
      item[j]=y;
      i++;
      j--;
    }
  } while (i<=j);
  if(left<j) qs(item,left,j);
  if(i<right) qs(item,i,right);
}
/* -------------------------------------------------------------
   HeadSlopeAspect
   This computes slope and aspect using the water table elevation. 

   Comment: rewritten to fill the sinks (Ning, 2013)
   ------------------------------------------------------------- */
void HeadSlopeAspect(MAPSIZE * Map, TOPOPIX ** TopoMap, SOILPIX ** SoilMap,
		     float **FlowGrad, unsigned char ***Dir, unsigned int **TotalDir, int MultiFlowDir, int x, int y)
{
  int n;
  float neighbor_elev[NNEIGHBORS];
  float slope, aspect;
  
  for (n = 0; n < NNEIGHBORS; n++) {
    int xn = x + xneighbor[n];
    int yn = y + yneighbor[n];			  
    if (valid_cell(Map, xn, yn))
      neighbor_elev[n] = ((TopoMap[yn][xn].Mask) ? SoilMap[yn][xn].WaterLevel : (float) OUTSIDEBASIN);
    else
      neighbor_elev[n] = (float) OUTSIDEBASIN;
  }
  slope_aspect(Map->DX, Map->DY, SoilMap[y][x].WaterLevel, neighbor_elev, &slope, &aspect);
  
  flow_fractions(Map->DX, Map->DY, slope, aspect, SoilMap[y][x].WaterLevel, neighbor_elev,
                 &(FlowGrad[y][x]), Dir[y][x], &(TotalDir[y][x]), MultiFlowDir);
  return;
}
