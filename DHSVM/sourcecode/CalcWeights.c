
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "constants.h"
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"

 /*****************************************************************************
   Function name: CalcWeights()

   Purpose      : Calculate interpolation weights to interpolate the
                  meteorological data from the various stations for every
                  individual pixel

   Required     :
     METLOCATION *Station - Location of meteorological stations
     int NStats           - Number of meteorological stations
     int NX               - Number of pixels in East - West direction
     int NY               - Number of pixels in North - South direction
     uchar ** BasinMask   - BasinMask
     uchar ***WeightArray - 3D array with interpolation weights

   Modifies     :
     The values stored at the addresses pointed to by WeightArray (i.e. it
     calculates the weights and stores them)

 *****************************************************************************/
void CalcWeights(METLOCATION * Station, int NStats, int NX, int NY,
  uchar ** BasinMask, uchar **** WeightArray,
  OPTIONSTRUCT * Options)
{
  double *Weights;		/* Array with weights for all stations */
  double *Distance;		/* Array with distances to all stations */
  double *InvDist2;		/* Array with inverse distance squared */
  double Denominator;		/* Sum of 1/Distance^2 */
  double mindistance;
  double tempdistance;
  double cr, crt;
  float uniformweight;
  int totalweight;
  int y;			/* Counter for rows */
  int x;			/* Counter for columns */
  int i, j;			/* Counter for stations */
  int CurrentStation;		/* Station at current location (if any) */
  int *stationid;		/* index array for sorted list of station distances */
  int *stat;
  int tempid;
  int closest = 0;
  int crstat;
  COORD Loc;			/* Location of current point */

  if (NStats > (int) MAXUCHAR && Options->Interpolation == UNIFORM) {
    printf("\nWARNING:\nCannot use more than 255 stations with UNIFORM interpolation.\n");
    printf("Setting interpolation method to INVDIST.\n\n");
    Options->Interpolation = INVDIST;
  }
  
  
  /* Allocate memory for a 3 dimensional array */

  if (DEBUG)
    printf("Calculating interpolation weights for %d stations\n", NStats);

  if (!((*WeightArray) = (uchar ***)calloc(NY, sizeof(uchar **))))
    ReportError("CalcWeights()", 1);

  for (y = 0; y < NY; y++)
    if (!((*WeightArray)[y] = (uchar **)calloc(NX, sizeof(uchar *))))
      ReportError("CalcWeights()", 1);

  for (y = 0; y < NY; y++)
    for (x = 0; x < NX; x++)
      if (!((*WeightArray)[y][x] = (uchar *)calloc(NStats, sizeof(uchar))))
        ReportError("CalcWeights()", 1);

  /* Allocate memory for the array that will contain weights, and the array for
     the distances to each of the towers, and the inverse distance squared */

  if (!(Weights = (double *) calloc(NStats, sizeof(double))))
    ReportError("CalcWeights()", 1);

  if (!(Distance = (double *)calloc(NStats, sizeof(double))))
    ReportError("CalcWeights()", 1);

  if (!(InvDist2 = (double *)calloc(NStats, sizeof(double))))
    ReportError("CalcWeights()", 1);

  if (!(stationid = (int *)calloc(NStats, sizeof(int))))
    ReportError("CalcWeights()", 1);

  if (!(stat = (int *)calloc(NStats + 1, sizeof(int))))
    ReportError("CalcWeights()", 1);

  /* Calculate the weights for each location that is inside the basin mask */
  /* note stations themselves can be outside the mask */
  /* this first scheme is an inverse distance squared scheme */

  if (Options->Interpolation == INVDIST) {
    for (y = 0; y < NY; y++) {
      Loc.N = y;
      for (x = 0; x < NX; x++) {
        Loc.E = x;
        if (INBASIN(BasinMask[y][x])) {
          if (IsStationLocation(&Loc, NStats, Station, &CurrentStation)) {
            for (i = 0; i < NStats; i++) {
              if (i == CurrentStation)
                (*WeightArray)[y][x][i] = MAXUCHAR;
              else
                (*WeightArray)[y][x][i] = 0;
            }
          }
          else {
            for (i = 0, Denominator = 0; i < NStats; i++) {
              /* N.B. these distances are in gridcell-units */
              Distance[i] = CalcDistance(&(Station[i].Loc), &Loc);
              if (Distance[i] < Options->MaxInterpDist)
                InvDist2[i] = 1 / (Distance[i] * Distance[i]);
              else
                InvDist2[i] = 0.0;
              Denominator += InvDist2[i];
            }
            if (Denominator > 0.0) {
              for (i = 0; i < NStats; i++) {
                (*WeightArray)[y][x][i] =
                  (uchar)Round(InvDist2[i] / Denominator * MAXUCHAR);
              }
            } else {
              /* Use nearest if none of the stations are within MaxInterpDist */
              mindistance = DHSVM_HUGE;
              for (i = 0; i < NStats; i++) {
                if (Distance[i] < mindistance) {
                  mindistance = Distance[i];
                  closest = i;
                }
              }
              for (i = 0; i < NStats; i++) {
                if (i == closest)
                  (*WeightArray)[y][x][i] = MAXUCHAR;
                else
                  (*WeightArray)[y][x][i] = 0;
              }
            } /* End handling nearest */
          }
        }
        else {
          for (i = 0; i < NStats; i++)
            (*WeightArray)[y][x][i] = 0;
        }
      }
    }
  }
  if (Options->Interpolation == NEAREST) {

    /* this next scheme is a nearest station */
    printf("Number of stations is %d \n", NStats);

    for (y = 0; y < NY; y++) {
      Loc.N = y;
      for (x = 0; x < NX; x++) {
        Loc.E = x;
        if (INBASIN(BasinMask[y][x])) {	/*we are inside the basin mask */
          /* Find the distance to nearest station */
          mindistance = DHSVM_HUGE;
          for (i = 0; i < NStats; i++) {
            Distance[i] = CalcDistance(&(Station[i].Loc), &Loc);
            if (Distance[i] < mindistance) {
              mindistance = Distance[i];
              closest = i;
            }
          }
          /* Got closest station */
          for (i = 0; i < NStats; i++) {
            if (i == closest)
              (*WeightArray)[y][x][i] = MAXUCHAR;
            else
              (*WeightArray)[y][x][i] = 0;
          }

        }			/* done in basin mask */
        else {
          for (i = 0; i < NStats; i++)
            (*WeightArray)[y][x][i] = 0;
        }
      }
    }
  }

  if (Options->Interpolation == VARCRESS) {

    /* this next scheme is a variable radius cressman */
    /* find the distance to the nearest station */
    /* make a decision based on the maximum allowable radius, cr */
    /* and the distance to the closest station */
    /* while limiting the number of interpolation stations to three */
    cr = (double)Options->CressRadius;
    if (cr < 2)
      ReportError("CalcWeights.c", 42);
    crstat = Options->CressStations;
    if (crstat < 2)
      ReportError("CalcWeights.c", 42);
    for (y = 0; y < NY; y++) {
      Loc.N = y;
      for (x = 0; x < NX; x++) {
        Loc.E = x;
        if (INBASIN(BasinMask[y][x])) {	/*we are inside the basin mask */
          /* find the distance to nearest station */
          for (i = 0; i < NStats; i++) {
            Distance[i] = CalcDistance(&(Station[i].Loc), &Loc);
            stationid[i] = i;
          }
          /* got distances for each station */
          /* now sort the list by distance */
          for (i = 0; i < NStats; i++) {
            for (j = 0; j < NStats; j++) {
              if (Distance[j] > Distance[i]) {
                tempdistance = Distance[i];
                tempid = stationid[i];
                Distance[i] = Distance[j];
                stationid[i] = stationid[j];
                Distance[j] = tempdistance;
                stationid[j] = tempid;
              }
            }
          }

          crt = Distance[0] * 2.0;
          if (crt < 1.0)
            crt = 1.0;
          for (i = 0, Denominator = 0; i < NStats; i++) {
            if (i < crstat && Distance[i] < crt) {
              InvDist2[i] =
                (crt * crt - Distance[i] * Distance[i]) /
                (crt * crt + Distance[i] * Distance[i]);
              Denominator += InvDist2[i];
            }
            else
              InvDist2[i] = 0.0;
          }

          for (i = 0; i < NStats; i++)
            (*WeightArray)[y][x][stationid[i]] =
            (uchar)Round(InvDist2[i] / Denominator * MAXUCHAR);

          /*at this point all weights have been assigned to one or more stations */

        }
      }
    }
  }
  if (Options->Interpolation == UNIFORM) {
    
    /* Simple uniform "interpolation," actually just the average of all stations */
    printf("Number of stations is %d, used as simple average for whole domain\n", NStats);
    printf("Note that lapse rates, precip multiplier, etc. will still apply\n");
    
    uniformweight = ((float) MAXUCHAR) / ((float) NStats);
    
    for (y = 0; y < NY; y++) {
      for (x = 0; x < NX; x++) {
        if (INBASIN(BasinMask[y][x])) {	/*we are inside the basin mask */
          
          for (i = 0; i < NStats; i++) {
            (*WeightArray)[y][x][i] = (uchar) uniformweight;
          }
          
        }			/* done in basin mask */
        else {
          for (i = 0; i < NStats; i++)
            (*WeightArray)[y][x][i] = 0;
        }
      }
    }
  }
  
  /*check that all weights add up to MAXUCHAR */
  /* and output some stats on the interpolation field */

  for (i = 0; i <= NStats; i++)
    stat[i] = 0;
  for (i = 0; i < NStats; i++)
    stationid[i] = 0;

  /*
  printf("\nChecking interpolation weights\n");
  printf("Sum should be 255 for all pixels \n");
  printf("Some error is expected due to roundoff \n");
  printf("Errors greater than +/- 2 Percent are: \n");
  */

  for (y = 0; y < NY; y++) {
    for (x = 0; x < NX; x++) {
      if (INBASIN(BasinMask[y][x])) {	/*we are inside the basin mask */
        tempid = 0;
        totalweight = 0;

        for (i = 0; i < NStats; i++) {
          totalweight += (int)(*WeightArray)[y][x][i];
          if ((*WeightArray)[y][x][i] > 0) {
            tempid += 1;
            stationid[i] = 1;
          }
        }

        /*if (totalweight < 250 || totalweight > 260)
          printf("error in interpolation weight at pixel y %d x %d : %d \n", y,
            x, totalweight); */
        stat[tempid] += 1;
      }

    }
  }
  
  free(Weights);
  free(Distance);
  free(InvDist2);
  free(stationid);
  free(stat);
}
