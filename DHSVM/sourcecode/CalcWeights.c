
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

   Purpose      : Calculate inverse distance squared interpolation weights to
                  interpolate the meteorological data from the various
                  stations for every individual pixel

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
  double *Distance;		/* Array with distances to all stations */
  double *InvDist2;		/* Array with inverse distance squared */
  double Denominator;		/* Sum of 1/Distance^2 */
  double mindistance;
  int y;			/* Counter for rows */
  int x;			/* Counter for columns */
  int i;			/* Counter for stations */
  int CurrentStation;		/* Station at current location (if any) */
  int closest = 0;
  COORD Loc;			/* Location of current point */

  /* Allocate memory for a 3 dimensional array */

  if (DEBUG)
    Message("Calculating interpolation weights for %d stations\n", NStats);

  if (!((*WeightArray) = (uchar ***)calloc(NY, sizeof(uchar **))))
    ReportError("CalcWeights()", 1);

  for (y = 0; y < NY; y++)
    if (!((*WeightArray)[y] = (uchar **)calloc(NX, sizeof(uchar *))))
      ReportError("CalcWeights()", 1);

  for (y = 0; y < NY; y++)
    for (x = 0; x < NX; x++)
      if (!((*WeightArray)[y][x] = (uchar *)calloc(NStats, sizeof(uchar))))
        ReportError("CalcWeights()", 1);

  /* Allocate memory for the array of the distances to each of the towers,
     and the inverse distance squared */

  if (!(Distance = (double *)calloc(NStats, sizeof(double))))
    ReportError("CalcWeights()", 1);

  if (!(InvDist2 = (double *)calloc(NStats, sizeof(double))))
    ReportError("CalcWeights()", 1);

  /* Calculate the weights for each location that is inside the basin mask */
  /* note stations themselves can be outside the mask */
  /* inverse distance squared scheme */

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

  free(Distance);
  free(InvDist2);
}
