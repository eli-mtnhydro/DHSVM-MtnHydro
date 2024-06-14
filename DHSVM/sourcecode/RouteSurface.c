/*
* SUMMARY:      RouteSurface.c - Route surface flow
* USAGE:        Part of DHSVM
*
* AUTHOR:       Bart Nijssen
* ORG:          University of Washington, Department of Civil Engineering
* E-MAIL:       nijssen@u.washington.edu
* ORIG-DATE:    Apr-96
* DESCRIPTION:  Route surface flow
* DESCRIP-END.
* FUNCTIONS:    RouteSurface()
* Modification: Changes are made to exclude the impervious channel cell (with
a non-zero impervious fraction) from surface routing. In the original
code, some impervious channel cells are routed to themselves causing
overestimated runoff in those cells (Ning, 2013).

* $Id: RouteSurface.c, v3.1.2  2013/3/21   Ning Exp $
*/
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "slopeaspect.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"
/*****************************************************************************
RouteSurface()
If the watertable calculated in WaterTableDepth() was negative, then water is
ponding on the surface.  At the moment no ponding of water is allowed in
DHSVM, and the "excess" water is routed to the outlet one pixel per time step
However, if the pixel contains an impervious fraction, then the surface water
is immediately routed to the nearest downslope pixel that contains a channel.
The net effect is that all pixels that have an impervious area are directly
connected (over the coarse of a single time step) to the channel network, this
assumption is likely to be true for small urban basins, and perhaps even for
large rural basins with some urban development
If Overland Routing = KINEMATIC, then "excess" water is routed to the outlet
using a infinite difference approximation to the kinematic wave solution of
the Saint-Venant equations.
*****************************************************************************/
void RouteSurface(MAPSIZE * Map, TIMESTRUCT * Time, TOPOPIX ** TopoMap,
  SOILPIX ** SoilMap, OPTIONSTRUCT *Options,
  UNITHYDR ** UnitHydrograph, UNITHYDRINFO * HydrographInfo, float *Hydrograph,
  DUMPSTRUCT *Dump, VEGPIX ** VegMap, VEGTABLE * VType,
  SOILTABLE *SType, CHANNEL *ChannelData, float Tair, float Rh)
{
  const char *Routine = "RouteSurface";
  int Lag;			/* Lag time for hydrograph */
  int Step;
  float StreamFlow;
  int TravelTime;
  int WaveLength;
  TIMESTRUCT NextTime;
  TIMESTRUCT VariableTime;
  int i, j, x, y, n, k;         /* Counters */
  float **Runon; /* (m3/s) */
  
  /* Kinematic wave routing */
  float knviscosity;           /* kinematic viscosity JSL */  
  double slope;                /* Slope is Manning's slope;  */
  double alpha;                /* alpha is channel parameter including wetted perimeter,  */
                               /* Manning's n, and Manning's slope; */
  double beta = 3./5.;         /* Beta is 3/5 */
  double outflow;              /* Outflow of water from a pixel during a sub-time step (m3/s) */
                               /* outflow is not entirely true for channel cells */
  float VariableDT;            /* Maximum stable time step (s) */
  
  if (Options->HasNetwork) {
    
    /* Option->Routing = false when routing = conventional */
    if(!Options->Routing) {
      for (y = 0; y < Map->NY; y++) {
        for (x = 0; x < Map->NX; x++) {
          if (INBASIN(TopoMap[y][x].Mask)) {
            SoilMap[y][x].Runoff = SoilMap[y][x].IExcess;
            SoilMap[y][x].IExcess = 0;
            SoilMap[y][x].DetentionIn = 0;
          }
        }
      }
      for (y = 0; y < Map->NY; y++) {
        for (x = 0; x < Map->NX; x++) {
          if (INBASIN(TopoMap[y][x].Mask)) {
            if (!channel_grid_has_channel(ChannelData->stream_map, x, y)) {
              if (VType[VegMap[y][x].Veg - 1].ImpervFrac > 0.0) {
                /* Calculate the outflow from impervious portion of urban cell straight to nearest channel cell */
                SoilMap[TopoMap[y][x].drains_y][TopoMap[y][x].drains_x].IExcess +=
                (1 - VType[VegMap[y][x].Veg - 1].DetentionFrac) *
                VType[VegMap[y][x].Veg - 1].ImpervFrac * SoilMap[y][x].Runoff;
                /* Retained water in detention storage */
                SoilMap[y][x].DetentionIn = VType[VegMap[y][x].Veg - 1].DetentionFrac *
                VType[VegMap[y][x].Veg - 1].ImpervFrac * SoilMap[y][x].Runoff;
                /* Retained water in detention storage routed to channel */
                SoilMap[y][x].DetentionStorage += SoilMap[y][x].DetentionIn;
                SoilMap[y][x].DetentionOut = SoilMap[y][x].DetentionStorage * VType[VegMap[y][x].Veg - 1].DetentionDecay;
                SoilMap[TopoMap[y][x].drains_y][TopoMap[y][x].drains_x].IExcess += SoilMap[y][x].DetentionOut;
                SoilMap[y][x].DetentionStorage -= SoilMap[y][x].DetentionOut;
                if (SoilMap[y][x].DetentionStorage < 0.0)
                  SoilMap[y][x].DetentionStorage = 0.0;
                /* Route the runoff from pervious portion of urban cell to the neighboring cell */
                for (n = 0; n < NDIRS; n++) {
                  int xn = x + xdirection[n];
                  int yn = y + ydirection[n];
                  if (valid_cell(Map, xn, yn)) {
                    SoilMap[yn][xn].IExcess += (1 - VType[VegMap[y][x].Veg - 1].ImpervFrac) * SoilMap[y][x].Runoff
                    * ((float) TopoMap[y][x].Dir[n] / (float) TopoMap[y][x].TotalDir);
                  }
                }
              } else {
                for (n = 0; n < NDIRS; n++) {
                  int xn = x + xdirection[n];
                  int yn = y + ydirection[n];
                  if (valid_cell(Map, xn, yn))
                    SoilMap[yn][xn].IExcess += SoilMap[y][x].Runoff * ((float) TopoMap[y][x].Dir[n] / (float) TopoMap[y][x].TotalDir);
                }
              }
            } else if (channel_grid_has_channel(ChannelData->stream_map, x, y)) {
              SoilMap[y][x].IExcess += SoilMap[y][x].Runoff;
            }
          }
        }
      }
      /* End if Options->routing = conventional */
    } else {
      /* Begin code for kinematic wave routing */
      
      NextTime = *Time;
      VariableTime = *Time;
      
      /* Holds the value of the next DHSVM time step. */ 
      IncreaseTime(&NextTime);
      
      /* Use the Courant condition to find the maximum stable time step (in seconds). */
      /* Must be an even increment of Dt. */
      VariableDT = FindDT(SoilMap, Map, Time, TopoMap, SType);
      
      /* Allocate memory for Runon Matrix */
      if ((Runon = (float **) calloc(Map->NY, sizeof(float *))) == NULL) 
        ReportError((char *) Routine, 1);
      for (y = 0; y < Map->NY; y++) {
        if ((Runon[y] = (float *) calloc(Map->NX, sizeof(float))) == NULL)
          ReportError((char *) Routine, 1);
      }
      
      /* Reset surface runoff and initialize runon */
      /* Initialize Runon variables */
      for (y = 0; y < Map->NY; y++) {
        for (x = 0; x < Map->NX; x++) {
          if (INBASIN(TopoMap[y][x].Mask)) {
            SoilMap[y][x].Runoff = 0.;
            Runon[y][x] = 0.;
            SoilMap[y][x].DetentionIn = 0;
          }
        }
      }
      
      /* Estimate kinematic viscosity through interpolation JSL */
      knviscosity = viscosity(Tair, Rh);
      /* Converting units to m2/sec */
      knviscosity /= 1000. * 1000.;
      
      /* Must loop through surface routing multiple times within one DHSVM  model time step */
      while (Before(&(VariableTime.Current), &(NextTime.Current))) {
        /* Loop thru all of the cells in descending order of elevation */
        for (k = (Map->NumCells - 1); k > -1;  k--) {
          y = Map->OrderedCells[k].y;
          x = Map->OrderedCells[k].x;
          
          /* Only compute kinematic routing parameters for cells with non-zero runoff */
          if (SoilMap[y][x].IExcess > 0.0 || Runon[y][x] > 0.0){
            
            // printf("%d k IExcess %6.6f Runon  %6.6f startRunoff %6.6f\n",k,SoilMap[y][x].IExcess,Runon[y][x],SoilMap[y][x].startRunoff);
            
            outflow = SoilMap[y][x].startRunoff;
            
            slope = TopoMap[y][x].Slope;
            if (slope == 0.0) slope=0.0001;
            assert(slope > 0.0);
            
            alpha = pow(SType[SoilMap[y][x].Soil - 1].Manning * pow((double) Map->DX, 2./3.) / sqrt(slope), beta);
            
            /* Calculate discharge (m3/s) from the grid cell using an explicit */
            /* Finite difference solution of the linear kinematic wave */
            if (Runon[y][x] > 0.0001 || outflow > 0.0001)
              outflow = ((VariableDT / Map->DX) * Runon[y][x] + alpha * beta * outflow *
                pow((outflow + Runon[y][x]) / 2.0, beta - 1.) +
                SoilMap[y][x].IExcess * Map->DX * VariableDT / Time->Dt) / ((VariableDT / Map->DX) + alpha * beta *
                pow((outflow + Runon[y][x]) / 2.0, beta - 1.));
            else if (SoilMap[y][x].IExcess > 0.0)
              outflow = SoilMap[y][x].IExcess * (Map->DX * Map->DY) / Time->Dt; 
            else
              outflow = 0.0; 
            if(outflow < 0.0)
              outflow = 0.0;
            
            /* Make sure calculated outflow doesn't exceed available water, and update surface water storage */
            if (outflow > (SoilMap[y][x].IExcess * (Map->DX * Map->DY) / Time->Dt + Runon[y][x])) 
              outflow = SoilMap[y][x].IExcess * (Map->DX * Map->DY) / Time->Dt + Runon[y][x];
            
            /* Save sub-timestep runoff for q(i)(t-1) and q(i-1)(t-1) of next time step. */
            SoilMap[y][x].startRunoff = outflow;
            
            /* Calculate total runoff in m per DHSVM timestep */
            /* This is serving the purpose of holding this value for Cournat condition calculation; */
            /* It does not truly represent the Runoff, because if there is a channel there is no outflow. */
            /* Instead, IExcess is updated based on Runon in the same manner as the original DHSVM */
            SoilMap[y][x].Runoff += outflow * VariableDT / (Map->DX * Map->DY); 
            
            if (channel_grid_has_channel(ChannelData->stream_map, x, y)
                  || (channel_grid_has_channel(ChannelData->road_map, x, y)
                        && !channel_grid_has_sink(ChannelData->road_map, x, y))) {
                        // SoilMap[y][x].IExcess += (Runon[y][x] - outflow) * VariableDT / (Map->DX * Map->DY);
                        outflow = 0.0;
            }
            
            SoilMap[y][x].IExcess += (Runon[y][x] - outflow) * VariableDT / (Map->DX * Map->DY);
            
            /* Redistribute surface water to downslope pixels */
            if(outflow > 0.) {
              for (n = 0; n < NDIRS; n++) {
                int xn = x + xdirection[n];
                int yn = y + ydirection[n];
                if (valid_cell(Map, xn, yn))
                  Runon[yn][xn] += outflow * ((float) TopoMap[y][x].Dir[n] / (float) TopoMap[y][x].TotalDir);
              } /* End loop thru possible flow directions */
            }
            
            /* Initialize runon for next timestep. */
            Runon[y][x] = 0.0;
            
          }
          
        } /* End loop thru ordered basin cells */
        IncreaseVariableTime(&VariableTime, VariableDT, &NextTime);
      } /* End of internal time step loop. */
      
      for (y = 0; y < Map->NY; y++) {
        free(Runon[y]);
      }
      free(Runon);
      
      /* Handle detention storage and impervious routing */
      for (y = 0; y < Map->NY; y++) {
        for (x = 0; x < Map->NX; x++) {
          if (INBASIN(TopoMap[y][x].Mask) &&
              !channel_grid_has_channel(ChannelData->stream_map, x, y) &&
              VType[VegMap[y][x].Veg - 1].ImpervFrac > 0.0) {
            /* Calculate the outflow from impervious portion of urban cell straight to nearest channel cell */
            SoilMap[TopoMap[y][x].drains_y][TopoMap[y][x].drains_x].IExcess +=
            (1 - VType[VegMap[y][x].Veg - 1].DetentionFrac) *
            VType[VegMap[y][x].Veg - 1].ImpervFrac * SoilMap[y][x].IExcess;
            /* Retained water in detention storage */
            SoilMap[y][x].DetentionIn = VType[VegMap[y][x].Veg - 1].DetentionFrac *
            VType[VegMap[y][x].Veg - 1].ImpervFrac * SoilMap[y][x].IExcess;
            /* Retained water in detention storage routed to channel */
            SoilMap[y][x].DetentionStorage += SoilMap[y][x].DetentionIn;
            SoilMap[y][x].DetentionOut = SoilMap[y][x].DetentionStorage * VType[VegMap[y][x].Veg - 1].DetentionDecay;
            SoilMap[TopoMap[y][x].drains_y][TopoMap[y][x].drains_x].IExcess += SoilMap[y][x].DetentionOut;
            SoilMap[y][x].DetentionStorage -= SoilMap[y][x].DetentionOut;
            if (SoilMap[y][x].DetentionStorage < 0.0)
              SoilMap[y][x].DetentionStorage = 0.0;
            /* Pervious IExcess remains in current cell */
            /* Because routing is already accomplished kinematically */
            SoilMap[y][x].IExcess = (1 - VType[VegMap[y][x].Veg - 1].ImpervFrac) * SoilMap[y][x].IExcess;
          }
        }
      }
      
    } /* End of code added for kinematic wave routing. */
  
  } else {
    /* No network, so use unit hydrograph method */
    /* MAKE SURE THIS WORKS WITH A TIMESTEP IN SECONDS */
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
        if (INBASIN(TopoMap[y][x].Mask)) {
          TravelTime = (int)TopoMap[y][x].Travel;
          if (TravelTime != 0) {
            WaveLength = HydrographInfo->WaveLength[TravelTime - 1];
            for (Step = 0; Step < WaveLength; Step++) {
              Lag = UnitHydrograph[TravelTime - 1][Step].TimeStep;
              Hydrograph[Lag] += SoilMap[y][x].Runoff * UnitHydrograph[TravelTime - 1][Step].Fraction;

            }
            SoilMap[y][x].Runoff = 0.0;
          }
        }
      }
    }

    StreamFlow = 0.0;
    for (i = 0; i < Time->Dt; i++)
      StreamFlow += (Hydrograph[i] * Map->DX * Map->DY) / Time->Dt;

    /* Advance Hydrograph */
    for (i = 0; i < Time->Dt; i++) {
      for (j = 0; j < HydrographInfo->TotalWaveLength - 1; j++) {
        Hydrograph[j] = Hydrograph[j + 1];
      }

    }

    /* Set the last elements of the hydrograph to zero */
    for (i = 0; i < Time->Dt; i++)
      Hydrograph[HydrographInfo->TotalWaveLength - (i + 1)] = 0.0;

    PrintDate(&(Time->Current), Dump->Stream.FilePtr);
    fprintf(Dump->Stream.FilePtr, " %g\n", StreamFlow);
  }
}

/*****************************************************************************
 FindDT()
 Find the variable time step that will satisfy the courant condition for stability 
 in overland flow routing.
 *****************************************************************************/
float FindDT(SOILPIX **SoilMap, MAPSIZE *Map, TIMESTRUCT *Time, 
             TOPOPIX **TopoMap, SOILTABLE *SType)
{
  int x, y;
  /* JSL: slope is manning's slope; alpha is channel parameter including wetted perimeter, 
   manning's n, and manning's slope.  Beta is 3/5 */
  float slope;
  double alpha;
  double beta = 3./5.;
  double Ck;
  float DT, minDT;
  float numinc;
  minDT = 36000.;
  
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
        if (SoilMap[y][x].Runoff > 0.0) {
          
          slope = TopoMap[y][x].Slope;
          
          if (slope <= 0) slope = 0.0001;
          alpha = pow((double) SType[SoilMap[y][x].Soil - 1].Manning *
            pow((double) Map->DX, (double) (2./3.)) / sqrt(slope), (double) beta);
          
          /* Calculate flow velocity from discharge  using Manning's equation */
          Ck = 1. / (alpha * beta * pow((double) SoilMap[y][x].Runoff, beta - 1.));
          
          if((Map->DX / Ck) < minDT)
            minDT = Map->DX / Ck;
        }
      }
    }
  }
  /* Find the time step that divides evenly into Time->DT */
  numinc = (float) ceil((double) Time->Dt/minDT);
  DT = ((float) Time->Dt) / numinc;
  
  if(DT > Time->Dt)
    DT = (float) Time->Dt;
  
  return DT;
}

