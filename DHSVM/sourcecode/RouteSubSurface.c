/*
 * SUMMARY:      RouteSubSurface.c - Route subsurface flow
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bill Perkins
 * ORG:          PNNL
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Route subsurface flow
 * DESCRIP-END.
 * FUNCTIONS:    RouteSubSurface()
 * COMMENTS:
 * $Id: RouteSubSurface.c,v3.1.2 2013/08/18 ning Exp $     
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"
#include "soilmoisture.h"
#include "slopeaspect.h"
#include "DHSVMChannel.h"

/*****************************************************************************
  RouteSubSurface()

  Sources: 
  Wigmosta, M. S., L. W. Vail, and D. P. Lettenmaier, A distributed 
      hydrology-vegetation model for complex terrain, Water Resour. Res.,
      30(6), 1665-1679, 1994.

  Quinn, P., K. Beven, P. Chevallier, and O. Planchon, The prediction of 
      hillslope flow paths for distributed hydrological modelling using 
      digital terrain models, Hydrological Processes, 5, 59-79, 1991.

  This routine follows Wigmosta et al. [1994] in calculating the subsurface
  flow.  The local gradient is based on the local hydraulic head, consisting 
  of the height of the pixel surface minus the depth of the water table 
  below the water surface.  This has the disadvantage that the local gradients
  have to be recalculated for each pixel for each timestep.  In Wigmosta et
  al. [1994] the local gradient is taken to be equal to the slope of the land
  surface, which is a reasonable assunption for mountainous areas.  For the 
  flat boreal forest landscape it is probably better to use the slope
  of the water surface.

  Set the gradient with pixels that are outside tha basin to zero.  This 
  ensures that there is no flux of water across the basin boundary.  In the 
  current implementation water can only leave the basin as surface flow.  
  This may not be entirely realistic, and should be analyzed further.  
  One consequence of this could be that the soil in the basin is more 
  saturated than it would be if subsurface flow out of the basin would
  be allowed.

  The surrounding grid cells are numbered in the following way

                |-----| DX

          0-----1-----2  ---
	  |\    |    /|   |
          | \   |   / |   |
          |  \  |  /  |   | DY
          |   \ | /   |   |
          |    \|/    |   |
          7-----*-----3  ---
          |    /|\    |
          |   / | \   |
          |  /  |  \  |
          | /   |   \ |
          |/    |    \|
          6-----5-----4

  For the current implementation it is assumed that the resolution is the 
  same in both the X and the Y direction.  If this is not the case an error
  message is generated and the program exits.  The reason is that the 
  formulation for the flow width in the diagonal direction changes if the
  grid is not square.  The method for defining the flow widths in the case
  of square grids is taken from Quinn et al [1991]

  Update Jan 2004 COD
  When Gradient = WATERTABLE, the watertable was used to route the
  surface water. This was because of the common use of TopoMap.Dir and 
  TopoMap.TotalDir. These are now for surface routing (always) and subsurface 
  routing (when Gradient = TOPOGRAPHY). Subsurface routing directions 
  and FlowGrad (SubDir, SubTotalDir, SubFlowGrad) for Gradient = WATERTABLE 
  are now determined locally here (in RouteSubsurface.c.)
 
  Update 2023-2024 Eli Boardman
  Initially the code had an if statement that prevented all
  lateral subsurface flow beneath stream channels. This version is considerably
  updated to enable hyporheic exchange and subsurface flow below channels.
 
  WORK IN PROGRESS
*****************************************************************************/
void RouteSubSurface(int Dt, MAPSIZE *Map, TOPOPIX **TopoMap,
		     VEGTABLE *VType, VEGPIX **VegMap,
		     ROADSTRUCT **Network, SOILTABLE *SType,
		     SOILPIX **SoilMap, CHANNEL *ChannelData,
		     TIMESTRUCT *Time, OPTIONSTRUCT *Options, 
		     char *DumpPath, int MaxStreamID, SNOWPIX **SnowMap)
{
  const char *Routine = "RouteSubSurface";
  int x, nx;			/* counter */
  int y, ny;			/* counter */
  int i, j;	        /* counters */
  int flag;
  float BankHeight;
  float ChannelWaterLevel;
  float EffThickness;
  float SoilDeficit;
  float *Adjust;
  float fract_used;
  float DeltaWaterLevel;
  float Depth;
  float OutFlow;
  float water_out_road;
  float water_out_stream;
  float water_in_stream;
  float Transmissivity;
  float TotalAvailableWater = 0.0; /* Including water that flows laterally and to channels/roads */
  float AvailableWater;
  float AdjTableDepth, AdjTableDepthK, AdjWaterLevel, AdjWaterLevelK;
  float PotentialSatFlow, ActualSatFlow, LayerContribWater, LayerStorageCap, DeltaTableDepth;
  float LayerContribWaterK, LayerStorageK, LayerStorageCapK, DeltaTableDepthK, LayerUseFrac;
  int k, q;
  float **SubFlowGrad;	        /* Magnitude of subsurface flow gradient slope * width */
  unsigned char ***SubDir;      /* Fraction of flux moving in each direction*/ 
  unsigned int **SubTotalDir;	/* Sum of Dir array */
  ITEM kOrdered[NDIRS];

  int count, totalcount;
  float mgrid, sat;
  char buffer[32];
  char satoutfile[100];         /* Character arrays to hold file name. */ 
  FILE *fs;                     /* File pointer. */
  /*****************************************************************************
   Allocate memory 
  ****************************************************************************/
  
  if (!(SubFlowGrad = (float **)calloc(Map->NY, sizeof(float *))))
    ReportError((char *) Routine, 1);
  for(i=0; i<Map->NY; i++) {
    if (!(SubFlowGrad[i] = (float *)calloc(Map->NX, sizeof(float))))
      ReportError((char *) Routine, 1);
  }
  
  if (!((SubDir) = (unsigned char ***) calloc(Map->NY, sizeof(unsigned char **))))
    ReportError((char *) Routine, 1);
  for (i=0; i<Map->NY; i++) {
    if (!((SubDir)[i] = (unsigned char **) calloc(Map->NX, sizeof(unsigned char*))))
	  ReportError((char *) Routine, 1);
	for (j=0; j<Map->NX; j++) {
      if (!(SubDir[i][j] = (unsigned char *)calloc(NDIRS, sizeof(unsigned char ))))
		ReportError((char *) Routine, 1);
    }
  }

  if (!(SubTotalDir = (unsigned int **)calloc(Map->NY, sizeof(unsigned int *))))
    ReportError((char *) Routine, 1);
  for (i=0; i<Map->NY; i++) {
    if (!(SubTotalDir[i] = (unsigned int *)calloc(Map->NX, sizeof(unsigned int))))
      ReportError((char *) Routine, 1);
  }
  
  /* Reset the saturated subsurface flow to zero 
     and update water table elevation */
  for (q = (Map->NumCells - 1); q > -1;  q--) {
    y = Map->OrderedCells[q].y;
    x = Map->OrderedCells[q].x;
    
    SoilMap[y][x].SatFlow = 0.0;
    SoilMap[y][x].RoadInt = 0.0;
    
    SoilMap[y][x].WaterLevel = TopoMap[y][x].Dem - SoilMap[y][x].TableDepth;
    
    if (Options->FlowGradient == WATERTABLE) {
      flag = 0;
      for (k = 0; k < NDIRS; k++) {
        nx = xdirection[k] + x;
        ny = ydirection[k] + y;
        if (SoilMap[y][x].WaterLevel > SoilMap[ny][nx].WaterLevel) {
          flag = 1;
          /* Use time-averaged WaterLevel if the change in local water level
           is greater than the hydraulic or topographic drop to down-gradient neighbors
           to mitigate water table oscillation in diffusion-dominated flat areas. */
          DeltaWaterLevel = ABSVAL(SoilMap[y][x].WaterLevel - SoilMap[y][x].WaterLevelLast);
          if (DeltaWaterLevel > (TopoMap[y][x].Dem - TopoMap[ny][nx].Dem) ||
              DeltaWaterLevel > (SoilMap[y][x].WaterLevel - SoilMap[ny][nx].WaterLevel)) {
            SoilMap[y][x].WaterLevel = (SoilMap[y][x].WaterLevel + SoilMap[y][x].WaterLevelLast) / 2.0;
            break;
          }
        }
      }
      /* Also use time-averaged WaterLevel if current cell is a sink */
      if (!flag)
        SoilMap[y][x].WaterLevel = (SoilMap[y][x].WaterLevel + SoilMap[y][x].WaterLevelLast) / 2.0;
      
      SoilMap[y][x].WaterLevelLast = SoilMap[y][x].WaterLevel;
    } /* End of water table adjustments */
  }
  
  /* Calculate flow directions and gradient */
  if (Options->FlowGradient == WATERTABLE) {
    for (q = (Map->NumCells - 1); q > -1;  q--)
      HeadSlopeAspect(Map, TopoMap, SoilMap, SubFlowGrad, SubDir, SubTotalDir, Options->MultiFlowDir,
                      Map->OrderedCells[q].x, Map->OrderedCells[q].y);
  }
  
  /* Next sweep through all the grid cells (by descending elevation),
     calculate the amount of flow in each direction,
     and divide the flow over the surrounding pixels */
  for (q = (Map->NumCells - 1); q > -1;  q--) {
    y = Map->OrderedCells[q].y;
    x = Map->OrderedCells[q].x;
    
    AdjTableDepth = TopoMap[y][x].Dem - SoilMap[y][x].WaterLevel;
    AdjWaterLevel = SoilMap[y][x].WaterLevel;
    BankHeight = (Network[y][x].BankHeight > SoilMap[y][x].Depth) ?
                 SoilMap[y][x].Depth : Network[y][x].BankHeight;
    Adjust = Network[y][x].Adjust;
    fract_used = 0.0f;
    water_out_road = 0.0;
    water_out_stream = 0.0;
    water_in_stream = 0.0;
    
    if (Options->FlowGradient == TOPOGRAPHY) {
      SubTotalDir[y][x] = TopoMap[y][x].TotalDir;
      SubFlowGrad[y][x] = TopoMap[y][x].FlowGrad;
      for (k = 0; k < NDIRS; k++) 
        SubDir[y][x][k] = TopoMap[y][x].Dir[k];
    }
    
    for (k = 0; k < NDIRS; k++) {
      fract_used += (float) SubDir[y][x][k];
    }
    if (SubTotalDir[y][x] > 0)
      fract_used /= (float) SubTotalDir[y][x];
    else
      fract_used = 0.;
    
    /* Only bother calculating subsurface flow if water table is above bedrock */
    if (AdjTableDepth < SoilMap[y][x].Depth) {
      Depth = ((AdjTableDepth > BankHeight) ?
                 AdjTableDepth : BankHeight);
      
      Transmissivity = CalcTransmissivity(SoilMap[y][x].Depth, Depth,
                                          SoilMap[y][x].KsLat,
                                          SoilMap[y][x].KsLatExp,
                                          SType[SoilMap[y][x].Soil - 1].DepthThresh);
      
      OutFlow = (Transmissivity * fract_used * SubFlowGrad[y][x] * Dt) / (Map->DX * Map->DY);
      
      /* Determine TOTAL amount of water available for redistribution */
      TotalAvailableWater =
      CalcAvailableWater(VType[VegMap[y][x].Veg - 1].NSoilLayers,
                         SoilMap[y][x].Depth, VType[VegMap[y][x].Veg - 1].RootDepth,
                         SoilMap[y][x].Porosity, SoilMap[y][x].FCap, SoilMap[y][x].Moist,
                         AdjTableDepth, Adjust);
    }
    else
      OutFlow = 0.0f;
    
    /* Compute road interception if water table is above road cut */
    if (AdjTableDepth < BankHeight &&
    channel_grid_has_channel(ChannelData->road_map, x, y)) {
      if (SubTotalDir[y][x] > 0)
        fract_used = ((float) Network[y][x].fraction /
          (float)SubTotalDir[y][x]);
      else
        fract_used = 0.;
      Transmissivity =
        CalcTransmissivity(BankHeight, AdjTableDepth,
                           SoilMap[y][x].KsLat,
                           SoilMap[y][x].KsLatExp,
                           SType[SoilMap[y][x].Soil - 1].DepthThresh);
      
      water_out_road = (Transmissivity * fract_used * SubFlowGrad[y][x] * Dt) / (Map->DX * Map->DY);
      
      AvailableWater =
        CalcAvailableWater(VType[VegMap[y][x].Veg - 1].NSoilLayers,
                           BankHeight, VType[VegMap[y][x].Veg - 1].RootDepth,
                           SoilMap[y][x].Porosity,
                           SoilMap[y][x].FCap, SoilMap[y][x].Moist,
                           AdjTableDepth, Adjust);
      
      water_out_road = (water_out_road > AvailableWater) ? AvailableWater : water_out_road;
      
      /* Increase lateral inflow to road channel */
      SoilMap[y][x].RoadInt = water_out_road;
      channel_grid_inc_inflow(ChannelData->road_map, x, y,
                              water_out_road * Map->DX * Map->DY);
    }
    
    /* Compute stream lateral inflow/outflow if water table is above channel cut */
    if (AdjTableDepth < BankHeight &&
    channel_grid_has_channel(ChannelData->stream_map, x, y)) {
      
      /* Also consider depth of water stored in channel */
      ChannelWaterLevel = BankHeight - channel_grid_cell_water_depth(ChannelData->stream_map, x, y);
      if (ChannelWaterLevel < 0.0)
        ChannelWaterLevel = 0.0;
      if (AdjTableDepth < ChannelWaterLevel) {
        /* Water table is above water level in channel,
           channel gains water laterally */
        Transmissivity =
          CalcTransmissivity(ChannelWaterLevel, AdjTableDepth,
                             SoilMap[y][x].KsLat,
                             SoilMap[y][x].KsLatExp,
                             SType[SoilMap[y][x].Soil - 1].DepthThresh);
        
        AvailableWater = 
          CalcAvailableWater(VType[VegMap[y][x].Veg - 1].NSoilLayers,
                             ChannelWaterLevel, VType[VegMap[y][x].Veg - 1].RootDepth,
                             SoilMap[y][x].Porosity,
                             SoilMap[y][x].FCap, SoilMap[y][x].Moist,
                             AdjTableDepth, Adjust);
        AvailableWater = AvailableWater - water_out_road;
        
        /* New method: contribute lateral inflow to each channel segment individually */
        water_out_stream = channel_grid_calc_satflow(ChannelData->stream_map, x, y,
                                                     AdjTableDepth,
                                                     Transmissivity, AvailableWater,
                                                     Map->DX, Map->DY, Dt);
        water_out_stream /= (Map->DX * Map->DY);
        
        SoilMap[y][x].ChannelInt += water_out_stream;
      } else {
        /* Water table is above bottom of channel but below channel water level,
           channel loses water laterally */
        Transmissivity =
          CalcTransmissivity(AdjTableDepth, ChannelWaterLevel,
                             SoilMap[y][x].KsLat,
                             SoilMap[y][x].KsLatExp,
                             SType[SoilMap[y][x].Soil - 1].DepthThresh);

        /* Find capacity (porosity - moist) of unsaturated soil below channel water level */
        SoilDeficit = 0.0;
        Depth = 0.0;
        for (i = 0; i < SType[SoilMap[y][x].Soil - 1].NLayers && Depth < AdjTableDepth; i++) {
          if (VType[VegMap[y][x].Veg - 1].RootDepth[i] < (SoilMap[y][x].Depth - Depth))
            Depth += VType[VegMap[y][x].Veg - 1].RootDepth[i];
          else
            Depth = SoilMap[y][x].Depth;

          if (Depth > ChannelWaterLevel) {
            if (Depth < AdjTableDepth)
              EffThickness = ((Depth - ChannelWaterLevel) < VType[VegMap[y][x].Veg - 1].RootDepth[i] ?
                                (Depth - ChannelWaterLevel) : VType[VegMap[y][x].Veg - 1].RootDepth[i]);
            else
              EffThickness = VType[VegMap[y][x].Veg - 1].RootDepth[i] - (Depth - AdjTableDepth);

            SoilDeficit += (SoilMap[y][x].Porosity[i] - SoilMap[y][x].Moist[i]) * EffThickness;
          }
        }
        /* Also add deep layer water capacity if water table is below root zone layers */
        if (AdjTableDepth > Depth) {
          i = SType[SoilMap[y][x].Soil - 1].NLayers;
          SoilDeficit += (SoilMap[y][x].Porosity[i] - SoilMap[y][x].Moist[i]) * (SoilMap[y][x].Depth - AdjTableDepth);
        }
        SoilDeficit /= (AdjTableDepth - ChannelWaterLevel); /* Convert from depth to percentage capacity */
        /* Note that we don't need to consider the Adjust value, since all layers considered here
           are intersected by the channel cut, so multiplying by depth
           and subsequently dividing by total depth cancels out the cut storage effect */
        
        water_in_stream = channel_grid_lateral_outflow(ChannelData->stream_map, x, y,
                                                       AdjTableDepth,
                                                       Transmissivity, SoilDeficit,
                                                       Map->DX, Map->DY, Dt);
        water_in_stream /= (Map->DX * Map->DY);
        water_out_stream -= water_in_stream;
        
        SoilMap[y][x].ChannelInfiltration += water_in_stream;
      }
    }
    
    /* Subsurface Component - decrease water change only by as much
     as possible (up to transmissivity) to not violate TotalAvailableWater */
    AvailableWater = TotalAvailableWater - water_out_road - water_out_stream;
    OutFlow = (OutFlow > AvailableWater) ? AvailableWater : OutFlow; 
    SoilMap[y][x].SatFlow -= (water_out_road + water_out_stream);
    
    /* Assign the water to appropriate surrounding pixels */
    if (SubTotalDir[y][x] > 0)
      OutFlow /= (float) SubTotalDir[y][x];
    else
      OutFlow = 0.;
    
    /* Sort flow directions in order of gradient */
    for (k = 0; k < NDIRS; k++) {
      nx = xdirection[k] + x;
      ny = ydirection[k] + y;
      kOrdered[k].x = nx;
      kOrdered[k].y = ny;
      if (valid_cell(Map, nx, ny))
        kOrdered[k].Rank = (float) SubDir[y][x][k];
      else
        kOrdered[k].Rank = 0.0;
    }
    
    quick(kOrdered, NDIRS);
    
    for (k = (NDIRS - 1); k > -1; k--) {
      nx = kOrdered[k].x;
      ny = kOrdered[k].y;
      
      if (valid_cell(Map, nx, ny) && INBASIN(TopoMap[ny][nx].Mask)) {
        
        PotentialSatFlow = OutFlow * kOrdered[k].Rank;
        
        if (Options->FlowGradient != WATERTABLE ||
           (TopoMap[y][x].Dem - SoilMap[y][x].Depth) > TopoMap[ny][nx].Dem) {
          ActualSatFlow = PotentialSatFlow;
        } else if (PotentialSatFlow > 0.0) {
          
          AdjTableDepthK = TopoMap[ny][nx].Dem - SoilMap[ny][nx].WaterLevel;
          AdjWaterLevelK = SoilMap[ny][nx].WaterLevel;
          
          /* Find layer containing initial water table in current cell */
          i = 0;
          Depth = 0.0;
          while (i < SType[SoilMap[y][x].Soil - 1].NLayers && Depth < AdjTableDepth) {
            if (VType[VegMap[y][x].Veg - 1].RootDepth[i] < (SoilMap[y][x].Depth - Depth))
              Depth += VType[VegMap[y][x].Veg - 1].RootDepth[i];
            else
              Depth = SoilMap[y][x].Depth;
            i++;
          }
          if (Depth > AdjTableDepth)
            i--; /* Water table in latest root zone layer */
          if (i < 0)
            i = 0;
          
          /* Extract water from top down until water level matches or outflow is satisfied */
          ActualSatFlow = 0.0;
          while (AdjWaterLevelK < AdjWaterLevel && PotentialSatFlow > 0.0
                && i <= SType[SoilMap[y][x].Soil - 1].NLayers) {
            
            if (i < SType[SoilMap[y][x].Soil - 1].NLayers) {
              LayerContribWater = (SoilMap[y][x].Moist[i] - SoilMap[y][x].FCap[i]) *
                                  Adjust[i] * VType[VegMap[y][x].Veg - 1].RootDepth[i];
              if (LayerContribWater > PotentialSatFlow)
                LayerContribWater = PotentialSatFlow;
              
              LayerStorageCap = (SoilMap[y][x].Porosity[i] - SoilMap[y][x].FCap[i]) *
                                Adjust[i] * VType[VegMap[y][x].Veg - 1].RootDepth[i];
              DeltaTableDepth = (-1.0 * LayerContribWater / LayerStorageCap) *
                                VType[VegMap[y][x].Veg - 1].RootDepth[i];
            } else { /* Water table in deep layer */
              LayerContribWater = (SoilMap[y][x].Moist[i] - SoilMap[y][x].FCap[i]) *
                                  Adjust[i] * (SoilMap[y][x].Depth - VType[VegMap[y][x].Veg - 1].TotalDepth);
              if (LayerContribWater > PotentialSatFlow)
                LayerContribWater = PotentialSatFlow;
              
              LayerStorageCap = (SoilMap[y][x].Porosity[i] - SoilMap[y][x].FCap[i]) *
                                Adjust[i] * (SoilMap[y][x].Depth - VType[VegMap[y][x].Veg - 1].TotalDepth);
              DeltaTableDepth = (-1.0 * LayerContribWater / LayerStorageCap) *
                                (SoilMap[y][x].Depth - VType[VegMap[y][x].Veg - 1].TotalDepth);
            }
            
            /* Find layer containing initial water table in downstream cell */
            j = 0;
            Depth = 0.0;
            while (j < SType[SoilMap[ny][nx].Soil - 1].NLayers && Depth < AdjTableDepthK) {
              if (VType[VegMap[ny][nx].Veg - 1].RootDepth[j] < (SoilMap[ny][nx].Depth - Depth))
                Depth += VType[VegMap[ny][nx].Veg - 1].RootDepth[j];
              else
                Depth = SoilMap[ny][nx].Depth;
              j++;
            }
            if (Depth > AdjTableDepthK)
              j--; /* Water table in last root zone layer */
            
            /* Contribute water to downstream cell as long as its water level remains lower */
            while (AdjWaterLevelK < AdjWaterLevel && LayerContribWater > 0.0 &&
                  j >= 0) {
              
              LayerContribWaterK = LayerContribWater;
              
              if (j < SType[SoilMap[ny][nx].Soil - 1].NLayers) {
                LayerStorageK = (SoilMap[ny][nx].Moist[j] - SoilMap[ny][nx].FCap[j]) *
                                Network[ny][nx].Adjust[j] * VType[VegMap[ny][nx].Veg - 1].RootDepth[j];
                LayerStorageCapK = (SoilMap[ny][nx].Porosity[j] - SoilMap[ny][nx].FCap[j]) *
                                   Network[ny][nx].Adjust[j] * VType[VegMap[ny][nx].Veg - 1].RootDepth[j];
                if (LayerContribWaterK > (LayerStorageCapK - LayerStorageK))
                  LayerContribWaterK = LayerStorageCapK - LayerStorageK;
                DeltaTableDepthK = (LayerContribWaterK / LayerStorageCapK) *
                                   VType[VegMap[ny][nx].Veg - 1].RootDepth[j];
              } else { /* Water table in deep layer */
                LayerStorageK = (SoilMap[ny][nx].Moist[j] - SoilMap[ny][nx].FCap[j]) *
                                Network[ny][nx].Adjust[j] * (SoilMap[ny][nx].Depth - VType[VegMap[ny][nx].Veg - 1].TotalDepth);
                LayerStorageCapK = (SoilMap[ny][nx].Porosity[j] - SoilMap[ny][nx].FCap[j]) *
                                   Network[ny][nx].Adjust[j] * (SoilMap[ny][nx].Depth - VType[VegMap[ny][nx].Veg - 1].TotalDepth);
                if (LayerContribWaterK > (LayerStorageCapK - LayerStorageK))
                  LayerContribWaterK = LayerStorageCapK - LayerStorageK;
                DeltaTableDepthK = (LayerContribWaterK / LayerStorageCapK) *
                                   (SoilMap[ny][nx].Depth - VType[VegMap[ny][nx].Veg - 1].TotalDepth);
              }
              
              LayerUseFrac = (AdjWaterLevelK - AdjWaterLevel) / (DeltaTableDepth - DeltaTableDepthK);
              if (LayerUseFrac > 1.0)
                LayerUseFrac = 1.0;
              if (LayerUseFrac < 0.0)
                LayerUseFrac = 0.0;
              
              AdjTableDepth -= DeltaTableDepth * LayerUseFrac;
              if (AdjTableDepth > SoilMap[y][x].Depth)
                AdjTableDepth = SoilMap[y][x].Depth;
              
              AdjTableDepthK -= DeltaTableDepthK * LayerUseFrac;
              if (AdjTableDepthK > SoilMap[ny][nx].Depth)
                AdjTableDepthK = SoilMap[ny][nx].Depth;
              
              AdjWaterLevel += DeltaTableDepth * LayerUseFrac;
              AdjWaterLevelK += DeltaTableDepthK * LayerUseFrac;
              
              LayerContribWater -= LayerContribWaterK * LayerUseFrac;
              ActualSatFlow += LayerContribWaterK * LayerUseFrac;
              PotentialSatFlow -= LayerContribWaterK * LayerUseFrac;
              j--; /* Move to next-highest layer in downhill cell */
            }
            
            /* Add excess water to downhill surface ponding */
            if (j < 0 && LayerContribWater > 0.0 && AdjWaterLevelK < AdjWaterLevel) {
              LayerUseFrac = (AdjWaterLevelK - AdjWaterLevel) / (DeltaTableDepth - LayerContribWater);
              if (LayerUseFrac > 1.0)
                LayerUseFrac = 1.0;
              if (LayerUseFrac < 0.0)
                LayerUseFrac = 0.0;
              
              AdjTableDepth -= DeltaTableDepth * LayerUseFrac;
              if (AdjTableDepth > SoilMap[y][x].Depth)
                AdjTableDepth = SoilMap[y][x].Depth;
              
              AdjTableDepthK -= LayerContribWater * LayerUseFrac;
              if (AdjTableDepthK > SoilMap[ny][nx].Depth)
                AdjTableDepthK = SoilMap[ny][nx].Depth;
              
              AdjWaterLevel += DeltaTableDepth * LayerUseFrac;
              AdjWaterLevelK += LayerContribWater * LayerUseFrac;
              
              ActualSatFlow += LayerContribWater * LayerUseFrac;
              PotentialSatFlow -= LayerContribWater * LayerUseFrac;
            }
            i++; /* Move to next-lowest layer in uphill cell */
          }
        } else {
          ActualSatFlow = 0.0;
        }
        SoilMap[ny][nx].SatFlow += ActualSatFlow;
        SoilMap[y][x].SatFlow -= ActualSatFlow;
      }
    }
  }
  
 for(i=0; i<Map->NY; i++) { 
    free(SubTotalDir[i]);
    free(SubFlowGrad[i]);
    for(j=0; j<Map->NX; j++){
      free(SubDir[i][j]);
    }
    free(SubDir[i]);
  }
  free(SubDir);
  free(SubTotalDir);
  free(SubFlowGrad);

  /**********************************************************************/
  /* Dump saturation extent file to screen.
     Saturation extent is based on the number of pixels with a water table 
     that is at least MTHRESH of soil depth. */ 
  
  count =0;
  totalcount = 0;
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
	     mgrid = (SoilMap[y][x].Depth - SoilMap[y][x].TableDepth)/SoilMap[y][x].Depth;
	     if (mgrid > MTHRESH) 
		   count += 1;
		 totalcount += 1;
      }
    }
  }
 
  sat = 100.*((float)count/(float)totalcount);
  
  sprintf(satoutfile, "%ssaturation_extent.txt", DumpPath);
  
  if((fs = fopen(satoutfile,"a")) == NULL){
    printf("Cannot open saturation extent output file.\n");
    exit(0);
  }
  
  SPrintDate(&(Time->Current), buffer);
  fprintf(fs, "%-20s %.4f \n", buffer, sat); 
  fclose(fs);    
}

