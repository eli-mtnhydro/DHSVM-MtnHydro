/* -------------------------------------------------------------
   file: DHSVMChannel.c
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Battelle Memorial Institute
   Pacific Northwest Laboratory
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Created August 30, 1996 by  William A Perkins
   $Id: DHSVMChannel.c, v3.1.2  2013/12/20   Ning Exp $
   ------------------------------------------------------------- */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "constants.h"
#include "getinit.h"
#include "DHSVMChannel.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "settings.h"
#include "errorhandler.h"
#include "fileio.h"

/* -----------------------------------------------------------------------------
   InitChannel
   Reads stream and road files and builds the networks.
   -------------------------------------------------------------------------- */
void
InitChannel(LISTPTR Input, MAPSIZE *Map, int deltat, CHANNEL *channel,
	    SOILTABLE *SType, SOILPIX ** SoilMap, VEGTABLE *VType, VEGPIX **VegMap,
	    int *MaxStreamID, int *MaxRoadID, OPTIONSTRUCT *Options)
{
  int i;
  STRINIENTRY StrEnv[] = {
    {"ROUTING", "STREAM NETWORK FILE", "", ""},
    {"ROUTING", "STREAM MAP FILE", "", ""},
    {"ROUTING", "STREAM CLASS FILE", "", ""},
    {"ROUTING", "RIPARIAN VEG FILE", "", ""},
    {"ROUTING", "ROAD NETWORK FILE", "", "none"},
    {"ROUTING", "ROAD MAP FILE", "", "none"},
    {"ROUTING", "ROAD CLASS FILE", "", "none"},
    {NULL, NULL, "", NULL}
  };

  printf("\nInitializing Road/Stream Networks\n");

  /* Read the key-entry pairs from the ROUTING section in the input file */
  for (i = 0; StrEnv[i].SectionName; i++) {
    GetInitString(StrEnv[i].SectionName, StrEnv[i].KeyName, StrEnv[i].Default,
      StrEnv[i].VarStr, (unsigned long)BUFSIZE, Input);
    if (!strncmp(StrEnv[i].KeyName, "RIPARIAN VEG FILE", 6)) {
      if (Options->StreamTemp) {
      if (IsEmptyStr(StrEnv[i].VarStr))
        ReportError(StrEnv[i].KeyName, 51);
      }
    }
    else {
      if (IsEmptyStr(StrEnv[i].VarStr))
      ReportError(StrEnv[i].KeyName, 51);
    }
  }

  channel->stream_class = NULL;
  channel->road_class = NULL;
  channel->streams = NULL;
  channel->roads = NULL;
  channel->stream_map = NULL;
  channel->road_map = NULL;

  channel_init();
  channel_grid_init(Map->NX, Map->NY);

  if (strncmp(StrEnv[stream_class].VarStr, "none", 4)) {

    printf("\tReading Stream data\n");

    if ((channel->stream_class =
	 channel_read_classes(StrEnv[stream_class].VarStr, stream_class)) == NULL) {
      ReportError(StrEnv[stream_class].VarStr, 5);
    }
    if ((channel->streams =
	 channel_read_network(StrEnv[stream_network].VarStr,
			      channel->stream_class, MaxStreamID)) == NULL) {
      ReportError(StrEnv[stream_network].VarStr, 5);
    }
    if ((channel->stream_map =
	 channel_grid_read_map(channel->streams,
			       StrEnv[stream_map].VarStr, SType, SoilMap, VType, VegMap)) == NULL) {
      ReportError(StrEnv[stream_map].VarStr, 5);
    }
    
    /* Associate each channel segment with its constituent map cells */
    channel_combine_map_network(channel->streams, channel->stream_map, Map);
    
    error_handler(ERRHDL_STATUS,
		  "InitChannel: computing stream network routing coefficients");
    channel_routing_parameters(channel->streams, (double) deltat);
  }

  if (Options->StreamTemp) {
	if (strncmp(StrEnv[riparian_veg].VarStr, "none", 4)) {
	  printf("\tReading channel riparian vegetation params\n");
	  channel_read_rveg_param(channel->streams, StrEnv[riparian_veg].VarStr, MaxStreamID);
	}
  }

  if (strncmp(StrEnv[road_class].VarStr, "none", 4)) {

    printf("\tReading Road data\n");

    if ((channel->road_class =
	 channel_read_classes(StrEnv[road_class].VarStr, road_class)) == NULL) {
      ReportError(StrEnv[road_class].VarStr, 5);
    }
    if ((channel->roads =
	 channel_read_network(StrEnv[road_network].VarStr,
			      channel->road_class, MaxRoadID)) == NULL) {
      ReportError(StrEnv[road_network].VarStr, 5);
    }
    if ((channel->road_map =
	 channel_grid_read_map(channel->roads,
			       StrEnv[road_map].VarStr, SType, SoilMap, VType, VegMap)) == NULL) {
      ReportError(StrEnv[road_map].VarStr, 5);
    }
    error_handler(ERRHDL_STATUS,
		  "InitChannel: computing road network routing coefficients");
    channel_routing_parameters(channel->roads, (double) deltat);
  }
}

/* -------------------------------------------------------------
   InitChannelDump
   ------------------------------------------------------------- */
void InitChannelDump(OPTIONSTRUCT *Options, CHANNEL * channel, 
					 char *DumpPath)
{
  char buffer[NAMESIZE];

  if (channel->streams != NULL) {
    sprintf(buffer, "%sStream.Flow", DumpPath);
    OpenFile(&(channel->streamout), buffer, "w", TRUE);
    sprintf(buffer, "%sStreamflow.Only", DumpPath);
    OpenFile(&(channel->streamflowout), buffer, "w", TRUE);
    /* output files for John's RBM model */
	if (Options->StreamTemp) {
      //inflow to segment
      sprintf(buffer, "%sInflow.Only", DumpPath);
      OpenFile(&(channel->streaminflow), buffer, "w", TRUE);
      // outflow ( redundant but it's a check
      sprintf(buffer, "%sOutflow.Only", DumpPath);
      OpenFile(&(channel->streamoutflow), buffer, "w", TRUE);
      //net incoming short wave
      sprintf(buffer, "%sNSW.Only", DumpPath);
      OpenFile(&(channel->streamNSW), buffer, "w", TRUE);
      // net incoming long wave
      sprintf(buffer, "%sNLW.Only", DumpPath);
      OpenFile(&(channel->streamNLW), buffer, "w", TRUE);
      //Vapor pressure
      sprintf(buffer, "%sVP.Only", DumpPath);
      OpenFile(&(channel->streamVP), buffer, "w", TRUE);
      //wind speed
      sprintf(buffer, "%sWND.Only", DumpPath);
      OpenFile(&(channel->streamWND), buffer, "w", TRUE);
      //air temperature
      sprintf(buffer, "%sATP.Only", DumpPath);
      OpenFile(&(channel->streamATP), buffer, "w", TRUE);
      //melt water in flow
      sprintf(buffer, "%sMelt.Only", DumpPath);
      OpenFile(&(channel->streamMelt), buffer, "w", TRUE);                      
	}
  }
  if (channel->roads != NULL) {
    sprintf(buffer, "%sRoad.Flow", DumpPath);
    OpenFile(&(channel->roadout), buffer, "w", TRUE);
    sprintf(buffer, "%sRoadflow.Only", DumpPath);
    OpenFile(&(channel->roadflowout), buffer, "w", TRUE);

  }
}


/* -------------------------------------------------------------
   ChannelCulvertFlow    
   computes outflow of channel/road network to a grid cell, if it
   contains a sink
   ------------------------------------------------------------- */
double ChannelCulvertFlow(int y, int x, CHANNEL * ChannelData)
{
  if (channel_grid_has_channel(ChannelData->road_map, x, y)) {
    return channel_grid_outflow(ChannelData->road_map, x, y);
  }
  else {
    return 0;
  }
}

/* -------------------------------------------------------------
   RouteChannel
   ------------------------------------------------------------- */
void
RouteChannel(CHANNEL *ChannelData, TIMESTRUCT *Time, MAPSIZE *Map,
	    TOPOPIX **TopoMap, SOILPIX **SoilMap, AGGREGATED *Total, 
	     OPTIONSTRUCT *Options, ROADSTRUCT **Network, SOILTABLE *SType,
	     VEGTABLE *VType, VEGPIX **VegMap,
	     PRECIPPIX **PrecipMap, float Tair, float Rh, SNOWPIX **SnowMap)
{
  int k, x, y;
  int flag;
  char buffer[32];
  float CulvertFlow;
  float temp;
  float max_bank_height = 0.0;
  float Depth;
  float EffThickness;
  float MaxInfiltrationCap;
  int i;
  float StreamInfiltration;
  
  SPrintDate(&(Time->Current), buffer);
  flag = IsEqualTime(&(Time->Current), &(Time->Start));
  
  /* Give any surface water to roads w/o sinks */
  if (ChannelData->roads != NULL) {
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
        if (INBASIN(TopoMap[y][x].Mask)) {
          if (channel_grid_has_channel(ChannelData->road_map, x, y) && 
              !channel_grid_has_sink(ChannelData->road_map, x, y)) {
            SoilMap[y][x].RoadInt += SoilMap[y][x].IExcess; /* road w/o sink */
            channel_grid_inc_inflow(ChannelData->road_map, x, y, SoilMap[y][x].IExcess * Map->DX * Map->DY);
            SoilMap[y][x].IExcess = 0.0f;
          }
        }
      }
    }
    
    /* Route the road network and save results */
    channel_route_network(ChannelData->roads, Time->Dt);
    channel_save_outflow_text(buffer, ChannelData->roads,
                              ChannelData->roadout, ChannelData->roadflowout, flag);
  }
  
  /* Add IExcess and ChannelInt to stream channels (by descending elevation) */
  /* Also add culvert outflow to surface water */
  Total->CulvertReturnFlow = 0.0;
  for (k = (Map->NumCells - 1); k > -1;  k--) {
    y = Map->OrderedCells[k].y;
    x = Map->OrderedCells[k].x;
    
    CulvertFlow = ChannelCulvertFlow(y, x, ChannelData);
    CulvertFlow /= Map->DX * Map->DY;
    /* CulvertFlow = (CulvertFlow > 0.0) ? CulvertFlow : 0.0; */
    
    if (channel_grid_has_channel(ChannelData->stream_map, x, y)) {
      channel_grid_inc_inflow(ChannelData->stream_map, x, y,
                              (SoilMap[y][x].IExcess + CulvertFlow) * Map->DX * Map->DY);
      
      channel_grid_satflow(ChannelData->stream_map, x, y);
      
      if (SnowMap[y][x].Outflow > SoilMap[y][x].IExcess)
        temp = SoilMap[y][x].IExcess;
      else
        temp = SnowMap[y][x].Outflow;
      channel_grid_inc_melt(ChannelData->stream_map, x, y, temp * Map->DX * Map->DY);
      
      SoilMap[y][x].ChannelInt += SoilMap[y][x].IExcess;
      Total->CulvertToChannel += CulvertFlow;
      SoilMap[y][x].IExcess = 0.0f;
    }
    else {
      SoilMap[y][x].IExcess += CulvertFlow;
      Total->CulvertReturnFlow += CulvertFlow;
    }
  }
  
  /* Route stream channels */
  if (ChannelData->streams != NULL) {
    
    /* Account for infiltration out of streams that are above the water table */
    /* Loop thru all of the cells in descending order of elevation */
    for (k = (Map->NumCells - 1); k > -1;  k--) {
      y = Map->OrderedCells[k].y;
      x = Map->OrderedCells[k].x;
      if (channel_grid_has_channel(ChannelData->stream_map, x, y)) {
        
        /* Only allow stream re-infiltration if water table is 1 mm below deepest channel */
        max_bank_height = channel_grid_cell_maxbankht(ChannelData->stream_map, x, y);
        if (SoilMap[y][x].TableDepth > (max_bank_height + 0.001)) {
          
          /* Find maximum amount of water that can be added to subsurface */
          /* So that water table just reaches bottom of lowest channel cut */
          MaxInfiltrationCap = 0.0;
          Depth = 0.0;
          for (i = 0; i < SType[SoilMap[y][x].Soil - 1].NLayers && Depth < SoilMap[y][x].Depth; i++) {
            if (VType[VegMap[y][x].Veg - 1].RootDepth[i] < (SoilMap[y][x].Depth - Depth))
              Depth += VType[VegMap[y][x].Veg - 1].RootDepth[i];
            else
              Depth = SoilMap[y][x].Depth;
            
            if (Depth > max_bank_height) {
              EffThickness = ((Depth - max_bank_height) < VType[VegMap[y][x].Veg - 1].RootDepth[i] ?
                                (Depth - max_bank_height) : VType[VegMap[y][x].Veg - 1].RootDepth[i]);
              MaxInfiltrationCap += (SoilMap[y][x].Porosity[i] - SoilMap[y][x].Moist[i]) * EffThickness;
            }
          }
          /* Also add deep layer water capacity if water table is below root zone layers */
          if (SoilMap[y][x].TableDepth > Depth) {
            i = SType[SoilMap[y][x].Soil - 1].NLayers;
            MaxInfiltrationCap += (SoilMap[y][x].Porosity[i] - SoilMap[y][x].Moist[i]) * (SoilMap[y][x].Depth - Depth);
          }
          MaxInfiltrationCap *= Map->DX * Map->DY;
        } else
          MaxInfiltrationCap = 0.0;
        
        channel_grid_calc_infiltration(ChannelData->stream_map, x, y, Time->Dt,
                                       SoilMap[y][x].TableDepth, MaxInfiltrationCap);
    }
    }
    
    channel_route_network(ChannelData->streams, Time->Dt);
    
    for (k = (Map->NumCells - 1); k > -1;  k--) {
      y = Map->OrderedCells[k].y;
      x = Map->OrderedCells[k].x;
      if (channel_grid_has_channel(ChannelData->stream_map, x, y)) {
        
        StreamInfiltration = channel_grid_infiltration(ChannelData->stream_map, x, y);
        StreamInfiltration /= Map->DX * Map->DY;
        
        SoilMap[y][x].SatFlow += StreamInfiltration;
        SoilMap[y][x].ChannelInfiltration += StreamInfiltration;
    }
    }
    
    channel_save_outflow_text(buffer, ChannelData->streams,
			      ChannelData->streamout,
			      ChannelData->streamflowout, flag);
	/* save parameters for John's RBM model */
	if (Options->StreamTemp)
	  channel_save_outflow_text_cplmt(Time, buffer,ChannelData->streams,ChannelData, flag);
  }
  
}

/* -------------------------------------------------------------
   ChannelCut
   computes necessary parameters for cell storage adjustment from
   channel/road dimensions
   ------------------------------------------------------------- */
void ChannelCut(int y, int x, CHANNEL * ChannelData, ROADSTRUCT * Network)
{
  float bank_height = 0.0;
  float cut_area = 0.0;

  if (channel_grid_has_channel(ChannelData->stream_map, x, y)) {
    bank_height = channel_grid_cell_bankht(ChannelData->stream_map, x, y);
    cut_area = channel_grid_cell_width(ChannelData->stream_map, x, y) * 
      channel_grid_cell_length(ChannelData->stream_map, x, y);
  }
  else if (channel_grid_has_channel(ChannelData->road_map, x, y)) {
    bank_height = channel_grid_cell_bankht(ChannelData->road_map, x, y);
    cut_area = channel_grid_cell_width(ChannelData->road_map, x, y) * 
      channel_grid_cell_length(ChannelData->road_map, x, y);
  }
  Network->Area = cut_area;
  Network->BankHeight = bank_height;
}

/* -------------------------------------------------------------
 ChannelLimitVegFC
 adjusts fractional cover (FC) if necessary
 so that FC does not exceed 1 - channel area / cell area
 ------------------------------------------------------------- */
void ChannelLimitVegFC(int y, int x, float DXDY, CHANNEL * ChannelData,
                       VEGTABLE *VType, VEGPIX *LocalVeg)
{
  float cut_area = 0.0;
  float FCmax = 0.0;
  
  /* Only understory FC is limited, while overstory FC can be > FCmax
   since riparian forest canopies may overhang the stream channel */
  
  if (channel_grid_has_channel(ChannelData->stream_map, x, y) &&
      VType->UnderStory == TRUE) {
    
    cut_area = channel_grid_cell_width(ChannelData->stream_map, x, y) * 
      channel_grid_cell_length(ChannelData->stream_map, x, y);
    
    FCmax = 1. - cut_area / DXDY;
    if (FCmax > 1.)
      FCmax = 1.;
    if (FCmax < 0.01) /* Require minimum of 1 percent FC if understory is present */
      FCmax = 0.01;
    
    if (VType->OverStory == TRUE && LocalVeg->Fract[1] > FCmax) {
      // printf("Warning: overriding understory FC %f with %f due to channel area\n",
      //        LocalVeg->Fract[1],FCmax);
      LocalVeg->Fract[1] = FCmax;
    } else if (LocalVeg->Fract[0] > FCmax) {
      // printf("Warning: overriding understory FC %f with %f due to channel area\n",
      //        LocalVeg->Fract[0],FCmax);
      LocalVeg->Fract[0] = FCmax;
    }
  }
}

/* -------------------------------------------------------------
   ChannelFraction
   This computes the (sub)surface flow fraction for a road
   ------------------------------------------------------------- */
uchar ChannelFraction(TOPOPIX * topo, ChannelMapRec * rds)
{
  float effective_width = 0;
  float total_width;
  ChannelMapRec *r;
  float fract = 0.0;

  if (rds == NULL) {
    return 0;
  }
  total_width = topo->FlowGrad / topo->Slope;
  effective_width = 0.0;

  for (r = rds; r != NULL; r = r->next) {
    effective_width += r->length * sin(fabs(topo->Aspect - r->aspect));
  }
  fract = effective_width / total_width * 255.0;
  fract = (fract > 255.0 ? 255.0 : floor(fract + 0.5));

  return (uchar) fract;
}
