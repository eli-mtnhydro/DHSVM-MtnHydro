
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
   Reads stream files and builds the network
   -------------------------------------------------------------------------- */
void
InitChannel(LISTPTR Input, MAPSIZE *Map, int deltat, CHANNEL *ChannelData,
	    SOILTABLE *SType, SOILPIX ** SoilMap, VEGTABLE *VType, VEGPIX **VegMap,
	    LAKETABLE *LType, TOPOPIX **TopoMap,
	    int *MaxStreamID,  OPTIONSTRUCT *Options)
{
  int i, x, y;
  ChannelMapPtr cell;
  STRINIENTRY StrEnv[] = {
    {"ROUTING", "STREAM NETWORK FILE", "", ""},
    {"ROUTING", "STREAM MAP FILE", "", ""},
    {"ROUTING", "STREAM CLASS FILE", "", ""},
    {NULL, NULL, "", NULL}
  };

  printf("\nInitializing Stream Networks\n");

  /* Read the key-entry pairs from the ROUTING section in the input file */
  for (i = 0; StrEnv[i].SectionName; i++) {
    GetInitString(StrEnv[i].SectionName, StrEnv[i].KeyName, StrEnv[i].Default,
      StrEnv[i].VarStr, (unsigned long)BUFSIZE, Input);
    if (IsEmptyStr(StrEnv[i].VarStr))
      ReportError(StrEnv[i].KeyName, 51);
  }

  ChannelData->stream_class = NULL;
  ChannelData->streams = NULL;
  ChannelData->stream_map = NULL;
  
  channel_grid_init(Map->NX, Map->NY);

  if (strncmp(StrEnv[stream_class].VarStr, "none", 4)) {
    
    printf("\tReading Stream data\n");

    if ((ChannelData->stream_class =
	 channel_read_classes(StrEnv[stream_class].VarStr, stream_class)) == NULL) {
      ReportError(StrEnv[stream_class].VarStr, 5);
    }
    if ((ChannelData->streams =
	 channel_read_network(StrEnv[stream_network].VarStr,
			      ChannelData->stream_class, MaxStreamID)) == NULL) {
      ReportError(StrEnv[stream_network].VarStr, 5);
    }
    if ((ChannelData->stream_map =
	 channel_grid_read_map(ChannelData->streams,
			       StrEnv[stream_map].VarStr, SType, SoilMap, VType, VegMap)) == NULL) {
      ReportError(StrEnv[stream_map].VarStr, 5);
    }
    
    /* Associate each channel segment with its constituent map cells */
    channel_combine_map_network(ChannelData->streams, ChannelData->stream_map, Map);
    
    if (Options->LakeDynamics) {
      
      /* Find channel segments that intersect lake cells */
      for (y = 0; y < Map->NY; y++) {
        for (x = 0; x < Map->NX; x++) {
          cell = ChannelData->stream_map[x][y];
          while (cell != NULL) {
            if (TopoMap[y][x].LakeID != 0) {
              cell->channel->IntersectsLake = TRUE;
              cell->channel->lake = &(LType[TopoMap[y][x].LakeID - 1]);
            }
            cell = cell->next;
      }
      }
      }
      
      /* Find outlets from lakes */
      for (i = 0; i < Map->NumLakes; i++) {
        if (LType[i].OutletID != 0) {
          LType[i].outlet = channel_find_segment(ChannelData->streams, LType[i].OutletID);
          if (LType[i].outlet == NULL) {
            printf("ERROR! cannot find outlet (%d) for lake %d", LType[i].OutletID, i);
          } else {
            LType[i].outlet->IntersectsLake = FALSE;
            LType[i].outlet->lake = NULL;
          }
        }
      }
    } /* End of lake dynamics */
    
    error_handler(ERRHDL_STATUS,
		  "InitChannel: computing stream network routing coefficients");
    channel_routing_parameters(ChannelData->streams, (double) deltat);
  }
}

/* -------------------------------------------------------------
   InitChannelDump
   ------------------------------------------------------------- */
void InitChannelDump(OPTIONSTRUCT *Options, CHANNEL * ChannelData, 
					 char *DumpPath)
{
  char buffer[NAMESIZE];

  if (ChannelData->streams != NULL) {
    sprintf(buffer, "%sStream.Flow", DumpPath);
    OpenFile(&(ChannelData->streamout), buffer, "w", TRUE);
    sprintf(buffer, "%sStreamflow.Only", DumpPath);
    OpenFile(&(ChannelData->streamflowout), buffer, "w", TRUE);
  }
}

/* -------------------------------------------------------------
   RouteChannel
   ------------------------------------------------------------- */
void
RouteChannel(CHANNEL *ChannelData, TIMESTRUCT *Time, MAPSIZE *Map,
	    TOPOPIX **TopoMap, SOILPIX **SoilMap, AGGREGATED *Total, 
	     OPTIONSTRUCT *Options, NETSTRUCT **Network, SOILTABLE *SType,
	     VEGTABLE *VType, VEGPIX **VegMap, EVAPPIX **Evap, LAKETABLE *LType)
{
  int k, x, y;
  int flag;
  char buffer[32];
  float max_bank_height = 0.0;
  float AdjTableDepth;
  float Transmissivity;
  float ChannelTableDepth;
  float Depth;
  float EffThickness;
  float MaxInfiltrationCap;
  int i;
  float StreamInfiltration, StreamEvap;
  
  SPrintDate(&(Time->Current), buffer);
  flag = IsEqualTime(&(Time->Current), &(Time->Start));
  
  /* Add IExcess and ChannelInt to stream channels (by descending elevation) */
  for (k = (Map->NumCells - 1); k > -1;  k--) {
    y = Map->OrderedCells[k].y;
    x = Map->OrderedCells[k].x;
    
    if (channel_grid_has_channel(ChannelData->stream_map, x, y) && TopoMap[y][x].LakeID == 0) {
      channel_grid_inc_inflow(ChannelData->stream_map, x, y,
                              SoilMap[y][x].IExcess * Map->DX * Map->DY);
      
      channel_grid_satflow(ChannelData->stream_map, x, y);
      
      SoilMap[y][x].ChannelInt += SoilMap[y][x].IExcess;
      SoilMap[y][x].IExcess = 0.0f;
    }
  }
  
  /* Route stream channels */
  /* Account for infiltration out of streams that are above the water table */
  /* Loop thru all of the cells in descending order of elevation */
  for (k = (Map->NumCells - 1); k > -1;  k--) {
    y = Map->OrderedCells[k].y;
    x = Map->OrderedCells[k].x;
    if (channel_grid_has_channel(ChannelData->stream_map, x, y) && TopoMap[y][x].LakeID == 0) {
      
      /* Only allow stream re-infiltration if local water table is 1 mm below deepest channel */
      /* Update water table depth directly below channel based on lateral diffusion on previous timestep */
      AdjTableDepth = TopoMap[y][x].Dem - SoilMap[y][x].WaterLevel;
      max_bank_height = channel_grid_cell_maxbankht(ChannelData->stream_map, x, y);
      Transmissivity = CalcTransmissivity(AdjTableDepth, max_bank_height,
                                          SoilMap[y][x].KsLat,
                                          SoilMap[y][x].KsLatExp,
                                          SType[SoilMap[y][x].Soil - 1].DepthThresh);
      ChannelTableDepth = channel_grid_table_depth(ChannelData->stream_map, x, y, Time->Dt,
                                                   AdjTableDepth, Transmissivity,
                                                   (SoilMap[y][x].Porosity[Network[y][x].CutBankZone] -
                                                     SoilMap[y][x].FCap[Network[y][x].CutBankZone]),
                                                     Map->DX);
      
      if (ChannelTableDepth > max_bank_height) {
        
        /* Find maximum amount of water that can be added to subsurface */
        /* So that water table immediately below channel just reaches bottom of lowest channel cut */
        MaxInfiltrationCap = 0.0;
        Depth = 0.0;
        for (i = 0; i < SType[SoilMap[y][x].Soil - 1].NLayers && Depth < ChannelTableDepth; i++) {
          if (VType[VegMap[y][x].Veg - 1].RootDepth[i] < (SoilMap[y][x].Depth - Depth))
            Depth += VType[VegMap[y][x].Veg - 1].RootDepth[i];
          else
            Depth = SoilMap[y][x].Depth;
          
          if (Depth > max_bank_height) {
            EffThickness = ((Depth - max_bank_height) < VType[VegMap[y][x].Veg - 1].RootDepth[i] ?
                              (Depth - max_bank_height) : VType[VegMap[y][x].Veg - 1].RootDepth[i]);
            if (Depth < ChannelTableDepth)
              MaxInfiltrationCap += (SoilMap[y][x].Porosity[i] - SoilMap[y][x].Moist[i]) * EffThickness;
            else {
              EffThickness -= (Depth - ChannelTableDepth);
              MaxInfiltrationCap += (SoilMap[y][x].Porosity[i] - SoilMap[y][x].FCap[i]) * EffThickness;
            }
          }
        }
        /* Also add deep layer water capacity if water table is below root zone layers */
        if (ChannelTableDepth > Depth) {
          i = SType[SoilMap[y][x].Soil - 1].NLayers;
          MaxInfiltrationCap += (SoilMap[y][x].Porosity[i] - SoilMap[y][x].FCap[i]) * (ChannelTableDepth - Depth);
        }
      } else
        MaxInfiltrationCap = 0.0;
      
      channel_grid_calc_infiltration(ChannelData->stream_map, x, y, Time->Dt,
                                     AdjTableDepth, MaxInfiltrationCap, Map->DX);
    }
  }
  
  /* Route lakes, then route streams */
  if (Options->LakeDynamics) {
    for (i = 0; i < Map->NumLakes; i++) {
      LType[i].outlet->lake_inflow += LType[i].Outflow * LType[i].Area;
    }
  }
  
  channel_route_network(ChannelData->streams, Time->Dt);
  
  for (k = (Map->NumCells - 1); k > -1;  k--) {
    y = Map->OrderedCells[k].y;
    x = Map->OrderedCells[k].x;
    if (channel_grid_has_channel(ChannelData->stream_map, x, y) && TopoMap[y][x].LakeID == 0) {
      
      StreamInfiltration = channel_grid_infiltration(ChannelData->stream_map, x, y);
      StreamInfiltration /= Map->DX * Map->DY;
      SoilMap[y][x].SatFlow += StreamInfiltration;
      SoilMap[y][x].ChannelInfiltration += StreamInfiltration;
      
      StreamEvap = channel_grid_evaporation(ChannelData->stream_map, x, y);
      StreamEvap /= Map->DX * Map->DY;
      VegMap[y][x].MoistureFlux += StreamEvap;
      Evap[y][x].ETot += StreamEvap;
      Evap[y][x].EvapChannel = StreamEvap;
    }
  }
  
  channel_save_outflow_text(buffer, ChannelData->streams,
                            ChannelData->streamout,
                            ChannelData->streamflowout, flag);
  
}

/* -------------------------------------------------------------
   ChannelCut
   computes necessary parameters for cell storage adjustment from
   channel dimensions
   ------------------------------------------------------------- */
void ChannelCut(int y, int x, CHANNEL * ChannelData, NETSTRUCT * Network)
{
  float bank_height = 0.0;
  float cut_area = 0.0;

  if (channel_grid_has_channel(ChannelData->stream_map, x, y)) {
    bank_height = channel_grid_cell_bankht(ChannelData->stream_map, x, y);
    cut_area = channel_grid_cell_width(ChannelData->stream_map, x, y) * 
      channel_grid_cell_length(ChannelData->stream_map, x, y);
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
