/* The Distributed Hydrology Soil Vegetation Model (DHSVM)
 * Custom version maintained by Eli Boardman | Mountain Hydrology LLC
 * https://mountainhydrology.com
 * Contact: eli.boardman@mountainhydrology.com
 * License: CC BY-NC-SA 4.0
*/

/******************************************************************************/
/* INCLUDES */
/******************************************************************************/
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "constants.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "fileio.h"
#include "getinit.h"
#include "DHSVMChannel.h"
#include "channel.h"

/******************************************************************************/
/* GLOBAL VARIABLES */
/******************************************************************************/

char *version = "Version X.2.1";    /* store version string */
char commandline[BUFSIZE + 1] = "";	/* store command line */
char fileext[BUFSIZ + 1] = "";			/* file extension */
char errorstr[BUFSIZ + 1] = "";			/* error message */

/******************************************************************************/
/* MAIN */
/******************************************************************************/

int main(int argc, char **argv) {
  float **PrismMap = NULL;
  float **SnowPatternMap = NULL;
  float **SnowPatternMapBase = NULL;
  unsigned char ***ShadowMap = NULL;
  float **SkyViewMap = NULL;
  float **PptMultiplierMap = NULL;                                  
  float **MeltMultiplierMap = NULL;                                  
  int MaxStreamID;
  clock_t start, finish1;
  double runtime = 0.0;
  int t = 0;
  int i, x, y;
  int NStats;
  uchar ***MetWeights = NULL;
  
  AGGREGATED Total = {			/* Total or average value of a  variable over the entire basin */
    {0.0, NULL, NULL, NULL, NULL, 0.0, 0.0},												/* EVAPPIX */
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, NULL, NULL, 0.0, 0},								/* PRECIPPIX */
    {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, 0.0, {0.0, 0.0}, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, /* PIXRAD */
    {0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
	  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0},     /* SNOWPIX */ 
    {0, 0.0, NULL, NULL, NULL, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, NULL, NULL, NULL}, /* SOILPIX */
    {0, 0.0, 0.0, 0.0, 0.0, NULL, NULL, NULL, NULL, NULL, 0.0, NULL},                             /* VEGPIX */
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0l, 0.0
  };
  CHANNEL ChannelData = {NULL, NULL, NULL, NULL, NULL};
  DUMPSTRUCT Dump;
  EVAPPIX **EvapMap = NULL;
  INPUTFILES InFiles;
  LAYER Soil;
  LAYER Veg;
  LISTPTR Input = NULL;
  MAPSIZE Map;
  METLOCATION *Stat = NULL;
  OPTIONSTRUCT Options;
  PIXMET LocalMet;
  PRECIPPIX **PrecipMap = NULL;
  PIXRAD **RadiationMap = NULL;
  NETSTRUCT **Network	= NULL;
  SNOWPIX **SnowMap = NULL;
  SOILPIX **SoilMap = NULL;
  SOILTABLE *SType = NULL;
  SOLARGEOMETRY SolarGeo;
  TIMESTRUCT Time;
  TOPOPIX **TopoMap = NULL;
  LAKETABLE *LType = NULL;
  VEGPIX **VegMap = NULL;
  VEGTABLE *VType = NULL;
  DYNAVEG DVeg;
  WATERBALANCE Mass =	{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  
/******************************************************************************/
/* INITIALIZATION */
/******************************************************************************/

  if (argc != 2) {
    fprintf(stderr, "\nUsage: %s inputfile\n\n", argv[0]);
    fprintf(stderr, "DHSVM uses two output streams: \n");
    fprintf(stderr, "Standard Out, for the majority of output \n");
    fprintf(stderr, "Standard Error, for the final mass balance \n");
    fprintf(stderr, "\nTo pipe output correctly to files: \n");
    fprintf(stderr, "(cmd > f1) >& f2 \n");
    fprintf(stderr, "where f1 is stdout_file and f2 is stderror_file\n");
    exit(EXIT_FAILURE);
  }
  
  sprintf(commandline, "%s %s", argv[0], argv[1]);
  printf("%s \n", commandline);
  fprintf(stderr, "%s \n", commandline);
  strcpy(InFiles.Const, argv[1]);
  
  printf("\nRunning DHSVM %s\n", version);
  printf("\nMountain Hydrology Version Copyright Eli Boardman\n");
  printf("\nLICENSE: CC BY-NC-SA 4.0 (Non-Commercial Use Only!)\n");
#ifdef SNOW_ONLY
  printf("----------------------------------\n");
  printf("WARNING: USING SNOW ONLY MODULES (prescribed in makefile)!\n");
  printf("----------------------------------\n");
#endif
  printf("\nSTARTING INITIALIZATION PROCEDURES\n\n");
  
  start = clock();
  
  ReadInitFile(InFiles.Const, &Input);
  InitConstants(Input, &Options, &Map, &SolarGeo, &Time);
  InitFileIO();
  InitTables(Time.NDaySteps, Input, &Options, &Map, &SType, &Soil, &VType, &Veg, &LType);
  InitTerrainMaps(Input, &Options, &Map, &Soil, &Veg, &TopoMap, SType, &SoilMap, VType, &VegMap, &DVeg, LType);
  InitSnowMap(&Map, &SnowMap, &Time);
  InitMappedConstants(Input, &Options, &Map, &SnowMap, VType, &VegMap);
  CheckOut(&Options, Veg, Soil, VType, SType, &Map, TopoMap, VegMap, SoilMap);
  if (Options.Extent != POINT)
    InitChannel(Input, &Map, Time.Dt, &ChannelData, SType, SoilMap, VType, VegMap,
                LType, TopoMap, &MaxStreamID, &Options);
  InitNetwork(Map.NY, Map.NX, Map.DX, Map.DY, TopoMap, SoilMap, 
	      VegMap, VType, &Network, &ChannelData, Veg, &Options);
  InitMetSources(Input, &Options, &Map, TopoMap, Soil.MaxLayers, &Time,
		 &InFiles, &NStats, &Stat);
  InitMetMaps(Input, Time.NDaySteps, &Map, &Options,
	      &PrismMap, &SnowPatternMap, &SnowPatternMapBase,
	      &ShadowMap, &SkyViewMap, &EvapMap, &PrecipMap, &PptMultiplierMap,
	      &MeltMultiplierMap, &RadiationMap, SoilMap, &Soil, VegMap, &Veg, TopoMap);
  InitInterpolationWeights(&Map, &Options, TopoMap, &MetWeights, Stat, NStats);
  InitDump(Input, &Options, &Map, Soil.MaxLayers, Veg.MaxLayers, Time.Dt,
	   TopoMap, &Dump);
  /* Done with initialization, delete the list with input strings */
  DeleteList(Input);
  
#ifndef SNOW_ONLY
  if (Options.Extent != POINT) {
    InitChannelDump(&Options, &ChannelData, Dump.Path);
    ReadChannelState(Dump.InitStatePath, &(Time.Start), ChannelData.streams);
  }
#endif
  
  InitAggregated(&Options, Veg.MaxLayers, Soil.MaxLayers, &Total);
  InitModelState(&(Time.Start), Time.NDaySteps, Time.Dt, &Map, &Options, PrecipMap, SnowMap, SoilMap,
		 Soil, SType, VegMap, Veg, VType, Dump.InitStatePath,
		 TopoMap, Network, &ChannelData);
  InitNewMonth(&Time, &Options, &Map, TopoMap, PrismMap, SnowPatternMap, SnowPatternMapBase, ShadowMap,
	       &InFiles, Veg.NTypes, VType, NStats, Stat, Dump.InitStatePath, &VegMap, SnowMap);
  InitNewDay(Time.Current.JDay, &SolarGeo);
  
  /* Setup for mass balance calculations */
  Aggregate(&Map, &Options, TopoMap, &Soil, &Veg, VegMap, EvapMap, PrecipMap,
	      RadiationMap, SnowMap, SoilMap, &Total, VType, Network, &ChannelData, Time.Dt, Time.NDaySteps);
  Mass.StartWaterStorage = Total.Soil.IExcess + Total.CanopyWater +
                           Total.SoilWater + Total.Snow.Swq + Total.Soil.SatFlow;
  Mass.OldWaterStorage = Mass.StartWaterStorage;
  
/******************************************************************************/
/* CALCULATIONS */
/******************************************************************************/

  while (Before(&(Time.Current), &(Time.End)) ||
  IsEqualTime(&(Time.Current), &(Time.End))) {
    
    ResetAggregate(&Soil, &Veg, &Total, &Options);
    
    if (Options.SnowSlide)
	    Avalanche(&Map, TopoMap, &Time, &Options, SnowMap);
    
    if (Options.DynamicVeg){
      if (IsVegDate(&(Time.Current), &DVeg))        
        UpdateVegMap(&(Time.Current), &Options, &Map, &Veg, &VegMap, VType, &DVeg);
    }
    
    if (IsNewWaterYear(&(Time.Current)))
      InitNewWaterYear(&Time, &Options, &Map, TopoMap, SnowMap, PrecipMap);
    if (IsNewMonth(&(Time.Current), Time.Dt))
      InitNewMonth(&Time, &Options, &Map, TopoMap, PrismMap, SnowPatternMap, SnowPatternMapBase, ShadowMap,
		   &InFiles, Veg.NTypes, VType, NStats, Stat, Dump.InitStatePath, &VegMap, SnowMap);
    if (IsNewDay(Time.DayStep)) {
      InitNewDay(Time.Current.JDay, &SolarGeo);
      PrintDate(&(Time.Current), stdout);
      printf("\n");
    }
    InitNewStep(&InFiles, &Map, &Time, Soil.MaxLayers, &Options, NStats, Stat,
                &SolarGeo, TopoMap, SoilMap);
    if (Options.Extent != POINT) {
      channel_step_initialize_network(ChannelData.streams);
    }
    
    for (y = 0; y < Map.NY; y++) {
      for (x = 0; x < Map.NX; x++) {
  	    if (INBASIN(TopoMap[y][x].Mask)) {
  	      LocalMet =
  	        MakeLocalMetData(y, x, &Map, Time.DayStep, Time.NDaySteps, &Options, NStats,
                            Stat, MetWeights[y][x], TopoMap[y][x].Dem,
                            &(RadiationMap[y][x]), &(PrecipMap[y][x]),
                            PrismMap, SnowPatternMap, &(SnowMap[y][x]),
                            &(VegMap[y][x].Type), &(VegMap[y][x]),
                            PptMultiplierMap[y][x], Time.Current.Month,
                            (Options.Shading ? SkyViewMap[y][x] : 0.0),
                            (Options.Shading ? ShadowMap[Time.DayStep][y][x] : 0.0),
                            SolarGeo.SunMax, SolarGeo.SineSolarAltitude);
  	      
  		    /* Get surface temperature of each soil layer */
    		  for (i = 0; i < Soil.MaxLayers; i++) {
    	        SoilMap[y][x].Temp[i] = LocalMet.Tair;
    		  }
    		  
    		  MassEnergyBalance(&Options, y, x, SolarGeo.SineSolarAltitude,
                          Map.DX, Map.DY, Time.Dt,
                          Options.HeatFlux, Options.CanopyRadAtt,
                          Options.Infiltration, Soil.MaxLayers,
                          Veg.MaxLayers, &LocalMet, &(Network[y][x]),
                          &(PrecipMap[y][x]), MeltMultiplierMap[y][x],
                          &(VType[VegMap[y][x].Veg - 1]), &(VegMap[y][x]),
                          &(SType[SoilMap[y][x].Soil - 1]), &(SoilMap[y][x]),
                          &(SnowMap[y][x]), &(RadiationMap[y][x]),
                          &(EvapMap[y][x]), &(Total.Rad), &ChannelData, SkyViewMap);
          
  		    PrecipMap[y][x].SumPrecip += PrecipMap[y][x].Precip;
  		    PrecipMap[y][x].SnowAccum += PrecipMap[y][x].SnowFall;
  		    PrecipMap[y][x].SnowMelt += SnowMap[y][x].Outflow;
		    }
	    }
    }
    
#ifndef SNOW_ONLY
    
    RouteSubSurface(Time.Dt, &Map, TopoMap, VType, VegMap, Network, 
		    SType, SoilMap, &ChannelData, &Time, &Options, Dump.Path);
    
    if (Options.Extent != POINT)
      RouteChannel(&ChannelData, &Time, &Map, TopoMap, SoilMap, &Total, 
		   &Options, Network, SType, VType, VegMap, EvapMap, LType);
    
    if (Options.Extent == BASIN)
      RouteSurface(&Map, &Time, TopoMap, SoilMap, &Options,
        &Dump, VegMap, VType, LType, SType, &ChannelData, LocalMet.Tair, LocalMet.Rh);
    
#endif
    
    Aggregate(&Map, &Options, TopoMap, &Soil, &Veg, VegMap, EvapMap, PrecipMap,
	      RadiationMap, SnowMap, SoilMap, &Total, VType, Network, &ChannelData, Time.Dt, Time.NDaySteps);
    
    if (Options.SnowStats)
      SnowStats(&(Time.Current), &Map, &Options, TopoMap, SnowMap, Time.Dt);
    
    MassBalance(&(Time.Current), &(Time.Start), &(Dump.Balance), &Total, &Mass);
    
    ExecDump(&Map, &(Time.Current), &(Time.Start), &Options, &Dump, TopoMap,
             EvapMap, RadiationMap, PrecipMap, SnowMap, VegMap, &Veg,
             SoilMap, Network, &ChannelData, &Soil, &Total);
    
    IncreaseTime(&Time);
	  t += 1;
  } /* End of calculation loop over time steps */
  
  ExecDump(&Map, &(Time.Current), &(Time.Start), &Options, &Dump, TopoMap,
	   EvapMap, RadiationMap, PrecipMap, SnowMap, VegMap, &Veg, SoilMap,
	   Network, &ChannelData, &Soil, &Total);
  
#ifndef SNOW_ONLY
  FinalMassBalance(&(Dump.FinalBalance), &Total, &Mass, &Options);
#endif
  
  printf("\nEND OF MODEL RUN\n\n");
  
  /* Record the total simulation run time */
  finish1 = clock();
  runtime = (finish1-start)/CLOCKS_PER_SEC;
  printf("***********************************************************************************");
  printf("\nRuntime Summary:\n");
  printf("%6.2f hours elapsed for the simulation period of %d hours (%.1f days) \n", 
	  runtime/3600, t*Time.Dt/3600, (float)t*Time.Dt/3600/24);
  
  return EXIT_SUCCESS;
}
