
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "fileio.h"
#include "functions.h"
#include "constants.h"
#include "sizeofnt.h"
#include "varid.h"
#include "channel_grid.h"

 /*****************************************************************************
   StoreModelState()

   Store the current state of the model.

   The state variables for DHSVM include the following variables:

     - Canopy interception for each vegetation layer

     - Snow pack conditions:
       - presence/absence
       - number of days since last snowfall (used in albedo calculation)
       - snow water equivalent
       - for each layer of the snow pack:
         - liquid water content
         - temperature
       - cold content

     - Soil conditions:
       - for each soil layer:
         - soil moisture (also for the layer below the deepest root zone)
         - temperature
       - surface temperature
       - ground heat storage

     - Restart memory, written after the maps above and read back by
       InitModelState() when present, so that a run restarted from a dump
       reproduces the continuous run (layout in InitModelState.c)
 *****************************************************************************/
void StoreModelState(char *Path, DATE * Current, MAPSIZE * Map,
  OPTIONSTRUCT * Options, TOPOPIX ** TopoMap,
  PRECIPPIX ** PrecipMap, SNOWPIX ** SnowMap,
  VEGPIX ** VegMap, 
  LAYER * Veg, SOILPIX ** SoilMap, LAYER * Soil, 
  NETSTRUCT ** Network, CHANNEL * ChannelData, LAKETABLE *LType)
{
  const char *Routine = "StoreModelState";
  char Str[NAMESIZE + 1];
  char FileLabel[MAXSTRING + 1];
  char FileName[BUFSIZE * 2 + 1];
  int i, j;			/* counter */
  int NRecords;			/* Most stream map records in one cell */
  ChannelMapPtr Cell;
  int x;			/* counter */
  int y;			/* counter */
  int NSoil;			/* Number of soil layers for current pixel */
  int NVeg;			/* Number of veg layers for current pixel */
  MAPDUMP DMap;			/* Dump Info */
  void *Array;
  
  printf("Storing model state\n");
  PrintDate(Current, stdout);
  printf("\n");
  
  /* Store the canopy interception */

  sprintf(Str, "%02d.%02d.%04d.%02d.%02d.%02d", Current->Month, Current->Day,
    Current->Year, Current->Hour, Current->Min, Current->Sec);
  snprintf(FileName, sizeof(FileName), "%sInterception.State.%s%s", Path, Str, fileext);
  strcpy(FileLabel, "Interception storage for each vegetation layer");

  CreateMapFile(FileName, FileLabel, Map);

  if (!(Array = (float *)calloc(Map->NY * Map->NX, sizeof(float))))
    ReportError((char *)Routine, 1);

  for (i = 0; i < Veg->MaxLayers; i++) {
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
        if (INBASIN(TopoMap[y][x].Mask)) {
          NVeg = Veg->NLayers[(VegMap[y][x].Veg - 1)];
          if (i < NVeg)
            ((float *)Array)[y * Map->NX + x] = PrecipMap[y][x].IntRain[i];
          else
            ((float *)Array)[y * Map->NX + x] = NA;
        }
        else
          ((float *)Array)[y * Map->NX + x] = NA;
      }
    }
    DMap.ID = 202;
    DMap.Layer = i;
    DMap.Resolution = MAP_OUTPUT;
    strcpy(DMap.FileName, "");
    GetVarAttr(&DMap);
    Write2DMatrix(FileName, Array, DMap.NumberType, Map, &DMap, 0);
  }

  for (i = 0; i < Veg->MaxLayers; i++) {
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
        if (INBASIN(TopoMap[y][x].Mask)) {
          NVeg = Veg->NLayers[(VegMap[y][x].Veg - 1)];
          if (i < NVeg)
            ((float *)Array)[y * Map->NX + x] = PrecipMap[y][x].IntSnow[i];
          else
            ((float *)Array)[y * Map->NX + x] = NA;
        }
        else
          ((float *)Array)[y * Map->NX + x] = NA;
      }
    }
    DMap.ID = 203;
    DMap.Layer = i;
    DMap.Resolution = MAP_OUTPUT;
    strcpy(DMap.FileName, "");
    GetVarAttr(&DMap);
    Write2DMatrix(FileName, Array, DMap.NumberType, Map, &DMap, 0);
  }

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
        ((float *)Array)[y * Map->NX + x] = PrecipMap[y][x].TempIntStorage;
      }
      else {
        ((float *)Array)[y * Map->NX + x] = NA;
      }
    }
  }
  DMap.ID = 204;
  DMap.Resolution = MAP_OUTPUT;
  strcpy(DMap.FileName, "");
  GetVarAttr(&DMap);
  Write2DMatrix(FileName, Array, DMap.NumberType, Map, &DMap, 0);

  /* Restart memory: canopy temperature (C) */
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask))
        ((float *)Array)[y * Map->NX + x] = VegMap[y][x].Tcanopy;
      else
        ((float *)Array)[y * Map->NX + x] = NA;
    }
  }
  Write2DMatrix(FileName, Array, NC_FLOAT, Map, &DMap, 0);

  free(Array);

  /* Store the snow pack conditions */

  snprintf(FileName, sizeof(FileName), "%sSnow.State.%s%s", Path, Str, fileext);
  strcpy(FileLabel, "Snow pack moisture and temperature state");
  CreateMapFile(FileName, FileLabel, Map);

  if (!(Array = (float *)calloc(Map->NY * Map->NX, sizeof(float))))
    ReportError((char *)Routine, 1);

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask))
        ((float *)Array)[y * Map->NX + x] = (float)SnowMap[y][x].HasSnow;
      else
        ((float *)Array)[y * Map->NX + x] = NA;
    }
  }
  DMap.ID = 401;
  DMap.Resolution = MAP_OUTPUT;
  strcpy(DMap.FileName, "");
  GetVarAttr(&DMap);
  Write2DMatrix(FileName, Array, DMap.NumberType, Map, &DMap, 0);

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask))
        ((float *)Array)[y * Map->NX + x] = SnowMap[y][x].LastSnow;
      else
        ((float *)Array)[y * Map->NX + x] = NA;
    }
  }
  DMap.ID = 403;
  DMap.Resolution = MAP_OUTPUT;
  strcpy(DMap.FileName, "");
  GetVarAttr(&DMap);
  Write2DMatrix(FileName, Array, DMap.NumberType, Map, &DMap, 0);

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask))
        ((float *)Array)[y * Map->NX + x] = SnowMap[y][x].Swq;
      else
        ((float *)Array)[y * Map->NX + x] = NA;
    }
  }
  DMap.ID = 404;
  DMap.Resolution = MAP_OUTPUT;
  strcpy(DMap.FileName, "");
  GetVarAttr(&DMap);
  Write2DMatrix(FileName, Array, DMap.NumberType, Map, &DMap, 0);

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask))
        ((float *)Array)[y * Map->NX + x] = SnowMap[y][x].PackWater;
      else
        ((float *)Array)[y * Map->NX + x] = NA;
    }
  }
  DMap.ID = 406;
  DMap.Resolution = MAP_OUTPUT;
  strcpy(DMap.FileName, "");
  GetVarAttr(&DMap);
  Write2DMatrix(FileName, Array, DMap.NumberType, Map, &DMap, 0);

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask))
        ((float *)Array)[y * Map->NX + x] = SnowMap[y][x].TPack;
      else
        ((float *)Array)[y * Map->NX + x] = NA;
    }
  }
  DMap.ID = 407;
  DMap.Resolution = MAP_OUTPUT;
  strcpy(DMap.FileName, "");
  GetVarAttr(&DMap);
  Write2DMatrix(FileName, Array, DMap.NumberType, Map, &DMap, 0);

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask))
        ((float *)Array)[y * Map->NX + x] = SnowMap[y][x].SurfWater;
      else
        ((float *)Array)[y * Map->NX + x] = NA;
    }
  }
  DMap.ID = 408;
  DMap.Resolution = MAP_OUTPUT;
  strcpy(DMap.FileName, "");
  GetVarAttr(&DMap);
  Write2DMatrix(FileName, Array, DMap.NumberType, Map, &DMap, 0);

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask))
        ((float *)Array)[y * Map->NX + x] = SnowMap[y][x].TSurf;
      else
        ((float *)Array)[y * Map->NX + x] = NA;
    }
  }
  DMap.ID = 409;
  DMap.Resolution = MAP_OUTPUT;
  strcpy(DMap.FileName, "");
  GetVarAttr(&DMap);
  Write2DMatrix(FileName, Array, DMap.NumberType, Map, &DMap, 0);

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask))
        ((float *)Array)[y * Map->NX + x] = SnowMap[y][x].ColdContent;
      else
        ((float *)Array)[y * Map->NX + x] = NA;
    }
  }
  DMap.ID = 410;
  DMap.Resolution = MAP_OUTPUT;
  strcpy(DMap.FileName, "");
  GetVarAttr(&DMap);
  Write2DMatrix(FileName, Array, DMap.NumberType, Map, &DMap, 0);

  /* Restart memory: accumulation season flag (0 or 1) and albedo of snow */
  for (i = 0; i < 2; i++) {
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
        if (INBASIN(TopoMap[y][x].Mask))
          ((float *)Array)[y * Map->NX + x] = (i == 0) ? (float)SnowMap[y][x].AccumSeason : SnowMap[y][x].Albedo;
        else
          ((float *)Array)[y * Map->NX + x] = NA;
      }
    }
    Write2DMatrix(FileName, Array, NC_FLOAT, Map, &DMap, 0);
  }

  free(Array);

  /* Store the soil conditions */

  snprintf(FileName, sizeof(FileName), "%sSoil.State.%s%s", Path, Str, fileext);
  strcpy(FileLabel, "Soil moisture and temperature state");
  CreateMapFile(FileName, FileLabel, Map);

  if (!(Array = (float *)calloc(Map->NY * Map->NX, sizeof(float))))
    ReportError((char *)Routine, 1);

  for (i = 0; i < Soil->MaxLayers + 1; i++) {
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
        if (INBASIN(TopoMap[y][x].Mask)) {
          NSoil = Soil->NLayers[(SoilMap[y][x].Soil - 1)];
          if (i <= NSoil)
            ((float *)Array)[y * Map->NX + x] = SoilMap[y][x].Moist[i];
          else
            ((float *)Array)[y * Map->NX + x] = NA;
        }
        else
          ((float *)Array)[y * Map->NX + x] = NA;
      }
    }
    DMap.ID = 501;
    DMap.Layer = i;
    DMap.Resolution = MAP_OUTPUT;
    strcpy(DMap.FileName, "");
    GetVarAttr(&DMap);
    Write2DMatrix(FileName, Array, DMap.NumberType, Map, &DMap, 0);
  }

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask))
        ((float *)Array)[y * Map->NX + x] = SoilMap[y][x].TSurf;
      else
        ((float *)Array)[y * Map->NX + x] = SoilMap[y][x].TSurf;
    }
  }
  DMap.ID = 505;
  DMap.Resolution = MAP_OUTPUT;
  strcpy(DMap.FileName, "");
  GetVarAttr(&DMap);
  Write2DMatrix(FileName, Array, DMap.NumberType, Map, &DMap, 0);

  for (i = 0; i < Soil->MaxLayers; i++) {
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
        if (INBASIN(TopoMap[y][x].Mask)) {
          NSoil = Soil->NLayers[SoilMap[y][x].Soil - 1];
          if (i < NSoil)
            ((float *)Array)[y * Map->NX + x] = SoilMap[y][x].Temp[i];
          else
            ((float *)Array)[y * Map->NX + x] = NA;
        }
        else
          ((float *)Array)[y * Map->NX + x] = NA;
      }
    }
    DMap.ID = 511;
    DMap.Layer = i;
    DMap.Resolution = MAP_OUTPUT;
    strcpy(DMap.FileName, "");
    GetVarAttr(&DMap);
    Write2DMatrix(FileName, Array, DMap.NumberType, Map, &DMap, 0);
  }

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask))
        ((float *)Array)[y * Map->NX + x] = SoilMap[y][x].Qst;
      else
        ((float *)Array)[y * Map->NX + x] = NA;
    }
  }
  DMap.ID = 510;
  DMap.Resolution = MAP_OUTPUT;
  strcpy(DMap.FileName, "");
  GetVarAttr(&DMap);
  Write2DMatrix(FileName, Array, DMap.NumberType, Map, &DMap, 0);

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
        ((float *)Array)[y * Map->NX + x] = SoilMap[y][x].IExcess;
      }
      else
        ((float *)Array)[y * Map->NX + x] = NA;
    }
  }
  DMap.ID = 512;
  DMap.Resolution = MAP_OUTPUT;
  strcpy(DMap.FileName, "");
  GetVarAttr(&DMap);
  Write2DMatrix(FileName, Array, DMap.NumberType, Map, &DMap, 0);

  /* Restart memory: percolation from each root layer in the last time step
     (averaged with the next one in UnsaturatedFlow) */
  for (i = 0; i < Soil->MaxLayers; i++) {
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
        ((float *)Array)[y * Map->NX + x] = NA;
        if (INBASIN(TopoMap[y][x].Mask)) {
          NSoil = Soil->NLayers[SoilMap[y][x].Soil - 1];
          if (i < NSoil)
            ((float *)Array)[y * Map->NX + x] = SoilMap[y][x].Perc[i];
        }
      }
    }
    Write2DMatrix(FileName, Array, NC_FLOAT, Map, &DMap, 0);
  }

  /* Restart memory: saturated flow still to be distributed, water table
     depth and level, dynamic infiltration, kinematic runoff, detention and
     lake outflow */
  for (j = 0; j < 10; j++) {
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
        ((float *)Array)[y * Map->NX + x] = NA;
        if (INBASIN(TopoMap[y][x].Mask)) {
          switch (j) {
          case 0: ((float *)Array)[y * Map->NX + x] = SoilMap[y][x].SatFlow; break;
          case 1: ((float *)Array)[y * Map->NX + x] = SoilMap[y][x].TableDepth; break;
          case 2: ((float *)Array)[y * Map->NX + x] = SoilMap[y][x].WaterLevel; break;
          case 3: ((float *)Array)[y * Map->NX + x] = SoilMap[y][x].InfiltAcc; break;
          case 4: ((float *)Array)[y * Map->NX + x] = SoilMap[y][x].MoistInit; break;
          case 5: ((float *)Array)[y * Map->NX + x] = (float)PrecipMap[y][x].PrecipStart; break;
          case 6: ((float *)Array)[y * Map->NX + x] = SoilMap[y][x].Runoff; break;
          case 7: ((float *)Array)[y * Map->NX + x] = SoilMap[y][x].startRunoff; break;
          case 8: ((float *)Array)[y * Map->NX + x] = SoilMap[y][x].DetentionStorage; break;
          case 9:
            if (TopoMap[y][x].LakeID > 0 && LType != NULL)
              ((float *)Array)[y * Map->NX + x] = LType[TopoMap[y][x].LakeID - 1].Outflow;
            break;
          }
        }
      }
    }
    Write2DMatrix(FileName, Array, NC_FLOAT, Map, &DMap, 0);
  }

  /* Restart memory: water table depth below the channel and water available
     for channel infiltration, for the k-th stream map record of each cell */
  NRecords = channel_grid_max_records(ChannelData->stream_map);
  for (j = 0; j < 2 * NRecords; j++) {
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
        ((float *)Array)[y * Map->NX + x] = NA;
        Cell = channel_grid_record(ChannelData->stream_map, x, y, j % NRecords);
        if (Cell != NULL)
          ((float *)Array)[y * Map->NX + x] = (j < NRecords) ? Cell->table_depth : Cell->avail_water;
      }
    }
    Write2DMatrix(FileName, Array, NC_FLOAT, Map, &DMap, 0);
  }

  free(Array);
}
