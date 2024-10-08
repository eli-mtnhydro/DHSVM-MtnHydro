/*
 * SUMMARY:      InitTerrainMaps() - Initialize terrain coverages
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Initialize terrain coverages
 * DESCRIP-END.
 * FUNCTIONS:    InitTerrainMaps()
 *               InitTopoMap()
 *               InitSoilMap()
 *               InitVegMap()
 * COMMENTS:
 * $Id: InitTerrainMaps.c,v 3.1 2013/2/3 00:08:33 Ning Exp $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "fileio.h"
#include "functions.h"
#include "constants.h"
#include "getinit.h"
#include "sizeofnt.h"
#include "slopeaspect.h"
#include "varid.h"
#include "assert.h"

 /*****************************************************************************
   InitTerrainMaps()
 *****************************************************************************/
void InitTerrainMaps(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map,
  LAYER *Soil, LAYER *Veg, TOPOPIX ***TopoMap, SOILTABLE *SType, SOILPIX ***SoilMap, 
  VEGTABLE *VType, VEGPIX ***VegMap, DYNAVEG *DVeg, WINDPIX ***WindMap)

{
  printf("\nInitializing terrain maps\n");

  InitTopoMap(Input, Options, Map, TopoMap);
  InitVegMap(Options, Input, Map, VegMap, VType, DVeg);
  InitSoilMap(Input, Options, Map, Soil, *TopoMap, SoilMap, SType, VegMap, VType);
  if (Options->CanopyGapping)
    InitCanopyGapMap(Options, Input, Map, Soil, Veg, VType, VegMap, SType, SoilMap);
  if (Options->WindDrift)
    InitWindMap(Options, Input, Map, *TopoMap, WindMap);
}

/*****************************************************************************
  InitTopoMap()
*****************************************************************************/
void InitTopoMap(LISTPTR Input, OPTIONSTRUCT * Options, MAPSIZE * Map,
  TOPOPIX *** TopoMap)
{
  const char *Routine = "InitTopoMap";
  char VarName[BUFSIZE + 1];	/* Variable name */
  int i;			/* Counter */
  int x;			/* Counter */
  int y;			/* Counter */
  int flag;         /* either or not reverse the matrix */
  int NumberType;		/* Number type of data set */
  unsigned char *Mask = NULL;	/* Basin mask */
  float *Elev;			/* Surface elevation */
  STRINIENTRY StrEnv[] = {
    {"TERRAIN", "DEM FILE", "", ""},
    {"TERRAIN", "BASIN MASK FILE", "", ""},
    {NULL, NULL, "", NULL}
  };

  /* Process the [TERRAIN] section in the input file */
  if (!(*TopoMap = (TOPOPIX **)calloc(Map->NY, sizeof(TOPOPIX *))))
    ReportError((char *)Routine, 1);
  for (y = 0; y < Map->NY; y++) {
    if (!((*TopoMap)[y] = (TOPOPIX *)calloc(Map->NX, sizeof(TOPOPIX))))
      ReportError((char *)Routine, 1);
  }

  /* Read the key-entry pairs from the input file */
  for (i = 0; StrEnv[i].SectionName; i++) {
    GetInitString(StrEnv[i].SectionName, StrEnv[i].KeyName, StrEnv[i].Default,
      StrEnv[i].VarStr, (unsigned long)BUFSIZE, Input);
    if (IsEmptyStr(StrEnv[i].VarStr))
      ReportError(StrEnv[i].KeyName, 51);
  }

  /* Read the elevation data from the DEM dataset */
  GetVarName(001, 0, VarName);
  GetVarNumberType(001, &NumberType);
  if (!(Elev = (float *)calloc(Map->NX * Map->NY,
    SizeOfNumberType(NumberType))))
    ReportError((char *)Routine, 1);

  flag = Read2DMatrix(StrEnv[demfile].VarStr, Elev, NumberType, Map, 0,
    VarName, 0);

  /* Assign the attributes to the map pixel */
  /* Reverse the matrix is flag = 1 & netcdf option is selected */
  if ((Options->FileFormat == NETCDF && flag == 0) || (Options->FileFormat == BIN)) {
    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*TopoMap)[y][x].Dem = Elev[i];
      }
    }
  }
  else if (Options->FileFormat == NETCDF && flag == 1) {
    for (y = Map->NY - 1, i = 0; y >= 0; y--) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*TopoMap)[y][x].Dem = Elev[i];
      }
    }
  }
  else ReportError((char *)Routine, 57);
  free(Elev);

  /* Read the mask */
  GetVarName(002, 0, VarName);
  GetVarNumberType(002, &NumberType);
  if (!(Mask = (unsigned char *)calloc(Map->NX * Map->NY,
    SizeOfNumberType(NumberType))))
    ReportError((char *)Routine, 1);
  flag = Read2DMatrix(StrEnv[maskfile].VarStr, Mask, NumberType, Map, 0,
    VarName, 0);

  if ((Options->FileFormat == NETCDF && flag == 0)
    || (Options->FileFormat == BIN))
  {
    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*TopoMap)[y][x].Mask = Mask[i];
      }
    }
  }
  else if (Options->FileFormat == NETCDF && flag == 1) {
    for (y = Map->NY - 1, i = 0; y >= 0; y--) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*TopoMap)[y][x].Mask = Mask[i];
      }
    }
  }
  else ReportError((char *)Routine, 57);
  free(Mask);


  /* Calculate slope, aspect, magnitude of subsurface flow gradient, and
     fraction of flow flowing in each direction based on the land surface
     slope. */
  ElevationSlopeAspect(Map, *TopoMap, Options->MultiFlowDir);

  /* After calculating the slopes and aspects for all the points, reset the
     mask if the model is to be run in point mode */
  if (Options->Extent == POINT) {
    for (y = 0; y < Map->NY; y++)
      for (x = 0; x < Map->NX; x++)
        (*TopoMap)[y][x].Mask = OUTSIDEBASIN;
    (*TopoMap)[Options->PointY][Options->PointX].Mask = (1 != OUTSIDEBASIN);
  }
  /* find out the minimum grid elevation of the basin */
  MINELEV = 9999;
  for (y = 0, i = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++, i++) {
      if (INBASIN((*TopoMap)[y][x].Mask)) {
      if ((*TopoMap)[y][x].Dem < MINELEV) {
        MINELEV = (*TopoMap)[y][x].Dem;
      }
      }
    }
  }
}

/*****************************************************************************
  InitSoilMap()
*****************************************************************************/
void InitSoilMap(LISTPTR Input, OPTIONSTRUCT * Options, MAPSIZE * Map,
  LAYER * Soil, TOPOPIX ** TopoMap, SOILPIX *** SoilMap, SOILTABLE * SType,
  VEGPIX *** VegMap, VEGTABLE * VType)
{
  const char *Routine = "InitSoilMap";
  char VarName[BUFSIZE + 1];	/* Variable name */
  int i;			/* counter */
  int x;			/* counter */
  int y;			/* counter */
  int NumberType;		/* number type */
  unsigned char *Type;		/* Soil Type */
  float *Depth;			/* Soil Depth */
  float *KsLat = NULL;		/* Soil Lateral Conductivity */
  float *KsLatExp = NULL;		/* Soil Exponential Decrease */
  float *Porosity = NULL;	/* Soil Porosity */
  float *FC = NULL;		/* Soil Field Capacity */
  int flag;
  int NSet, VSet;
  int sidx;
  float LayerDepth, Transmissivity, KsVertCalc;
  
  STRINIENTRY StrEnv[] = {
    {"SOILS", "SOIL MAP FILE", "", ""},
    {"SOILS", "SOIL DEPTH FILE", "", ""},
    {"SOILS", "SOIL CONDUCTIVITY MAP FILE", "", "none"},
    {"SOILS", "SOIL EXPONENTIAL DECREASE MAP FILE", "", "none"},
    {"SOILS", "SOIL POROSITY MAP FILE", "", "none"},
    {"SOILS", "SOIL FIELD CAPACITY FILE", "", "none"},
    {NULL, NULL, "", NULL}
  };

  /* Process the filenames in the [SOILS] section in the input file */
  /* Assign the attributes to the correct map pixel */
  if (!(*SoilMap = (SOILPIX **)calloc(Map->NY, sizeof(SOILPIX *))))
    ReportError((char *)Routine, 1);
  for (y = 0; y < Map->NY; y++) {
    if (!((*SoilMap)[y] = (SOILPIX *)calloc(Map->NX, sizeof(SOILPIX))))
      ReportError((char *)Routine, 1);
  }

  /* Read the key-entry pairs from the input file */
  for (i = 0; StrEnv[i].SectionName; i++) {
    GetInitString(StrEnv[i].SectionName, StrEnv[i].KeyName, StrEnv[i].Default,
      StrEnv[i].VarStr, (unsigned long)BUFSIZE, Input);
    if (IsEmptyStr(StrEnv[i].VarStr))
      ReportError(StrEnv[i].KeyName, 51);
  }

  /* Read the soil type */
  GetVarName(003, 0, VarName);
  GetVarNumberType(003, &NumberType);
  if (!(Type = (unsigned char *)calloc(Map->NX * Map->NY,
    SizeOfNumberType(NumberType))))
    ReportError((char *)Routine, 1);
  flag = Read2DMatrix(StrEnv[soiltype_file].VarStr, Type, NumberType, 
	Map, 0, VarName, 0);

  if ((Options->FileFormat == NETCDF && flag == 0)
    || (Options->FileFormat == BIN))
  {
    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        if (((int)Type[i]) > Soil->NTypes)
          ReportError(StrEnv[soiltype_file].VarStr, 32);
        (*SoilMap)[y][x].Soil = Type[i];
      }
    }
  }
  else if (Options->FileFormat == NETCDF && flag == 1) {
    for (y = Map->NY - 1, i = 0; y >= 0; y--) {
      for (x = 0; x < Map->NX; x++, i++) {
        if (((int)Type[i]) > Soil->NTypes)
          ReportError(StrEnv[soiltype_file].VarStr, 32);
        (*SoilMap)[y][x].Soil = Type[i];
      }
    }
  }
  else ReportError((char *)Routine, 57);
  
  /******************************************************************/
  
  /* Read the total soil depth  */
  GetVarName(004, 0, VarName);
  GetVarNumberType(004, &NumberType);
  if (!(Depth = (float *)calloc(Map->NX * Map->NY,
    SizeOfNumberType(NumberType))))
    ReportError((char *)Routine, 1);
  flag = Read2DMatrix(StrEnv[soildepth_file].VarStr, Depth, NumberType, 
	Map, 0, VarName, 0);

  /* Assign the attributes to the correct map pixel */
  if ((Options->FileFormat == NETCDF && flag == 0)
    || (Options->FileFormat == BIN))
  {
    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*SoilMap)[y][x].Depth = Depth[i];
      }
    }
  }
  else if (Options->FileFormat == NETCDF && flag == 1) {
    for (y = Map->NY - 1, i = 0; y >= 0; y--) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*SoilMap)[y][x].Depth = Depth[i];
      }
    }
  }
  else ReportError((char *)Routine, 57);

  /******************************************************************/

  /* Read the spatial Lateral Conductivity map */
  GetVarName(012, 0, VarName);
  GetVarNumberType(012, &NumberType);

  if (strncmp(StrEnv[kslat_file].VarStr, "none", 4)) {
      printf("Spatial lateral conductivity map provided, reading map\n");
      if (!(KsLat = (float *)calloc(Map->NX * Map->NY,
        SizeOfNumberType(NumberType))))
        ReportError((char *)Routine, 1);
      flag = Read2DMatrix(StrEnv[kslat_file].VarStr, KsLat, NumberType, 
      Map, 0, VarName, 0);

      if ((Options->FileFormat == NETCDF && flag == 0)
        || (Options->FileFormat == BIN))
      {
        for (y = 0, i = 0; y < Map->NY; y++) {
          for (x = 0; x < Map->NX; x++, i++) {
            if (KsLat[i] > 0.0)
              (*SoilMap)[y][x].KsLat = KsLat[i]/1000.0;
            else
              (*SoilMap)[y][x].KsLat = SType[(*SoilMap)[y][x].Soil - 1].KsLat;
          }
        }
      }
      else if (Options->FileFormat == NETCDF && flag == 1) {
        for (y = Map->NY - 1, i = 0; y >= 0; y--) {
          for (x = 0; x < Map->NX; x++, i++) {
            if (KsLat[i] > 0.0)
              (*SoilMap)[y][x].KsLat = KsLat[i]/1000.0;
            else
              (*SoilMap)[y][x].KsLat = SType[(*SoilMap)[y][x].Soil - 1].KsLat;
          }
        }
      }
      else ReportError((char *)Routine, 57);
      free(KsLat);
      KsLat = NULL;
  }
  else{
    printf("Spatial lateral conductivity map not provided, generating map\n");
    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
          (*SoilMap)[y][x].KsLat = SType[(*SoilMap)[y][x].Soil - 1].KsLat;
      }
    }
  }
  
  /******************************************************************/
  
  /* Read the spatial exponential decrease map */
  GetVarName(016, 0, VarName);
  GetVarNumberType(016, &NumberType);
  
  if (strncmp(StrEnv[expdec_file].VarStr, "none", 4)) {
    printf("Spatial exponential decrease map provided, reading map\n");
    if (!(KsLatExp = (float *)calloc(Map->NX * Map->NY,
                      SizeOfNumberType(NumberType))))
      ReportError((char *)Routine, 1);
    flag = Read2DMatrix(StrEnv[expdec_file].VarStr, KsLatExp, NumberType, 
                        Map, 0, VarName, 0);
    
    if ((Options->FileFormat == NETCDF && flag == 0)
          || (Options->FileFormat == BIN))
    {
      for (y = 0, i = 0; y < Map->NY; y++) {
        for (x = 0; x < Map->NX; x++, i++) {
          (*SoilMap)[y][x].KsLatExp = KsLatExp[i];
        }
      }
    }
    else if (Options->FileFormat == NETCDF && flag == 1) {
      for (y = Map->NY - 1, i = 0; y >= 0; y--) {
        for (x = 0; x < Map->NX; x++, i++) {
          (*SoilMap)[y][x].KsLatExp = KsLatExp[i];
        }
      }
    }
    else ReportError((char *)Routine, 57);
    free(KsLatExp);
    KsLatExp = NULL;
  }
  else{
    printf("Spatial exponential decrease map not provided, generating map\n");
    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*SoilMap)[y][x].KsLatExp = SType[(*SoilMap)[y][x].Soil - 1].KsLatExp;
      }
    }
  }
  
  /******************************************************************/
  
  /* Allocate memory for vertical conductivity */
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (!((*SoilMap)[y][x].KsVert =
          (float *)calloc(Soil->MaxLayers + 1, sizeof(float *))))
        ReportError((char *)Routine, 1);
    }
  }
  
  /* Creating spatial layered vertical conductivity */
  for (NSet = 0; NSet < Soil->MaxLayers; NSet++) {
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
        if (INBASIN((TopoMap)[y][x].Mask)) {
          sidx = (*SoilMap)[y][x].Soil - 1;
          if (NSet < Soil->NLayers[sidx]) {
            if (Options->UseKsatAnisotropy &&
                VType[(*VegMap)[y][x].Veg - 1].NSoilLayers > NSet &&
                VType[(*VegMap)[y][x].Veg - 1].RootDepth[NSet] > 0.001) {
              
              /* Get total depth of current soil layer bottom */
              LayerDepth = 0.0;
              for (VSet = 0; VSet <= NSet; VSet++)
                LayerDepth += VType[(*VegMap)[y][x].Veg - 1].RootDepth[VSet];
              
              /* Calculate effective lateral conductivity of current soil layer */
              Transmissivity =
                CalcTransmissivity(LayerDepth,
                                   LayerDepth - VType[(*VegMap)[y][x].Veg - 1].RootDepth[NSet],
                                   (*SoilMap)[y][x].KsLat,
                                   (*SoilMap)[y][x].KsLatExp,
                                   SType[(*SoilMap)[y][x].Soil - 1].DepthThresh);
              KsVertCalc = Transmissivity / VType[(*VegMap)[y][x].Veg - 1].RootDepth[NSet];
              
              /* Account for user-supplied vertical anisotropy */
              KsVertCalc /= SType[(*SoilMap)[y][x].Soil - 1].KsAnisotropy;
              
              (*SoilMap)[y][x].KsVert[NSet] = KsVertCalc;
            }
            else
              (*SoilMap)[y][x].KsVert[NSet] = SType[sidx].Ks[NSet];
          }
        }            
      }
    }
  }
  
  /* Calculate deep layer vertical conductivity */
  NSet = Soil->MaxLayers;
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN((TopoMap)[y][x].Mask)) {
        sidx = (*SoilMap)[y][x].Soil - 1;
        if (NSet == Soil->NLayers[sidx]) {
          if (Options->UseKsatAnisotropy &&
              VType[(*VegMap)[y][x].Veg - 1].NSoilLayers > NSet &&
              ((*SoilMap)[y][x].Depth - VType[(*VegMap)[y][x].Veg - 1].TotalDepth) > 0.001) {
            
            /* Calculate effective lateral conductivity of current soil layer */
            Transmissivity =
              CalcTransmissivity((*SoilMap)[y][x].Depth,
                                 (*SoilMap)[y][x].Depth - VType[(*VegMap)[y][x].Veg - 1].TotalDepth,
                                 (*SoilMap)[y][x].KsLat,
                                 (*SoilMap)[y][x].KsLatExp,
                                 SType[(*SoilMap)[y][x].Soil - 1].DepthThresh);
            KsVertCalc = Transmissivity / ((*SoilMap)[y][x].Depth - VType[(*VegMap)[y][x].Veg - 1].TotalDepth);
            
            /* Account for user-supplied vertical anisotropy */
            KsVertCalc /= SType[(*SoilMap)[y][x].Soil - 1].KsAnisotropy;
            
            (*SoilMap)[y][x].KsVert[NSet] = KsVertCalc;
          }
          else
            (*SoilMap)[y][x].KsVert[NSet] = SType[sidx].Ks[NSet - 1];
        }
      }
    }
  }
  
  /* Copy surface layer KsLat to MaxInfiltrationRate map */
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN((TopoMap)[y][x].Mask)) {
        sidx = (*SoilMap)[y][x].Soil - 1;
        if (Options->UseKsatAnisotropy) {
          (*SoilMap)[y][x].MaxInfiltrationRate = (*SoilMap)[y][x].KsLat;
          (*SoilMap)[y][x].MaxInfiltrationRate /= SType[(*SoilMap)[y][x].Soil - 1].KsAnisotropy;
        }
        else
          (*SoilMap)[y][x].MaxInfiltrationRate = SType[sidx].MaxInfiltrationRate;
      }
    }
  }
  
  /******************************************************************/
  
  /* Read the spatial field capacity map */
  GetVarNumberType(015, &NumberType);

  /*Allocate memory*/  
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (!((*SoilMap)[y][x].FCap =
            (float *)calloc(Soil->MaxLayers + 1, sizeof(float *))))
        ReportError((char *)Routine, 1);
    }
  }

  /*Creating spatial layered field capacity*/
  if (strncmp(StrEnv[fc_file].VarStr, "none", 4)) {
    printf("Spatial soil field capacity provided, reading map\n"); 
    /*Read data layer by layer*/
    for (NSet = 0; NSet < Soil->MaxLayers; NSet++) {
      GetVarName(015, NSet, VarName);
      if (!(FC = (float *)calloc(Map->NX * Map->NY,
        SizeOfNumberType(NumberType))))
        ReportError((char *)Routine, 1);
      flag = Read2DMatrix(StrEnv[fc_file].VarStr, FC, NumberType, Map, NSet, VarName, 0);

      if ((Options->FileFormat == NETCDF && flag == 0)
        || (Options->FileFormat == BIN))
      {
        for (y = 0, i = 0; y < Map->NY; y++) {
          for (x = 0; x < Map->NX; x++, i++) {
            if (INBASIN((TopoMap)[y][x].Mask)) {
              sidx = (*SoilMap)[y][x].Soil - 1;
              if (NSet < Soil->NLayers[sidx]) {
                if (FC[i] > 0.0)
                  (*SoilMap)[y][x].FCap[NSet] = FC[i];
                else
                  (*SoilMap)[y][x].FCap[NSet] = SType[sidx].FCap[NSet];
                /*Make sure FCap larger than WP*/
                if (((*SoilMap)[y][x].FCap[NSet] < SType[sidx].WP[NSet]))
                  ReportError(SType[sidx].Desc, 11);
              }
            }            
          }
        }
      } else if (Options->FileFormat == NETCDF && flag == 1) {
        for (y = Map->NY - 1, i = 0; y >= 0; y--) {
          for (x = 0; x < Map->NX; x++, i++) {
            if (INBASIN((TopoMap)[y][x].Mask)) {
              sidx = (*SoilMap)[y][x].Soil - 1;
              if (NSet < Soil->NLayers[sidx]) {
                if (FC[i] > 0.0)
                  (*SoilMap)[y][x].FCap[NSet] = FC[i];
                else
                  (*SoilMap)[y][x].FCap[NSet] = SType[sidx].FCap[NSet];
                if (((*SoilMap)[y][x].FCap[NSet] <SType[sidx].WP[NSet]))
                  ReportError(SType[sidx].Desc, 11);
              }
            }
          }
        }
      }
      else ReportError((char *)Routine, 57);
    }
    free(FC);
    FC = NULL;
  }
  else{
    printf("Spatial soil field capacity map not provided, generating map\n"); 
    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        if (INBASIN((TopoMap)[y][x].Mask)) {
          /* FIXME: this assumes a valid soil type index */
          sidx = (*SoilMap)[y][x].Soil - 1;
          for (NSet = 0; NSet < Soil->NLayers[sidx]; NSet++) {
            (*SoilMap)[y][x].FCap[NSet] = SType[sidx].FCap[NSet];
            if (((*SoilMap)[y][x].FCap[NSet] <SType[sidx].WP[NSet]))
              ReportError(SType[sidx].Desc, 11);
          }
        }
      }
    }
  }
  
  /* Copy deep layer field capacity from deepest root layer */
  NSet = Soil->MaxLayers;
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN((TopoMap)[y][x].Mask)) {
        sidx = (*SoilMap)[y][x].Soil - 1;
        if (NSet == Soil->NLayers[sidx])
          (*SoilMap)[y][x].FCap[NSet] = (*SoilMap)[y][x].FCap[NSet - 1];
      }
    }
  }
  
  /******************************************************************/
  
  /* Read the spatial porosity map */
  GetVarNumberType(013, &NumberType);

  /*Allocate memory for porosity*/
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (!((*SoilMap)[y][x].Porosity =
            (float *)calloc(Soil->MaxLayers + 1, sizeof(float *))))
        ReportError((char *)Routine, 1);
    }
  }

  /*Creating spatial layered porosity*/
  if (strncmp(StrEnv[porosity_file].VarStr, "none", 4)) {
    printf("Spatial soil porosity map provided, reading map\n");   
    /*Read data layer by layer*/
    for (NSet = 0; NSet < Soil->MaxLayers; NSet++) {
      GetVarName(013, NSet, VarName);
      if (!(Porosity = (float *)calloc(Map->NX * Map->NY,
        SizeOfNumberType(NumberType))))
        ReportError((char *)Routine, 1);
      flag = Read2DMatrix(StrEnv[porosity_file].VarStr, Porosity, NumberType, Map, NSet, VarName, 0);

      if ((Options->FileFormat == NETCDF && flag == 0)
        || (Options->FileFormat == BIN))
      {
        for (y = 0, i = 0; y < Map->NY; y++) {
          for (x = 0; x < Map->NX; x++, i++) {
            if (INBASIN((TopoMap)[y][x].Mask)) {
              sidx = (*SoilMap)[y][x].Soil - 1;
              if (NSet < Soil->NLayers[sidx]) {
                if (Porosity[i] > 0.0)
                  (*SoilMap)[y][x].Porosity[NSet] = Porosity[i];
                else
                  (*SoilMap)[y][x].Porosity[NSet] = SType[sidx].Porosity[NSet];
                /*Make sure porosity larger than FCap and WP*/
                if (((*SoilMap)[y][x].Porosity[NSet] < (*SoilMap)[y][x].FCap[NSet])
                    || ((*SoilMap)[y][x].Porosity[NSet] < SType[sidx].WP[NSet]))
                  ReportError(SType[sidx].Desc, 11);
              }
            }            
          }
        }
      } else if (Options->FileFormat == NETCDF && flag == 1) {
        for (y = Map->NY - 1, i = 0; y >= 0; y--) {
          for (x = 0; x < Map->NX; x++, i++) {
            if (INBASIN((TopoMap)[y][x].Mask)) {
              sidx = (*SoilMap)[y][x].Soil - 1;
              if (NSet < Soil->NLayers[sidx]) {
                if (Porosity[i] > 0.0)
                  (*SoilMap)[y][x].Porosity[NSet] = Porosity[i];
                else
                  (*SoilMap)[y][x].Porosity[NSet] = SType[sidx].Porosity[NSet];
                /*Make sure porosity larger than FCap and WP*/
                if (((*SoilMap)[y][x].Porosity[NSet] < (*SoilMap)[y][x].FCap[NSet])
                    || ((*SoilMap)[y][x].Porosity[NSet] < SType[sidx].WP[NSet]))
                  ReportError(SType[sidx].Desc, 11);
              }
            }
          }
        }
      }
      else ReportError((char *)Routine, 57);
    }
    free(Porosity);
    Porosity = NULL;
  }
  else{
    printf("Spatial soil porosity map not provided, generating map\n"); 
    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        if (INBASIN((TopoMap)[y][x].Mask)) {
          /* FIXME: this assumes a valid soil type index */
          sidx = (*SoilMap)[y][x].Soil - 1;
          for (NSet = 0; NSet < Soil->NLayers[sidx]; NSet++) {
            (*SoilMap)[y][x].Porosity[NSet] = SType[sidx].Porosity[NSet];
            /*Make sure porosity larger than FCap and WP*/
            if (((*SoilMap)[y][x].Porosity[NSet] < (*SoilMap)[y][x].FCap[NSet])
                || ((*SoilMap)[y][x].Porosity[NSet] < SType[sidx].WP[NSet]))
              ReportError(SType[sidx].Desc, 11);
          }
        }
      }
    }
  }
  
  /* Copy deep layer porosity from deepest root layer */
  NSet = Soil->MaxLayers;
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN((TopoMap)[y][x].Mask)) {
        sidx = (*SoilMap)[y][x].Soil - 1;
        if (NSet == Soil->NLayers[sidx])
          (*SoilMap)[y][x].Porosity[NSet] = (*SoilMap)[y][x].Porosity[NSet - 1];
      }
    }
  }

   /******************************************************************/
   /******************************************************************/

  for (y = 0, i = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++, i++) {
      if (Options->Infiltration == DYNAMIC)
        (*SoilMap)[y][x].InfiltAcc = 0.;
      (*SoilMap)[y][x].MoistInit = 0.;

      /* allocate memory for the number of root layers, plus an additional
       layer below the deepest root layer */
      if (INBASIN(TopoMap[y][x].Mask)) {
        if (!((*SoilMap)[y][x].Moist =
          (float *)calloc((Soil->NLayers[Type[i] - 1] + 1), sizeof(float))))
          ReportError((char *)Routine, 1);
        if (!((*SoilMap)[y][x].Perc =
          (float *)calloc(Soil->NLayers[Type[i] - 1], sizeof(float))))
          ReportError((char *)Routine, 1);
        if (!((*SoilMap)[y][x].Temp =
          (float *)calloc(Soil->NLayers[Type[i] - 1], sizeof(float))))
          ReportError((char *)Routine, 1);
      }
      else {
        (*SoilMap)[y][x].Moist = NULL;
        (*SoilMap)[y][x].Perc = NULL;
        (*SoilMap)[y][x].Temp = NULL;
      }
    }
  }
  free(Type);
  free(Depth);
}

/*****************************************************************************
  InitVegMap()
*****************************************************************************/
void InitVegMap(OPTIONSTRUCT * Options, LISTPTR Input, MAPSIZE * Map, VEGPIX *** VegMap,
                VEGTABLE *VType, DYNAVEG *DVeg)
{
  const char *Routine = "InitVegMap";
  char VarName[BUFSIZE + 1];
  //char VegMapFileName[BUFSIZE + 1];
  int i;			/* counter */
  int x;			/* counter */
  int y;			/* counter */
  int j;			/* counter */
  int flag;
  int NumberType;		/* number type */
  unsigned char *Type;		/* Vegetation type */
  float *FC = NULL;		/* Vegetation Fractional Coverage */
  float *LAIMonthly= NULL; /* Vegetation Leaf Area Index, monthly */
  float *Height = NULL; /*Vegetaion Height*/
  int NSet; /*Counter for LAI map month*/

  /* Get the map filename from the [VEGETATION] section */
  STRINIENTRY StrEnv[] = {
    {"VEGETATION", "VEGETATION MAP FILE", "", ""},
    {"VEGETATION", "VEGETATION FC MAP FILE", "", "none"},
    {"VEGETATION", "VEGETATION LAI MAP FILE", "", "none"},
    {"VEGETATION", "VEGETATION HEIGHT MAP FILE", "", "none"},
    {"VEGETATION", "DYNAMIC VEGETATION MAP PATH", "", "none"},
    {"VEGETATION", "NUMBER OF DYNAMIC VEGETATION MAPS", "","none"},
    {NULL, NULL, "", NULL}
  };

  /* Assign the attributes to the correct map pixel */
  if (!(*VegMap = (VEGPIX **)calloc(Map->NY, sizeof(VEGPIX *))))
    ReportError((char *)Routine, 1);
  for (y = 0; y < Map->NY; y++) {
    if (!((*VegMap)[y] = (VEGPIX *)calloc(Map->NX, sizeof(VEGPIX))))
      ReportError((char *)Routine, 1);
  }

  /* Read the key-entry pairs from the input file */
  for (i = 0; StrEnv[i].SectionName; i++) {
    GetInitString(StrEnv[i].SectionName, StrEnv[i].KeyName, StrEnv[i].Default,
      StrEnv[i].VarStr, (unsigned long)BUFSIZE, Input);
    if (IsEmptyStr(StrEnv[i].VarStr))
      ReportError(StrEnv[i].KeyName, 51);
  }

/**********************************/
  if (Options->DynamicVeg){
    printf("Warning: Dynamic Vegetation Mode, not compataible with SNOTEL option\n");

    /*Read from Config of dynamic veg path and dates*/
    if (IsEmptyStr(StrEnv[dynaveg_num].VarStr))  
      ReportError(StrEnv[dynaveg_num].KeyName, 51);
    else if (!CopyInt(&(DVeg->NUpdate), StrEnv[dynaveg_num].VarStr, 1))
      ReportError(StrEnv[dynaveg_num].KeyName, 51);

    //printf("number of dyna vegs %d",&(DVeg->NUpdate));

    if (IsEmptyStr(StrEnv[dynaveg_path].VarStr)) 
      ReportError(StrEnv[dynaveg_path].KeyName, 51);
    strcpy(DVeg->DynaVegPath, StrEnv[dynaveg_path].VarStr);

    /*Initiate the information*/
    InitVegUpdate(Input, DVeg->NUpdate, &(DVeg->DUpdate));
  }

    /******************************/
  
  /* Read the vegetation type */
  GetVarName(005, 0, VarName);
  GetVarNumberType(005, &NumberType);
  if (!(Type = (unsigned char *)calloc(Map->NX * Map->NY,
    SizeOfNumberType(NumberType))))
    ReportError((char *)Routine, 1);
  flag = Read2DMatrix(StrEnv[vegtype_file].VarStr, Type, NumberType, Map, 0, VarName, 0);
  
  if ((Options->FileFormat == NETCDF && flag == 0)
    || (Options->FileFormat == BIN))
  {
    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*VegMap)[y][x].Veg = Type[i];
        (*VegMap)[y][x].Tcanopy = 0.0;
      }
    }
  }
  else if (Options->FileFormat == NETCDF && flag == 1) {
    for (y = Map->NY - 1, i = 0; y >= 0; y--) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*VegMap)[y][x].Veg = Type[i];
        (*VegMap)[y][x].Tcanopy = 0.0;
      }
    }
  }
  else ReportError((char *)Routine, 57);

  free(Type);

  /* Read the vegetation fractional coverage map */
  GetVarName(010, 0, VarName);
  GetVarNumberType(010, &NumberType);

  if (strncmp(StrEnv[vegfc_file].VarStr, "none", 4)) {
    printf("Spatial fractional cover map provided, reading FC from map\n");
    if (!(FC = (float *)calloc(Map->NX * Map->NY,
      SizeOfNumberType(NumberType))))
      ReportError((char *)Routine, 1);
    flag = Read2DMatrix(StrEnv[vegfc_file].VarStr, FC, NumberType, Map, 0, VarName, 0);

    if ((Options->FileFormat == NETCDF && flag == 0)
      || (Options->FileFormat == BIN))
    {
      for (y = 0, i = 0; y < Map->NY; y++) {
        for (x = 0; x < Map->NX; x++, i++) {
          /*Allocate Memory*/
          if (!((*VegMap)[y][x].Fract = (float *)calloc(VType[(*VegMap)[y][x].Veg - 1].NVegLayers, sizeof(float))))
            ReportError((char *)Routine, 1);
          if ( VType[(*VegMap)[y][x].Veg - 1].OverStory == TRUE) {
            if (FC[i] > 0.0)
              (*VegMap)[y][x].Fract[0] = FC[i];
            else
              (*VegMap)[y][x].Fract[0] = VType[(*VegMap)[y][x].Veg - 1].Fract[0];
            /*If understory exists, set default understory FC=1.0*/
            if (VType[(*VegMap)[y][x].Veg - 1].UnderStory == TRUE)
              (*VegMap)[y][x].Fract[1] = 1.0;
          }
          else{
            if (VType[(*VegMap)[y][x].Veg - 1].UnderStory == TRUE)
              (*VegMap)[y][x].Fract[0] = 1.0;
          }

        }
      }
    }
    else if (Options->FileFormat == NETCDF && flag == 1) {
      for (y = Map->NY - 1, i = 0; y >= 0; y--) {
        for (x = 0; x < Map->NX; x++, i++) {
          /*Allocate memory*/
          if (!((*VegMap)[y][x].Fract = (float *)calloc(VType[(*VegMap)[y][x].Veg - 1].NVegLayers, sizeof(float))))
            ReportError((char *)Routine, 1);

          if ( VType[(*VegMap)[y][x].Veg - 1].OverStory == TRUE) {  
            if (FC[i] > 0.0)
              (*VegMap)[y][x].Fract[0] = FC[i];
            else
            /* If value from the fractional cover map is NaN, then read value from attribute table*/
              (*VegMap)[y][x].Fract[0] = VType[(*VegMap)[y][x].Veg - 1].Fract[0];
            /*If understory exists, set default understory FC=1.0*/
            if (VType[(*VegMap)[y][x].Veg - 1].UnderStory == TRUE)
              (*VegMap)[y][x].Fract[1] = 1.0;
          }
          else{
            if ( VType[(*VegMap)[y][x].Veg - 1].UnderStory == TRUE)
              (*VegMap)[y][x].Fract[0] = 1.0;	   
          }
        }
      }
    }
    else ReportError((char *)Routine, 57);
    free(FC);
  }
  else{
    printf("Vegetation fractional coverage created from vegetation table\n");
    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
          /*Allocate Memory*/
          if (!((*VegMap)[y][x].Fract = (float *)calloc(VType[(*VegMap)[y][x].Veg - 1].NVegLayers, sizeof(float))))
            ReportError((char *)Routine, 1);

          if ( VType[(*VegMap)[y][x].Veg - 1].OverStory == TRUE) {
              (*VegMap)[y][x].Fract[0] = VType[(*VegMap)[y][x].Veg - 1].Fract[0];
            /*If understory exists, set default understory FC=1.0*/
            if (VType[(*VegMap)[y][x].Veg - 1].UnderStory == TRUE)
              (*VegMap)[y][x].Fract[1] = 1.0;
          }
          else{
            if (VType[(*VegMap)[y][x].Veg - 1].UnderStory == TRUE)
              (*VegMap)[y][x].Fract[0] = 1.0;
          }
        }
      }
  }

  /*Calculate Vf */
  if (Options->ImprovRadiation == TRUE) {
    for (y = 0, i = 0; y < Map->NY; y++) {
        for (x = 0; x < Map->NX; x++, i++) {
          if ( VType[(*VegMap)[y][x].Veg - 1].NVegLayers >0) 
            (*VegMap)[y][x].Vf = (*VegMap)[y][x].Fract[0] * VType[(*VegMap)[y][x].Veg - 1].VfAdjust;
        }
    }
  }

  /* Read the vegetation LAI map */
  /*Instead of reading LAI data by month, read data all together as Bill suggested*/
    
  GetVarName(011, 0, VarName);
  GetVarNumberType(011, &NumberType);
 
  if (strncmp(StrEnv[veglai_file].VarStr, "none", 4)) {
    printf("Spatial LAI provided, reading LAI from map\n");
    /*Allocate Memory: if FC file avaiable, assume max 2 layers of vegtation*/  
      for (y = 0; y < Map->NY; y++) {
        for (x = 0; x < Map->NX; x++) {

          if (!((*VegMap)[y][x].LAIMonthly = (float **)calloc(VType[(*VegMap)[y][x].Veg - 1].NVegLayers , sizeof(float *))))
            ReportError((char *)Routine, 1);
          for (j = 0; j < VType[(*VegMap)[y][x].Veg - 1].NVegLayers; j++) {
              if (!((*VegMap)[y][x].LAIMonthly[j] = (float *)calloc(12, sizeof(float))))
              ReportError((char *)Routine, 1);
              }
        }
      }
   
  /*Read data monthy by month*/
  for (NSet = 0; NSet < 12; NSet++) {   
    if (!(LAIMonthly = (float *)calloc(Map->NX * Map->NY,
      SizeOfNumberType(NumberType))))
      ReportError((char *)Routine, 1);
    flag = Read2DMatrix(StrEnv[veglai_file].VarStr, LAIMonthly, NumberType, Map, NSet, VarName, 0);
    
    printf("begining month %d\n",NSet);
    
    if ((Options->FileFormat == NETCDF && flag == 0)
      || (Options->FileFormat == BIN))
    {
      for (y = 0, i = 0; y < Map->NY; y++) {
        for (x = 0; x < Map->NX; x++, i++) {
          if ( VType[(*VegMap)[y][x].Veg - 1].OverStory == TRUE) {
           if (LAIMonthly[i] > 0.0)
              (*VegMap)[y][x].LAIMonthly[0][NSet] = LAIMonthly[i];
            else
              (*VegMap)[y][x].LAIMonthly[0][NSet] = VType[(*VegMap)[y][x].Veg - 1].LAIMonthly[0][NSet];
           
            if ( VType[(*VegMap)[y][x].Veg - 1].UnderStory  == TRUE )
              (*VegMap)[y][x].LAIMonthly[1][NSet] = VType[(*VegMap)[y][x].Veg - 1].LAIMonthly[1][NSet];
          }
          else{
            if ( VType[(*VegMap)[y][x].Veg - 1].UnderStory == TRUE)
              (*VegMap)[y][x].LAIMonthly[0][NSet] = VType[(*VegMap)[y][x].Veg - 1].LAIMonthly[0][NSet];
          }
        }
      }
    }
    else if (Options->FileFormat == NETCDF && flag == 1) {
      for (y = Map->NY - 1, i = 0; y >= 0; y--) {
        for (x = 0; x < Map->NX; x++, i++) {
        
          if ( VType[(*VegMap)[y][x].Veg - 1].OverStory == TRUE) {
            if (LAIMonthly[i] > 0.0)
              (*VegMap)[y][x].LAIMonthly[0][NSet] = LAIMonthly[i];
            else
              (*VegMap)[y][x].LAIMonthly[0][NSet] = VType[(*VegMap)[y][x].Veg - 1].LAIMonthly[0][NSet];
            /*If understory exists, set default understory FC=1.0*/
            if (VType[(*VegMap)[y][x].Veg - 1].UnderStory == TRUE)
              (*VegMap)[y][x].LAIMonthly[1][NSet] = VType[(*VegMap)[y][x].Veg - 1].LAIMonthly[1][NSet];
          }
          else{
            if (VType[(*VegMap)[y][x].Veg - 1].UnderStory == TRUE)
              (*VegMap)[y][x].LAIMonthly[0][NSet] = VType[(*VegMap)[y][x].Veg - 1].LAIMonthly[0][NSet];
          }
        }
      }
    }
    else ReportError((char *)Routine, 57);  

    free(LAIMonthly);     
  }   
  }
  else{
    printf("No spatial LAI provided, generating from vegetation table\n");

      for (y = 0; y < Map->NY; y++) {
        for (x = 0; x < Map->NX; x++) {

          if (!((*VegMap)[y][x].LAIMonthly = (float **)calloc(VType[(*VegMap)[y][x].Veg - 1].NVegLayers , sizeof(float *))))
            ReportError((char *)Routine, 1);
          for (j = 0; j < VType[(*VegMap)[y][x].Veg - 1].NVegLayers; j++) {
              if (!((*VegMap)[y][x].LAIMonthly[j] = (float *)calloc(12, sizeof(float))))
              ReportError((char *)Routine, 1);
              }
        }
      }

    for (NSet = 0; NSet < 12; NSet++) {
      for (y = 0, i = 0; y < Map->NY; y++) {
          for (x = 0; x < Map->NX; x++, i++) {

            if ( VType[(*VegMap)[y][x].Veg - 1].OverStory == TRUE) {        
              (*VegMap)[y][x].LAIMonthly[0][NSet] = VType[(*VegMap)[y][x].Veg - 1].LAIMonthly[0][NSet];
              if ( VType[(*VegMap)[y][x].Veg - 1].UnderStory == TRUE ){
                (*VegMap)[y][x].LAIMonthly[1][NSet] = VType[(*VegMap)[y][x].Veg - 1].LAIMonthly[1][NSet]; 
              }
            }
            else{
              if ( VType[(*VegMap)[y][x].Veg - 1].UnderStory == TRUE){
                (*VegMap)[y][x].LAIMonthly[0][NSet] = VType[(*VegMap)[y][x].Veg - 1].LAIMonthly[0][NSet];
              }
            }
          }
        }
    }
  }
  /*Need to be careful here about the MaxInt*/
  for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
      /*Allocate memory to LAI values*/
      if (!((*VegMap)[y][x].LAI = (float *)calloc(VType[(*VegMap)[y][x].Veg - 1].NVegLayers, sizeof(float))))
        ReportError((char *)Routine, 1);
      if (!((*VegMap)[y][x].MaxInt = (float *)calloc(VType[(*VegMap)[y][x].Veg - 1].NVegLayers, sizeof(float))))
        ReportError((char *)Routine, 1);
    }
  }

  /*Read Tree Height Map*/
   /* Read the vegetation fractional coverage map */
  GetVarName(014, 0, VarName);
  GetVarNumberType(014, &NumberType);

  /*Allcate Memory*/
  for (y = 0, i = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++, i++) {
      if (!((*VegMap)[y][x].Height = (float *)calloc(VType[(*VegMap)[y][x].Veg - 1].NVegLayers, sizeof(float))))
      ReportError((char *)Routine, 1);
    }
  }

  if (strncmp(StrEnv[vegheight_file].VarStr, "none", 4)) {
    printf("Spatial tree height map provided, reading height from map\n");

    /*Looks like both ovrestory height and understory height would need to be updated?*/
    for (NSet = 0; NSet < 2; NSet++) {
      if (!(Height = (float *)calloc(Map->NX * Map->NY, SizeOfNumberType(NumberType))))
        ReportError((char *)Routine, 1);
      flag = Read2DMatrix(StrEnv[vegheight_file].VarStr, Height, NumberType, Map, NSet, VarName, 0);

      if ((Options->FileFormat == NETCDF && flag == 0)
        || (Options->FileFormat == BIN))
      {
        for (y = 0, i = 0; y < Map->NY; y++) {
          for (x = 0; x < Map->NX; x++, i++) {
            if (NSet <= VType[(*VegMap)[y][x].Veg - 1].NVegLayers -1 ){
              if (Height[i] > 0.0)
                (*VegMap)[y][x].Height[NSet] = Height[i];
              else
                (*VegMap)[y][x].Height[NSet] = VType[(*VegMap)[y][x].Veg - 1].Height[NSet];
            }
          }
        }
      }
      else if (Options->FileFormat == NETCDF && flag == 1) {
        for (y = Map->NY - 1, i = 0; y >= 0; y--) {
          for (x = 0; x < Map->NX; x++, i++) {
            if (NSet <= VType[(*VegMap)[y][x].Veg - 1].NVegLayers -1){
              if (Height[i] > 0.0)
                (*VegMap)[y][x].Height[NSet] = Height[i];
              else
                (*VegMap)[y][x].Height[NSet] = VType[(*VegMap)[y][x].Veg - 1].Height[NSet];
            }
          }
        }
      }
      else ReportError((char *)Routine, 57);

      free(Height);
    }
    
  }
  else{
    printf("Vegetation tree height created from vegetation table\n");
    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
          if ( VType[(*VegMap)[y][x].Veg - 1].OverStory == TRUE) {
            (*VegMap)[y][x].Height[0] = VType[(*VegMap)[y][x].Veg - 1].Height[0];
            /*If understory exists, set default understory FC=1.0*/
            if (VType[(*VegMap)[y][x].Veg - 1].UnderStory == TRUE)
              (*VegMap)[y][x].Height[1] = VType[(*VegMap)[y][x].Veg - 1].Height[1];
          }
          else{
            if (VType[(*VegMap)[y][x].Veg - 1].UnderStory == TRUE)
              (*VegMap)[y][x].Height[0] = VType[(*VegMap)[y][x].Veg - 1].Height[0];
          }
        }
      }
  }
}


/*****************************************************************************
InitCanopyGapMap()
*****************************************************************************/
void InitCanopyGapMap(OPTIONSTRUCT *Options, LISTPTR Input, MAPSIZE *Map,
  LAYER *Soil, LAYER *Veg, VEGTABLE *VType, VEGPIX ***VegMap, 
  SOILTABLE *SType, SOILPIX ***SoilMap)
{
  const char *Routine = "InitCanopyGapMap";
  char VarName[BUFSIZE + 1];
  char CanopyMapFileName[BUFSIZE + 1];
  int i, j;			/* counter */
  int x;			/* counter */
  int y;			/* counter */
  int flag;
  int NVeg;
  int NSoil;
  int NumberType;		/* number type */
  float *Gap;		/* gap diameter */

  /* Get the canopy gap map filename from the [VEGETATION] section */
  GetInitString("VEGETATION", "CANOPY GAP MAP FILE", "", CanopyMapFileName,
    (unsigned long)BUFSIZE, Input);
  if (IsEmptyStr(CanopyMapFileName))
    ReportError("CANOPY GAP MAP FILE", 51);

  /* Read the vegetation type */
  GetVarName(007, 0, VarName);
  GetVarNumberType(007, &NumberType);
  if (!(Gap = (float *)calloc(Map->NX * Map->NY,
    SizeOfNumberType(NumberType))))
    ReportError((char *)Routine, 1);
  flag = Read2DMatrix(CanopyMapFileName, Gap, NumberType, Map, 0, VarName, 0);

  /* if NetCDF, may need to reverse the matrix */
  if ((Options->FileFormat == NETCDF && flag == 0)
    || (Options->FileFormat == BIN))
  {
    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*VegMap)[y][x].Gapping = Gap[i];
        /* set gapping to false for cells with no overstory */
        if (VType[(*VegMap)[y][x].Veg - 1].OverStory == FALSE)
          (*VegMap)[y][x].Gapping = 0.0;
      }
    }
  }
  else if (Options->FileFormat == NETCDF && flag == 1) {
    for (y = Map->NY - 1, i = 0; y >= 0; y--) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*VegMap)[y][x].Gapping = Gap[i];
        /* set gapping to false for cells with no overstory */
        if (VType[(*VegMap)[y][x].Veg - 1].OverStory == FALSE)
          (*VegMap)[y][x].Gapping = 0.0;
        /* set gapping to false given glacier cell */
        if (VType[(*VegMap)[y][x].Veg - 1].Index == GLACIER)
          (*VegMap)[y][x].Gapping = 0.0;
      }
    }
  }
  else ReportError((char *)Routine, 57);

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      NVeg = Veg->MaxLayers;
      NSoil = Soil->MaxLayers;
      if (Options->CanopyGapping) {
        if (!((*VegMap)[y][x].Type = (CanopyGapStruct *)calloc(2, sizeof(CanopyGapStruct))))
          ReportError((char *)Routine, 1);
        for (i = 0; i < CELL_PARTITION; i++) {
          if (!((*VegMap)[y][x].Type[i].IntRain = (float *)calloc(NVeg, sizeof(float))))
            ReportError((char *)Routine, 1);
          if (!((*VegMap)[y][x].Type[i].IntSnow = (float *)calloc(NVeg, sizeof(float))))
            ReportError((char *)Routine, 1);
          if (!((*VegMap)[y][x].Type[i].Moist = (float *)calloc(NSoil+1, sizeof(float))))
            ReportError((char *)Routine, 1);
          if (!((*VegMap)[y][x].Type[i].EPot = (float *)calloc(NVeg+1, sizeof(float))))
            ReportError((char *)Routine, 1);
          if (!((*VegMap)[y][x].Type[i].EAct = (float *)calloc(NVeg+1, sizeof(float))))
            ReportError((char *)Routine, 1);
          if (!((*VegMap)[y][x].Type[i].EInt = (float *)calloc(NVeg, sizeof(float))))
            ReportError((char *)Routine, 1);
          if (!((*VegMap)[y][x].Type[i].ESoil = (float **)calloc(NVeg, sizeof(float *))))
            ReportError((char *)Routine, 1);

          for (j = 0; j < NVeg; j++) {
            if (!((*VegMap)[y][x].Type[i].ESoil[j] = (float *)calloc(NSoil, sizeof(float))))
              ReportError((char *)Routine, 1);
          }
        }
      }
    }
  }
  free(Gap);
}

/*****************************************************************************
  InitWindMap()
*****************************************************************************/
void InitWindMap(OPTIONSTRUCT * Options, LISTPTR Input, MAPSIZE * Map,
                 TOPOPIX ** TopoMap, WINDPIX *** WindMap)
{
  const char *Routine = "InitWindMap";
  char VarName[BUFSIZE + 1];
  char VarStr[BUFSIZE + 1];
  int i, x, y;        /* counters */
  int flag;
  int NumberType;		  /* number type */
  float *FetchDist;		/* Fetch distance */
  float *Ux = NULL;		/* West-to-east wind speed */
  float *Uy = NULL;		/* North-to-south wind speed */
  float *Uz = NULL;		/* Down-to-up wind speed */
  float *Ke = NULL;		/* Turbulent kinetic energy */
  int NSet;           /* Counter for wind layer */
  float *WindHeight = NULL;
  float TotalXYabs;
  float CloudMinBase, LocalCloudBase, CloudExpDec;
  
  /* Get the map filename from the [TERRAIN] section */
  STRINIENTRY StrEnv[] = {
    {"TERRAIN", "DEM FILE", "", ""},
    {"TERRAIN", "BASIN MASK FILE", "", ""},
    {"TERRAIN", "SNOW CLOUD BASE", "",""},
    {"TERRAIN", "MINIMUM CLOUD BASE", "",""},
    {"TERRAIN", "CLOUD SPEED EXP DEC", "",""},
    {"TERRAIN", "SNOW SETTLING SPEED", "",""},
    {"TERRAIN", "NUMBER OF WIND LAYERS", "",""},
    {"TERRAIN", "FETCH DISTANCE FILE", "", ""},
    {"TERRAIN", "WIND UX MAP FILE", "", ""},
    {"TERRAIN", "WIND UY MAP FILE", "", ""},
    {"TERRAIN", "WIND UZ MAP FILE", "", ""},
    {"TERRAIN", "WIND KE MAP FILE", "", ""},
    {NULL, NULL, "", NULL}
  };
  
  /* Read the key-entry pairs from the input file */
  for (i = 0; StrEnv[i].SectionName; i++) {
    GetInitString(StrEnv[i].SectionName, StrEnv[i].KeyName, StrEnv[i].Default,
      StrEnv[i].VarStr, (unsigned long)BUFSIZE, Input);
    if (IsEmptyStr(StrEnv[i].VarStr))
      ReportError(StrEnv[i].KeyName, 51);
  }
  
  /* Read number of wind layers from config */
  if (IsEmptyStr(StrEnv[windlayers_num].VarStr))  
    ReportError(StrEnv[windlayers_num].KeyName, 51);
  else if (!CopyInt(&NWINDLAYERS, StrEnv[windlayers_num].VarStr, 1))
    ReportError(StrEnv[windlayers_num].KeyName, 51);
  printf("\nReading %d wind layers\n",NWINDLAYERS);
  
  /* Read other data */
  if (IsEmptyStr(StrEnv[snowcloudbase].VarStr))
    ReportError(StrEnv[snowcloudbase].KeyName, 51);
  else if (!CopyFloat(&CloudBaseElev, StrEnv[snowcloudbase].VarStr, 1))
    ReportError(StrEnv[snowcloudbase].KeyName, 51);
  printf("Snowfall will be instantiated at %.0f m elevation,\n",CloudBaseElev);
  
  if (IsEmptyStr(StrEnv[mincloudheight].VarStr))
    ReportError(StrEnv[mincloudheight].KeyName, 51);
  else if (!CopyFloat(&CloudMinBase, StrEnv[mincloudheight].VarStr, 1))
    ReportError(StrEnv[mincloudheight].KeyName, 51);
  printf("but no less than %.0f m above ground level\n\n",CloudMinBase);
  
  if (IsEmptyStr(StrEnv[snowexpdec].VarStr))
    ReportError(StrEnv[snowexpdec].KeyName, 51);
  else if (!CopyFloat(&CloudExpDec, StrEnv[snowexpdec].VarStr, 1))
    ReportError(StrEnv[snowexpdec].KeyName, 51);
  
  if (IsEmptyStr(StrEnv[snowfallvel].VarStr))
    ReportError(StrEnv[snowfallvel].KeyName, 51);
  else if (!CopyFloat(&SNOWFALLVEL, StrEnv[snowfallvel].VarStr, 1))
    ReportError(StrEnv[snowfallvel].KeyName, 51);
  
  /* Allocate memory for horizontal map */
  if (!(*WindMap = (WINDPIX **)calloc(Map->NY, sizeof(WINDPIX *))))
    ReportError((char *)Routine, 1);
  for (y = 0; y < Map->NY; y++) {
    if (!((*WindMap)[y] = (WINDPIX *)calloc(Map->NX, sizeof(WINDPIX))))
      ReportError((char *)Routine, 1);
  }
  /* Allocate memory for vertical layers */
  for (y = 0, i = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (!((*WindMap)[y][x].LayerElevUpper = (float *)calloc(NWINDLAYERS, sizeof(float))))
        ReportError((char *)Routine, 1);
      if (!((*WindMap)[y][x].LayerElevLower = (float *)calloc(NWINDLAYERS, sizeof(float))))
        ReportError((char *)Routine, 1);
      if (!((*WindMap)[y][x].WindSpeedXY = (float *)calloc(NWINDLAYERS, sizeof(float))))
        ReportError((char *)Routine, 1);
      if (!((*WindMap)[y][x].WindSpeedZ = (float *)calloc(NWINDLAYERS, sizeof(float))))
        ReportError((char *)Routine, 1);
      if (!((*WindMap)[y][x].TurbulenceK = (float *)calloc(NWINDLAYERS, sizeof(float))))
        ReportError((char *)Routine, 1);
      if (!((*WindMap)[y][x].SnowingScale = (float *)calloc(NWINDLAYERS, sizeof(float))))
        ReportError((char *)Routine, 1);
      if (!((*WindMap)[y][x].Qsusp = (float *)calloc(NWINDLAYERS, sizeof(float))))
        ReportError((char *)Routine, 1);
      if (!((*WindMap)[y][x].QsuspLastIt = (float *)calloc(NWINDLAYERS, sizeof(float))))
        ReportError((char *)Routine, 1);
      if (!((*WindMap)[y][x].SublimationFrac = (float *)calloc(NWINDLAYERS, sizeof(float))))
        ReportError((char *)Routine, 1);
      if (!((*WindMap)[y][x].WindDirFrac = (float **)calloc(NWINDLAYERS, sizeof(float *))))
        ReportError((char *)Routine, 1);
      for (NSet = 0; NSet < NWINDLAYERS; NSet++) {
        if (!((*WindMap)[y][x].WindDirFrac[NSet] = (float *)calloc(4, sizeof(float))))
          ReportError((char *)Routine, 1);
      }
    }
  }
  
  /* Read the fetch distance */
  GetVarName(900, 0, VarName);
  GetVarNumberType(900, &NumberType);
  if (!(FetchDist = (float *)calloc(Map->NX * Map->NY,
                     SizeOfNumberType(NumberType))))
    ReportError((char *)Routine, 1);
  flag = Read2DMatrix(StrEnv[fetchfile].VarStr, FetchDist, NumberType, Map, 0, VarName, 0);

  if ((Options->FileFormat == NETCDF && flag == 0)
        || (Options->FileFormat == BIN)) {
    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*WindMap)[y][x].FetchDist = FetchDist[i];
      }
    }
  }
  else if (Options->FileFormat == NETCDF && flag == 1) {
    for (y = Map->NY - 1, i = 0; y >= 0; y--) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*WindMap)[y][x].FetchDist = FetchDist[i];
      }
    }
  }
  else ReportError((char *)Routine, 57);
  free(FetchDist);
  
  /* Read the Ux, Uy, Uz wind speeds and Ke layer-by-layer */
  GetVarName(903, 0, VarName);
  GetVarNumberType(903, &NumberType);
  if (!(Ux = (float *)calloc(Map->NX * Map->NY,
              SizeOfNumberType(NumberType))))
    ReportError((char *)Routine, 1);
  if (!(Uy = (float *)calloc(Map->NX * Map->NY,
              SizeOfNumberType(NumberType))))
    ReportError((char *)Routine, 1);
  if (!(Uz = (float *)calloc(Map->NX * Map->NY,
              SizeOfNumberType(NumberType))))
    ReportError((char *)Routine, 1);
  if (!(Ke = (float *)calloc(Map->NX * Map->NY,
              SizeOfNumberType(NumberType))))
    ReportError((char *)Routine, 1);
  
  for (NSet = 0; NSet < NWINDLAYERS; NSet++) {
    /* Read current layer */
    flag = Read2DMatrix(StrEnv[winduxfile].VarStr, Ux, NumberType, Map, NSet, VarName, 0);
    flag = Read2DMatrix(StrEnv[winduyfile].VarStr, Uy, NumberType, Map, NSet, VarName, 0);
    flag = Read2DMatrix(StrEnv[winduzfile].VarStr, Uz, NumberType, Map, NSet, VarName, 0);
    flag = Read2DMatrix(StrEnv[windkefile].VarStr, Ke, NumberType, Map, NSet, VarName, 0);
    
    if ((Options->FileFormat == NETCDF && flag == 0)
          || (Options->FileFormat == BIN)) {
      for (y = 0, i = 0; y < Map->NY; y++) {
        for (x = 0; x < Map->NX; x++, i++) {
          
          (*WindMap)[y][x].TurbulenceK[NSet] = Ke[i];
          
          // Ux[i] = 1.0;
          // Uy[i] = 0.0;
          // Uz[i] = 0.0;
          
          (*WindMap)[y][x].WindSpeedXY[NSet] = pow(Ux[i]*Ux[i] + Uy[i]*Uy[i], 0.5);
          (*WindMap)[y][x].WindSpeedZ[NSet] = Uz[i];
          
          /* Calculate fractional directions */
          TotalXYabs = ABSVAL(Ux[i]) + ABSVAL(Uy[i]);
          (*WindMap)[y][x].WindDirFrac[NSet][0] = (Uy[i] < 0.0 ? -1.0 * Uy[i] / TotalXYabs : 0.0);
          (*WindMap)[y][x].WindDirFrac[NSet][1] = (Ux[i] > 0.0 ? Ux[i] / TotalXYabs : 0.0);
          (*WindMap)[y][x].WindDirFrac[NSet][2] = (Uy[i] > 0.0 ? Uy[i] / TotalXYabs : 0.0);
          (*WindMap)[y][x].WindDirFrac[NSet][3] = (Ux[i] < 0.0 ? -1.0 * Ux[i] / TotalXYabs : 0.0);
        }
      }
    }
    else if (Options->FileFormat == NETCDF && flag == 1) {
      for (y = Map->NY - 1, i = 0; y >= 0; y--) {
        for (x = 0; x < Map->NX; x++, i++) {
          (*WindMap)[y][x].WindSpeedXY[NSet] = pow(Ux[i]*Ux[i] + Uy[i]*Uy[i], 0.5);
          (*WindMap)[y][x].WindSpeedZ[NSet] = Uz[i];
          
          /* Calculate fractional directions */
          TotalXYabs = ABSVAL(Ux[i]) + ABSVAL(Uy[i]);
          (*WindMap)[y][x].WindDirFrac[NSet][0] = (Uy[i] < 0.0 ? -1.0 * Uy[i] / TotalXYabs : 0.0);
          (*WindMap)[y][x].WindDirFrac[NSet][1] = (Ux[i] > 0.0 ? Ux[i] / TotalXYabs : 0.0);
          (*WindMap)[y][x].WindDirFrac[NSet][2] = (Uy[i] > 0.0 ? Uy[i] / TotalXYabs : 0.0);
          (*WindMap)[y][x].WindDirFrac[NSet][3] = (Ux[i] < 0.0 ? -1.0 * Ux[i] / TotalXYabs : 0.0);
        }
      }
    }
    else ReportError((char *)Routine, 57);
  }
  free(Ux);
  free(Uy);
  free(Uz);
  free(Ke);
  
  /* Construct layer elevations from max height above ground and DEM */
  if (!(WindHeight = (float *)calloc(NWINDLAYERS, sizeof(float))))
    ReportError((char *)Routine, 1);
  /* Read max layer height above surface */
  for (NSet = 0; NSet < NWINDLAYERS; NSet++) {
    sprintf(VarName, "WIND LAYER HEIGHT %d", NSet + 1);
    GetInitString("TERRAIN", VarName, "", VarStr,
      (unsigned long)BUFSIZE, Input);
    if (!CopyFloat(&(WindHeight[NSet]), VarStr, 1))
      ReportError((char *)VarName, 51);
  }
  /* Add to DEM elevation for each layer and grid cell */
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      for (NSet = 0; NSet < NWINDLAYERS; NSet++) {
        (*WindMap)[y][x].LayerElevUpper[NSet] = WindHeight[NSet] + TopoMap[y][x].Dem;
        if (NSet == 0)
          (*WindMap)[y][x].LayerElevLower[NSet] = TopoMap[y][x].Dem;
        else
          (*WindMap)[y][x].LayerElevLower[NSet] = (*WindMap)[y][x].LayerElevUpper[NSet - 1];
      }
    }
  }
  free(WindHeight);
  
  /* Identify which layer will be used to initialize snowfall */
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      LocalCloudBase = MAX(CloudBaseElev, TopoMap[y][x].Dem + CloudMinBase);
      
      /* Also consider horizontal distance to nearest terrain (i.e., cliffs) */
      
      for (NSet = 0; NSet < NWINDLAYERS; NSet++) {
        if ((*WindMap)[y][x].LayerElevUpper[NSet] >= LocalCloudBase ||
            NSet == NWINDLAYERS - 1) {
          (*WindMap)[y][x].SnowfallLayer = (int) NSet;
          // printf("Elev: %.0f Layer: %d\n",TopoMap[y][x].Dem,(*WindMap)[y][x].SnowfallLayer);
          break;
        }
      }
    }
  }
  
  /* Calculate scale factors to adjust wind speed profiles during snowfall */
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      LocalCloudBase = MAX(CloudBaseElev, TopoMap[y][x].Dem + CloudMinBase);
      for (NSet = 0; NSet < NWINDLAYERS; NSet++) {
        if (CloudExpDec > 0.0) {
          if ((*WindMap)[y][x].LayerElevUpper[NSet] >= LocalCloudBase) {
            (*WindMap)[y][x].SnowingScale[NSet] = 0.0;
            
          } else if ((*WindMap)[y][x].LayerElevUpper[NSet] >= (LocalCloudBase - CloudMinBase)) {
            (*WindMap)[y][x].SnowingScale[NSet] = (exp(-1.0 * ((*WindMap)[y][x].LayerElevUpper[NSet] - (LocalCloudBase - CloudMinBase)) *
              CloudExpDec / CloudMinBase) - exp(-1.0 * CloudExpDec)) /
                (1.0 - exp(-1.0 * CloudExpDec));
            
          } else {
            (*WindMap)[y][x].SnowingScale[NSet] = 1.0;
          }
        } else {
          (*WindMap)[y][x].SnowingScale[NSet] = 1.0;
        }
        if (x==30 && y==175)
          printf("Elev: %.0f Layer: %d CloudBase: %.0f WindScale: %.2f\n",
                 TopoMap[y][x].Dem,NSet,LocalCloudBase,(*WindMap)[y][x].SnowingScale[NSet]);
      }
    }
  }
  
  /* Initialize fluxes to zero */
  NwindIters = 1;
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      for (NSet = 0; NSet < NWINDLAYERS; NSet++) {
        (*WindMap)[y][x].Qsusp[NSet] = 0.0;
        (*WindMap)[y][x].QsuspLastIt[NSet] = 0.0;
      }
      (*WindMap)[y][x].Qsalt = 0.0;
      (*WindMap)[y][x].WindDeposition = 0.0;
      (*WindMap)[y][x].nIters = NwindIters;
    }
  }
}





