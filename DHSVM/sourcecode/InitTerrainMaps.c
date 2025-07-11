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

 /*****************************************************************************
   InitTerrainMaps()
 *****************************************************************************/
void InitTerrainMaps(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map,
  LAYER *Soil, LAYER *Veg, TOPOPIX ***TopoMap, SOILTABLE *SType, SOILPIX ***SoilMap, 
  VEGTABLE *VType, VEGPIX ***VegMap, DYNAVEG *DVeg)

{
  printf("\nInitializing terrain maps\n");

  InitTopoMap(Input, Options, Map, TopoMap);
  InitVegMap(Options, Input, Map, VegMap, VType, DVeg);
  InitSoilMap(Input, Options, Map, Soil, *TopoMap, SoilMap, SType, VegMap, VType);
  if (Options->CanopyGapping)
    InitCanopyGapMap(Options, Input, Map, Soil, Veg, VType, VegMap, SType, SoilMap);
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
  float TotalRootDepth, LayerDepth, Transmissivity, KsVertCalc;
  
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
  
  /* Apply calibration parameter adjustment */
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      (*SoilMap)[y][x].Depth += SOIL_DEPTH_ADJ;
    }
  }
  
  /* Ensure soil is at least 1 cm deeper than rooting depth */
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN((TopoMap)[y][x].Mask)) {
        TotalRootDepth = 0.0;
        sidx = (*SoilMap)[y][x].Soil - 1;
        for (NSet = 0; NSet < Soil->NLayers[sidx]; NSet++) {
          TotalRootDepth += VType[(*VegMap)[y][x].Veg - 1].RootDepth[NSet];
        }
        if ((*SoilMap)[y][x].Depth < (TotalRootDepth + 0.01))
          (*SoilMap)[y][x].Depth = (TotalRootDepth + 0.01);
      }
    }
  }
  
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
  
  /* Apply calibration parameter adjustment */
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      (*SoilMap)[y][x].KsLat *= SOIL_KSAT_ADJ;
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
  
  /* Apply calibration parameter adjustment */
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if ((*SoilMap)[y][x].KsLatExp < 0.001 &&
          SOIL_EXPDEC_ADJ != 1.0) {
        (*SoilMap)[y][x].KsLatExp = 0.001;
      }
      (*SoilMap)[y][x].KsLatExp *= SOIL_EXPDEC_ADJ;
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
                  (*SoilMap)[y][x].FCap[NSet] = FC[i] * SOIL_FIELDCAP_ADJ;
                else
                  (*SoilMap)[y][x].FCap[NSet] = SType[sidx].FCap[NSet] * SOIL_FIELDCAP_ADJ;
                /* Make sure FCap is larger than WP */
                if ((*SoilMap)[y][x].FCap[NSet] < (SType[sidx].WP[NSet] + 0.01))
                  (*SoilMap)[y][x].FCap[NSet] = SType[sidx].WP[NSet] + 0.01;
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
                  (*SoilMap)[y][x].FCap[NSet] = FC[i] * SOIL_FIELDCAP_ADJ;
                else
                  (*SoilMap)[y][x].FCap[NSet] = SType[sidx].FCap[NSet] * SOIL_FIELDCAP_ADJ;
                /* Make sure FCap is larger than WP */
                if ((*SoilMap)[y][x].FCap[NSet] < (SType[sidx].WP[NSet] + 0.01))
                  (*SoilMap)[y][x].FCap[NSet] = SType[sidx].WP[NSet] + 0.01;
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
            (*SoilMap)[y][x].FCap[NSet] = SType[sidx].FCap[NSet] * SOIL_FIELDCAP_ADJ;
            /* Make sure FCap is larger than WP */
            if ((*SoilMap)[y][x].FCap[NSet] < (SType[sidx].WP[NSet] + 0.01))
              (*SoilMap)[y][x].FCap[NSet] = SType[sidx].WP[NSet] + 0.01;
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
                  (*SoilMap)[y][x].Porosity[NSet] = Porosity[i] * SOIL_POROSITY_ADJ;
                else
                  (*SoilMap)[y][x].Porosity[NSet] = SType[sidx].Porosity[NSet] * SOIL_POROSITY_ADJ;
                /* Make sure porosity is larger than FCap, which is already larger than WP */
                if ((*SoilMap)[y][x].Porosity[NSet] < ((*SoilMap)[y][x].FCap[NSet] + 0.01))
                  (*SoilMap)[y][x].Porosity[NSet] = (*SoilMap)[y][x].FCap[NSet] + 0.01;
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
                  (*SoilMap)[y][x].Porosity[NSet] = Porosity[i] * SOIL_POROSITY_ADJ;
                else
                  (*SoilMap)[y][x].Porosity[NSet] = SType[sidx].Porosity[NSet] * SOIL_POROSITY_ADJ;
                /* Make sure porosity is larger than FCap, which is already larger than WP */
                if ((*SoilMap)[y][x].Porosity[NSet] < ((*SoilMap)[y][x].FCap[NSet] + 0.01))
                  (*SoilMap)[y][x].Porosity[NSet] = (*SoilMap)[y][x].FCap[NSet] + 0.01;
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
            (*SoilMap)[y][x].Porosity[NSet] = SType[sidx].Porosity[NSet] * SOIL_POROSITY_ADJ;
            /* Make sure porosity is larger than FCap, which is already larger than WP */
            if ((*SoilMap)[y][x].Porosity[NSet] < ((*SoilMap)[y][x].FCap[NSet] + 0.01))
              (*SoilMap)[y][x].Porosity[NSet] = (*SoilMap)[y][x].FCap[NSet] + 0.01;
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
          if (VType[(*VegMap)[y][x].Veg - 1].OverStory == TRUE) {
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

          if (VType[(*VegMap)[y][x].Veg - 1].OverStory == TRUE) {  
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
            if (VType[(*VegMap)[y][x].Veg - 1].UnderStory == TRUE)
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
  
  /* Apply calibration parameter adjustment - overstory only */
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (VType[(*VegMap)[y][x].Veg - 1].OverStory == TRUE) {
        (*VegMap)[y][x].Fract[0] *= VEG_FC_ADJ;
        /* Limit to maximum of 100% cover */
        if ((*VegMap)[y][x].Fract[0] > 1.0)
          (*VegMap)[y][x].Fract[0] = 1.0;
      }
    }
  }
  
  /* Calculate Vf */
  if (Options->ImprovRadiation == TRUE) {
    for (y = 0; y < Map->NY; y++) {
        for (x = 0; x < Map->NX; x++) {
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
          if (VType[(*VegMap)[y][x].Veg - 1].OverStory == TRUE) {
           if (LAIMonthly[i] > 0.0)
              (*VegMap)[y][x].LAIMonthly[0][NSet] = LAIMonthly[i];
            else
              (*VegMap)[y][x].LAIMonthly[0][NSet] = VType[(*VegMap)[y][x].Veg - 1].LAIMonthly[0][NSet];
           
            if (VType[(*VegMap)[y][x].Veg - 1].UnderStory  == TRUE )
              (*VegMap)[y][x].LAIMonthly[1][NSet] = VType[(*VegMap)[y][x].Veg - 1].LAIMonthly[1][NSet];
          }
          else{
            if (VType[(*VegMap)[y][x].Veg - 1].UnderStory == TRUE)
              (*VegMap)[y][x].LAIMonthly[0][NSet] = VType[(*VegMap)[y][x].Veg - 1].LAIMonthly[0][NSet];
          }
        }
      }
    }
    else if (Options->FileFormat == NETCDF && flag == 1) {
      for (y = Map->NY - 1, i = 0; y >= 0; y--) {
        for (x = 0; x < Map->NX; x++, i++) {
        
          if (VType[(*VegMap)[y][x].Veg - 1].OverStory == TRUE) {
            if (LAIMonthly[i] > 0.0)
              (*VegMap)[y][x].LAIMonthly[0][NSet] = LAIMonthly[i];
            else
              (*VegMap)[y][x].LAIMonthly[0][NSet] = VType[(*VegMap)[y][x].Veg - 1].LAIMonthly[0][NSet];
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

            if (VType[(*VegMap)[y][x].Veg - 1].OverStory == TRUE) {        
              (*VegMap)[y][x].LAIMonthly[0][NSet] = VType[(*VegMap)[y][x].Veg - 1].LAIMonthly[0][NSet];
              if ( VType[(*VegMap)[y][x].Veg - 1].UnderStory == TRUE ){
                (*VegMap)[y][x].LAIMonthly[1][NSet] = VType[(*VegMap)[y][x].Veg - 1].LAIMonthly[1][NSet]; 
              }
            }
            else{
              if (VType[(*VegMap)[y][x].Veg - 1].UnderStory == TRUE){
                (*VegMap)[y][x].LAIMonthly[0][NSet] = VType[(*VegMap)[y][x].Veg - 1].LAIMonthly[0][NSet];
              }
            }
          }
        }
    }
  }
  
  /* Apply calibration parameter adjustment - overstory only */
  for (NSet = 0; NSet < 12; NSet++) {
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
        if (VType[(*VegMap)[y][x].Veg - 1].OverStory == TRUE) {
          (*VegMap)[y][x].LAIMonthly[0][NSet] *= VEG_LAI_ADJ;
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







