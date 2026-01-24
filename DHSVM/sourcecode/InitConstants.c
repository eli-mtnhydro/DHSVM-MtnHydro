
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "data.h"
#include "Calendar.h"
#include "fileio.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "getinit.h"
#include "constants.h"
#include "rad.h"

/*****************************************************************************
  Function name: InitConstants()

  Purpose      : Initialize constants and settings for DHSVM run
                 Processes the following sections in InFile:
                 [OPTIONS]
                 [AREA]
                 [TIME]
                 [CONSTANTS}

  Required     :
    LISTPTR Input          - Linked list with input strings
    OPTIONSTRUCT *Options   - Structure with different program options
    MAPSIZE *Map            - Coverage and resolution of model area
    SOLARGEOMETRY *SolarGeo - Solar geometry information
    TIMESTRUCT *Time        - Begin and end times, model timestep

*****************************************************************************/
void InitConstants(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map,
		   SOLARGEOMETRY *SolarGeo, TIMESTRUCT *Time)
{
  int i;			/* counter */
  double PointModelX;		/* X-coordinate for POINT model mode */
  double PointModelY;		/* Y-coordinate for POINT model mode */
  float TimeStep;		/* Timestep in hours */
  DATE End;			/* End of run */
  DATE Start;			/* Start of run */

  STRINIENTRY StrEnv[] = {
    {"OPTIONS", "EXTENT", "", ""},
    {"OPTIONS", "GRADIENT", "", ""},
    {"OPTIONS", "ROUTING NEIGHBORS", "", "8"},
    {"OPTIONS", "MULTIPLE FLOW DIRECTIONS", "", "TRUE"},
    {"OPTIONS", "SENSIBLE HEAT FLUX", "", ""},
    {"OPTIONS", "OVERLAND ROUTING", "", ""},
    {"OPTIONS", "LAKE DYNAMICS", "", "FALSE"},
    {"OPTIONS", "VERTICAL KSAT SOURCE", "", "TABLE"},
    {"OPTIONS", "INFILTRATION", "", ""},
    {"OPTIONS", "INTERPOLATION", "", ""},
    {"OPTIONS", "MAX INTERPOLATION DISTANCE", "", ""},
    {"OPTIONS", "PRISM", "", ""},
    {"OPTIONS", "SNOW PATTERN", "", ""},
    {"OPTIONS", "CANOPY RADIATION ATTENUATION MODE", "", ""},
    {"OPTIONS", "SHADING", "", ""},
    {"OPTIONS", "OUTSIDE", "", ""},
    {"OPTIONS", "RHOVERRIDE", "", ""},
    {"OPTIONS", "TEMPERATURE LAPSE RATE", "", ""},
    {"OPTIONS", "CRESSMAN RADIUS", "", ""},
    {"OPTIONS", "CRESSMAN STATIONS", "", ""},
    {"OPTIONS", "PRISM DATA PATH", "", ""},
    {"OPTIONS", "PRISM DATA EXTENSION", "", ""},
    {"OPTIONS", "SNOW PATTERN DATA PATH", "", ""},
    {"OPTIONS", "SHADING DATA PATH", "", ""},
    {"OPTIONS", "SHADING DATA EXTENSION", "", ""},
    {"OPTIONS", "SKYVIEW DATA PATH", "", ""},
	  {"OPTIONS", "VARIABLE LIGHT TRANSMITTANCE", "", "" },
    {"OPTIONS", "CANOPY GAPPING", "", "" },
    {"OPTIONS", "SNOW SLIDING", "", "" },
    {"OPTIONS", "PRECIPITATION SEPARATION", "", "FALSE" },
    {"OPTIONS", "SNOW STATISTICS", "", "FALSE" },
    {"OPTIONS", "DYNAMIC VEGETATION", "", "FALSE" },
    {"OPTIONS", "EXTRA STREAM STATE DATA", "", "FALSE" },
    {"OPTIONS", "GROUNDWATER SPINUP", "", "FALSE" },
    {"OPTIONS", "GROUNDWATER SPINUP YEARS", "", "0" },
    {"OPTIONS", "GROUNDWATER SPINUP RECHARGE", "", "0.0" },
    {"AREA", "COORDINATE SYSTEM", "", ""},
    {"AREA", "EXTREME NORTH", "", ""},
    {"AREA", "EXTREME WEST", "", ""},
    {"AREA", "CENTER LATITUDE", "", ""},
    {"AREA", "CENTER LONGITUDE", "", ""},
    {"AREA", "TIME ZONE MERIDIAN", "", ""},
    {"AREA", "NUMBER OF ROWS", "", ""},
    {"AREA", "NUMBER OF COLUMNS", "", ""},
    {"AREA", "GRID SPACING", "", ""},
    {"AREA", "POINT NORTH", "", ""},
    {"AREA", "POINT EAST", "", ""},
    {"TIME", "TIME STEP", "", ""},
    {"TIME", "MODEL START", "", ""},
    {"TIME", "MODEL END", "", ""},
    {"CONSTANTS", "GROUND ROUGHNESS", "", ""},
    {"CONSTANTS", "SNOW ROUGHNESS", "", ""},
    {"CONSTANTS", "SNOW WATER CAPACITY", "", ""},
    {"CONSTANTS", "REFERENCE HEIGHT", "", ""},
    {"CONSTANTS", "RAIN LAI MULTIPLIER", "", ""},
    {"CONSTANTS", "SNOW LAI MULTIPLIER", "", ""},
    {"CONSTANTS", "MIN INTERCEPTED SNOW", "", ""},
    {"CONSTANTS", "MIN ALBEDO RESET SNOWFALL", "", "" },
    {"CONSTANTS", "ALBEDO LAMBDA FOREST OFFSET", "", "0.0" },
    {"CONSTANTS", "OUTSIDE BASIN VALUE", "", ""},
    {"CONSTANTS", "TEMPERATURE LAPSE RATE", "", ""},
    {"CONSTANTS", "MAX SURFACE SNOW LAYER DEPTH", "", "0.125" },
    {"CONSTANTS", "SNOWSLIDE PARAMETER1", "", "" },
    {"CONSTANTS", "SNOWSLIDE PARAMETER2", "", "" },
    {"CONSTANTS", "GAP WIND ADJ FACTOR", "", "" },
    {"CONSTANTS", "SNOW PATTERN WEIGHT", "", "1.0" },
    {"CONSTANTS", "TEMPERATURE OFFSET", "", "0.0" },
    {"CONSTANTS", "LAPSE RATE BIAS", "", "0.0" },
    {"CONSTANTS", "LAPSE RATE BIAS ELEVATION", "", "0.0" },
    {"CONSTANTS", "SOIL DEPTH ADJUST", "", "0.0" },
    {"CONSTANTS", "SOIL CONDUCTIVITY ADJUST", "", "1.0" },
    {"CONSTANTS", "SOIL EXP DEC ADJUST", "", "1.0" },
    {"CONSTANTS", "SOIL POROSITY ADJUST", "", "1.0" },
    {"CONSTANTS", "SOIL FIELD CAP ADJUST", "", "1.0" },
    {"CONSTANTS", "VEG LAI ADJUST", "", "1.0" },
    {"CONSTANTS", "VEG TREE COVER ADJUST", "", "1.0" },
    {NULL, NULL, "", NULL}
  };

  /* Read the key-entry pairs from the input file */
  for (i = 0; StrEnv[i].SectionName; i++)
    GetInitString(StrEnv[i].SectionName, StrEnv[i].KeyName, StrEnv[i].Default,
		  StrEnv[i].VarStr, (unsigned long) BUFSIZE, Input);

  /**************** Determine model options ****************/
  
  /* Determine whether the model should be run in POINT mode or in BASIN mode.
     If in POINT mode also read which pixel to model */
  if (strncmp(StrEnv[extent].VarStr, "POINT", 5) == 0) {
    Options->Extent = POINT;
  }
  else if (strncmp(StrEnv[extent].VarStr, "BASIN", 5) == 0) {
    Options->Extent = BASIN;
  }
  else
    ReportError(StrEnv[extent].KeyName, 51);
  
  /* Determine how the flow gradient should be calculated */
  if (Options->Extent != POINT) {
    if (strncmp(StrEnv[gradient].VarStr, "TOPO", 4) == 0)
      Options->FlowGradient = TOPOGRAPHY;
    else if (strncmp(StrEnv[gradient].VarStr, "WATER", 5) == 0)
      Options->FlowGradient = WATERTABLE;
    else
      ReportError(StrEnv[gradient].KeyName, 51);
  }
  else
    Options->FlowGradient = NOT_APPLICABLE;
  
  /* Determine how many neighbors are used for surface/subsurface routing */
  if (!CopyInt(&(NDIRS), StrEnv[routing_neighbors].VarStr, 1))
    ReportError(StrEnv[routing_neighbors].KeyName, 51);
  switch (NDIRS) {
  case (4):
    xdirection = &xdirection4[0];
    ydirection = &ydirection4[0];
    break;
  case (8):
    xdirection = &xdirection8[0];
    ydirection = &ydirection8[0];
    break;
  default:
    ReportError(StrEnv[routing_neighbors].KeyName, 51);
  }
  /* Determine whether to use MFD-8 or D8 routing */
  if (strncmp(StrEnv[routing_mfd].VarStr, "TRUE", 4) == 0) {
    Options->MultiFlowDir = TRUE;
    printf("Using %d neighbors with multiple flow directions for surface/subsurface routing\n", NDIRS);
  }
  else if (strncmp(StrEnv[routing_mfd].VarStr, "FALSE", 5) == 0) {
    Options->MultiFlowDir = FALSE;
    if (NDIRS == 8)
      printf("Using %d neighbors with steepest descent for surface/subsurface routing\n", NDIRS);
    else
      printf("Using %d neighbors with multiple flow directions for surface/subsurface routing\n"
               "(Steepest descent is only available when ROUTING NEIGHBORS = 8)\n", NDIRS);
  }
  else
    ReportError(StrEnv[routing_mfd].KeyName, 51);
  
  /* Determine what meterological interpolation to use */
  if (strncmp(StrEnv[interpolation].VarStr, "INVDIST", 7) == 0)
    Options->Interpolation = INVDIST;
  else if (strncmp(StrEnv[interpolation].VarStr, "NEAREST", 7) == 0)
    Options->Interpolation = NEAREST;
  else if (strncmp(StrEnv[interpolation].VarStr, "VARCRESS", 8) == 0)
    Options->Interpolation = VARCRESS;
  else if (strncmp(StrEnv[interpolation].VarStr, "UNIFORM", 7) == 0)
    Options->Interpolation = UNIFORM;
  else
    ReportError(StrEnv[interpolation].KeyName, 51);
  
  /* If INVDIST interpolation then get maximum distance */
  if (Options->Interpolation == INVDIST) {
    if (!CopyFloat(&(Options->MaxInterpDist), StrEnv[max_interp_dist].VarStr, 1)) {
      printf("\nMaximum meteorological station interpolation distance not supplied\n");
      printf("Setting no limit on interpolation distance\n\n");
      Options->MaxInterpDist = DHSVM_HUGE;
    }
  }
  
  /* If VARIABLE CRESTMAN interpolation then get parameters */
  if (Options->Interpolation == VARCRESS) {
    if (!CopyInt(&(Options->CressRadius), StrEnv[cressman_radius].VarStr, 1))
      ReportError(StrEnv[cressman_radius].KeyName, 51);
    if (!CopyInt
	(&(Options->CressStations), StrEnv[cressman_stations].VarStr, 1))
      ReportError(StrEnv[cressman_stations].KeyName, 51);
  }
  
  /* Determine whether a sensible heat flux should be calculated */
  if (strncmp(StrEnv[sensible_heat_flux].VarStr, "TRUE", 4) == 0)
    Options->HeatFlux = TRUE;
  else if (strncmp(StrEnv[sensible_heat_flux].VarStr, "FALSE", 5) == 0)
    Options->HeatFlux = FALSE;
  else
    ReportError(StrEnv[sensible_heat_flux].KeyName, 51);
 
  /* Determine overland flow routing method to use */
  if (strncmp(StrEnv[routing].VarStr, "KINEMATIC", 9) == 0)
    Options->Routing = TRUE;
  else if (strncmp(StrEnv[routing].VarStr, "CONVENTIONAL", 12) == 0)
    Options->Routing = FALSE;
  else
    ReportError(StrEnv[routing].KeyName, 51);
  
  /* Determine whether to simulate lake dynamics */
  if (strncmp(StrEnv[lakedyna].VarStr, "TRUE", 4) == 0)
    Options->LakeDynamics = TRUE;
  else if (strncmp(StrEnv[lakedyna].VarStr, "FALSE", 5) == 0)
    Options->LakeDynamics = FALSE;
  else
    ReportError(StrEnv[lakedyna].KeyName, 51);
  
  /* Determine if vertical conductivity is calculated from map anisotropy
   or read from default soil table */
  if (strncmp(StrEnv[vertksatsource].VarStr, "ANISOTROPY", 6) == 0)
    Options->UseKsatAnisotropy = TRUE;
  else if (strncmp(StrEnv[vertksatsource].VarStr, "TABLE", 12) == 0)
    Options->UseKsatAnisotropy = FALSE;
  else
    ReportError(StrEnv[vertksatsource].KeyName, 51);
  
  /* Determine if the maximum infiltration rate is static or dynamic */
  if (strncmp(StrEnv[infiltration].VarStr, "STATIC", 6) == 0)
    Options->Infiltration = STATIC;
  else if (strncmp(StrEnv[infiltration].VarStr, "DYNAMIC", 7) == 0) {
    Options->Infiltration = DYNAMIC ;
    printf("WARNING: Dynamic maximum infiltration capacity has\n");
    printf("not been fully tested. It is a work in progress.\n\n");
  }
  else
    ReportError(StrEnv[infiltration].KeyName, 51);
  
  /* Determine if PRISM maps will be used to interpolate precip fields */
  if (strncmp(StrEnv[prism].VarStr, "TRUE", 4) == 0)
    Options->Prism = TRUE;
  else if (strncmp(StrEnv[prism].VarStr, "FALSE", 5) == 0)
    Options->Prism = FALSE;
  else
    ReportError(StrEnv[prism].KeyName, 51);
  
  /* Determine if snow pattern map will be used to redistribute snow */
  if (strncmp(StrEnv[snowpattern].VarStr, "TRUE", 4) == 0)
    Options->SnowPattern = TRUE;
  else if (strncmp(StrEnv[snowpattern].VarStr, "FALSE", 5) == 0)
    Options->SnowPattern = FALSE;
  else
    ReportError(StrEnv[snowpattern].KeyName, 51);
  if (Options->Prism == FALSE && Options->SnowPattern == TRUE) {
    Options->SnowPattern = FALSE;
    printf("\nError: SnowPattern must be FALSE if PRISM is FALSE; Setting SnowPattern = FALSE\n");
  }
  
  /* Determine the kind of canopy radiation attenuation to be used */
  if (strncmp(StrEnv[canopy_radatt].VarStr, "FIXED", 3) == 0)
    Options->CanopyRadAtt = FIXED;
  else if (strncmp(StrEnv[canopy_radatt].VarStr, "VARIABLE", 3) == 0)
    Options->CanopyRadAtt = VARIABLE;
  else
    ReportError(StrEnv[canopy_radatt].KeyName, 51);

  /* Determine if solar shading maps will be used */
  if (strncmp(StrEnv[shading].VarStr, "TRUE", 4) == 0)
    Options->Shading = TRUE;
  else if (strncmp(StrEnv[shading].VarStr, "FALSE", 5) == 0)
    Options->Shading = FALSE;
  else
    ReportError(StrEnv[shading].KeyName, 51);
  
  /* Determine if the improved radiation scheme will be used */
  if (strncmp(StrEnv[improv_radiation].VarStr, "TRUE", 4) == 0)
    Options->ImprovRadiation = TRUE;
  else if (strncmp(StrEnv[improv_radiation].VarStr, "FALSE", 5) == 0)
    Options->ImprovRadiation = FALSE;
  else
    ReportError(StrEnv[improv_radiation].KeyName, 51);

  /* Determine if canopy gapping will be modeled */
  if (strncmp(StrEnv[gapping].VarStr, "TRUE", 4) == 0)
    Options->CanopyGapping = TRUE;
  else if (strncmp(StrEnv[gapping].VarStr, "FALSE", 5) == 0)
    Options->CanopyGapping = FALSE;
  else
    ReportError(StrEnv[gapping].KeyName, 51);

  /* Determine if snow sliding will be modeled */
  if (strncmp(StrEnv[snowslide].VarStr, "TRUE", 4) == 0)
    Options->SnowSlide = TRUE;
  else if (strncmp(StrEnv[snowslide].VarStr, "FALSE", 5) == 0)
    Options->SnowSlide = FALSE;
  else
    ReportError(StrEnv[snowslide].KeyName, 51);
  
  /* Determine if dumps snow stats */
  if (strncmp(StrEnv[snowstats].VarStr, "TRUE", 4) == 0)
    Options->SnowStats = TRUE;
  else if (strncmp(StrEnv[snowstats].VarStr, "FALSE", 5) == 0)
    Options->SnowStats = FALSE;
  else
    ReportError(StrEnv[snowstats].KeyName, 51);
  
  /* Determine if uses dynamic vegetation */
  if (strncmp(StrEnv[dynaveg].VarStr, "TRUE", 4) == 0)
    Options->DynamicVeg = TRUE;
  else if (strncmp(StrEnv[dynaveg].VarStr, "FALSE", 5) == 0)
    Options->DynamicVeg = FALSE;
  else
    ReportError(StrEnv[dynaveg].KeyName, 51);
  
  /* Determine if use separate input of rain and snow */
  if (strncmp(StrEnv[sepr].VarStr, "TRUE", 4) == 0)
    Options->PrecipSepr = TRUE;
  else if (strncmp(StrEnv[sepr].VarStr, "FALSE", 5) == 0)
    Options->PrecipSepr = FALSE;
  else
    ReportError(StrEnv[sepr].KeyName, 51);
  
  /* Determine whether to save extra information in stream states */
  if (strncmp(StrEnv[streamdata].VarStr, "TRUE", 4) == 0)
    Options->DumpExtraStream = TRUE;
  else if (strncmp(StrEnv[streamdata].VarStr, "FALSE", 5) == 0)
    Options->DumpExtraStream = FALSE;
  else
    ReportError(StrEnv[streamdata].KeyName, 51);
  
  /* Determine whether to spinup groundwater state before starting model run */
  if (strncmp(StrEnv[gw_spinup].VarStr, "TRUE", 4) == 0)
    Options->GW_SPINUP = TRUE;
  else if (strncmp(StrEnv[gw_spinup].VarStr, "FALSE", 5) == 0)
    Options->GW_SPINUP = FALSE;
  else
    ReportError(StrEnv[gw_spinup].KeyName, 51);
  
  if (Options->GW_SPINUP == TRUE) {
    if (!CopyInt(&(Options->GW_SPINUP_YRS), StrEnv[gw_spinup_yrs].VarStr, 1))
      ReportError(StrEnv[gw_spinup_yrs].KeyName, 51);
    if (!CopyFloat(&(Options->GW_SPINUP_RECHARGE), StrEnv[gw_spinup_recharge].VarStr, 1))
      ReportError(StrEnv[gw_spinup_recharge].KeyName, 51);
  }
  
  /* If canopy gapping option is true, the improved radiation scheme must be true */
  if (Options->CanopyGapping == TRUE && Options->ImprovRadiation == FALSE) {
    ReportError(StrEnv[gapping].KeyName, 71);
  }
  /* Determine if listed met stations outside bounding box are used */
  if (strncmp(StrEnv[outside].VarStr, "TRUE", 4) == 0)
    Options->Outside = TRUE;
  else if (strncmp(StrEnv[outside].VarStr, "FALSE", 5) == 0)
    Options->Outside = FALSE;
  else
    ReportError(StrEnv[outside].KeyName, 51);

  /* The file path to PRISM files */
  if (Options->Prism == TRUE) {
    if (IsEmptyStr(StrEnv[prism_data_path].VarStr))
      ReportError(StrEnv[prism_data_path].KeyName, 51);
    strcpy(Options->PrismDataPath, StrEnv[prism_data_path].VarStr);
    if (IsEmptyStr(StrEnv[prism_data_ext].VarStr))
      ReportError(StrEnv[prism_data_ext].KeyName, 51);
    strcpy(Options->PrismDataExt, StrEnv[prism_data_ext].VarStr);
  }

  /* The file path to snow pattern map */
  if (Options->SnowPattern == TRUE) {
    if (IsEmptyStr(StrEnv[snowpattern_data_path].VarStr))
      ReportError(StrEnv[snowpattern_data_path].KeyName, 51);
    strcpy(Options->SnowPatternDataPath, StrEnv[snowpattern_data_path].VarStr);
  }
  
  if (Options->Shading == TRUE) {
    if (IsEmptyStr(StrEnv[shading_data_path].VarStr))
      ReportError(StrEnv[shading_data_path].KeyName, 51);
    strcpy(Options->ShadingDataPath, StrEnv[shading_data_path].VarStr);
    if (IsEmptyStr(StrEnv[shading_data_ext].VarStr))
      ReportError(StrEnv[shading_data_ext].KeyName, 51);
    strcpy(Options->ShadingDataExt, StrEnv[shading_data_ext].VarStr);
    if (IsEmptyStr(StrEnv[skyview_data_path].VarStr))
      ReportError(StrEnv[skyview_data_path].KeyName, 51);
    strcpy(Options->SkyViewDataPath, StrEnv[skyview_data_path].VarStr);
  }

  /* Determine if rh override is used */
  if (strncmp(StrEnv[rhoverride].VarStr, "TRUE", 4) == 0)
    Options->Rhoverride = TRUE;
  else if (strncmp(StrEnv[rhoverride].VarStr, "FALSE", 5) == 0)
    Options->Rhoverride = FALSE;
  else
    ReportError(StrEnv[rhoverride].KeyName, 51);

  /* Determine the type of temperature lapse rate */
  if (strncmp(StrEnv[temp_lapse].VarStr, "CONSTANT", 8) == 0)
    Options->TempLapse = CONSTANT;
  else if (strncmp(StrEnv[temp_lapse].VarStr, "VARIABLE", 8) == 0)
    Options->TempLapse = VARIABLE;
  else
    ReportError(StrEnv[temp_lapse].KeyName, 51);
  
  /**************** Determine areal extent ****************/

  if (IsEmptyStr(StrEnv[coordinate_system].VarStr))
    ReportError(StrEnv[coordinate_system].KeyName, 51);
  strcpy(Map->System, StrEnv[coordinate_system].VarStr);

  if (!CopyDouble(&(Map->Yorig), StrEnv[extreme_north].VarStr, 1))
    ReportError(StrEnv[extreme_north].KeyName, 51);

  if (!CopyDouble(&(Map->Xorig), StrEnv[extreme_west].VarStr, 1))
    ReportError(StrEnv[extreme_west].KeyName, 51);

  if (!CopyFloat(&(SolarGeo->Latitude), StrEnv[center_latitude].VarStr, 1))
    ReportError(StrEnv[center_latitude].KeyName, 51);
  SolarGeo->Latitude *= (float) RADPDEG;

  if (!CopyFloat(&(SolarGeo->Longitude), StrEnv[center_longitude].VarStr, 1))
    ReportError(StrEnv[center_longitude].KeyName, 51);
  SolarGeo->Longitude *= (float) RADPDEG;

  if (!CopyFloat(&(SolarGeo->StandardMeridian),
		 StrEnv[time_zone_meridian].VarStr, 1))
    ReportError(StrEnv[time_zone_meridian].KeyName, 51);
  SolarGeo->StandardMeridian *= (float) RADPDEG;

  if (!CopyInt(&(Map->NY), StrEnv[number_of_rows].VarStr, 1))
    ReportError(StrEnv[number_of_rows].KeyName, 51);
  if (Map->NY < 1 || Map->NY > 1000000)
    ReportError(StrEnv[number_of_rows].KeyName, 51);
  
  if (!CopyInt(&(Map->NX), StrEnv[number_of_columns].VarStr, 1))
    ReportError(StrEnv[number_of_columns].KeyName, 51);
  if (Map->NX < 1 || Map->NX > 1000000)
    ReportError(StrEnv[number_of_rows].KeyName, 51);

  if (!CopyFloat(&(Map->DY), StrEnv[grid_spacing].VarStr, 1))
    ReportError(StrEnv[grid_spacing].KeyName, 51);
  
  Options->MaxInterpDist /= Map->DY;
  
  Map->DX = Map->DY;
  Map->DXY = (float) sqrt(Map->DX * Map->DX + Map->DY * Map->DY);
  Map->X = 0;
  Map->Y = 0;
  Map->OffsetX = 0;
  Map->OffsetY = 0;
  Map->NumCells = 0;

  if (Options->Extent == POINT) {
    if (!CopyDouble(&PointModelY, StrEnv[point_north].VarStr, 1))
      ReportError(StrEnv[point_north].KeyName, 51);
    if (!CopyDouble(&PointModelX, StrEnv[point_east].VarStr, 1))
      ReportError(StrEnv[point_east].KeyName, 51);
    Options->PointY =
      Round(((Map->Yorig - 0.5 * Map->DY) - PointModelY) / Map->DY);
    Options->PointX =
      Round((PointModelX - (Map->Xorig + 0.5 * Map->DX)) / Map->DX);
  }
  else {
    Options->PointY = 0;
    Options->PointX = 0;
  }

  /**************** Determine model period ****************/

  if (!CopyFloat(&(TimeStep), StrEnv[time_step].VarStr, 1))
    ReportError(StrEnv[time_step].KeyName, 51);
  TimeStep *= SECPHOUR;

  if (!SScanDate(StrEnv[model_start].VarStr, &(Start)))
    ReportError(StrEnv[model_start].KeyName, 51);

  if (!SScanDate(StrEnv[model_end].VarStr, &(End)))
    ReportError(StrEnv[model_end].KeyName, 51);

  InitTime(Time, &Start, &End, (int) TimeStep);

   /**************** Determine model constants ****************/

  if (!CopyFloat(&Z0_GROUND, StrEnv[ground_roughness].VarStr, 1))
    ReportError(StrEnv[ground_roughness].KeyName, 51);

  if (!CopyFloat(&Z0_SNOW, StrEnv[snow_roughness].VarStr, 1))
    ReportError(StrEnv[snow_roughness].KeyName, 51);

  if (!CopyFloat(&LIQUID_WATER_CAPACITY, StrEnv[snow_water_capacity].VarStr, 1))
    ReportError(StrEnv[snow_water_capacity].KeyName, 51);

  if (!CopyFloat(&Zref, StrEnv[reference_height].VarStr, 1))
    ReportError(StrEnv[reference_height].KeyName, 51);

  if (!CopyFloat(&LAI_WATER_MULTIPLIER, StrEnv[rain_lai_multiplier].VarStr, 1))
    ReportError(StrEnv[rain_lai_multiplier].KeyName, 51);

  if (!CopyFloat(&LAI_SNOW_MULTIPLIER, StrEnv[snow_lai_multiplier].VarStr, 1))
    ReportError(StrEnv[snow_lai_multiplier].KeyName, 51);

  if (!CopyFloat(&MIN_INTERCEPTION_STORAGE,
		 StrEnv[min_intercepted_snow].VarStr, 1))
    ReportError(StrEnv[min_intercepted_snow].KeyName, 51);

  if (!CopyFloat(&MIN_SNOW_RESET_ALBEDO,
		 StrEnv[min_snow_reset_albedo].VarStr, 1))
    ReportError(StrEnv[min_snow_reset_albedo].KeyName, 51);
  
  if (!CopyFloat(&LAMBDA_FOREST_OFFSET,
		 StrEnv[albedo_lambda_forest].VarStr, 1))
    ReportError(StrEnv[albedo_lambda_forest].KeyName, 51);
  
  if (!CopyUChar(&OUTSIDEBASIN, StrEnv[outside_basin].VarStr, 1))
    ReportError(StrEnv[outside_basin].KeyName, 51);

  if (Options->TempLapse == CONSTANT) {
    if (!CopyFloat(&TEMPLAPSE, StrEnv[temp_lapse_rate].VarStr, 1))
      ReportError(StrEnv[temp_lapse_rate].KeyName, 51);
  }
  else
    TEMPLAPSE = NOT_APPLICABLE;

  if (Options->SnowPattern == TRUE) {
    if (!CopyFloat(&SNOWPAT_WEIGHT, StrEnv[snowpattern_weight].VarStr, 1))
      ReportError(StrEnv[snowpattern_weight].KeyName, 51);
  }
  else
    SNOWPAT_WEIGHT = NOT_APPLICABLE;
  
  if (!CopyFloat(&TEMPERATURE_OFFSET, StrEnv[temperature_offset].VarStr, 1))
    ReportError(StrEnv[temperature_offset].KeyName, 51);
  
  if (!CopyFloat(&LAPSE_RATE_BIAS, StrEnv[lapse_bias].VarStr, 1))
    ReportError(StrEnv[lapse_bias].KeyName, 51);
  
  if (!CopyFloat(&LAPSE_BIAS_ELEV, StrEnv[lapse_elev].VarStr, 1))
    ReportError(StrEnv[lapse_elev].KeyName, 51);
  
  printf("\nTemperature inputs are adjusted by %.3f deg. C\n",TEMPERATURE_OFFSET);
  printf("Temperature also adjusted by %.3e deg. C / m rel. to %.0f m\n",
         LAPSE_RATE_BIAS,LAPSE_BIAS_ELEV);
  
  if (!CopyFloat(&SOIL_DEPTH_ADJ, StrEnv[soil_depth_adj].VarStr, 1))
    ReportError(StrEnv[soil_depth_adj].KeyName, 51);
  
  if (!CopyFloat(&SOIL_KSAT_ADJ, StrEnv[soil_ksat_adj].VarStr, 1))
    ReportError(StrEnv[soil_ksat_adj].KeyName, 51);
  
  if (!CopyFloat(&SOIL_EXPDEC_ADJ, StrEnv[soil_expdec_adj].VarStr, 1))
    ReportError(StrEnv[soil_expdec_adj].KeyName, 51);
  
  if (!CopyFloat(&SOIL_POROSITY_ADJ, StrEnv[soil_porosity_adj].VarStr, 1))
    ReportError(StrEnv[soil_porosity_adj].KeyName, 51);
  
  if (!CopyFloat(&SOIL_FIELDCAP_ADJ, StrEnv[soil_fieldcap_adj].VarStr, 1))
    ReportError(StrEnv[soil_fieldcap_adj].KeyName, 51);
  
  if (!CopyFloat(&VEG_LAI_ADJ, StrEnv[veg_lai_adj].VarStr, 1))
    ReportError(StrEnv[veg_lai_adj].KeyName, 51);
  
  if (!CopyFloat(&VEG_FC_ADJ, StrEnv[veg_fc_adj].VarStr, 1))
    ReportError(StrEnv[veg_fc_adj].KeyName, 51);
  
  /* maximum depth of the surface layer in snow water equivalent (m) */
  if (!CopyFloat(&MAX_SURFACE_SWE,
    StrEnv[max_swe].VarStr, 1))
    ReportError(StrEnv[max_swe].KeyName, 51);

  /* if turn on canopy gap module */
  if (Options->CanopyGapping) {
    if (!CopyFloat(&GAPWIND_FACTOR, StrEnv[gapwind_adj].VarStr, 1))
      ReportError(StrEnv[gapwind_adj].KeyName, 51);
    if (GAPWIND_FACTOR <= 0 || GAPWIND_FACTOR>1)
      ReportError(StrEnv[gapwind_adj].KeyName, 74);
  }
  if (Options->SnowSlide) {
    if (!CopyFloat(&SNOWSLIDE1, StrEnv[snowslide_parameter1].VarStr, 1))
      ReportError(StrEnv[snowslide_parameter1].KeyName, 51);

    if (!CopyFloat(&SNOWSLIDE2, StrEnv[snowslide_parameter2].VarStr, 1))
      ReportError(StrEnv[snowslide_parameter2].KeyName, 51);
  }
}


/******************************************************************************
 InitMappedConstants
 ******************************************************************************/
void
InitMappedConstants(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map,
                    SNOWPIX ***SnowMap, VEGTABLE *VType, VEGPIX ***VegMap)
{
  STRINIENTRY StrEnv[] =
    {
     {"CONSTANTS", "RAIN THRESHOLD", "", ""},
     {"CONSTANTS", "SNOW THRESHOLD", "", ""},
     {"CONSTANTS", "FRESH SNOW ALBEDO", "", "0.85" },                                         
     {"CONSTANTS", "ALBEDO ACCUMULATION LAMBDA", "", "" },
     {"CONSTANTS", "ALBEDO MELTING LAMBDA", "", "" },
     {"CONSTANTS", "ALBEDO ACCUMULATION MIN", "", "" },
     {"CONSTANTS", "ALBEDO MELTING MIN", "", "" },
     {"CONSTANTS", "PRECIPITATION MULTIPLIER MAP", "", "" },
     {NULL, NULL, "", NULL}
    };
  int i, x, y;
  int MapId;
  int ParamType;
  char FileName[BUFSIZE + 1];	      /* Variable name */
  
  /* Read the key-entry pairs from the input file */
  for (i = 0; StrEnv[i].SectionName; i++)
    GetInitString(StrEnv[i].SectionName, StrEnv[i].KeyName, StrEnv[i].Default,
		  StrEnv[i].VarStr, (unsigned long) BUFSIZE, Input);

  /****** snow parameters (take either spatial input or constants) *****/
  if (IsEmptyStr(StrEnv[rain_threshold].VarStr)) {
    ReportError(StrEnv[rain_threshold].KeyName, 51);
  }
  else {
	MapId = 801;
	if (!CopyFloat(&MIN_RAIN_TEMP, StrEnv[rain_threshold].VarStr, 1)) {
	  printf("%s: spatial parameters are used\n", StrEnv[rain_threshold].KeyName);
	  ParamType = MAP;
	  strcpy(FileName, StrEnv[rain_threshold].VarStr);
	}
	else
	  ParamType = CONSTANT;

	/* Initiate spatial input of parameters */
	InitParameterMaps(Options, Map, MapId, FileName, SnowMap, ParamType, MIN_RAIN_TEMP);
  }

  if (IsEmptyStr(StrEnv[snow_threshold].VarStr)) {
    ReportError(StrEnv[snow_threshold].KeyName, 51);
  }
  else {
	MapId = 800;
    if (!CopyFloat(&MAX_SNOW_TEMP, StrEnv[snow_threshold].VarStr, 1)) {
      printf("%s: spatial parameters are used\n", StrEnv[snow_threshold].KeyName);
	  ParamType = MAP;
      strcpy(FileName, StrEnv[snow_threshold].VarStr);      
    }
	else 
	  ParamType = CONSTANT;
	InitParameterMaps(Options, Map, MapId, FileName, SnowMap, ParamType, MAX_SNOW_TEMP);
  }

  if (IsEmptyStr(StrEnv[alb_acc_lambda].VarStr)) {
    ReportError(StrEnv[alb_acc_lambda].KeyName, 51);
  }
  else {
	MapId = 803;
    if (!CopyFloat(&ALB_ACC_LAMBDA, StrEnv[alb_acc_lambda].VarStr, 1)) {
      printf("%s: spatial parameters are used\n", StrEnv[alb_acc_lambda].KeyName);
	  ParamType = MAP;
      strcpy(FileName, StrEnv[alb_acc_lambda].VarStr);
    }
	else
	  ParamType = CONSTANT;
	InitParameterMaps(Options, Map, MapId, FileName, SnowMap, ParamType, ALB_ACC_LAMBDA);
  }


  if (IsEmptyStr(StrEnv[alb_melt_lambda].VarStr)) {
    ReportError(StrEnv[alb_melt_lambda].KeyName, 51);
  }
  else {
	MapId = 804;
	if (!CopyFloat(&ALB_MELT_LAMBDA, StrEnv[alb_melt_lambda].VarStr, 1)) {
	  printf("%s: spatial parameters are used\n", StrEnv[alb_melt_lambda].KeyName);
	  ParamType = MAP;
	  strcpy(FileName, StrEnv[alb_melt_lambda].VarStr);
	}
	else
	  ParamType = CONSTANT;
	InitParameterMaps(Options, Map, MapId, FileName, SnowMap, ParamType, ALB_MELT_LAMBDA);
  }
  
  /* Adjust albedo melt/acc lambda in forested regions */
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (VType[(*VegMap)[y][x].Veg - 1].OverStory == TRUE) {
        (*SnowMap)[y][x].LamdaAcc += (LAMBDA_FOREST_OFFSET * (*VegMap)[y][x].Fract[0]);
        (*SnowMap)[y][x].LamdaMelt += (LAMBDA_FOREST_OFFSET * (*VegMap)[y][x].Fract[0]);
      }
    }
  }
  
  if (IsEmptyStr(StrEnv[alb_acc_min].VarStr)) {
    ReportError(StrEnv[alb_acc_min].KeyName, 51);
  }
  else {
	MapId = 805;
    if (!CopyFloat(&ALB_ACC_MIN, StrEnv[alb_acc_min].VarStr, 1)) {
      printf("%s: spatial parameters are used\n", StrEnv[alb_acc_min].KeyName);
	  ParamType = MAP;
      strcpy(FileName, StrEnv[alb_acc_min].VarStr);
	}
	else
	  ParamType = CONSTANT;
	InitParameterMaps(Options, Map, MapId, FileName, SnowMap, ParamType, ALB_ACC_MIN);
  }

  if (IsEmptyStr(StrEnv[alb_melt_min].VarStr)) {
    ReportError(StrEnv[alb_melt_min].KeyName, 51);
  }
  else {
	MapId = 806;
	if (!CopyFloat(&ALB_MELT_MIN, StrEnv[alb_melt_min].VarStr, 1)) {
	  printf("%s: spatial parameters are used\n", StrEnv[alb_melt_min].KeyName);
	  ParamType = MAP;
	  strcpy(FileName, StrEnv[alb_melt_min].VarStr);
	}
	else
	  ParamType = CONSTANT;
	InitParameterMaps(Options, Map, MapId, FileName, SnowMap, ParamType, ALB_MELT_MIN);
  }

  /* fresh albedo - this was made a constant 0.85 in previous versions */
  if (IsEmptyStr(StrEnv[fresh_alb].VarStr)) {
    ReportError(StrEnv[fresh_alb].KeyName, 51);
  }
  else {
	MapId = 802;
    if (!CopyFloat(&ALB_MAX, StrEnv[fresh_alb].VarStr, 1)) {
      printf("%s: spatial parameters are used\n", StrEnv[fresh_alb].KeyName);

	  ParamType = MAP;
      strcpy(FileName, StrEnv[fresh_alb].VarStr);
	}
	else
	  ParamType = CONSTANT;
	InitParameterMaps(Options, Map, MapId, FileName, SnowMap, ParamType, ALB_MAX);
  }

  /* precipitation multiplier that bias correct the precipitation */
  strcpy(Options->PrecipMultiplierMapPath, "");
  if (IsEmptyStr(StrEnv[multiplier].VarStr)) {
    PRECIP_MULTIPLIER = 0;
    printf("No input of precipitation multiplier map - no correction is made\n");
  }
  else {
    if (!CopyFloat(&PRECIP_MULTIPLIER, StrEnv[multiplier].VarStr, 1)) {
      printf("%s: spatial parameters are used\n", StrEnv[multiplier].KeyName);
      PRECIP_MULTIPLIER = NA;
      strcpy(Options->PrecipMultiplierMapPath, StrEnv[multiplier].VarStr);
    }
    else
      printf("Precipitation inputs are rescaled by a factor of %.3f\n",PRECIP_MULTIPLIER);
  }
}
