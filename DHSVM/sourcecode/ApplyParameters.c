/*
 * DESCRIPTION:  The [PARAMETERS] values applied to the class and lake tables
 *               as read (ApplyTableParameters), then to the input maps, once,
 *               at initialization, each by its one formula, in the order
 *               numbered here:
 *                 T1 DeepFlux = (soil table DeepFlux + offset) / (365 x steps
 *                   per day), every soil class
 *                 T2 RsMin = vegetation table minimum resistance x scale,
 *                   every class and layer
 *                 T3 SnowIntEff = min(vegetation table SnowInterceptionEff
 *                   x scale, 1), every overstory class
 *                 T4 PowLawScale = lake table OutflowScale x scale, every
 *                   lake (Lake Dynamics)
 *                 T5 PowLawExponent = lake table OutflowExponent x scale,
 *                   every lake (Lake Dynamics)
 *                 1 Depth = max(Soil.Depth + offset, root zone depth + 0.01)
 *                 2 KsLat = Soil.LateralConductivity x scale
 *                 3 KsLatExp = Soil.ExponentialDecrease x scale
 *                   then the vertical conductivity, the maximum infiltration
 *                   rate and the deep flux of each cell
 *                 4 FCap[l] = max(Soil.FieldCapacity[l] x scale, WP[l] + 0.01)
 *                 5 Porosity[l] = max(Soil.Porosity[l] x scale, FCap[l] + 0.01)
 *                 6 overstory cover = min(Veg.OverstoryCover x scale, 1)
 *                 7 overstory LAI of a month = Veg.OverstoryLAI x the class
 *                   fraction of the month x scale
 *                 8-15 the snow parameters of every cell, with the forest
 *                   offset of the albedo decay rates
 *               and the melt scale of every cell from Snow.MeltScale (1
 *               without it).
 *               Only cells inside the basin mask are set. The input maps
 *               were checked before (Validate.c), so the bounds of rows 1,
 *               4, 5 and 6 act only through the parameters.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "Calendar.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"

static void ApplySoilParameters(OPTIONSTRUCT *Options, MAPSIZE *Map, INPUTMAPS *Raw,
  TOPOPIX **TopoMap, LAYER *Soil, SOILTABLE *SType, VEGTABLE *VType, SOILPIX ***SoilMap,
  int Dt);
static void ApplyVegParameters(OPTIONSTRUCT *Options, MAPSIZE *Map, INPUTMAPS *Raw,
  TOPOPIX **TopoMap, VEGTABLE *VType, VEGPIX ***VegMap);
static void ApplySnowParameters(MAPSIZE *Map, INPUTMAPS *Raw, TOPOPIX **TopoMap,
  VEGTABLE *VType, VEGPIX **VegMap, SNOWPIX **SnowMap);

/*****************************************************************************
  ApplyTableParameters()

  Rows T1-T5, on the class and lake tables as read, before the maps: the
  deep flux offset of every soil class (m/yr, then m per time step), the
  minimum stomatal resistance of every vegetation layer, the snow
  interception efficiency of every overstory class (at most 1), the outflow
  scale and exponent of every lake
*****************************************************************************/
void ApplyTableParameters(OPTIONSTRUCT *Options, int NDaySteps, SOILTABLE *SType, int NSoils,
  VEGTABLE *VType, int NVegs, LAKETABLE *LType, int NLakes)
{
  int i, j;

  /* T1 Deep Flux Offset, then m/yr to m per time step */
  for (i = 0; i < NSoils; i++) {
    SType[i].DeepFlux += SOIL_DEEPFLUX_ADJ;
    SType[i].DeepFlux /= (DAYPYEAR * NDaySteps);
  }

  for (i = 0; i < NVegs; i++) {
    /* T2 Min Stomatal Resistance Scale: every layer (the maximum unchanged) */
    for (j = 0; j < VType[i].NVegLayers; j++)
      VType[i].RsMin[j] *= VEG_RSMIN_ADJ;

    /* T3 Snow Interception Efficiency Scale: at most 1 */
    if (VType[i].OverStory == TRUE) {
      VType[i].SnowIntEff *= VEG_SNOWINT_ADJ;
      if (VType[i].SnowIntEff > 1.0)
        VType[i].SnowIntEff = 1.0;
    }
  }

  /* T4 Lake Outflow Rate Scale, T5 Lake Outflow Exponent Scale */
  if (Options->LakeDynamics) {
    for (i = 0; i < NLakes; i++) {
      LType[i].PowLawScale *= LAKE_OUTFLOW_ADJ;
      LType[i].PowLawExponent *= LAKE_EXPONENT_ADJ;
    }
  }
}

/*****************************************************************************
  ApplyParameters()

  Sets up the soil and vegetation maps of every cell from the input maps and
  the [PARAMETERS] values, and the snow parameters of every cell. Frees the
  input maps it uses.
*****************************************************************************/
void ApplyParameters(OPTIONSTRUCT *Options, MAPSIZE *Map, INPUTMAPS *Raw,
  TOPOPIX **TopoMap, LAYER *Soil, SOILTABLE *SType, VEGTABLE *VType,
  SOILPIX ***SoilMap, VEGPIX ***VegMap, SNOWPIX **SnowMap, int Dt)
{
  Message("\nApplying parameters\n");

  ApplySoilParameters(Options, Map, Raw, TopoMap, Soil, SType, VType, SoilMap, Dt);
  ApplyVegParameters(Options, Map, Raw, TopoMap, VType, VegMap);
  ApplySnowParameters(Map, Raw, TopoMap, VType, *VegMap, SnowMap);
}

/*****************************************************************************
  ApplySoilParameters()

  Rows 1-5 and the values derived from them
*****************************************************************************/
static void ApplySoilParameters(OPTIONSTRUCT *Options, MAPSIZE *Map, INPUTMAPS *Raw,
  TOPOPIX **TopoMap, LAYER *Soil, SOILTABLE *SType, VEGTABLE *VType, SOILPIX ***SoilMap,
  int Dt)
{
  const char *Routine = "ApplySoilParameters";
  int i;			/* counter */
  int x;			/* counter */
  int y;			/* counter */
  int NCells = Map->NY * Map->NX;
  int NLayers = Map->SoilLayers;	/* root zone layers; the deep layer is index NLayers */
  int NSet, VSet;
  int sidx;
  SOILPIX *Cell;
  VEGTABLE *Veg;
  float TotalRootDepth, LayerDepth, Transmissivity, KsVertCalc;
  double MaxSeepage;		/* Largest seepage loss per time step (negative) */

  if (!(*SoilMap = (SOILPIX **)calloc(Map->NY, sizeof(SOILPIX *))))
    ReportError((char *)Routine, 1);
  for (y = 0; y < Map->NY; y++) {
    if (!((*SoilMap)[y] = (SOILPIX *)calloc(Map->NX, sizeof(SOILPIX))))
      ReportError((char *)Routine, 1);
  }

  for (y = 0, i = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++, i++) {
      Cell = &((*SoilMap)[y][x]);
      Cell->Soil = Raw->SoilClass[i];
      if (!INBASIN(TopoMap[y][x].Mask))
        continue;
      sidx = Cell->Soil - 1;
      Veg = &(VType[Raw->VegClass[i] - 1]);

      /* allocate memory for the root zone layers and the deep layer */
      if (!(Cell->KsVert = (float *)calloc(NLayers + 1, sizeof(float))) ||
          !(Cell->FCap = (float *)calloc(NLayers + 1, sizeof(float))) ||
          !(Cell->Porosity = (float *)calloc(NLayers + 1, sizeof(float))) ||
          !(Cell->Moist = (float *)calloc(NLayers + 1, sizeof(float))) ||
          !(Cell->Perc = (float *)calloc(NLayers, sizeof(float))) ||
          !(Cell->InterFlow = (float *)calloc(NLayers + 1, sizeof(float))) ||
          !(Cell->Temp = (float *)calloc(NLayers, sizeof(float))))
        ReportError((char *)Routine, 1);
      if (Options->Infiltration == DYNAMIC)
        Cell->InfiltAcc = 0.;
      Cell->MoistInit = 0.;

      /* 1 Soil Depth Offset: at least 0.01 m deeper than the root zone */
      Cell->Depth = Raw->Depth[i];
      Cell->Depth += SOIL_DEPTH_ADJ;
      TotalRootDepth = 0.0;
      for (NSet = 0; NSet < NLayers; NSet++)
        TotalRootDepth += Veg->RootDepth[NSet];
      if (Cell->Depth < (TotalRootDepth + 0.01))
        Cell->Depth = (TotalRootDepth + 0.01);

      /* 2 Lateral Conductivity Scale */
      Cell->KsLat = Raw->LateralConductivity[i];
      Cell->KsLat *= SOIL_KSAT_ADJ;

      /* 3 Exponential Decrease Scale */
      Cell->KsLatExp = Raw->ExponentialDecrease[i];
      Cell->KsLatExp *= SOIL_EXPDEC_ADJ;

      /* Vertical conductivity of each root zone layer: the transmissivity of
         the layer over its thickness and the anisotropy, or the soil table */
      for (NSet = 0; NSet < NLayers; NSet++) {
        if (Options->UseKsatAnisotropy) {
          /* Get total depth of current soil layer bottom */
          LayerDepth = 0.0;
          for (VSet = 0; VSet <= NSet; VSet++)
            LayerDepth += Veg->RootDepth[VSet];

          /* Calculate effective lateral conductivity of current soil layer */
          Transmissivity =
            CalcTransmissivity(LayerDepth, LayerDepth - Veg->RootDepth[NSet],
                               Cell->KsLat, Cell->KsLatExp, SType[sidx].DepthThresh);
          KsVertCalc = Transmissivity / Veg->RootDepth[NSet];

          /* Account for user-supplied vertical anisotropy */
          KsVertCalc /= SType[sidx].KsAnisotropy;

          Cell->KsVert[NSet] = KsVertCalc;
        }
        else
          Cell->KsVert[NSet] = SType[sidx].Ks[NSet];
      }

      /* Vertical conductivity of the deep layer */
      NSet = NLayers;
      if (Options->UseKsatAnisotropy) {
        /* Calculate effective lateral conductivity of current soil layer */
        Transmissivity =
          CalcTransmissivity(Cell->Depth, Cell->Depth - Veg->TotalDepth,
                             Cell->KsLat, Cell->KsLatExp, SType[sidx].DepthThresh);
        KsVertCalc = Transmissivity / (Cell->Depth - Veg->TotalDepth);

        /* Account for user-supplied vertical anisotropy */
        KsVertCalc /= SType[sidx].KsAnisotropy;

        Cell->KsVert[NSet] = KsVertCalc;
      }
      else
        Cell->KsVert[NSet] = SType[sidx].Ks[NSet];

      /* Maximum infiltration rate: surface layer KsLat over the anisotropy,
         or the soil table */
      if (Options->UseKsatAnisotropy) {
        Cell->MaxInfiltrationRate = Cell->KsLat;
        Cell->MaxInfiltrationRate /= SType[sidx].KsAnisotropy;
      }
      else
        Cell->MaxInfiltrationRate = SType[sidx].MaxInfiltrationRate;

      /* Deep flux of the soil class; a seepage loss (negative deep flux) is
         limited to the deep layer's vertical conductivity over one time step */
      Cell->MaxDeepFlux = SType[sidx].DeepFlux;
      if (Cell->MaxDeepFlux < 0.0) {
        MaxSeepage = Cell->KsVert[NLayers] * -1.0;
        MaxSeepage *= Dt;
        if (Cell->MaxDeepFlux < MaxSeepage)
          Cell->MaxDeepFlux = MaxSeepage;
      }

      /* 4 Field Capacity Scale: at least 0.01 above the wilting point */
      for (NSet = 0; NSet <= NLayers; NSet++) {
        Cell->FCap[NSet] = Raw->FieldCapacity[NSet * NCells + i] * SOIL_FIELDCAP_ADJ;
        if (Cell->FCap[NSet] < (SType[sidx].WP[NSet] + 0.01))
          Cell->FCap[NSet] = SType[sidx].WP[NSet] + 0.01;
      }

      /* 5 Porosity Scale: at least 0.01 above the field capacity */
      for (NSet = 0; NSet <= NLayers; NSet++) {
        Cell->Porosity[NSet] = Raw->Porosity[NSet * NCells + i] * SOIL_POROSITY_ADJ;
        if (Cell->Porosity[NSet] < (Cell->FCap[NSet] + 0.01))
          Cell->Porosity[NSet] = Cell->FCap[NSet] + 0.01;
      }
    }
  }

  free(Raw->SoilClass);
  free(Raw->Depth);
  free(Raw->LateralConductivity);
  free(Raw->ExponentialDecrease);
  free(Raw->FieldCapacity);
  free(Raw->Porosity);
  Raw->SoilClass = NULL;
  Raw->Depth = Raw->LateralConductivity = Raw->ExponentialDecrease = NULL;
  Raw->FieldCapacity = Raw->Porosity = NULL;
}

/*****************************************************************************
  ApplyVegParameters()

  Rows 6 and 7: the overstory cover and LAI of overstory cells from the maps,
  the understory cover 1.0 (lowered later on channel cells), the understory
  LAI from the vegetation table; the heights from the height maps and the
  wind profile and aerodynamic resistances of every mask cell
*****************************************************************************/
static void ApplyVegParameters(OPTIONSTRUCT *Options, MAPSIZE *Map, INPUTMAPS *Raw,
  TOPOPIX **TopoMap, VEGTABLE *VType, VEGPIX ***VegMap)
{
  const char *Routine = "ApplyVegParameters";
  int i;			/* counter */
  int x;			/* counter */
  int y;			/* counter */
  int j;			/* counter */
  int NSet; /*Counter for LAI month*/
  VEGPIX *Cell;
  VEGTABLE *Veg;

  if (!(*VegMap = (VEGPIX **)calloc(Map->NY, sizeof(VEGPIX *))))
    ReportError((char *)Routine, 1);
  for (y = 0; y < Map->NY; y++) {
    if (!((*VegMap)[y] = (VEGPIX *)calloc(Map->NX, sizeof(VEGPIX))))
      ReportError((char *)Routine, 1);
  }

  for (y = 0, i = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++, i++) {
      Cell = &((*VegMap)[y][x]);
      Cell->Veg = Raw->VegClass[i];
      Cell->Tcanopy = 0.0;
      if (!INBASIN(TopoMap[y][x].Mask))
        continue;
      Veg = &(VType[Cell->Veg - 1]);

      if (!(Cell->Fract = (float *)calloc(Veg->NVegLayers, sizeof(float))) ||
          !(Cell->LAIMonthly = (float **)calloc(Veg->NVegLayers, sizeof(float *))) ||
          !(Cell->LAI = (float *)calloc(Veg->NVegLayers, sizeof(float))) ||
          !(Cell->MaxInt = (float *)calloc(Veg->NVegLayers, sizeof(float))) ||
          !(Cell->Height = (float *)calloc(Veg->NVegLayers, sizeof(float))))
        ReportError((char *)Routine, 1);
      for (j = 0; j < Veg->NVegLayers; j++) {
        if (!(Cell->LAIMonthly[j] = (float *)calloc(12, sizeof(float))))
          ReportError((char *)Routine, 1);
      }

      /* 6 Overstory Cover Scale: at most 1; the understory cover is 1 */
      if (Veg->OverStory == TRUE) {
        Cell->Fract[0] = Raw->OverstoryCover[i];
        if (Veg->UnderStory == TRUE)
          Cell->Fract[1] = 1.0;
        Cell->Fract[0] *= VEG_FC_ADJ;
        /* Limit to maximum of 100% cover */
        if (Cell->Fract[0] > 1.0)
          Cell->Fract[0] = 1.0;
      }
      else {
        if (Veg->UnderStory == TRUE)
          Cell->Fract[0] = 1.0;
      }

      /* Canopy view factor */
      if (Options->ImprovRadiation == TRUE && Veg->NVegLayers > 0)
        Cell->Vf = Cell->Fract[0] * Veg->VfAdjust;

      /* 7 Overstory LAI Scale, January to December: the map times the
         class fraction of the month, times the scale; the understory LAI
         from the vegetation table */
      for (NSet = 0; NSet < 12; NSet++) {
        if (Veg->OverStory == TRUE) {
          Cell->LAIMonthly[0][NSet] = Raw->OverstoryLAI[i] * Veg->LAIFraction[NSet];
          if (Veg->UnderStory == TRUE)
            Cell->LAIMonthly[1][NSet] = Veg->LAIMonthly[1][NSet];
          Cell->LAIMonthly[0][NSet] *= VEG_LAI_ADJ;
        }
        else {
          if (Veg->UnderStory == TRUE)
            Cell->LAIMonthly[0][NSet] = Veg->LAIMonthly[0][NSet];
        }
      }

      /* Vegetation heights of the cell, layer order (overstory first), from
         the overstory and understory height maps */
      if (Veg->OverStory == TRUE) {
        Cell->Height[0] = Raw->OverstoryHeight[i];
        if (Veg->UnderStory == TRUE)
          Cell->Height[1] = Raw->UnderstoryHeight[i];
      }
      else {
        if (Veg->UnderStory == TRUE)
          Cell->Height[0] = Raw->UnderstoryHeight[i];
      }

      /* Wind profile and aerodynamic resistances of the cell's layers,
         normalized for a reference height wind speed of 1 m/s (scaled by
         the actual wind each time step) */
      CalcAerodynamic(y, x, Veg->NVegLayers, Veg->OverStory, Veg->Cn, Cell->Height, Veg->Trunk,
                      Cell->U, &(Cell->USnow), Cell->Ra, &(Cell->RaSnow));
    }
  }

  free(Raw->VegClass);
  free(Raw->OverstoryCover);
  free(Raw->OverstoryLAI);
  free(Raw->OverstoryHeight);
  free(Raw->UnderstoryHeight);
  Raw->VegClass = NULL;
  Raw->OverstoryCover = Raw->OverstoryLAI = Raw->OverstoryHeight = Raw->UnderstoryHeight = NULL;
}

/*****************************************************************************
  ApplySnowParameters()

  Rows 8-15: the rain and snow thresholds, the fresh snow albedo, the albedo
  decay rates (with the forest offset on overstory cells, after row 6) and
  minima of every cell; the melt scale of every mask cell from its map, or 1
*****************************************************************************/
static void ApplySnowParameters(MAPSIZE *Map, INPUTMAPS *Raw, TOPOPIX **TopoMap,
  VEGTABLE *VType, VEGPIX **VegMap, SNOWPIX **SnowMap)
{
  int x, y, i;

  for (y = 0, i = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++, i++) {
      if (INBASIN(TopoMap[y][x].Mask))
        SnowMap[y][x].MeltScale = (Raw->MeltScale != NULL) ? Raw->MeltScale[i] : 1.0;
      SnowMap[y][x].Tr = MIN_RAIN_TEMP;
      SnowMap[y][x].Ts = MAX_SNOW_TEMP;
      SnowMap[y][x].amax = ALB_MAX;
      SnowMap[y][x].LamdaAcc = ALB_ACC_LAMBDA;
      SnowMap[y][x].LamdaMelt = ALB_MELT_LAMBDA;
      SnowMap[y][x].AccMin = ALB_ACC_MIN;
      SnowMap[y][x].MeltMin = ALB_MELT_MIN;
    }
  }

  /* Adjust albedo melt/acc lambda in forested regions */
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask) && VType[VegMap[y][x].Veg - 1].OverStory == TRUE) {
        SnowMap[y][x].LamdaAcc += (LAMBDA_FOREST_OFFSET * VegMap[y][x].Fract[0]);
        SnowMap[y][x].LamdaMelt += (LAMBDA_FOREST_OFFSET * VegMap[y][x].Fract[0]);
      }
    }
  }

  free(Raw->MeltScale);
  Raw->MeltScale = NULL;
}
