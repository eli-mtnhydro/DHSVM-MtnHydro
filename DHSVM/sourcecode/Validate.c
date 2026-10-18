/*
 * DESCRIPTION:  Checks of the inputs before the run starts. The run stops
 *               at the first failure with one line on stderr and the exit
 *               code of the check:
 *               104 an input file or directory cannot be opened
 *               108 cell values: the input maps and the vegetation
 *                   updates as given, before any parameter applies; the
 *                   porosity after the parameters; the station weights
 *               109 the channel tables and the lakes
 *               Cell messages name the map, the condition, the number of
 *               cells and the first one (rows and columns from 1, row after
 *               row from north to south, each row from west to east).
 *               The summary of the classes in the basin is printed here.
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
#include "mapio.h"
#include "DHSVMChannel.h"

static const char *MonthNames[12] = {
  "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
};

/* Cells failing one condition: how many, and the first one */
typedef struct {
  long N;
  int Row;
  int Column;
  double Value[2];
} FAILURES;

static void Fail(FAILURES *F, int y, int x, double Value1, double Value2)
{
  if (F->N == 0) {
    F->Row = y + 1;
    F->Column = x + 1;
    F->Value[0] = Value1;
    F->Value[1] = Value2;
  }
  F->N++;
}

/* Stops the run if a condition failed: "<What>: <Condition> in <n> cells;
   first at row <r>, column <c> (<Name1> = <v>[, <Name2> = <v>])" */
static void ReportCells(FAILURES *F, const char *What, const char *Condition,
                        const char *Name1, const char *Name2)
{
  if (F->N == 0)
    return;
  if (Name2 == NULL)
    ReportErrorCode(108, "%s: %s in %ld cells; first at row %d, column %d (%s = %g)",
                    What, Condition, F->N, F->Row, F->Column, Name1, F->Value[0]);
  else
    ReportErrorCode(108, "%s: %s in %ld cells; first at row %d, column %d (%s = %g, %s = %g)",
                    What, Condition, F->N, F->Row, F->Column, Name1, F->Value[0], Name2,
                    F->Value[1]);
}

/* Name of a soil layer of a map: <Name>.1 .. <Name>.L, <Name>.Deep */
static void SoilLayerName(char *Buffer, size_t Size, const char *Name, int Layer,
                          int SoilLayers)
{
  if (Layer < SoilLayers)
    snprintf(Buffer, Size, "%s.%d", Name, Layer + 1);
  else
    snprintf(Buffer, Size, "%s.Deep", Name);
}

/*****************************************************************************
  CheckInputFiles()

  Every input file the run uses can be opened, and the station data
  directory exists (the station data files are opened with the station
  table). PointTable: TRUE when [OUTPUT] Point names a variable (a missing
  Point Table key is an output request error)
*****************************************************************************/
void CheckInputFiles(LISTPTR Input, OPTIONSTRUCT *Options, int PointTable)
{
  const char *Always[] = {
    "Terrain.Elevation", "Terrain.Mask", "Soil.Class", "Soil.Depth",
    "Soil.LateralConductivity", "Soil.ExponentialDecrease", "Soil.Porosity",
    "Soil.FieldCapacity", "Veg.Class", "Veg.OverstoryCover", "Veg.OverstoryLAI",
    "Veg.OverstoryHeight", "Veg.UnderstoryHeight", "Soil Table",
    "Vegetation Table", "Channel Class Table", "Channel Segment Table", "Channel Cell Table",
    "Station Table"
  };
  const char *Used[10];
  char Path[BUFSIZE + 1];
  FILE *File;
  int n = 0;
  int i;

  if (Options->Shading) {
    Used[n++] = "Terrain.HorizonAngle";
    Used[n++] = "Terrain.SkyView";
  }
  if (Options->Prism)
    Used[n++] = "Met.PrecipPattern";
  if (Options->SnowPattern)
    Used[n++] = "Met.SnowPattern";
  if (Options->MeltScale)
    Used[n++] = "Snow.MeltScale";
  if (Options->WindScale)
    Used[n++] = "Met.WindScale";
  if (Options->LakeDynamics) {
    Used[n++] = "Lake.ID";
    Used[n++] = "Lake Table";
  }
  if (PointTable && GetConfigValue(Input, "INPUTS", "Point Table", Path) && !IsEmptyStr(Path))
    Used[n++] = "Point Table";

  for (i = 0; i < (int) (sizeof(Always) / sizeof(Always[0])) + n; i++) {
    RequireConfigValue(Input, "INPUTS",
                       (i < (int) (sizeof(Always) / sizeof(Always[0]))) ? Always[i] :
                       Used[i - (int) (sizeof(Always) / sizeof(Always[0]))], Path);
    if (!(File = fopen(Path, "rb")))
      ReportErrorCode(104, "%s: cannot open", Path);
    fclose(File);
  }
  if (!IsEmptyStr(Options->InitialSWEPath)) {
    if (!(File = fopen(Options->InitialSWEPath, "rb")))
      ReportErrorCode(104, "%s: cannot open", Options->InitialSWEPath);
    fclose(File);
  }

  RequireConfigValue(Input, "INPUTS", "Station Data Directory", Path);
  if (!IsDirectory(Path))
    ReportErrorCode(104, "%s: cannot open", Path);
}

/*****************************************************************************
  ValidateVegMaps()

  The vegetation maps as given: the class of every mask cell in the
  vegetation table; no no-data in the overstory maps on overstory cells and
  in the understory height on understory cells; the overstory cover in
  (0, 1], the overstory LAI and the heights positive. With BaseClass (the
  classes of Veg.Class, for a vegetation update) a cell's class keeps the
  root depths of its class there. Suffix follows the map names in messages
*****************************************************************************/
void ValidateVegMaps(MAPSIZE *Map, unsigned char *Mask, unsigned char *VegClass, float *Cover,
  float *LAI, float *OverstoryHeight, float *UnderstoryHeight, unsigned char *BaseClass,
  LAYER *Veg, VEGTABLE *VType, const char *Suffix)
{
  int x, y, i, k;
  char Name[BUFSIZE + 1];
  char What[BUFSIZE + 1];
  FAILURES F;

#define VMASK(i) (Mask[i] == 1)
#define VCLASSOK(i) (VegClass[i] >= 1 && VegClass[i] <= Veg->NTypes)
#define VOVER(i) (VType[VegClass[i] - 1].OverStory == TRUE)
#define VUNDER(i) (VType[VegClass[i] - 1].UnderStory == TRUE)
#define VCHECK(Array, Where, MapName, Condition, Text, ValueName)    \
  memset(&F, 0, sizeof(F));                                         \
  for (y = 0, i = 0; y < Map->NY; y++)                              \
    for (x = 0; x < Map->NX; x++, i++)                              \
      if (VMASK(i) && (Where) && (Condition))                       \
        Fail(&F, y, x, (Array)[i], 0);                              \
  snprintf(What, sizeof(What), "%s%s", MapName, Suffix);            \
  ReportCells(&F, What, Text, ValueName, NULL);

  /* Classes: no no-data, in the table */
  VCHECK(VegClass, TRUE, "Veg.Class", VegClass[i] == MAP_NODATA_BYTE, "no-data", "Veg.Class");
  snprintf(Name, sizeof(Name), "class not in 1..%d (the vegetation table rows)", Veg->NTypes);
  VCHECK(VegClass, TRUE, "Veg.Class", !VCLASSOK(i), Name, "Veg.Class");

  /* No no-data where a value is used */
  VCHECK(Cover, VOVER(i), "Veg.OverstoryCover", Cover[i] == MAP_NODATA_FLOAT, "no-data",
         "Veg.OverstoryCover");
  VCHECK(LAI, VOVER(i), "Veg.OverstoryLAI", LAI[i] == MAP_NODATA_FLOAT, "no-data",
         "Veg.OverstoryLAI");
  VCHECK(OverstoryHeight, VOVER(i), "Veg.OverstoryHeight",
         OverstoryHeight[i] == MAP_NODATA_FLOAT, "no-data", "Veg.OverstoryHeight");
  VCHECK(UnderstoryHeight, VUNDER(i), "Veg.UnderstoryHeight",
         UnderstoryHeight[i] == MAP_NODATA_FLOAT, "no-data", "Veg.UnderstoryHeight");

  /* Overstory cover in (0, 1], overstory LAI and heights positive */
  VCHECK(Cover, VOVER(i), "Veg.OverstoryCover", !(Cover[i] > 0.0 && Cover[i] <= 1.0),
         "cover not in (0, 1] on overstory cells", "cover");
  VCHECK(LAI, VOVER(i), "Veg.OverstoryLAI", !(LAI[i] > 0.0), "LAI <= 0 on overstory cells",
         "LAI");
  VCHECK(OverstoryHeight, VOVER(i), "Veg.OverstoryHeight", !(OverstoryHeight[i] > 0.0),
         "height <= 0 on overstory cells", "height");
  VCHECK(UnderstoryHeight, VUNDER(i), "Veg.UnderstoryHeight", !(UnderstoryHeight[i] > 0.0),
         "height <= 0 on understory cells", "height");

  /* The soil layers of a cell follow its class's root depths, which an
     update keeps */
  if (BaseClass != NULL) {
    memset(&F, 0, sizeof(F));
    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        if (VMASK(i) && VegClass[i] != BaseClass[i]) {
          for (k = 0; k < Map->SoilLayers; k++)
            if (VType[VegClass[i] - 1].RootDepth[k] != VType[BaseClass[i] - 1].RootDepth[k])
              break;
          if (k < Map->SoilLayers)
            Fail(&F, y, x, VegClass[i], BaseClass[i]);
        }
      }
    }
    snprintf(What, sizeof(What), "Veg.Class%s", Suffix);
    ReportCells(&F, What, "root depths differ from those of the cell's class in Veg.Class",
                "class", "Veg.Class");
  }

#undef VMASK
#undef VCLASSOK
#undef VOVER
#undef VUNDER
#undef VCHECK
}

/*****************************************************************************
  ValidateInputs()

  The input maps as given (before the parameters), the class tables' root
  depths and the station normals
*****************************************************************************/
void ValidateInputs(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map, INPUTMAPS *Raw,
  LAYER *Soil, SOILTABLE *SType, LAYER *Veg, VEGTABLE *VType, int NStats,
  METLOCATION *Stat)
{
  int NCells = Map->NY * Map->NX;
  int L = Map->SoilLayers;
  int x, y, i, k, Ones;

  float TotalRootDepth, Fc, Wp, Por;
  float *Layer;
  char Name[BUFSIZE + 1];
  char Path[BUFSIZE + 1];
  FAILURES F;

  Message("\nChecking the inputs\n");

#define INMASK(i) (Raw->Mask[i] == 1)
#define NODATA(v) ((v) == MAP_NODATA_FLOAT)
#define SOILOK(i) (Raw->SoilClass[i] >= 1 && Raw->SoilClass[i] <= Soil->NTypes)

  /* The mask is 0 or 1 everywhere, with at least one 1 */
  memset(&F, 0, sizeof(F));
  for (y = 0, i = 0, Ones = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++, i++) {
      if (Raw->Mask[i] != 0 && Raw->Mask[i] != 1)
        Fail(&F, y, x, Raw->Mask[i], 0);
      if (Raw->Mask[i] == 1)
        Ones++;
    }
  }
  ReportCells(&F, "Terrain.Mask", "value other than 0 or 1", "Terrain.Mask", NULL);
  if (Ones == 0)
    ReportErrorCode(108, "Terrain.Mask: no cell has the value 1");

  /* The elevation has a value everywhere */
  memset(&F, 0, sizeof(F));
  for (y = 0, i = 0; y < Map->NY; y++)
    for (x = 0; x < Map->NX; x++, i++)
      if (NODATA(Raw->Elevation[i]))
        Fail(&F, y, x, Raw->Elevation[i], 0);
  ReportCells(&F, "Terrain.Elevation", "no-data", "Terrain.Elevation", NULL);

  /* No no-data where a value is used: every mask cell, the overstory
     maps on overstory cells */
#define CHECK_NODATA(Array, Where, MapName)                         \
  memset(&F, 0, sizeof(F));                                         \
  for (y = 0, i = 0; y < Map->NY; y++)                              \
    for (x = 0; x < Map->NX; x++, i++)                              \
      if ((Where) && NODATA((Array)[i]))                            \
        Fail(&F, y, x, (Array)[i], 0);                              \
  ReportCells(&F, MapName, "no-data", MapName, NULL);
  if (Options->HorizonAngle.Ptr != NULL) {
    if (!(Layer = calloc(NCells, sizeof(float))))
      ReportError("ValidateInputs", 1);
    for (k = 0; k < MAP_NAZIMUTHS; k++) {
      ReadMapLayer(&(Options->HorizonAngle), k, Layer, Map);
      snprintf(Name, sizeof(Name), "Terrain.HorizonAngle.Az%03d", k * 360 / MAP_NAZIMUTHS);
      CHECK_NODATA(Layer, INMASK(i), Name);
    }
    free(Layer);
  }
  if (Raw->SkyView != NULL) {
    CHECK_NODATA(Raw->SkyView, INMASK(i), "Terrain.SkyView");
  }
  memset(&F, 0, sizeof(F));
  for (y = 0, i = 0; y < Map->NY; y++)
    for (x = 0; x < Map->NX; x++, i++)
      if (INMASK(i) && Raw->SoilClass[i] == MAP_NODATA_BYTE)
        Fail(&F, y, x, Raw->SoilClass[i], 0);
  ReportCells(&F, "Soil.Class", "no-data", "Soil.Class", NULL);
  CHECK_NODATA(Raw->Depth, INMASK(i), "Soil.Depth");
  CHECK_NODATA(Raw->LateralConductivity, INMASK(i), "Soil.LateralConductivity");
  CHECK_NODATA(Raw->ExponentialDecrease, INMASK(i), "Soil.ExponentialDecrease");
  for (k = 0; k <= L; k++) {
    Layer = Raw->Porosity + k * NCells;
    SoilLayerName(Name, sizeof(Name), "Soil.Porosity", k, L);
    CHECK_NODATA(Layer, INMASK(i), Name);
  }
  for (k = 0; k <= L; k++) {
    Layer = Raw->FieldCapacity + k * NCells;
    SoilLayerName(Name, sizeof(Name), "Soil.FieldCapacity", k, L);
    CHECK_NODATA(Layer, INMASK(i), Name);
  }
  if (Options->Prism) {
    if (!(Layer = calloc(NCells, sizeof(float))))
      ReportError("ValidateInputs", 1);
    for (k = 0; k < 12; k++) {
      ReadMapLayer(&(Options->PrecipPattern), k, Layer, Map);
      snprintf(Name, sizeof(Name), "Met.PrecipPattern.%s", MonthNames[k]);
      CHECK_NODATA(Layer, INMASK(i), Name);
    }
    free(Layer);
  }
  if (Options->SnowPattern) {
    CHECK_NODATA(Raw->SnowPattern, INMASK(i), "Met.SnowPattern");
  }
  if (Options->MeltScale) {
    CHECK_NODATA(Raw->MeltScale, INMASK(i), "Snow.MeltScale");
  }
  if (Options->WindScale) {
    CHECK_NODATA(Raw->WindScale, INMASK(i), "Met.WindScale");
  }
  if (Options->LakeDynamics) {
    memset(&F, 0, sizeof(F));
    for (y = 0, i = 0; y < Map->NY; y++)
      for (x = 0; x < Map->NX; x++, i++)
        if (INMASK(i) && Raw->LakeID[i] == MAP_NODATA_BYTE)
          Fail(&F, y, x, Raw->LakeID[i], 0);
    ReportCells(&F, "Lake.ID", "no-data", "Lake.ID", NULL);
  }
  if (Raw->InitialSWE != NULL) {
    CHECK_NODATA(Raw->InitialSWE, INMASK(i), "Snow.InitialSWE");
    free(Raw->InitialSWE);
    Raw->InitialSWE = NULL;
  }

  /* Classes in the tables */
  memset(&F, 0, sizeof(F));
  for (y = 0, i = 0; y < Map->NY; y++)
    for (x = 0; x < Map->NX; x++, i++)
      if (INMASK(i) && !SOILOK(i))
        Fail(&F, y, x, Raw->SoilClass[i], 0);
  snprintf(Name, sizeof(Name), "class not in 1..%d (the soil table rows)", Soil->NTypes);
  ReportCells(&F, "Soil.Class", Name, "Soil.Class", NULL);

  /* The vegetation maps */
  ValidateVegMaps(Map, Raw->Mask, Raw->VegClass, Raw->OverstoryCover, Raw->OverstoryLAI,
                  Raw->OverstoryHeight, Raw->UnderstoryHeight, NULL, Veg, VType, "");

  /* The depth map at least 0.01 m below the root zone */
  memset(&F, 0, sizeof(F));
  for (y = 0, i = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++, i++) {
      if (INMASK(i)) {
        TotalRootDepth = 0.0;
        for (k = 0; k < L; k++)
          TotalRootDepth += VType[Raw->VegClass[i] - 1].RootDepth[k];
        if (Raw->Depth[i] < (TotalRootDepth + 0.01))
          Fail(&F, y, x, Raw->Depth[i], TotalRootDepth);
      }
    }
  }
  ReportCells(&F, "Soil.Depth", "depth < root zone depth + 0.01", "depth", "root zone depth");

  /* Conductivity and exponential decrease, and their scales */
  memset(&F, 0, sizeof(F));
  for (y = 0, i = 0; y < Map->NY; y++)
    for (x = 0; x < Map->NX; x++, i++)
      if (INMASK(i) && !(Raw->LateralConductivity[i] > 0.0))
        Fail(&F, y, x, Raw->LateralConductivity[i], 0);
  ReportCells(&F, "Soil.LateralConductivity", "value <= 0", "Soil.LateralConductivity", NULL);
  if (!(SOIL_KSAT_ADJ > 0.0))
    ReportErrorCode(108, "[PARAMETERS] Lateral Conductivity Scale = %g: must be > 0",
                    SOIL_KSAT_ADJ);
  memset(&F, 0, sizeof(F));
  for (y = 0, i = 0; y < Map->NY; y++)
    for (x = 0; x < Map->NX; x++, i++)
      if (INMASK(i) && !(Raw->ExponentialDecrease[i] >= 0.0))
        Fail(&F, y, x, Raw->ExponentialDecrease[i], 0);
  ReportCells(&F, "Soil.ExponentialDecrease", "value < 0", "Soil.ExponentialDecrease", NULL);
  if (!(SOIL_EXPDEC_ADJ >= 0.0))
    ReportErrorCode(108, "[PARAMETERS] Exponential Decrease Scale = %g: must be >= 0",
                    SOIL_EXPDEC_ADJ);

  /* Field capacity map at least 0.01 above the wilting point */
  for (k = 0; k <= L; k++) {
    memset(&F, 0, sizeof(F));
    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        if (INMASK(i)) {
          Fc = Raw->FieldCapacity[k * NCells + i];
          Wp = SType[Raw->SoilClass[i] - 1].WP[k];
          if (Fc < (Wp + 0.01))
            Fail(&F, y, x, Fc, Wp);
        }
      }
    }
    SoilLayerName(Name, sizeof(Name), "Soil.FieldCapacity", k, L);
    ReportCells(&F, Name, "field capacity < wilting point + 0.01", "field capacity",
                "wilting point");
  }

  /* Porosity map at least 0.01 above the field capacity map */
  for (k = 0; k <= L; k++) {
    memset(&F, 0, sizeof(F));
    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        if (INMASK(i)) {
          Por = Raw->Porosity[k * NCells + i];
          Fc = Raw->FieldCapacity[k * NCells + i];
          if (Por < (Fc + 0.01))
            Fail(&F, y, x, Por, Fc);
        }
      }
    }
    SoilLayerName(Name, sizeof(Name), "Soil.Porosity", k, L);
    ReportCells(&F, Name, "porosity < field capacity + 0.01", "porosity", "field capacity");
  }

  /* Porosity map at most 1 */
  for (k = 0; k <= L; k++) {
    memset(&F, 0, sizeof(F));
    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        if (INMASK(i)) {
          Por = Raw->Porosity[k * NCells + i];
          if (Por > 1.0)
            Fail(&F, y, x, Por, 0);
        }
      }
    }
    SoilLayerName(Name, sizeof(Name), "Soil.Porosity", k, L);
    ReportCells(&F, Name, "porosity > 1", "porosity", NULL);
  }

  /* The scale of the overstory cover positive, the monthly LAI fractions of
     every class with an overstory positive */
  if (!(VEG_FC_ADJ > 0.0))
    ReportErrorCode(108, "[PARAMETERS] Overstory Cover Scale = %g: must be > 0", VEG_FC_ADJ);
  RequireConfigValue(Input, "INPUTS", "Vegetation Table", Path);
  for (i = 0; i < Veg->NTypes; i++) {
    for (k = 0; VType[i].OverStory == TRUE && k < 12; k++) {
      if (!(VType[i].LAIFraction[k] > 0.0))
        ReportErrorCode(108, "%s: line %d: OverstoryLAIFraction_%s = %g must be > 0", Path, i + 3,
                        MonthNames[k], VType[i].LAIFraction[k]);
    }
  }

  /* Lake IDs in the lake table (Lake Dynamics) */
  if (Options->LakeDynamics) {
    memset(&F, 0, sizeof(F));
    for (y = 0, i = 0; y < Map->NY; y++)
      for (x = 0; x < Map->NX; x++, i++)
        if (INMASK(i) && Raw->LakeID[i] > Map->NumLakes)
          Fail(&F, y, x, Raw->LakeID[i], 0);
    snprintf(Name, sizeof(Name), "lake not in 0..%d (the lake table rows)", Map->NumLakes);
    ReportCells(&F, "Lake.ID", Name, "Lake.ID", NULL);
  }

  /* Station normals positive */
  RequireConfigValue(Input, "INPUTS", "Station Table", Path);
  for (i = 0; i < NStats; i++) {
    for (k = 0; Options->Prism && k < 12; k++) {
      if (!(Stat[i].PrismPrecip[k] > 0.0))
        ReportErrorCode(108, "%s: line %d: PrecipNormal_%s = %g must be > 0", Path, i + 3,
                        MonthNames[k], Stat[i].PrismPrecip[k]);
    }
    if (Options->SnowPattern && !(Stat[i].SnowPatternBase > 0.0))
      ReportErrorCode(108, "%s: line %d: SnowNormal = %g must be > 0", Path, i + 3,
                      Stat[i].SnowPatternBase);
  }

  /* Root zone layers thicker than 0.001 m */
  RequireConfigValue(Input, "INPUTS", "Vegetation Table", Path);
  for (i = 0; i < Veg->NTypes; i++) {
    for (k = 0; k < L; k++) {
      if (!(VType[i].RootDepth[k] > 0.001))
        ReportErrorCode(108, "%s: line %d: RootDepth_%d = %g must be > 0.001", Path, i + 3,
                        k + 1, VType[i].RootDepth[k]);
    }
  }

#undef CHECK_NODATA
#undef SOILOK
#undef NODATA
#undef INMASK
}

/*****************************************************************************
  ValidateParameters()

  The porosity after the field capacity and porosity parameters at
  most 1
*****************************************************************************/
void ValidateParameters(MAPSIZE *Map, TOPOPIX **TopoMap, SOILPIX **SoilMap)
{
  int x, y, k;
  int L = Map->SoilLayers;
  char Name[BUFSIZE + 1];
  FAILURES F;

  for (k = 0; k <= L; k++) {
    memset(&F, 0, sizeof(F));
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
        if (INBASIN(TopoMap[y][x].Mask) && (double) SoilMap[y][x].Porosity[k] > 1.0)
          Fail(&F, y, x, SoilMap[y][x].Porosity[k], SoilMap[y][x].FCap[k]);
      }
    }
    SoilLayerName(Name, sizeof(Name), "Soil.Porosity", k, L);
    strcat(Name, " x Porosity Scale");
    ReportCells(&F, Name, "porosity > 1", "porosity", "field capacity");
  }
}

/*****************************************************************************
  ValidateWeights()

  The station weights of every mask cell (bytes) sum to more than 0
*****************************************************************************/
void ValidateWeights(MAPSIZE *Map, TOPOPIX **TopoMap, uchar ***MetWeights, int NStats)
{
  int x, y, i;
  long Sum;
  FAILURES F;

  memset(&F, 0, sizeof(F));
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
        for (i = 0, Sum = 0; i < NStats; i++)
          Sum += MetWeights[y][x][i];
        if (Sum <= 0)
          Fail(&F, y, x, Sum, 0);
      }
    }
  }
  ReportCells(&F, "station weights", "weights sum to 0", "sum", NULL);
}

/*****************************************************************************
  PrintClassSummary()

  The vegetation and soil classes in the basin and their fractions
*****************************************************************************/
void PrintClassSummary(OPTIONSTRUCT *Options, LAYER Veg, LAYER Soil, VEGTABLE *VType,
  SOILTABLE *SType, MAPSIZE *Map, TOPOPIX **TopoMap, VEGPIX **VegMap, SOILPIX **SoilMap)
{
  int y, x, i, k;
  int *count = NULL, *scount = NULL;
  double *usnow = NULL;
  int npixels;

  if (!(count = calloc(Veg.NTypes, sizeof(int))) || !(scount = calloc(Soil.NTypes, sizeof(int))) ||
      !(usnow = calloc(Veg.NTypes, sizeof(double))))
    ReportError("PrintClassSummary", 1);

  npixels = 0;
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
        count[VegMap[y][x].Veg - 1]++;
        usnow[VegMap[y][x].Veg - 1] += VegMap[y][x].USnow;
        scount[SoilMap[y][x].Soil - 1]++;
        npixels++;
      }
    }
  }
  Message("\nBasin has %d active pixels \n", npixels);

  Message("\nThe following VEG types are in the current basin \n");
  for (i = 0; i < Veg.NTypes; i++) {
    if (count[i] > 0)
      Message("Class # %d of Type: %s has fraction basin area: %5.3f\n",
        i + 1, VType[i].Desc, (float)count[i] / (float)npixels);
  }

  Message("\nThe following SOIL types are in the current basin \n");
  for (i = 0; i < Soil.NTypes; i++)
    if (scount[i] > 0)
      Message("Class # %d of Type: %s has fraction basin area: %5.3f\n",
        i + 1, SType[i].Desc, (float)scount[i] / (float)npixels);

  Message("\nSome estimates for current vegetation specification\n");
  for (i = 0; i < Veg.NTypes; i++) {
    if (count[i] > 0) {
      Message("\nVegetation Type: %s\n", VType[i].Desc);
      Message("2meter    wind speed fraction of ref level %1.3f (mean of the cells)\n",
              usnow[i] / count[i]);
    }
  }

  if (Options->ImprovRadiation) {
    Message("\nMonthly attenuation coeff. for the VEG Type with an overstory: \n\n");
    for (i = 0; i < Veg.NTypes; i++) {
      if (VType[i].OverStory == TRUE) {
        Message("\n%s:\n", VType[i].Desc);
        for (k = 0; k < 12; k++) {
          Message("%.3f ", VType[i].MonthlyExtnCoeff[k]);
        }
        Message("\n");
      }
    }
    Message("\n");
  }

  free(count);
  free(scount);
  free(usnow);
}

/*****************************************************************************
  ValidateChannels()

  The channel class, segment and cell tables and the lakes, after the
  parameters (the channel cuts are compared with the soil depth)
*****************************************************************************/
void ValidateChannels(OPTIONSTRUCT *Options, MAPSIZE *Map, TOPOPIX **TopoMap,
                      SOILPIX **SoilMap, CHANNELTABLES *Tables, LAKETABLE *LType)
{
  int i, j, k, MaxOrder, Found, InMask;
  int x, y;
  const char *c;

  /* Segment and channel class IDs positive and unique */
  for (i = 0; i < Tables->NSegments; i++) {
    if (Tables->SegmentID[i] < 1 || Tables->SegmentID[i] > 65535)
      ReportErrorCode(109, "%s: line %d: Segment %d is not in 1..65535", Tables->SegmentPath,
                      i + 3, Tables->SegmentID[i]);
    for (j = 0; j < i; j++) {
      if (Tables->SegmentID[j] == Tables->SegmentID[i])
        ReportErrorCode(109, "%s: line %d: Segment %d repeated", Tables->SegmentPath, i + 3,
                        Tables->SegmentID[i]);
    }
  }
  for (i = 0; i < Tables->NClasses; i++) {
    if (Tables->ClassID[i] < 1 || Tables->ClassID[i] > 65535)
      ReportErrorCode(109, "%s: line %d: Class %d is not in 1..65535", Tables->ClassPath,
                      i + 3, Tables->ClassID[i]);
    for (j = 0; j < i; j++) {
      if (Tables->ClassID[j] == Tables->ClassID[i])
        ReportErrorCode(109, "%s: line %d: Class %d repeated", Tables->ClassPath, i + 3,
                        Tables->ClassID[i]);
    }
  }

  /* The class of every segment is in the class table */
  for (i = 0; i < Tables->NSegments; i++) {
    for (j = 0, Found = FALSE; j < Tables->NClasses && !Found; j++)
      Found = (Tables->ClassID[j] == Tables->Class[i]);
    if (!Found)
      ReportErrorCode(109, "%s: line %d: Class %d not in %s", Tables->SegmentPath, i + 3,
                      Tables->Class[i], Tables->ClassPath);
  }

  /* Every outlet 0 or a segment */
  for (i = 0; i < Tables->NSegments; i++) {
    for (j = 0, Found = (Tables->Outlet[i] == 0); j < Tables->NSegments && !Found; j++)
      Found = (Tables->SegmentID[j] == Tables->Outlet[i]);
    if (!Found)
      ReportErrorCode(109, "%s: line %d: Outlet %d is not a segment", Tables->SegmentPath,
                      i + 3, Tables->Outlet[i]);
  }

  /* Orders 1..max all present, every segment drains to a higher order */
  for (i = 0, MaxOrder = 0; i < Tables->NSegments; i++) {
    if (Tables->Order[i] < 1)
      ReportErrorCode(109, "%s: line %d: Order %d must be >= 1", Tables->SegmentPath, i + 3,
                      Tables->Order[i]);
    if (Tables->Order[i] > MaxOrder)
      MaxOrder = Tables->Order[i];
  }
  for (k = 1; k <= MaxOrder; k++) {
    for (i = 0, Found = FALSE; i < Tables->NSegments && !Found; i++)
      Found = (Tables->Order[i] == k);
    if (!Found)
      ReportErrorCode(109, "%s: no segment of order %d", Tables->SegmentPath, k);
  }
  for (i = 0; i < Tables->NSegments; i++) {
    for (j = 0; j < Tables->NSegments && Tables->Outlet[i] != 0; j++) {
      if (Tables->SegmentID[j] == Tables->Outlet[i] && Tables->Order[j] <= Tables->Order[i])
        ReportErrorCode(109, "%s: line %d: segment %d (order %d) drains to segment %d (order %d)",
                        Tables->SegmentPath, i + 3, Tables->SegmentID[i], Tables->Order[i],
                        Tables->SegmentID[j], Tables->Order[j]);
    }
  }

  /* Positive slopes, lengths, widths, bank heights, roughness */
  for (i = 0; i < Tables->NSegments; i++) {
    if (!(Tables->Slope[i] > 0.0))
      ReportErrorCode(109, "%s: line %d: Slope = %g must be > 0", Tables->SegmentPath, i + 3,
                      Tables->Slope[i]);
    if (!(Tables->Length[i] > 0.0))
      ReportErrorCode(109, "%s: line %d: Length = %g must be > 0", Tables->SegmentPath, i + 3,
                      Tables->Length[i]);
  }
  for (i = 0; i < Tables->NClasses; i++) {
    if (!(Tables->Width[i] > 0.0))
      ReportErrorCode(109, "%s: line %d: Width = %g must be > 0", Tables->ClassPath, i + 3,
                      Tables->Width[i]);
    if (!(Tables->BankHeight[i] > 0.0))
      ReportErrorCode(109, "%s: line %d: BankHeight = %g must be > 0", Tables->ClassPath, i + 3,
                      Tables->BankHeight[i]);
    if (!(Tables->Friction[i] > 0.0))
      ReportErrorCode(109, "%s: line %d: ManningsN = %g must be > 0", Tables->ClassPath, i + 3,
                      Tables->Friction[i]);
  }

  /* Channel cells inside the grid and the mask */
  for (i = 0; i < Tables->NCells; i++) {
    y = Tables->Row[i] - 1;
    x = Tables->Column[i] - 1;
    InMask = (y >= 0 && y < Map->NY && x >= 0 && x < Map->NX && INBASIN(TopoMap[y][x].Mask));
    if (!InMask)
      ReportErrorCode(109, "%s: line %d: row %d, column %d is outside the mask", Tables->CellPath,
                      i + 3, Tables->Row[i], Tables->Column[i]);
  }

  /* The segments of the channel cells exist, lengths positive (the lateral
     inflow of a cell is split by length), cuts not negative */
  for (i = 0; i < Tables->NCells; i++) {
    for (j = 0, Found = FALSE; j < Tables->NSegments && !Found; j++)
      Found = (Tables->SegmentID[j] == Tables->Segment[i]);
    if (!Found)
      ReportErrorCode(109, "%s: line %d: Segment %d is not in %s", Tables->CellPath, i + 3,
                      Tables->Segment[i], Tables->SegmentPath);
    if (!(Tables->CellLength[i] > 0.0))
      ReportErrorCode(109, "%s: line %d: Length = %g must be > 0", Tables->CellPath, i + 3,
                      Tables->CellLength[i]);
    if (!(Tables->CutWidth[i] >= 0.0))
      ReportErrorCode(109, "%s: line %d: CutWidth = %g must be >= 0", Tables->CellPath, i + 3,
                      Tables->CutWidth[i]);
    if (!(Tables->CutHeight[i] >= 0.0))
      ReportErrorCode(109, "%s: line %d: CutHeight = %g must be >= 0", Tables->CellPath, i + 3,
                      Tables->CutHeight[i]);
  }

  /* Channel cuts not deeper than the soil (after the depth offset) */
  for (i = 0; i < Tables->NCells; i++) {
    y = Tables->Row[i] - 1;
    x = Tables->Column[i] - 1;
    if (Tables->CutHeight[i] > SoilMap[y][x].Depth)
      ReportErrorCode(109, "%s: line %d: CutHeight %g is deeper than the soil (%g m) at row %d, column %d",
                      Tables->CellPath, i + 3, Tables->CutHeight[i], SoilMap[y][x].Depth,
                      Tables->Row[i], Tables->Column[i]);
  }

  /* Save names of letters, digits, _ . - and unique */
  for (i = 0; i < Tables->NSegments; i++) {
    if (Tables->SaveName[i][0] == '\0')
      continue;
    for (c = Tables->SaveName[i]; *c != '\0'; c++) {
      if (!((*c >= 'A' && *c <= 'Z') || (*c >= 'a' && *c <= 'z') || (*c >= '0' && *c <= '9') ||
            *c == '_' || *c == '.' || *c == '-'))
        ReportErrorCode(109, "%s: line %d: SaveName \"%s\" has a character outside A-Z a-z 0-9 _ . -",
                        Tables->SegmentPath, i + 3, Tables->SaveName[i]);
    }
    for (j = 0; j < i; j++) {
      if (strcmp(Tables->SaveName[j], Tables->SaveName[i]) == 0)
        ReportErrorCode(109, "%s: line %d: SaveName \"%s\" repeated", Tables->SegmentPath, i + 3,
                        Tables->SaveName[i]);
    }
  }

  /* Lake outlets are segments, every lake has a cell in the mask */
  if (Options->LakeDynamics) {
    for (i = 0; i < Map->NumLakes; i++) {
      for (j = 0, Found = FALSE; j < Tables->NSegments && !Found; j++)
        Found = (Tables->SegmentID[j] == LType[i].Outlet);
      if (!Found)
        ReportErrorCode(109, "%s: line %d: OutletSegment %d is not a segment",
                        Options->LakeTablePath, i + 3, LType[i].Outlet);
    }
    for (i = 0; i < Map->NumLakes; i++) {
      for (y = 0, Found = FALSE; y < Map->NY && !Found; y++)
        for (x = 0; x < Map->NX && !Found; x++)
          Found = (INBASIN(TopoMap[y][x].Mask) && TopoMap[y][x].LakeID == i + 1);
      if (!Found)
        ReportErrorCode(109, "Lake.ID: lake %d has no cell in the mask", i + 1);
    }
  }
}
