/*
 * DESCRIPTION:  Vegetation updates: at each [INPUTS] Vegetation Update <n>
 *               instant the five vegetation maps of every mask cell are
 *               replaced by <Vegetation Update Directory>/<Name>.<stamp>.map
 *               (Veg.Class, Veg.OverstoryCover, Veg.OverstoryLAI,
 *               Veg.OverstoryHeight, Veg.UnderstoryHeight), before the time
 *               step starting then. The latest update before Start gives
 *               the vegetation at Start. A cell's class keeps the root
 *               depths of its class in Veg.Class, so its soil layers stay.
 *               The water of a canopy the new class lacks, and the
 *               detention storage of a class without impervious fraction,
 *               becomes surface water of the cell.
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
#include "mapio.h"
#include "DHSVMChannel.h"

/* The maps of one update */
typedef struct {
  unsigned char *Class;
  float *Cover;
  float *LAI;
  float *OverstoryHeight;
  float *UnderstoryHeight;
} VEGMAPS;

/*****************************************************************************
  InitVegUpdateTimes()

  The update instants [INPUTS] Vegetation Update 1, 2, ... (the list ends at
  the first missing number), on the step grid and increasing, with
  Vegetation Update Directory; any instant before, during or after the run
*****************************************************************************/
void InitVegUpdateTimes(LISTPTR Input, OPTIONSTRUCT *Options, TIMESTRUCT *Time)
{
  VEGUPDATES *U = &(Options->VegUpdates);
  char KeyName[BUFSIZE + 1];
  char Value[BUFSIZE + 1];
  long long Seconds;
  int n, k;

  for (n = 0; ; n++) {
    snprintf(KeyName, sizeof(KeyName), "Vegetation Update %d", n + 1);
    if (!GetConfigValue(Input, "INPUTS", KeyName, Value))
      break;
  }
  U->N = n;
  U->Next = 0;
  U->Time = NULL;
  U->Step = NULL;
  U->BaseClass = NULL;
  U->Present = NULL;
  if (n > 0 && IsEmptyStr(U->Directory))
    ReportErrorCode(103, "[INPUTS] Vegetation Update 1 needs Vegetation Update Directory");
  if (n == 0 && !IsEmptyStr(U->Directory))
    ReportErrorCode(103, "[INPUTS] Vegetation Update Directory needs Vegetation Update 1");
  if (n == 0)
    return;

  if (!(U->Time = (DATE *) calloc(n, sizeof(DATE))) ||
      !(U->Step = (long long *) calloc(n, sizeof(long long))))
    ReportError("InitVegUpdateTimes", 1);

  for (k = 0; k < n; k++) {
    snprintf(KeyName, sizeof(KeyName), "Vegetation Update %d", k + 1);
    RequireConfigValue(Input, "INPUTS", KeyName, Value);
    if (!SScanIsoDate(Value, &(U->Time[k])))
      ReportErrorCode(103, "[INPUTS] %s = %s: expected YYYY-MM-DD HH:MM", KeyName, Value);
    Seconds = (MinutesOf(&(U->Time[k])) - MinutesOf(&(Time->Start))) * SECPMIN;
    if (Seconds % Time->Dt != 0)
      ReportErrorCode(110, "[INPUTS] %s = %s: not on the time step grid", KeyName, Value);
    U->Step[k] = Seconds / Time->Dt;
    if (k > 0 && U->Step[k] <= U->Step[k - 1])
      ReportErrorCode(110, "[INPUTS] %s = %s: not after Vegetation Update %d", KeyName, Value, k);
  }
}

/* Reads one map of update k into Values */
static void ReadUpdateMap(VEGUPDATES *U, int k, const char *Name, int NumberType,
                          const char *Units, MAPSIZE *Map, void *Values)
{
  char Path[BUFSIZE * 2 + 1];
  char Stamp[NAMESIZE + 1];
  MAPFILE File;

  SPrintIsoStamp(&(U->Time[k]), Stamp);
  if (snprintf(Path, sizeof(Path), "%s%s.%s.map", U->Directory, Name, Stamp) >= (int) sizeof(Path))
    ReportErrorCode(104, "%s%s.%s.map: path too long", U->Directory, Name, Stamp);
  OpenInputMap(&File, Path, Name, NumberType, Units, MAP_ONE, Map);
  ReadMapLayer(&File, 0, Values, Map);
  CloseMap(&File);
}

/*****************************************************************************
  ReadVegUpdate()

  Reads and checks the five maps of update k (the checks of the input
  vegetation maps, and the root depths of the classes in Veg.Class) and
  marks the classes they give to mask cells
*****************************************************************************/
static void ReadVegUpdate(VEGUPDATES *U, int k, MAPSIZE *Map, unsigned char *Mask, LAYER *Veg,
                          VEGTABLE *VType, VEGMAPS *V)
{
  size_t NCells = (size_t) Map->NY * Map->NX;
  char Suffix[NAMESIZE + 2];
  size_t i;

  if (!(V->Class = (unsigned char *) calloc(NCells, sizeof(unsigned char))) ||
      !(V->Cover = (float *) calloc(NCells, sizeof(float))) ||
      !(V->LAI = (float *) calloc(NCells, sizeof(float))) ||
      !(V->OverstoryHeight = (float *) calloc(NCells, sizeof(float))) ||
      !(V->UnderstoryHeight = (float *) calloc(NCells, sizeof(float))))
    ReportError("ReadVegUpdate", 1);
  ReadUpdateMap(U, k, "Veg.Class", NC_BYTE, "-", Map, V->Class);
  ReadUpdateMap(U, k, "Veg.OverstoryCover", NC_FLOAT, "fraction", Map, V->Cover);
  ReadUpdateMap(U, k, "Veg.OverstoryLAI", NC_FLOAT, "m2/m2", Map, V->LAI);
  ReadUpdateMap(U, k, "Veg.OverstoryHeight", NC_FLOAT, "m", Map, V->OverstoryHeight);
  ReadUpdateMap(U, k, "Veg.UnderstoryHeight", NC_FLOAT, "m", Map, V->UnderstoryHeight);

  Suffix[0] = '.';
  SPrintIsoStamp(&(U->Time[k]), Suffix + 1);
  ValidateVegMaps(Map, Mask, V->Class, V->Cover, V->LAI, V->OverstoryHeight,
                  V->UnderstoryHeight, U->BaseClass, Veg, VType, Suffix);

  for (i = 0; i < NCells; i++)
    if (Mask[i] == 1)
      U->Present[V->Class[i] - 1] = TRUE;
}

static void FreeVegUpdate(VEGMAPS *V)
{
  free(V->Class);
  free(V->Cover);
  free(V->LAI);
  free(V->OverstoryHeight);
  free(V->UnderstoryHeight);
}

/*****************************************************************************
  InitVegUpdates()

  Reads and checks every update before End (so that a bad map stops the run
  before it starts); the latest one before Start replaces the vegetation
  maps as read, which then give the vegetation at Start
*****************************************************************************/
void InitVegUpdates(OPTIONSTRUCT *Options, MAPSIZE *Map, TIMESTRUCT *Time, INPUTMAPS *Raw,
                    LAYER *Veg, VEGTABLE *VType)
{
  VEGUPDATES *U = &(Options->VegUpdates);
  size_t NCells = (size_t) Map->NY * Map->NX;
  char Stamp[NAMESIZE + 1];
  VEGMAPS V;
  int Latest = -1;
  int k;

  if (U->N == 0)
    return;

  Message("\nChecking the vegetation updates\n");

  if (!(U->BaseClass = (unsigned char *) malloc(NCells * sizeof(unsigned char))) ||
      !(U->Present = (unsigned char *) calloc(Veg->NTypes, sizeof(unsigned char))))
    ReportError("InitVegUpdates", 1);
  memcpy(U->BaseClass, Raw->VegClass, NCells * sizeof(unsigned char));

  /* Classes in Veg.Class are checked first: the root depth check reads them */
  ValidateVegMaps(Map, Raw->Mask, Raw->VegClass, Raw->OverstoryCover, Raw->OverstoryLAI,
                  Raw->OverstoryHeight, Raw->UnderstoryHeight, NULL, Veg, VType, "");

  U->Next = U->N;
  for (k = 0; k < U->N && U->Step[k] < Time->NSteps; k++) {
    if (U->Step[k] < 0)
      Latest = k;
    else if (U->Next == U->N)
      U->Next = k;
    ReadVegUpdate(U, k, Map, Raw->Mask, Veg, VType, &V);
    FreeVegUpdate(&V);
  }

  if (Latest >= 0) {
    ReadVegUpdate(U, Latest, Map, Raw->Mask, Veg, VType, &V);
    free(Raw->VegClass);
    free(Raw->OverstoryCover);
    free(Raw->OverstoryLAI);
    free(Raw->OverstoryHeight);
    free(Raw->UnderstoryHeight);
    Raw->VegClass = V.Class;
    Raw->OverstoryCover = V.Cover;
    Raw->OverstoryLAI = V.LAI;
    Raw->OverstoryHeight = V.OverstoryHeight;
    Raw->UnderstoryHeight = V.UnderstoryHeight;
    SPrintIsoStamp(&(U->Time[Latest]), Stamp);
    Message("Vegetation at Start from the update of %s\n", Stamp);
  }
}

/*****************************************************************************
  ApplyVegUpdate()

  At the start of a time step with an update: the vegetation of every mask
  cell from the update's maps (rows 6 and 7, the channel limit on the
  understory cover, the albedo decay rates of the snow, the LAI of the
  month and the ground albedo under the snow); intercepted water moves with
  its canopy, a canopy the new class lacks leaves its water on the ground
*****************************************************************************/
void ApplyVegUpdate(OPTIONSTRUCT *Options, TIMESTRUCT *Time, MAPSIZE *Map, TOPOPIX **TopoMap,
                    LAYER *Veg, VEGTABLE *VType, VEGPIX **VegMap, SOILTABLE *SType,
                    SOILPIX **SoilMap, SNOWPIX **SnowMap, PRECIPPIX **PrecipMap,
                    CHANNEL *ChannelData)
{
  VEGUPDATES *U = &(Options->VegUpdates);
  unsigned char *Mask;
  char Stamp[NAMESIZE + 1];
  VEGMAPS V;
  VEGTABLE *Old, *New;
  PRECIPPIX *Precip;
  VEGPIX *Cell;
  float Rain[2], Snow[2];
  float Surface;
  int x, y, i, c, j;

  if (U->Next >= U->N || U->Step[U->Next] != (long long) Time->Step)
    return;

  SPrintIsoStamp(&(U->Time[U->Next]), Stamp);
  Message("Updating the vegetation (%s)\n", Stamp);

  if (!(Mask = (unsigned char *) calloc((size_t) Map->NY * Map->NX, sizeof(unsigned char))))
    ReportError("ApplyVegUpdate", 1);
  for (y = 0, i = 0; y < Map->NY; y++)
    for (x = 0; x < Map->NX; x++, i++)
      Mask[i] = INBASIN(TopoMap[y][x].Mask) ? 1 : 0;
  ReadVegUpdate(U, U->Next, Map, Mask, Veg, VType, &V);

  for (y = 0, i = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++, i++) {
      if (!Mask[i])
        continue;
      Cell = &(VegMap[y][x]);
      Precip = &(PrecipMap[y][x]);
      Old = &(VType[Cell->Veg - 1]);
      New = &(VType[V.Class[i] - 1]);

      /* Intercepted water by canopy (overstory, understory) */
      for (c = 0; c < 2; c++) {
        j = CanopyIndex(Old, c);
        Rain[c] = (j >= 0) ? Precip->IntRain[j] : 0.0;
        Snow[c] = (j >= 0) ? Precip->IntSnow[j] : 0.0;
      }
      Surface = 0.0;
      if (New->OverStory != TRUE) {
        Surface += Rain[0] + Snow[0] + Precip->TempIntStorage;
        Precip->TempIntStorage = 0.0;
      }
      if (New->UnderStory != TRUE)
        Surface += Rain[1] + Snow[1];
      for (j = 0; j < Veg->MaxLayers; j++) {
        Precip->IntRain[j] = 0.0;
        Precip->IntSnow[j] = 0.0;
      }
      for (c = 0; c < 2; c++) {
        j = CanopyIndex(New, c);
        if (j >= 0) {
          Precip->IntRain[j] = Rain[c];
          Precip->IntSnow[j] = Snow[c];
        }
      }

      /* Detention storage only drains from an impervious fraction */
      if (!(New->ImpervFrac > 0.0)) {
        Surface += SoilMap[y][x].DetentionStorage;
        SoilMap[y][x].DetentionStorage = 0.0;
      }
      SoilMap[y][x].IExcess += Surface;

      /* The cell's vegetation */
      Cell->Veg = V.Class[i];
      SetCellVegetation(Options, New, Cell, V.Cover[i], V.LAI[i], V.OverstoryHeight[i],
                        V.UnderstoryHeight[i], y, x);
      ChannelLimitVegFC(y, x, (Map->DX * Map->DY), ChannelData, New, Cell);
      SetCellMonth(Options, New, Cell, Time->Current.Month);

      /* Albedo decay rates of the snow (the forest offset on overstory
         cells) and the ground albedo under the snow, as at Start */
      SnowMap[y][x].LamdaAcc = ALB_ACC_LAMBDA;
      SnowMap[y][x].LamdaMelt = ALB_MELT_LAMBDA;
      if (New->OverStory == TRUE) {
        SnowMap[y][x].LamdaAcc += (LAMBDA_FOREST_OFFSET * Cell->Fract[0]);
        SnowMap[y][x].LamdaMelt += (LAMBDA_FOREST_OFFSET * Cell->Fract[0]);
      }
      if (New->UnderStory == TRUE)
        SnowMap[y][x].AlbedoGround =
          New->AlbedoMonthly[(New->OverStory == TRUE) ? 1 : 0][Time->Current.Month - 1];
      else
        SnowMap[y][x].AlbedoGround = SType[SoilMap[y][x].Soil - 1].Albedo;
    }
  }

  FreeVegUpdate(&V);
  free(Mask);
  U->Next++;
}
