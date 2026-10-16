/*
 * DESCRIPTION:  The basin water balance from the basin means: storage,
 *               error of each step and of the run, water input, rain or
 *               melt; the channel network sums; the final balance of the
 *               run on stdout.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"
#include "output.h"

static float Storage(AGGREGATE *Basin);
static float WaterInput(AGGREGATE *Aggregate, float *DeltaSWE);
static double NetworkStorage(MODELVIEW *Model, double Area);

/* Water in the basin (m): surface water, canopy storage, soil water, snow
   water equivalent, saturated flow and detention storage, in this order */
static float Storage(AGGREGATE *Basin)
{
  return AggregateValue(Basin, V_SOIL_SURFACEWATER, 0) +
    AggregateValue(Basin, V_CANOPY_STORAGE, 0) + AggregateValue(Basin, V_SOIL_WATER, 0) +
    AggregateValue(Basin, V_SNOW_SWE, 0) + AggregateValue(Basin, V_SOIL_SATFLOW, 0) +
    AggregateValue(Basin, V_SURFACE_DETENTION, 0);
}

/* Water input of the step: precipitation plus the decrease of the snow water
   equivalent plus the vapour flux to the pack, 0 when its absolute value is
   at most 1e-12; sets the mean snow water equivalent at the start of the
   next step */
static float WaterInput(AGGREGATE *Aggregate, float *DeltaSWE)
{
  float Swq = AggregateValue(Aggregate, V_SNOW_SWE, 0);
  float Input;

  *DeltaSWE = Aggregate->OldSwq - Swq;
  Input = AggregateValue(Aggregate, V_MET_PRECIPITATION, 0) + *DeltaSWE +
    AggregateValue(Aggregate, V_SNOW_VAPOR, 0);
  if (fabs(Input) <= 1.e-12)
    Input = 0.;
  SetAggregateValue(Aggregate, V_BALANCE_WATERINPUT, 0, Input);
  Aggregate->OldSwq = Swq;
  return Input;
}

/*****************************************************************************
  StartBalance()

  The state at Start, after the aggregation before the first step: storage
  and the mean snow water equivalent of the basin and of each class
*****************************************************************************/
void StartBalance(TOTALS *Totals)
{
  int i;

  Totals->Storage = Storage(&(Totals->Basin));
  Totals->OldStorage = Totals->Storage;
  Totals->StartStorage = Totals->Storage;
  Totals->CumulativeError = 0.0;
  Totals->Basin.OldSwq = AggregateValue(&(Totals->Basin), V_SNOW_SWE, 0);
  for (i = 0; i < Totals->NClasses; i++)
    Totals->Class[i].OldSwq = AggregateValue(&(Totals->Class[i]), V_SNOW_SWE, 0);
}

/*****************************************************************************
  MassBalance()

  Water balance of the time step from the basin means; water input of the
  basin and of each class
*****************************************************************************/
void MassBalance(TOTALS *Totals)
{
  AGGREGATE *Basin = &(Totals->Basin);
  float Precip = AggregateValue(Basin, V_MET_PRECIPITATION, 0);
  float Swq = AggregateValue(Basin, V_SNOW_SWE, 0);
  float Vapor = AggregateValue(Basin, V_SNOW_VAPOR, 0);
  float CanopyVapor = AggregateValue(Basin, V_CANOPY_SNOWVAPOR, 0);
  float DeepFlux = AggregateValue(Basin, V_SOIL_DEEPFLUX, 0);
  float ET = AggregateValue(Basin, V_EVAP_TOTAL, 0);
  float ChannelInt = AggregateValue(Basin, V_STREAM_INTERCEPTION, 0);
  float ChannelInfiltration = AggregateValue(Basin, V_STREAM_INFILTRATION, 0);
  float ChannelEvap = AggregateValue(Basin, V_EVAP_CHANNEL, 0);
  float DeltaSWE;
  float Output;
  float Input;
  int i;

  WaterInput(Basin, &DeltaSWE);
  for (i = 0; i < Totals->NClasses; i++)
    WaterInput(&(Totals->Class[i]), &Input);

  /* Rain or melt */
  if (Swq > 0 || (Swq == 0 && DeltaSWE > 0))
    Totals->RainOrMelt = AggregateValue(Basin, V_SNOW_OUTFLOW, 0);
  else
    Totals->RainOrMelt = Precip - AggregateValue(Basin, V_MET_SNOWFALL, 0);

  /* Storage change, output and input of the step */
  Totals->Storage = Storage(Basin);
  Output = (ChannelInt - ChannelInfiltration - ChannelEvap) + ET;
  Input = Precip + Vapor + CanopyVapor + DeepFlux;
  Totals->Error = (Totals->Storage - Totals->OldStorage) + Output - Input;
  Totals->CumulativeError += Totals->Error;
  Totals->OldStorage = Totals->Storage;

  /* Sums of the run */
  Totals->CumPrecip += Precip;
  Totals->CumChannelInt += ChannelInt;
  Totals->CumChannelInfiltration += ChannelInfiltration;
  Totals->CumChannelEvap += ChannelEvap;
  Totals->CumET += ET;
  Totals->CumSnowVapor += Vapor + CanopyVapor;
  Totals->CumDeepFlux += DeepFlux;
}

/* Water in the channel segments (segment table order) over Area (m) */
static double NetworkStorage(MODELVIEW *Model, double Area)
{
  Channel *Segment;
  float Storage = 0.0;

  for (Segment = Model->Streams; Segment != NULL; Segment = Segment->next)
    Storage += Segment->storage;
  return Storage / Area;
}

/*****************************************************************************
  StartNetworkStorage()

  The water in the channel network at Start, over the basin area (m)
*****************************************************************************/
void StartNetworkStorage(MODELVIEW *Model, int NMask, TOTALS *Totals)
{
  double Area = (double) NMask * Model->Map->DX * Model->Map->DY;

  Totals->Network[NET_STORAGE] = NetworkStorage(Model, Area);
}

/*****************************************************************************
  NetworkSums()

  Sums over the channel segments (segment table order), divided by the
  basin area: lateral inflow (from cells and lakes), outflow of the
  segments without an outlet, storage change, infiltration, evaporation,
  error, storage
*****************************************************************************/
void NetworkSums(MODELVIEW *Model, int NMask, TOTALS *Totals)
{
  Channel *Segment;
  float LateralInflow = 0.0;
  float Outflow = 0.0;
  float StorageChange = 0.0;
  float Infiltration = 0.0;
  float Evaporation = 0.0;
  float Error;
  double Area = (double) NMask * Model->Map->DX * Model->Map->DY;

  for (Segment = Model->Streams; Segment != NULL; Segment = Segment->next) {
    LateralInflow += (Segment->lateral_inflow + Segment->lake_inflow);
    if (Segment->outlet == NULL)
      Outflow += Segment->outflow;
    StorageChange += Segment->storage - Segment->last_storage;
    Infiltration += Segment->infiltration;
    Evaporation += Segment->evaporation;
  }
  Error = StorageChange - LateralInflow + Outflow + Infiltration + Evaporation;

  Totals->Network[NET_LATERALINFLOW] = LateralInflow / Area;
  Totals->Network[NET_OUTFLOW] = Outflow / Area;
  Totals->Network[NET_STORAGECHANGE] = StorageChange / Area;
  Totals->Network[NET_INFILTRATION] = Infiltration / Area;
  Totals->Network[NET_EVAPORATION] = Evaporation / Area;
  Totals->Network[NET_ERROR] = Error / Area;
  Totals->Network[NET_STORAGE] = NetworkStorage(Model, Area);
  Totals->CumOutlet += Totals->Network[NET_OUTFLOW];
  Totals->CumNetworkError += Totals->Network[NET_ERROR];
}

/*****************************************************************************
  BasinValue()

  Value of a Basin.csv column: a basin mean, or a basin-only variable
*****************************************************************************/
double BasinValue(TOTALS *Totals, int Code, int Layer)
{
  switch (Code) {
  case V_BALANCE_STORAGE: return Totals->Storage;
  case V_BALANCE_ERROR: return Totals->Error;
  case V_BALANCE_CUMULATIVEERROR: return Totals->CumulativeError;
  case V_BALANCE_RAINORMELT: return Totals->RainOrMelt;
  case V_SOIL_SATURATEDFRACTION: return Totals->SaturatedFraction;
  case V_CHANNEL_NETWORKLATERALINFLOW: return Totals->Network[NET_LATERALINFLOW];
  case V_CHANNEL_NETWORKOUTFLOW: return Totals->Network[NET_OUTFLOW];
  case V_CHANNEL_NETWORKSTORAGECHANGE: return Totals->Network[NET_STORAGECHANGE];
  case V_CHANNEL_NETWORKINFILTRATION: return Totals->Network[NET_INFILTRATION];
  case V_CHANNEL_NETWORKEVAPORATION: return Totals->Network[NET_EVAPORATION];
  case V_CHANNEL_NETWORKERROR: return Totals->Network[NET_ERROR];
  case V_CHANNEL_NETWORKSTORAGE: return Totals->Network[NET_STORAGE];
  default:
    return AggregateValue(&(Totals->Basin), Code, Layer);
  }
}

/*****************************************************************************
  FinalMassBalance()

  The water balance of the whole run on stdout, in mm over the basin; then
  the two mass errors: the sum of the step errors, and the final minus
  initial balance (storage change plus output minus input, from the values
  printed); last the sum of the step errors of the channel network (water
  the segments create or lose, which the basin balance, ending at the
  channel inflow, does not see)
*****************************************************************************/
void FinalMassBalance(TOTALS *Totals)
{
  AGGREGATE *Basin = &(Totals->Basin);
  float Input;
  float Output;
  float ChannelNet;

  Input = Totals->CumPrecip + Totals->CumSnowVapor + Totals->CumDeepFlux;
  ChannelNet = Totals->CumChannelInt - Totals->CumChannelInfiltration - Totals->CumChannelEvap;
  Output = ChannelNet + Totals->CumET;

  fflush(stderr);
  printf("Final mass balance (mm over the basin)\n");
  printf("Total input                      %.6f\n", Input * 1000);
  printf("  Precipitation                  %.6f\n", Totals->CumPrecip * 1000);
  printf("  Snow vapour flux               %.6f\n", Totals->CumSnowVapor * 1000);
  printf("  Deep flux                      %.6f\n", Totals->CumDeepFlux * 1000);
  printf("Total output                     %.6f\n", Output * 1000);
  printf("  Evapotranspiration             %.6f\n", Totals->CumET * 1000);
  printf("  Channel inflow, net            %.6f\n", ChannelNet * 1000);
  printf("    Channel inflow               %.6f\n", Totals->CumChannelInt * 1000);
  printf("    Channel re-infiltration      %.6f\n", -Totals->CumChannelInfiltration * 1000);
  printf("    Channel evaporation          %.6f\n", -Totals->CumChannelEvap * 1000);
  printf("Outlet outflow                   %.6f\n", Totals->CumOutlet * 1000);
  printf("Storage change                   %.6f\n", (Totals->Storage - Totals->StartStorage) * 1000);
  printf("  Initial storage                %.6f\n", Totals->StartStorage * 1000);
  printf("  Final storage                  %.6f\n", Totals->Storage * 1000);
  printf("  Final snow water equivalent    %.6f\n", AggregateValue(Basin, V_SNOW_SWE, 0) * 1000);
  printf("  Final soil water               %.6f\n",
         (AggregateValue(Basin, V_SOIL_WATER, 0) + AggregateValue(Basin, V_SOIL_SATFLOW, 0)) * 1000);
  printf("  Final surface water            %.6f\n",
         (AggregateValue(Basin, V_SOIL_SURFACEWATER, 0) + AggregateValue(Basin, V_CANOPY_STORAGE, 0) +
          AggregateValue(Basin, V_SURFACE_DETENTION, 0)) * 1000);
  printf("Cumulative mass error            %.6f\n", Totals->CumulativeError * 1000);
  printf("Final minus initial mass error   %.6f\n",
         (double) ((Totals->Storage - Totals->StartStorage) * 1000) + (double) (Output * 1000) -
         (double) (Input * 1000));
  printf("Channel network error            %.6f\n", Totals->CumNetworkError * 1000);
  fflush(stdout);
}
