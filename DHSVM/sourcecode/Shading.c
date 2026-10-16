/*
 * DESCRIPTION:  Beam shading computed each time step from the horizon
 *               angles of every cell (36 azimuths, every 10 degrees
 *               clockwise from north). The time step is divided into equal
 *               intervals of at most 10 minutes, with a sun position at the
 *               middle of each. At a position with the sun above the
 *               horizon a cell is lit when its horizon angle toward the sun
 *               (interpolated between the two neighbouring azimuths) is
 *               below the sun altitude, and then receives the beam at the
 *               incidence angle of its slope and aspect. The factor of a
 *               cell is the sum of the cosine of incidence over its lit
 *               positions divided by the sum of the sine of the sun
 *               altitude over the positions above the horizon: the beam on
 *               the cell relative to the beam on a flat, unshaded surface.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"
#include "mapio.h"

#define SUN_INTERVAL 600			/* Longest interval of a sun position (s) */

/*****************************************************************************
  InitShading()

  Reads the horizon angle map (opened and checked by ReadInputMaps(), its
  values checked by ValidateInputs()) into hundredths of a degree for the
  mask cells, closes it, keeps the sky view map read by ReadInputMaps(), and
  allocates the sun positions of a time step of Dt seconds
*****************************************************************************/
void InitShading(OPTIONSTRUCT *Options, MAPSIZE *Map, TOPOPIX **TopoMap, INPUTMAPS *Raw,
                 int Dt, SHADE *Shade)
{
  const char *Routine = "InitShading";
  size_t NCells = (size_t) Map->NY * Map->NX;
  size_t i;
  float *Layer;
  int x, y, k, n;

  Message("Initializing the terrain shading\n");

  Shade->NX = Map->NX;
  if (!(Shade->Horizon = (unsigned short *) calloc(NCells * MAP_NAZIMUTHS,
                                                    sizeof(unsigned short))))
    ReportError((char *) Routine, 1);
  if (!(Layer = (float *) calloc(NCells, sizeof(float))))
    ReportError((char *) Routine, 1);
  for (k = 0; k < MAP_NAZIMUTHS; k++) {
    ReadMapLayer(&(Options->HorizonAngle), k, Layer, Map);
    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        if (INBASIN(TopoMap[y][x].Mask))
          Shade->Horizon[i * MAP_NAZIMUTHS + k] =
            (unsigned short) floor((double) Layer[i] * 100.0 + 0.5);
      }
    }
  }
  free(Layer);
  CloseMap(&(Options->HorizonAngle));
  Shade->SkyView = Raw->SkyView;
  Raw->SkyView = NULL;

  n = (Dt + SUN_INTERVAL - 1) / SUN_INTERVAL;
  Shade->NSun = n;
  Shade->NUp = 0;
  Shade->SumSinAlt = 0.0;
  if (!(Shade->SinAlt = (double *) calloc(n, sizeof(double))) ||
      !(Shade->CosAlt = (double *) calloc(n, sizeof(double))) ||
      !(Shade->AltDeg = (double *) calloc(n, sizeof(double))) ||
      !(Shade->CosAz = (double *) calloc(n, sizeof(double))) ||
      !(Shade->SinAz = (double *) calloc(n, sizeof(double))) ||
      !(Shade->Az0 = (int *) calloc(n, sizeof(int))) ||
      !(Shade->Az1 = (int *) calloc(n, sizeof(int))) ||
      !(Shade->AzWeight = (double *) calloc(n, sizeof(double))))
    ReportError((char *) Routine, 1);
}

/*****************************************************************************
  ShadeStep()

  The sun positions of the time step DayStep (from 0 at midnight) of Dt
  seconds that are above the horizon, from the day's latitude, declination
  and time adjustment (SolarDay()). Local standard hour of position j
  (1..n): (DayStep x Dt + (j - 0.5) x Dt / n) / 3600.
*****************************************************************************/
void ShadeStep(SHADE *Shade, SOLARGEOMETRY *SolarGeo, int DayStep, int Dt)
{
  double Lat = SolarGeo->Latitude;
  double Decl = SolarGeo->Declination;
  double Hour, w, SinAlt, Alt, Az, a;
  int j, k, n = Shade->NSun;

  Shade->NUp = 0;
  Shade->SumSinAlt = 0.0;
  for (j = 1; j <= n; j++) {
    /* Hour angle of the solar hour */
    Hour = ((double) DayStep * Dt + (j - 0.5) * Dt / n) / 3600.0;
    w = (Hour + SolarGeo->TimeAdjustment - 12.0) * PI / 12.0;

    SinAlt = sin(Lat) * sin(Decl) + cos(Lat) * cos(Decl) * cos(w);
    if (!(SinAlt > 0.0))
      continue;

    /* Azimuth clockwise from north in [0, 2 pi) */
    Az = atan2(-cos(Decl) * sin(w), sin(Decl) * cos(Lat) - cos(Decl) * cos(w) * sin(Lat));
    if (Az < 0.0)
      Az += 2.0 * PI;
    Alt = asin(SinAlt);

    k = Shade->NUp;
    Shade->SinAlt[k] = SinAlt;
    Shade->CosAlt[k] = cos(Alt);
    Shade->AltDeg[k] = Alt * 180.0 / PI;
    Shade->CosAz[k] = cos(Az);
    Shade->SinAz[k] = sin(Az);

    /* The two horizon azimuths around the sun azimuth */
    a = Az * 180.0 / PI / (360.0 / MAP_NAZIMUTHS);
    Shade->Az0[k] = ((int) floor(a)) % MAP_NAZIMUTHS;
    Shade->Az1[k] = (Shade->Az0[k] + 1) % MAP_NAZIMUTHS;
    Shade->AzWeight[k] = a - floor(a);

    Shade->SumSinAlt += SinAlt;
    Shade->NUp++;
  }
}

/*****************************************************************************
  ShadeFactor()

  The beam factor of the time step for the mask cell (y, x): 0 when the sun
  is below the horizon during the whole step
*****************************************************************************/
double ShadeFactor(SHADE *Shade, TOPOPIX *Topo, int y, int x)
{
  unsigned short *H = Shade->Horizon + ((size_t) y * Shade->NX + x) * MAP_NAZIMUTHS;
  double Slope, CosSlope, SinSlope, CosAspect, SinAspect;
  double H0, H1, Horizon, CosI;
  double Lit = 0.0;
  int j;

  if (Shade->NUp == 0)
    return 0.0;

  /* Slope angle from the slope (a tangent); aspect clockwise from north */
  Slope = atan((double) Topo->Slope);
  CosSlope = cos(Slope);
  SinSlope = sin(Slope);
  CosAspect = cos((double) Topo->Aspect);
  SinAspect = sin((double) Topo->Aspect);

  for (j = 0; j < Shade->NUp; j++) {
    H0 = H[Shade->Az0[j]] / 100.0;
    H1 = H[Shade->Az1[j]] / 100.0;
    Horizon = H0 + Shade->AzWeight[j] * (H1 - H0);
    if (Horizon < Shade->AltDeg[j]) {
      CosI = CosSlope * Shade->SinAlt[j] + Shade->CosAlt[j] * SinSlope *
        (Shade->CosAz[j] * CosAspect + Shade->SinAz[j] * SinAspect);
      if (CosI > 0.0)
        Lit += CosI;
    }
  }

  return Lit / Shade->SumSinAlt;
}
