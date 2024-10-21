/*
 * SUMMARY:      RedistributeSnow.c - Wind redistribution of snowpack
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Eli Boardman
 * ORG:          Mountain Hydrology LLC
 * E-MAIL:       eli.boardman@mountainhydrology.com
 * ORIG-DATE:    2024
 * DESCRIPTION:  Wind redistribution of snowpack based on CHM / PBSM
 *               The novel part here is that a 3D wind field
 *               from the Wind Ninja fluid dynamics model (OpenFOAM)
 *               is used to derive longer-distance, higher-altitude
 *               wind suspension trajectories, which makes wind redistribution
 *               independent of the surface saltation/suspension equilibrium.
 * DESCRIP-END.
 * FUNCTIONS:    RedistributeSnow()
 * COMMENTS:
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "functions.h"
#include "constants.h"
#include "slopeaspect.h"
#include "DHSVMerror.h"

/*****************************************************************************
RedistributeSnow()
*****************************************************************************/
void RedistributeSnow(int y, int x, float DX, float DY, int Dt, int WindIter,
                      SNOWPIX *LocalSnow, WINDPIX **WindMap,
                      TOPOPIX **TopoMap, MAPSIZE *Map)
{
  const char *Routine = "RedistributeSnow";
  int j, k, l, nx, ny;
  float RefWindSpeed, RefWindSpeedOriginal;
  float Qsalt, Qsusp;
  float QsaltIncoming, QsuspIncoming;
  float TravelTime, DeltaZ, DeltaXY, UpperZ, LowerZ, LandingFrac, ExportFrac;
  float OverlapWidth, LayerFrac, LiftRatio;
  float MaxErode, AvgDownwindElev, DownwindSlope, SaltScale, QsaltToSusp;
  float *QsuspLocal;
  float DepthConversion = (WATER_DENSITY*DX*DX);
  float ThreshSaltSlope = 0.5;
  float MeanLayerHeight, SettlingVel, KinEnergy;
  int flag;
  float AvyFrac, AvySlope, SnowHoldingDepth;
  
  if (!(QsuspLocal = (float *) calloc(NWINDLAYERS, sizeof(float))))
    ReportError((char *) Routine, 1);
  
  /* Use adaptive number of wind redistribution iterations for each grid cell */
  if (WindIter == (WindMap[y][x].nIters - 2)) {
    for (l = 0; l < NWINDLAYERS; l++) {
      WindMap[y][x].QsuspLastIt[l] = WindMap[y][x].Qsusp[l];
    }
  }
  flag = 0;
  if (WindIter == (WindMap[y][x].nIters - 1)) {
    /* Check if any layers changed by at least 50% on the last iteration */
    for (l = 0; l < NWINDLAYERS; l++) {
      if (ABSVAL(WindMap[y][x].Qsusp[l] - WindMap[y][x].QsuspLastIt[l]) /
          WindMap[y][x].Qsusp[l] > 0.5 &&
            WindMap[y][x].Qsusp[l] > (1e-6 * DepthConversion)) {
        if (WindMap[y][x].nIters < NwindIters)
          WindMap[y][x].nIters++;
        flag = 1;
        break;
      }
    }
    if (flag == 0 && WindMap[y][x].nIters > 2)
      WindMap[y][x].nIters--;
  }
  
  /* Reference wind speed from meteorology data */
  RefWindSpeed = WindMap[y][x].Wind;
  RefWindSpeedOriginal = RefWindSpeed;
  
  /***************************************************************************
   * Local equilibrium ground conditions: saltation and suspension
   **************************************************************************/
  
  Qsalt = WindMap[y][x].QsaltLocal / (float) WindMap[y][x].nIters;
  Qsusp = WindMap[y][x].QsuspLocal / (float) WindMap[y][x].nIters;
  
  /***************************************************************************
   * Subtract sublimation losses from each suspended layer
   **************************************************************************/
  
  for (l = 0; l < NWINDLAYERS; l++) {
    WindMap[y][x].Qsusp[l] -= WindMap[y][x].SublimationFrac[l] * WindMap[y][x].Qsusp[l];
  }
  
  /***************************************************************************
   * Local cell ground mass balance
   **************************************************************************/
  
  /* Keep track of upwind contributions to calculate flux divergence */
  QsaltIncoming = WindMap[y][x].Qsalt;
  QsuspIncoming = WindMap[y][x].Qsusp[0];
  WindMap[y][x].Qsalt = Qsalt;
  WindMap[y][x].Qsusp[0] = Qsusp;
  
  /* Only subtract up to the available amount */
  MaxErode = (LocalSnow->Swq - MINSNOWBLOW) * DepthConversion;
  if (WindMap[y][x].Qsalt > QsaltIncoming + MaxErode) {
    WindMap[y][x].Qsalt = QsaltIncoming + MaxErode;
    WindMap[y][x].Qsusp[0] = 0.0;
  } else if (WindMap[y][x].Qsalt + WindMap[y][x].Qsusp[0] >
               QsaltIncoming + QsuspIncoming + MaxErode) {
    WindMap[y][x].Qsusp[0] = QsaltIncoming - WindMap[y][x].Qsalt + QsuspIncoming + MaxErode;
  }
  
  /* Add or subtract flux divergence from local ground snow storage */
  /* Also convert kg to water depth equivalent in each cell */
  WindMap[y][x].WindDeposition += (QsaltIncoming - WindMap[y][x].Qsalt) / DepthConversion;
  WindMap[y][x].WindDeposition += (QsuspIncoming - WindMap[y][x].Qsusp[0]) / DepthConversion;
  
  /***************************************************************************
   * Spatial propagation of saltation and suspension
   **************************************************************************/
  
  for (l = 0; l < NWINDLAYERS; l++) {
    QsuspLocal[l] = 0.0;
  }
  QsaltToSusp = 0.0;
  
  /* Saltation propagation */
  AvgDownwindElev = 0.0;
  for (k = 0; k < 4; k++) {
    if (WindMap[y][x].WindDirFrac[0][k] > 0.0) {
      nx = xdirection4[k] + x;
      ny = ydirection4[k] + y;
      if (valid_cell(Map, nx, ny) && INBASIN(TopoMap[ny][nx].Mask)) {
        
        DownwindSlope = (TopoMap[ny][nx].Dem - TopoMap[y][x].Dem) / DX;
        if (DownwindSlope >= 0.0) {
          SaltScale = 1.0;
        } else {
          /* Steep downhill saltation becomes partially airborne and is propagated with suspension */
          SaltScale = MAX(0.0, (ThreshSaltSlope + DownwindSlope) / ThreshSaltSlope);
          QsaltToSusp += (1.0 - SaltScale) * WindMap[y][x].Qsalt * WindMap[y][x].WindDirFrac[0][k];
        }
        WindMap[ny][nx].Qsalt += SaltScale * WindMap[y][x].Qsalt * WindMap[y][x].WindDirFrac[0][k];
        
        /* Also find the average downwind cell elevation for use later */
        AvgDownwindElev += TopoMap[ny][nx].Dem * WindMap[y][x].WindDirFrac[0][k];
      }
    }
  }
  WindMap[y][x].Qsalt = 0.0;
  
  /* Determine suspended snow travel path for each vertical layer */
  /* Note that each layer can have a different set of downwind directions */
  for (l = 0; l < NWINDLAYERS; l++) {
    if (WindMap[y][x].Qsusp[l] > 0.0 || (l == 0 && QsaltToSusp > 0.0)) {
      
      MeanLayerHeight = (WindMap[y][x].LayerElevUpper[l] +
                         WindMap[y][x].LayerElevLower[l]) / 2.0
                        - TopoMap[y][x].Dem;
      
      /* If it is snowing, wind speed is reduced in proportion to distance from cloud base,
         assuming that the clouds are effectively stationary (i.e., orographically driven) */
      RefWindSpeed = RefWindSpeedOriginal * 
                     (WindMap[y][x].IsSnowing == TRUE ? WindMap[y][x].SnowingScale[l] : 1.0);
      
      /* Lehning et al. 2008 eq. 10 */
      SettlingVel = SNOWFALLVEL -
                    (RefWindSpeed * WindMap[y][x].WindSpeedXY[l] * VON_KARMAN /
                     log(MeanLayerHeight / Z0_SNOW));
      SettlingVel = MAX(SettlingVel, 0.01 * SNOWFALLVEL);
      
      /* Calculate vertical trajectory given wind speeds and current position */
      if (RefWindSpeed * WindMap[y][x].WindSpeedXY[l] < SettlingVel)
        TravelTime = (WindMap[y][x].LayerElevUpper[l] - WindMap[y][x].LayerElevLower[l]) / SettlingVel;
      else
        TravelTime = DX / (RefWindSpeed * WindMap[y][x].WindSpeedXY[l]);
      
      DeltaXY = TravelTime * RefWindSpeed * WindMap[y][x].WindSpeedXY[l];
      
      // if (WindMap[y][x].WindSpeedZ[l] > WindMap[y][x].WindSpeedXY[l])
      //   LiftRatio = 0.0;
      // else if (WindMap[y][x].WindSpeedZ[l] > 0.0)
      //   LiftRatio = 1.0 - WindMap[y][x].WindSpeedZ[l] / WindMap[y][x].WindSpeedXY[l];
      // else
      //   LiftRatio = 1.0;
      LiftRatio = 1.0;
      DeltaZ = TravelTime * (RefWindSpeed * WindMap[y][x].WindSpeedZ[l] * LiftRatio - SettlingVel);
      
      UpperZ = WindMap[y][x].LayerElevUpper[l] + DeltaZ;
      LowerZ = WindMap[y][x].LayerElevLower[l] + DeltaZ;
      
      /* Account for amount potentially landing in current grid cell */
      if (UpperZ < AvgDownwindElev) {
        LandingFrac = 1.0;
        ExportFrac = 0.0;
      }
      else if (LowerZ < AvgDownwindElev) {
        LandingFrac = (AvgDownwindElev - LowerZ) / (UpperZ - LowerZ);
        ExportFrac = MIN(1.0 - LandingFrac, DeltaXY / (ABSVAL(DeltaZ) + DeltaXY));
      }
      else {
        LandingFrac = 0.0;
        ExportFrac = DeltaXY / (ABSVAL(DeltaZ) + DeltaXY);
      }
      
      if (LandingFrac > 0.0)
        WindMap[y][x].WindDeposition += WindMap[y][x].Qsusp[l] * LandingFrac / DepthConversion;
      
      if (DeltaZ > 0.0) {
        if (l < NWINDLAYERS - 1)
          QsuspLocal[l+1] += WindMap[y][x].Qsusp[l] * (1.0 - ExportFrac - LandingFrac);
        else
          QsuspLocal[l] += WindMap[y][x].Qsusp[l] * (1.0 - ExportFrac - LandingFrac);
      } else if (DeltaZ < 0.0) {
        if (l > 0)
          QsuspLocal[l-1] += WindMap[y][x].Qsusp[l] * (1.0 - ExportFrac - LandingFrac);
        else
          WindMap[y][x].Qsalt += WindMap[y][x].Qsusp[l] * (1.0 - ExportFrac - LandingFrac);
      } else {
        QsuspLocal[l] += WindMap[y][x].Qsusp[l] * (1.0 - ExportFrac - LandingFrac);
      }
      
      /* Loop through pre-determined downwind directions */
      if (ExportFrac > 0.0) {
        for (k = 0; k < 4; k++) {
          if (WindMap[y][x].WindDirFrac[l][k] > 0.0) {
            nx = xdirection4[k] + x;
            ny = ydirection4[k] + y;
            if (valid_cell(Map, nx, ny) && INBASIN(TopoMap[ny][nx].Mask)) {
              
              /* Add suspended snow that intersects the downwind ground surface to the saltation flux */
              if (LowerZ < TopoMap[ny][nx].Dem) {
                OverlapWidth = MIN(TopoMap[ny][nx].Dem,UpperZ) - LowerZ;
                LayerFrac = OverlapWidth / (UpperZ - LowerZ);
                
                WindMap[ny][nx].WindDeposition += WindMap[y][x].Qsusp[l] * WindMap[y][x].WindDirFrac[l][k] *
                                                  ExportFrac * LayerFrac / DepthConversion;
              } else {
                LayerFrac = 0.0;
              }
              
              /* Redistribute snow from current layer to appropriate layer(s) of downwind cell */
              if (LayerFrac < 1.0) {
                for (j = 0; j < NWINDLAYERS; j++) {
                  if (j == (NWINDLAYERS - 1))
                    OverlapWidth = UpperZ - MAX(WindMap[ny][nx].LayerElevLower[j],LowerZ);
                  else
                    OverlapWidth = MIN(WindMap[ny][nx].LayerElevUpper[j],UpperZ) -
                                   MAX(WindMap[ny][nx].LayerElevLower[j],LowerZ);
                  LayerFrac = (OverlapWidth > 0.0 ? OverlapWidth / (UpperZ - LowerZ) : 0.0);
                  WindMap[ny][nx].Qsusp[j] += WindMap[y][x].Qsusp[l] * WindMap[y][x].WindDirFrac[l][k] * ExportFrac * LayerFrac;
                  if (l == 0)
                    WindMap[ny][nx].Qsusp[j] += QsaltToSusp * WindMap[y][x].WindDirFrac[l][k] * LayerFrac;
                }
              } /* End of loop over vertical layers in downwind cell */
            } /* End of checking if in basin */
      }
      } /* End of loop over directions from current cell */
      }
  }
  /* Reset snow transport through current layer of current cell */
  WindMap[y][x].Qsusp[l] = 0.0;
  } /* End of loop over vertical layers in current cell */
  
  for (l = 0; l < NWINDLAYERS; l++) {
    WindMap[y][x].Qsusp[l] = QsuspLocal[l];
  }
  
  /* Simple avalanche routine */
  // for (k = 0; k < 8; k++) {
  //   nx = xdirection8[k] + x;
  //   ny = ydirection8[k] + y;
  //   if (valid_cell(Map, nx, ny) && INBASIN(TopoMap[ny][nx].Mask)) {
  // 
  //     AvySlope = (TopoMap[y][x].Dem - TopoMap[ny][nx].Dem) / DX;
  //     SnowHoldingDepth = 10.0 * exp(-5.0 * AvySlope);
  // 
  //     if (AvySlope > 0 && LocalSnow->Swq > SnowHoldingDepth) {
  //       // AvyFrac = MIN(1.0, (AvySlope - 0.5) / 0.5);
  //       WindMap[ny][nx].WindDeposition += (LocalSnow->Swq - SnowHoldingDepth);
  //       WindMap[y][x].WindDeposition -= (LocalSnow->Swq - SnowHoldingDepth);
  //     }
  //   }
  // }
}

/*****************************************************************************
 BlowingSnowConditions()
*****************************************************************************/
void BlowingSnowConditions(int y, int x, float DX, float DY, int Dt,
                           VEGTABLE *VType, VEGPIX *LocalVeg, SNOWPIX *LocalSnow,
                           WINDPIX **WindMap, TOPOPIX **TopoMap)
{
  float RefWindSpeed;
  float Z0eff; /* Effective roughness considering turbulence */
  float uStar; /* Friction velocity */
  float uStarT; /* Min friction velocity for transport */
  float LAI, HeightDiff;
  float cSalt, uSalt, hSalt, Qsalt;
  float cSusp, uSusp, hSusp, hSuspMid, Qsusp;
  
  /* Reference wind speed from meteorology data */
  RefWindSpeed = WindMap[y][x].Wind;
  
  /***************************************************************************
   * Local equilibrium ground conditions: saltation and suspension
   **************************************************************************/
  
  /* Consider a single uniform ground-equilibrium suspension layer */
  uSusp = RefWindSpeed * WindMap[y][x].WindSpeedXY[0];
  hSusp = WindMap[y][x].LayerElevUpper[0] - TopoMap[y][x].Dem;
  
  /* Friction velocity */
  /* Note that hSusp / 2.0 because uSusp is for the middle of layer 0 */
  uStar = uSusp * VON_KARMAN / log((hSusp / 2.0) / Z0_SNOW);
  
  /* Marsh et al. 2020 eq. 17, ref. Li & Pomeroy 1997 */
  uStarT = 0.35 + WindMap[y][x].Tair * (1.0 / 150.0 + WindMap[y][x].Tair / 8200.0);
  
  /* Pomeroy & Male 1992 */
  hSalt = 0.08436 * pow(uStar, 1.27);
  
  /* Pomeroy & Gray 1990 */
  uSalt = 2.8 * uStarT;
  
  /* Saltation calculations */
  if (LocalSnow->HasSnow == TRUE &&
      LocalSnow->Swq > MINSNOWBLOW &&
      LocalSnow->TSurf < -1.0 &&
      uStar > uStarT) {
    
    if (VType->OverStory == TRUE) {
      LAI = LocalVeg->LAI[0];
      HeightDiff = LocalVeg->Height[0] - LocalSnow->Swq / CONST_SNOW_DENSITY;
    }
    else {
      LAI = 0.0;
      HeightDiff = 0.0;
    }
    
    cSalt = SaltationConcentration(WindMap[y][x].AirDens, WindMap[y][x].FetchDist,
                                   LAI, HeightDiff, uStar, uStarT);
    Qsalt = cSalt * uSalt * hSalt * Dt * DX; /* kg */
    
  } else {
    Qsalt = 0.0;
  }
  
  /* Average with upwind incoming saltation flux to smooth out the gradient */
  Qsalt = (Qsalt + WindMap[y][x].Qsalt) / 2.0;
  
  /* Suspension calculations */
  if (Qsalt > 0.0) {
    
    cSalt = Qsalt / (uSalt * hSalt * Dt * DX);
    
    if (hSusp < (hSalt + 1.0))
      hSusp = hSalt + 1.0;
    hSuspMid = (hSusp - hSalt) / 2.0 + hSalt;
    cSusp = SuspensionConcentration(cSalt, uStar, hSuspMid);
    Qsusp = cSusp * uSusp * (hSusp - hSalt) * Dt * DX; /* kg */
    
    if (Qsusp < 0.0)
      Qsusp = 0.0;
  } else {
    Qsalt = 0.0;
    Qsusp = 0.0;
  }
  
  /* Average with upwind incoming suspension flux to smooth out the gradient */
  Qsusp = (Qsusp + WindMap[y][x].Qsusp[0]) / 2.0;
  
  WindMap[y][x].QsaltLocal = Qsalt;
  WindMap[y][x].QsuspLocal = Qsusp;
}

/*****************************************************************************
 SaltationConcentration()
 *****************************************************************************/
float SaltationConcentration(float AirDens, float FetchDist,
                             float LAI, float HeightDiff,
                             float uStar, float uStarT) {
  float cSalt; /* Saltation concentration, kg/m^3 */
  float lambda, ShearPart; /* For vegetation roughness */
  
  if (LAI > 0.0){
    /* From CHM code: 
       LAI/2.0 suggestion from Raupach 1994 (DOI:10.1007/BF00709229) Section 3(a) */
    lambda = 0.5 * LAI * HeightDiff;
    
    /* Marsh eq. 18-19, ref. Raupach et al. 1993) */
    ShearPart = 0.16 * 202.0 * lambda / (1 + 0.16 * 202.0 * lambda);
  } else {
    ShearPart = 0.0;
  }
  
  /* Marsh eq. 16, ref. Pomeroy & Gray 1990, Pomeroy & Male 1992 */
  cSalt = (AirDens / (3.29 * uStar)) * (1.0 - ShearPart - (uStarT * uStarT) / (uStar * uStar));
  
  /* Marsh eq. 22, ref. Pomeroy & Male 1986) */
  /* Note that the 0.5 must be added before multiplying
     (there is a missing parentheses in the Marsh paper eq.) */
  if (FetchDist < 300.0)
    cSalt *= (tanh(4.0 * FetchDist / 300.0 - 2.0) / 2.0 + 0.5);
  
  return cSalt;
}

/*****************************************************************************
 SuspensionConcentration()
 *****************************************************************************/
float SuspensionConcentration(float cSalt, float uStar, float hSusp) {
  float cSusp; /* Suspension concentration, kg/m^3 */
  
  /* Pomeroy et al. 1993 eq. 7 */
  cSusp = cSalt * exp(-1.55 * (pow(0.05628 * uStar, -0.544) - pow(hSusp, -0.544)));
  
  return cSusp;
}


/***************************************************************************
 * WindSublimation()
 **************************************************************************/
void WindSublimation(int y, int x, float DX, float DY, int Dt,
                     PIXMET *LocalMet, WINDPIX **WindMap, TOPOPIX **TopoMap)
{
  int l;
  float RefWindSpeed, RefWindSpeedOriginal;
  float MeanLayerHeight, TravelTime;
  float cSublimation;
  float DepthConversion = (WATER_DENSITY*DX*DX);
  
  /* Reference wind speed from meteorology data */
  RefWindSpeed = WindMap[y][x].Wind;
  RefWindSpeedOriginal = RefWindSpeed;
  
  for (l = 0; l < NWINDLAYERS; l++) {
    if (WindMap[y][x].Qsusp[l] > (1e-9 * DepthConversion)) {
      
      MeanLayerHeight = (WindMap[y][x].LayerElevUpper[l] +
                         WindMap[y][x].LayerElevLower[l]) / 2.0
                        - TopoMap[y][x].Dem;
      
      if (WindMap[y][x].IsSnowing == TRUE &&
          l >= WindMap[y][x].SnowfallLayer) {
        cSublimation = 0.0;
      } else {
        
        /* If it is snowing, wind speed is reduced in proportion to distance from cloud base,
         assuming that the clouds are effectively stationary (i.e., orographically driven) */
        RefWindSpeed = RefWindSpeedOriginal *
          (WindMap[y][x].IsSnowing == TRUE ? WindMap[y][x].SnowingScale[l] : 1.0);
        
        /* Subtract sublimation from current layer */
        cSublimation = CalcSublimation(MeanLayerHeight,
                                       RefWindSpeed * WindMap[y][x].WindSpeedXY[l],
                                       (WindMap[y][x].IsSnowing == TRUE ? 1.0 : LocalMet->Rh / 100.0),
                                       LocalMet->Tair + 273.15, LocalMet->Es);
        
        /* Sublimation is NOT calculated over the full model time step,
         but rather just over the inter-grid-cell effective travel time */
        if (RefWindSpeed * WindMap[y][x].WindSpeedXY[l] <  SNOWFALLVEL)
          TravelTime = (WindMap[y][x].LayerElevUpper[l] - WindMap[y][x].LayerElevLower[l]) / SNOWFALLVEL;
        else
          TravelTime = DX / (RefWindSpeed * WindMap[y][x].WindSpeedXY[l]);
        
        cSublimation *= TravelTime;
        
        if (cSublimation < 0.0)
          cSublimation = 0.0;
        if (cSublimation > 1.0)
          cSublimation = 1.0;
      }
    } else {
      cSublimation = 1.0;
    }
    WindMap[y][x].SublimationFrac[l] = cSublimation;
  }
}

/*****************************************************************************
 CalcSublimation()
 *****************************************************************************/
float CalcSublimation(float zHeight, float xySpeed,
                      float RH, float TairK, float SatVapPress) {
  float cSublimation; /* Sublimation coeff, 1/s */
  float dmdt;
  float sigma, lambdaT, SatDensity, rm, Qr, expr1;
  float mm_alpha, mm, rmm, xrz, Vr, Re, NuSh, Diffusivity;
  
  /* Marsh eq. 13, ref. Pomeroy & Li 2000 */
  sigma = (RH - 1.0) * (1.019 + 0.27 * log(zHeight));
  
  /* From CHM code:
     Thermal conductivity, use Harder 2013 A.9, Pomeroy's
     is off by an order of magnitude, this matches this
     https://www.engineeringtoolbox.com/air-properties-d_156.html */
  lambdaT = 0.000063 * TairK + 0.00673;
  
  /* CHM code */
  SatDensity = (MOLWEIGHTH20 * SatVapPress) / (GASR * TairK);
  
  /* From CHM code:
     Eq. 18, mean particle radius
     This is 'r_r' in Pomeroy and Gray 1995, eq. 53 */
   rm = 4.6e-5 * pow(zHeight, -0.258);
  
  /* From CHM code:
     radiative energy absorbed by the particle -- take from CRHM's
     PBSM implementation
     120.0 = PBSM_constants::Qstar (Solar Radiation Input)
     0.9 comes from Schmidt (1972) assuming a snow particle albedo of 0.5
     and a snow surface albedo of 0.8
     rm is used here as in Liston and Sturm (1998) */
  Qr = 0.9 * PI * rm * rm * 120.0;
  
  /* CHM code, mainly from Pomeroy et al. 1993 */
  /* Mean particle mass */
  mm_alpha = 4.08 + 12.6 * zHeight;
  mm = 4.0 / 3.0 * PI * ICE_DENSITY * rm * rm * rm *
       (1.0 + 3.0 / mm_alpha + 2.0 / (mm_alpha * mm_alpha));
  /* Mean radius of mean-mass particle */
  rmm = pow((3.0 * mm) / (4 * PI * ICE_DENSITY), 0.333);
  xrz = 0.005 * pow(xySpeed, 1.36);
  Vr = SNOWFALLVEL + 3.0 * xrz * cos(M_PI / 4.0);
  Re = 2.0 * rmm * Vr / AIRKINVISC;
  NuSh = 1.79 + 0.606 * pow(Re, 0.5);
  
  /* From CHM code:
     Diffusivity of water vapour in air, TairK in K,
     Eq. A-7 in Liston 1998 or Harder 2013 A.6 */
  Diffusivity = 2.06e-5 * pow(TairK / 273.15, 1.75);
  
  /* Pomeroy et al. 1993 eq. 11 */
  expr1 = (LSUB * MOLWEIGHTH20 / (GASR * TairK) - 1.0) / lambdaT * TairK * NuSh;
  dmdt = 2 * PI * rm * sigma - Qr * expr1;
  dmdt /= LSUB * expr1 + 1.0 / (Diffusivity * SatDensity * NuSh);
  
  cSublimation = -1.0 * dmdt / mm;
  
  return cSublimation;
}
