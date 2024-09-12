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

/*****************************************************************************
RedistributeSnow()
*****************************************************************************/
void RedistributeSnow(OPTIONSTRUCT *Options, int y, int x,
                      float DX, float DY, int Dt,
                      PIXMET *LocalMet, VEGTABLE *VType, VEGPIX *LocalVeg,
                      SNOWPIX *LocalSnow, WINDPIX **WindMap, TOPOPIX **TopoMap, MAPSIZE *Map)
{
  // const char *Routine = "RedistributeSnow";
  int j, k, l, nx, ny;
  float RefWindSpeed;
  float uStar; /* Friction velocity */
  float uStarT; /* Min friction velocity for transport */
  float LAI, HeightDiff;
  float cSalt, uSalt, hSalt, Qsalt;
  float cSusp, uSusp, hSusp, hSuspMid, Qsusp;
  float QsaltIncoming, QsuspIncoming;
  float cSublimation;
  float TravelTime, DeltaZ, UpperZ, LowerZ;
  float OverlapWidth, LayerFrac, MeanLayerHeight;
  float MaxErode;
  float DepthConversion = (WATER_DENSITY*DX*DX);
  
  if (LocalSnow->HasSnow == TRUE && LocalSnow->Swq > MINSNOWBLOW && LocalSnow->TPack < -1.0) {
    
    /* 2 m reference wind speed from meteorology data */
    RefWindSpeed = VType->USnow * LocalMet->Wind;
    
    /***************************************************************************
     * Local equilibrium ground conditions: saltation and suspension
     **************************************************************************/
    
    /* Consider a single uniform ground-equilibrium suspension layer */
    uSusp = RefWindSpeed * WindMap[y][x].WindSpeedXY[0];
    hSusp = WindMap[y][x].LayerElevUpper[0] - TopoMap[y][x].Dem;
    
    /* Friction velocity */
    uStar = uSusp * VON_KARMAN / log(hSusp / Z0_SNOW);
    
    /* Marsh et al. 2020 eq. 17, ref. Li & Pomeroy 1997) */
    uStarT = 0.35 + LocalMet->Tair * (1.0 / 150.0 + LocalMet->Tair / 8200.0);
    
    if (uStar > uStarT) {
      
      if (VType->OverStory == TRUE) {
        LAI = LocalVeg->LAI[0];
        HeightDiff = LocalVeg->Height[0] - LocalSnow->Swq / CONST_SNOW_DENSITY;
      }
      else {
        LAI = 0.0;
        HeightDiff = 0.0;
      }
      
      cSalt = SaltationConcentration(LocalMet->AirDens, WindMap[y][x].FetchDist,
                                     LAI, HeightDiff, uStar, uStarT);
      uSalt = 2.8 * uStarT; /* Pomeroy & Gray 1990 */
      hSalt = 0.08436 * pow(uStar, 1.27); /* Pomeroy & Male 1992 */
      Qsalt = cSalt * uSalt * hSalt * Dt * DX; /* kg */
      
      if (hSusp < (hSalt + 1.0))
        hSusp = hSalt + 1.0;
      hSuspMid = (hSusp - hSalt) / 2.0 + hSalt;
      cSusp = SuspensionConcentration(cSalt, uStar, hSuspMid);
      Qsusp = cSusp * uSusp * (hSusp - hSalt) * Dt * DX; /* kg */
      
      if (Qsalt < 0.0)
        Qsalt = 0.0;
      if (Qsusp < 0.0)
        Qsusp = 0.0;
      
    } else {
      Qsalt = 0.0;
      Qsusp = 0.0;
    }
  } else {
    Qsalt = 0.0;
    Qsusp = 0.0;
  }
  
  /***************************************************************************
   * Subtract sublimation losses from each suspended layer
   **************************************************************************/
  
  for (l = 0; l < NWINDLAYERS; l++) {
    if (WindMap[y][x].Qsusp[l] > 0.0) {

      MeanLayerHeight = (WindMap[y][x].LayerElevUpper[l] +
                         WindMap[y][x].LayerElevLower[l]) / 2.0
                        - TopoMap[y][x].Dem;

      /* Subtract sublimation from current layer */
      cSublimation = CalcSublimation(MeanLayerHeight,
                                     RefWindSpeed * WindMap[y][x].WindSpeedXY[l],
                                     LocalMet->Rh, LocalMet->Tair + 273.15, LocalMet->Es);

      if (cSublimation < 0.0)
        cSublimation = 0.0;
      if (cSublimation > 1.0)
        cSublimation = 1.0;

      WindMap[y][x].Qsusp[l] -= cSublimation * Dt * WindMap[y][x].Qsusp[l];
  }
  }
  
  /***************************************************************************
   * Local cell ground mass balance
   **************************************************************************/
  
  /* Keep track of upwind contributions to calculate flux divergence */
  QsaltIncoming = WindMap[y][x].Qsalt;
  QsuspIncoming = WindMap[y][x].Qsusp[0];
  WindMap[y][x].Qsalt = Qsalt;
  WindMap[y][x].Qsusp[0] = Qsusp;
  
  /* Only subtract the available amount */
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
  LocalSnow->Swq += (QsaltIncoming - WindMap[y][x].Qsalt) / DepthConversion;
  LocalSnow->Swq += (QsuspIncoming - WindMap[y][x].Qsusp[0]) / DepthConversion;
  
  if (isnan(QsaltIncoming))
    printf("\n\nNan SaltIn %d %d\n",x,y);
  if (isnan(QsuspIncoming))
    printf("\n\nNan SuspIn %d %d\n",x,y);
  if (isnan(WindMap[y][x].Qsalt))
    printf("\n\nNan Salt %d %d\n",x,y);
  if (isnan(WindMap[y][x].Qsusp[0]))
    printf("\n\nNan Susp %d %d\n",x,y);
  
  /***************************************************************************
   * Spatial propagation of saltation and suspension
   **************************************************************************/
  
  /* Saltation propagation */
  for (k = 0; k < 4; k++) {
    if (WindMap[y][x].WindDirFrac[0][k] > 0.0) {
      nx = xdirection4[k] + x;
      ny = ydirection4[k] + y;
      
      if (valid_cell(Map, nx, ny) && INBASIN(TopoMap[ny][nx].Mask)) {
        WindMap[ny][nx].Qsalt += WindMap[y][x].Qsalt * WindMap[y][x].WindDirFrac[0][k];
      }
  }
  }
  WindMap[y][x].Qsalt = 0.0;
  
  /* Determine suspended snow travel path for each vertical layer */
  /* Note that each layer can have a different set of downwind directions */
  for (l = 0; l < NWINDLAYERS; l++) {
    if (WindMap[y][x].Qsusp[l] > 0.0) {
      
      /* Calculate vertical trajectory given wind speeds and current position */
      if (RefWindSpeed * WindMap[y][x].WindSpeedXY[l] < 0.1)
        TravelTime = Dt;
      else
        TravelTime = DX / (RefWindSpeed * WindMap[y][x].WindSpeedXY[l]);
      DeltaZ = TravelTime * (RefWindSpeed * WindMap[y][x].WindSpeedZ[l] - SNOWFALLVEL);
      UpperZ = WindMap[y][x].LayerElevUpper[l] + DeltaZ;
      LowerZ = WindMap[y][x].LayerElevLower[l] + DeltaZ;
      
      // if (x==83 && y==75)
      //   printf("DX: %f RefWindSpeed: %f SpeedXY: %f\n",
      //          DX,RefWindSpeed,WindMap[y][x].WindSpeedXY[l]);
      // if (x==83 && y==75)
      //   printf("TravelTime: %f DeltaZ: %f DEM: %f Upper: %f Lower: %f\n",
      //          TravelTime,DeltaZ,WindMap[y][x].LayerElevUpper[l],WindMap[y][x].LayerElevLower[l]);
      
      /* Loop through pre-determined downwind directions */
      for (k = 0; k < 4; k++) {
        if (WindMap[y][x].WindDirFrac[l][k] > 0.0) {
          nx = xdirection4[k] + x;
          ny = ydirection4[k] + y;
          if (valid_cell(Map, nx, ny) && INBASIN(TopoMap[ny][nx].Mask)) {
            
            /* Add suspended snow that intersects the downwind ground surface to the saltation flux */
            if (LowerZ < TopoMap[ny][nx].Dem) {
              OverlapWidth = MIN(TopoMap[ny][nx].Dem,UpperZ) - LowerZ;
              LayerFrac = OverlapWidth / (UpperZ - LowerZ);
              
              // if (x==83 && y==75)
              //   printf("LowerZ: %f UpperZ: %f DEM: %f\n",
              //          LowerZ,UpperZ,TopoMap[ny][nx].Dem);
              // if (x==83 && y==75)
              //   printf("Overlap: %f LayerFrac: %f QsaltDownwind: %f Qsusp: %f WindDirFrac: %f",
              //          OverlapWidth,LayerFrac,WindMap[ny][nx].Qsalt,WindMap[y][x].Qsusp[l],WindMap[y][x].WindDirFrac[l][k]);
              
              WindMap[ny][nx].Qsalt += WindMap[y][x].Qsusp[l] * WindMap[y][x].WindDirFrac[l][k] * LayerFrac;
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
                WindMap[ny][nx].Qsusp[j] += WindMap[y][x].Qsusp[l] * WindMap[y][x].WindDirFrac[l][k] * LayerFrac;
              }
            } /* End of loop over vertical layers in downwind cell */
          } /* End of checking if in basin */
      }
      } /* End of loop over directions from current cell */
  }
  /* Reset snow transport through current layer of current cell */
  WindMap[y][x].Qsusp[l] = 0.0;
  } /* End of loop over vertical layers in current cell */
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
  
  cSublimation = dmdt / mm;
  
  return cSublimation;
}
