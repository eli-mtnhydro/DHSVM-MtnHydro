
#ifndef MASSENERGY_H
#define MASSENERGY_H

#include "data.h"
#include <stdarg.h>
#include "DHSVMChannel.h"

float CanopyResistance(float LAI, float RsMin, float RsMax, float Rpc,
		       float VpdThres, float MoistThres, float WP,
		       float TSoil, float SoilMoisture, float Vpd, float Rp);

float Desorption(int Dt, float Moisture, float Porosity, float Ks, 
			   float Press, float m);

void EvapoTranspiration(int Layer, int Dt, PIXMET *Met,
              float NetRad, float Rp, STORYTABLE *Story, SOILTABLE *SType,
              float MoistureFlux, float *Moist, float *Temp, float *Int,
              float *EPot, float *EInt, float **ESoil, float *EAct,
              float *ETot, float *Adjust, float *LayerThickness, float Ra,
              VEGPIX *LocalVeg);

void InterceptionStorage(int NAct, float *MaxInt, float *Fract, float *Int,
               float *Precip);

void LongwaveBalance(unsigned char OverStory,
			   float F, float Ld, float Tcanopy, float Tsurf,
               float Cu, float Tair, PIXRAD *LocalRad);

void NoSensibleHeatFlux(int Dt, PIXMET *LocalMet, float ETot, SOILPIX *LocalSoil);

void RadiationBalance(OPTIONSTRUCT *Options,
              int Overstory, float FuR, float Cu,
              float SineSolarAltitude, float VICRs,
              float Rs, float Rsb, float Rsd, float Ld, float Tair, float Tcanopy,
              float SoilAlbedo, OVERTABLE *Over, UNDERTABLE *Under, SNOWPIX *LocalSnow,
              PIXRAD *LocalRad, VEGPIX *LocalVeg);

void ShortwaveBalance(unsigned char OverStory,
              float F, float Rs, float Rsb, float Rsd, float Tau,
              float Taud, float *Albedo, float FuR, float AlbedoU, float AlbedoG,
              float TauU, PIXRAD * LocalRad);

float PondEvaporation(int Dt, float DXDY, float ChannelArea,
                      float Temp, float Slope, float Gamma,
                      float Lv, float AirDens, float Vpd, float NetRad, float LowerRa,
                      float Evapotranspiration, float *IExcess);

void ChannelEvaporation(int Dt, float DXDY,
                        float Temp, float Slope, float Gamma,
                        float Lv, float AirDens, float Vpd, float NetRad, float LowerRa,
                        float Evapotranspiration, int x, int y, CHANNEL *ChannelData);

float ChannelSoilEvaporation(int Dt, float DXDY,
                             float Temp, float Slope, float Gamma, float Lv,
                             float AirDens, float Vpd, float NetRad, float RaSoil,
                             float Evapotranspiration, float *Porosity, float *FCap, float *Ks,
                             float *Press, float *m, float LayerThickness,
                             float *MoistContent, float *Adjust,
                             int x, int y, CHANNEL *ChannelData, int CutBankZone);

float SoilEvaporation(int Dt,
                      float Temp, float Slope, float Gamma, float Lv,
                      float AirDens, float Vpd, float NetRad, float RaSoil,
                      float Evapotranspiration, float Porosity, float FCap, float Ks,
                      float Press, float m, float RootDepth,
                      float *MoistContent, float Adjust, float *PotEvap,
                      float *DesorptionAcc);

float StabilityCorrection(float Z, float d, float Tsurf, float Tair,
			  float Wind, float Z0);

#endif
