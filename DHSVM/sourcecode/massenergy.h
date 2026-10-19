
#ifndef MASSENERGY_H
#define MASSENERGY_H

#include "data.h"
#include <stdarg.h>
#include "DHSVMChannel.h"
#include "photosynthesis.h"
#include "planthydraulics.h"
#include "roothydraulics.h"
#include "stomatalscheme.h"

void AggregateRadiation(int MaxVegLayers, int NVegL, PIXRAD * Rad,
			PIXRAD * TotalRad);

float CanopyResistance(float LAI, float RsMin, float RsMax, float Rpc,
		       float VpdThres, float MoistThres, float WP,
		       float TSoil, float SoilMoisture, float Vpd, float Rp);




  
float Desorption(int Dt, float Moisture, float Porosity, float Ks, 
			   float Press, float m);

void EvapoTranspiration(int Layer, int impvRad, int Dt, PIXMET *Met,
              float NetRad, float Rp, VEGTABLE *VType, SOILTABLE *SType,
              float MoistureFlux, float *Moist, float *Temp, float *Int,
              float *EPot, float *EInt, float **ESoil, float *EAct,
              float *ETot, float *Adjust, float Ra, VEGPIX *LocalVeg, SOILPIX *LocalSoil, OPTIONSTRUCT *Options);

void InitLocalRad(int HeatFluxOption, float Rs, float Ld, float Tair, 
               float Tcanopy, float Tsoil, VEGTABLE *VType, 
               SNOWPIX *LocalSnow, PIXRAD *LocalRad);

void InterceptionStorage(int NAct, float *MaxInt, float *Fract, float *Int,
               float *Precip);

void LongwaveBalance(OPTIONSTRUCT *Options, unsigned char OverStory, 
			   float F, float Vf, float Ld, float Tcanopy, float Tsurf, 
               PIXRAD *LocalRad);

void NoSensibleHeatFlux(int Dt, PIXMET *LocalMet, float ETot, SOILPIX *LocalSoil);

void RadiationBalance(OPTIONSTRUCT *Options, int HeatFluxOption,
              int CanopyRadAttOption, int Overstory, int Understory,
              float SineSolarAltitude, float VICRs,
              float Rs, float Rsd, float Rsb, float Ld, float Tair, float Tcanopy,
              float Tsoil, float SoilAlbedo, VEGTABLE *VType, SNOWPIX *LocalSnow,
              PIXRAD *LocalRad, VEGPIX *LocalVeg);

void SensibleHeatFlux(int y, int x, int Dt, float Ra, float ZRef,
		      float Displacement, float Z0, PIXMET *LocalMet,
		      float NetShort, float LongIn, float ETot, int NSoilLayers, 
		      float *SoilDepth, SOILTABLE *SoilType, float MeltEnergy, 
              SOILPIX *LocalSoil);

void ShortwaveBalance(OPTIONSTRUCT *Options, unsigned char OverStory,
              float F, float Rs, float Rsb, float Rsd, float Tau,
              float Taud, float *Albedo, PIXRAD * LocalRad);

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
                      float *MoistContent, float Adjust);

float StabilityCorrection(float Z, float d, float Tsurf, float Tair,
			  float Wind, float Z0);

float SurfaceEnergyBalance(float TSurf, va_list ap);

/* -------------------------------------------------------------------------
   The single seam between DHSVM and the physiology modules.  Everything the
   hydraulic schemes report comes back through here; unit conversion happens
   in CanopyResistance.c and nowhere else.
   ------------------------------------------------------------------------- */
typedef struct {
  float Rc;                        /* canopy resistance (s/m)              */
  float An;                        /* umol/m2 ground/s                     */
  float Escheme;                   /* m/s, what the scheme predicts        */
  float Ecrit;                     /* m/s, hydraulic failure threshold     */
  float Tsupply;                   /* m/s; < 0 means no supply limit       */
  float PsiLeaf, PsiRoot, Heff;    /* MPa                                  */
  float Tleaf;                     /* degC                                 */
  float VpdLeaf;                   /* kPa                                  */
  float Beta;                      /* empirical factor; 1 with hydraulics  */
  float PLC;                       /* percent loss of conductivity (%)     */
  float SafetyMargin;              /* PsiLeaf - P50 (MPa)                  */
  float Efrac[ROOT_MAXLAYERS];     /* per-layer uptake weights, sum to 1   */
  int   NLayers;
  int   EConsistent;               /* 1: the scheme solved its own leaf
                                      energy balance and Escheme is the
                                      transpiration Penman-Monteith must
                                      reproduce (ProfitMax, ProfitMax2, SOX);
                                      0: Rc is the state (Medlyn paths)   */
} CANOPYHYD;

/* Root-fraction-weighted soil water potential (MPa).  Defined for EVERY
   scheme, including Jarvis: it describes the soil, not the stomata. */
float CanopyRootZonePotential(int Layer, VEGTABLE *VType, SOILTABLE *SType,
                              SOILPIX *LocalSoil, float *Moist);

float CanopyResistanceScheme(int StomScheme, int Hydraulics,
  float Lai, int Layer, float Rp, float NetRad,
  PIXMET *Met, VEGTABLE *VType, SOILTABLE *SType, SOILPIX *LocalSoil,
  float *Moist, float Dormancy, CANOPYHYD *Out);

/* The canopy resistance that makes DHSVM's Penman-Monteith expression return
   a prescribed transpiration Eflux (m/s) given the layer potential rate EPot
   (m/s): the inverse of the factor (Slope+Gamma)/(Slope+Gamma(1+Rc/Ra)).
   Bounded by [RcMin, RcMax].  This is the leaf-basis -> air-basis conversion
   of docs/provenance.md departure #1. */
float CanopyResistanceFromFlux(float Eflux, float EPot, float Slope,
                               float Gamma, float Ra, float RcMin, float RcMax);

#endif
