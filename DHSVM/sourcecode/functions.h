
#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "data.h"
#include "channel.h"
#include "DHSVMChannel.h"

void Aggregate(MAPSIZE *Map, OPTIONSTRUCT *Options, TOPOPIX **TopoMap,
	       LAYER *Soil, LAYER *Veg, VEGPIX **VegMap, EVAPPIX **Evap,
	       PRECIPPIX **Precip, PIXRAD **RadMap, SNOWPIX **Snow,
	       SOILPIX **SoilMap, AGGREGATED *Total, VEGTABLE *VType,
	       NETSTRUCT **Network, CHANNEL *ChannelData, int Dt, int NDaySteps);

void Avalanche(MAPSIZE *Map, TOPOPIX **TopoMap, TIMESTRUCT *Time, OPTIONSTRUCT *Options,
  SNOWPIX **SnowMap);

void CalcAerodynamic(int NVegLayers, unsigned char OverStory,
		     float n, float *Height, float Trunk, float *U,
		     float *U2mSnow, float *Ra, float *RaSnow);

float CalcBagnold(float DS, TIMESTRUCT *Time, float, float, float, float);

double CalcDistance(COORD *LocA, COORD *LocB);

float CalcEffectiveKh(int NSoilLayers, float Top, float Bottom,
		      float *SoilDepth, float *KhDry, float *KhSol,
		      float *Moisture, float *Porosity, float *TSoil);

float CalcKhDry(float Density);

float CalcSnowAlbedo(SNOWPIX *LocalSnow, int StepsPerDay);

float CalcTransmissivity(float SoilDepth, float WaterTable, float LateralKs,
			 float KsExponent, float DepthThresh);

void CalcWeights(METLOCATION *Station, int NStats, int NX, int NY,
		 uchar **BasinMask, uchar ****WeightArray,
		 OPTIONSTRUCT *Options);

double ChannelCulvertSedFlow(int y, int x, CHANNEL * ChannelData, int i);

void CheckOut(OPTIONSTRUCT *Options, LAYER Veg, LAYER Soil,
	      VEGTABLE *VType, SOILTABLE *SType, MAPSIZE *Map, 
	      TOPOPIX **TopoMap, VEGPIX **VegMap, SOILPIX **SoilMap);

unsigned char dequal(double a, double b);

void DumpMap(MAPSIZE *Map, DATE *Current, MAPDUMP *DMap, TOPOPIX **TopoMap,
	     EVAPPIX **EvapMap, PRECIPPIX **PrecipMap, PIXRAD **RadMap,
	     SNOWPIX **Snowap, SOILPIX **SoilMap, LAYER *Soil, VEGPIX **VegMap, 
         LAYER *Veg, NETSTRUCT **Network, OPTIONSTRUCT *Options);

void DumpPix(DATE *Current, int first, FILES *OutFile, EVAPPIX *Evap,
        PRECIPPIX *Precip, PIXRAD *Rad, SNOWPIX *Snow, SOILPIX *Soil,
        VEGPIX *Veg, int NSoil, int NVeg, OPTIONSTRUCT *Options, int flag);

#ifdef TOPO_DUMP
void DumpTopo(MAPSIZE *Map, TOPOPIX **TopoMap);
#endif

void ExecDump(MAPSIZE *Map, DATE *Current, DATE *Start, OPTIONSTRUCT *Options,
	      DUMPSTRUCT *Dump, TOPOPIX **TopoMap, EVAPPIX **EvapMap, PIXRAD **RadiMap,
	      PRECIPPIX ** PrecipMap, SNOWPIX **SnowMap, 
          VEGPIX **VegMap, LAYER *Veg, SOILPIX **SoilMap, NETSTRUCT **Network, 
          CHANNEL *ChannelData, LAYER *Soil, AGGREGATED *Total);

unsigned char fequal(float a, float b);

void FinalMassBalance(FILES *Out, AGGREGATED *Total, WATERBALANCE *Mass, OPTIONSTRUCT *Options);

float FindDT(SOILPIX **SoilMap, MAPSIZE *Map, TIMESTRUCT *Time, 
             TOPOPIX **TopoMap, SOILTABLE *SType); 

void GenerateScales(MAPSIZE *Map, int NumberType, void **XScale,
		    void **YScale);

void GetMetData(OPTIONSTRUCT *Options, TIMESTRUCT *Time, int NSoilLayers,
		int NStats, float SunMax, METLOCATION *Stat);

uchar InArea(MAPSIZE *Map, COORD *Loc);

void InitAggregated(OPTIONSTRUCT *Options, int MaxVegLayers, int MaxSoilLayers,
  AGGREGATED *Total);

void InitCharArray(char *Array, int Size);

void InitConstants(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map,
	SOLARGEOMETRY *SolarGeo, TIMESTRUCT *Time);  

void InitDump(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map,
	      int MaxSoilLayers, int MaxVegLayers, int Dt,
	      TOPOPIX **TopoMap, DUMPSTRUCT *Dump);

void InitVegUpdate(LISTPTR Input, int NUpdate, DATE ** DUpdate);

void InitEvapMap(MAPSIZE *Map, EVAPPIX ***EvapMap, SOILPIX **SoilMap,
		 LAYER *Soil, VEGPIX **VegMap, LAYER *Veg, TOPOPIX **TopoMap);

void InitImageDump(LISTPTR Input, int Dt, MAPSIZE *Map, int MaxSoilLayers,
		   int MaxVegLayers, char *Path, int NMaps, int NImages, MAPDUMP **DMap);

void InitInFiles(INPUTFILES *InFiles);

void InitInterpolationWeights(MAPSIZE *Map, OPTIONSTRUCT *Options,
			      TOPOPIX **TopoMap, uchar ****MetWeights,
			      METLOCATION *Stats, int NStats);

void InitMapDump(LISTPTR Input, MAPSIZE *Map, int MaxSoilLayers, int MaxVegLayers,
		 char *Path, int NMaps, MAPDUMP **DMap);

void InitMappedConstants(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map,
                         SNOWPIX ***SnowMap, VEGTABLE *VType, VEGPIX ***VegMap);

void InitMassWaste(LISTPTR Input, TIMESTRUCT *Time);

void InitGridMet(OPTIONSTRUCT *Options, LISTPTR Input, MAPSIZE *Map, TOPOPIX **TopoMap,  
         GRID *Grid, METLOCATION **Stat, int *NStats);

void InitMetMaps(LISTPTR Input, int NDaySteps, MAPSIZE *Map, 
                 OPTIONSTRUCT *Options,
                 float ***PrismMap, float ***SnowPatternMap, float ***SnowPatternMapBase,
                 unsigned char ****ShadowMap, float ***SkyViewMap,
                 EVAPPIX ***EvapMap, PRECIPPIX ***PrecipMap, float ***PptMultiplierMap,
                 float ***MeltMultiplierMap, PIXRAD ***RadMap,
                 SOILPIX **SoilMap, LAYER *Soil, VEGPIX **VegMap,
                 LAYER *Veg, TOPOPIX **TopoMap);

void InitMetSources(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map,
            TOPOPIX **TopoMap, int NSoilLayers, TIMESTRUCT *Time, 
            INPUTFILES *InFiles, int *NStats, METLOCATION **Stat);

void InitModelState(DATE *Start, int StepsPerDay, int Dt,
        MAPSIZE *Map, OPTIONSTRUCT *Options,
		    PRECIPPIX **PrecipMap, SNOWPIX **SnowMap,
		    SOILPIX **SoilMap, LAYER Soil, SOILTABLE *SType,
		    VEGPIX **VegMap, LAYER Veg, VEGTABLE *VType, char *Path,
		    TOPOPIX **TopoMap,
		    NETSTRUCT **Network, CHANNEL *ChannelData);

void InitNetwork(int NY, int NX, float DX, float DY, TOPOPIX **TopoMap, 
		 SOILPIX **SoilMap, VEGPIX **VegMap, VEGTABLE *VType, 
		 NETSTRUCT ***Network, CHANNEL *ChannelData, 
		 LAYER Veg, OPTIONSTRUCT *Options);

void InitNewDay(int DayOfYear, SOLARGEOMETRY *SolarGeo);

void InitNewMonth(TIMESTRUCT *Time, OPTIONSTRUCT *Options, MAPSIZE *Map,
		  TOPOPIX **TopoMap, float **PrismMap, float **SnowPatternMap, float **SnowPatternMapBase, unsigned char ***ShadowMap, 
		  INPUTFILES *InFiles, int NVegs, VEGTABLE *VType, int NStats,
		  METLOCATION *Stat, char *Path, VEGPIX ***VegMap, SNOWPIX **SnowMap);

void InitNewStep(INPUTFILES *InFiles, MAPSIZE *Map, TIMESTRUCT *Time,
		 int NSoilLayers, OPTIONSTRUCT *Options, int NStats,
		 METLOCATION *Stat, SOLARGEOMETRY *SolarGeo, 
		 TOPOPIX **TopoMap, SOILPIX **SoilMap);

void InitNewWaterYear(TIMESTRUCT *Time, OPTIONSTRUCT *Options, MAPSIZE *Map,
                TOPOPIX **TopoMap, SNOWPIX **SnowMap, PRECIPPIX **PrecipMap);

void InitParameterMaps(OPTIONSTRUCT *Options, MAPSIZE *Map, int Id,
  char *FileName, SNOWPIX ***SnowMap, int ParamType, float temp);

int InitPixDump(LISTPTR Input, MAPSIZE *Map, uchar **BasinMask, char *Path,
		int NPix, PIXDUMP **Pix, OPTIONSTRUCT *Options);
    
void InitMultiplierMaps(OPTIONSTRUCT *Options, MAPSIZE *Map,
                        float ***PptMultiplierMap, float ***MeltMultiplierMap);

void InitPrismMap(int NY, int NX, float ***PrismMap);

void InitSnowPatternMap(float ***SnowPatternMap, float ***SnowPatternMapBase,
                        MAPSIZE *Map, OPTIONSTRUCT *Options);

void InitShadeMap(OPTIONSTRUCT *Options, int NDaySteps, MAPSIZE *Map,
		  unsigned char ****ShadowMap, float ***SkyViewMap);

void InitPrecipMap(MAPSIZE *Map, PRECIPPIX ***PrecipMap, VEGPIX **VegMap,
		   LAYER *Veg, TOPOPIX **TopoMap);

void InitRadMap(MAPSIZE *Map, PIXRAD ***RadMap);

void InitSatVaporTable(void);

void InitSnowMap(MAPSIZE *Map, SNOWPIX ***SnowMap, TIMESTRUCT *Time);

void InitSoilMap(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map,
		 LAYER *Soil, TOPOPIX **TopoMap, SOILPIX ***SoilMap, SOILTABLE * SType,
		 VEGPIX *** VegMap, VEGTABLE * VType);

int InitSoilTable(OPTIONSTRUCT *Options, SOILTABLE **SType, 
			LISTPTR Input, LAYER *Soil, int InfiltOption);

void InitStateDump(LISTPTR Input, int NStates, DATE **DState);

void InitStations(LISTPTR Input, MAPSIZE *Map, int NDaySteps,
		  OPTIONSTRUCT *Options, int *NStats, METLOCATION **Stat);

void InitTables(int StepsPerDay, LISTPTR Input, OPTIONSTRUCT *Options, 
  MAPSIZE *Map, SOILTABLE **SType, LAYER *Soil, VEGTABLE **VType, LAYER *Veg,
  LAKETABLE **LType);
    
void InitTerrainMaps(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map,
  LAYER *Soil, LAYER *Veg, TOPOPIX ***TopoMap, SOILTABLE *SType, SOILPIX ***SoilMap, 
  VEGTABLE *VType, VEGPIX ***VegMap, DYNAVEG *DVeg, LAKETABLE *LType);

void InitTopoMap(LISTPTR Input, OPTIONSTRUCT * Options, MAPSIZE * Map,
  TOPOPIX *** TopoMap, LAKETABLE *LType);

int InitLakeTable(LAKETABLE **LType, 
                  LISTPTR Input, OPTIONSTRUCT *Options);

void InitVegMap( OPTIONSTRUCT *Options, LISTPTR Input, MAPSIZE *Map, VEGPIX ***VegMap,
				VEGTABLE *VType, DYNAVEG *DVeg);

uchar IsVegDate(DATE *Current, DYNAVEG *DVeg);

int InitVegTable(VEGTABLE **VType, LISTPTR Input, OPTIONSTRUCT *Options, LAYER *Veg);

float evalexpint(int n, float x);

uchar IsStationLocation(COORD *Loc, int NStats, METLOCATION *Station,
			int *WhichStation);

float LapseT(float Temp, float FromElev, float ToElev, float LapseRate);
 
PIXMET MakeLocalMetData(int y, int x, MAPSIZE *Map, int DayStep, int NDaySteps,
			OPTIONSTRUCT *Options, int NStats, METLOCATION *Stat, 
      uchar *MetWeights, float LocalElev, PIXRAD *RadMap,
			PRECIPPIX *PrecipMap, 
			float **PrismMap, float **SnowPatternMap, SNOWPIX *LocalSnow, 
      CanopyGapStruct **Gap, VEGPIX *VegMap,
			float precipMultiplier, int Month, float skyview,
			unsigned char shadow, float SunMax, float SineSolarAltitude);

void MassBalance(DATE *Current, DATE *Start, FILES *Out, AGGREGATED *Total, WATERBALANCE *Mass);

void MassEnergyBalance(OPTIONSTRUCT *Options, int y, int x, float SineSolarAltitude,
            float DX, float DY, int Dt, int HeatFluxOption, int CanopyRadAttOption,
            int InfiltOption, int MaxSoilLayer, int MaxVegLayers, PIXMET *LocalMet,
            NETSTRUCT *LocalNetwork, PRECIPPIX *LocalPrecip,  float SnowMeltMultiplier, VEGTABLE *VType,
            VEGPIX *LocalVeg, SOILTABLE *SType, SOILPIX *LocalSoil,
            SNOWPIX *LocalSnow, PIXRAD *LocalRad, EVAPPIX *LocalEvap, PIXRAD *TotalRad,
            CHANNEL *ChannelData, float **skyview);

double pow (double a, double b);

void quick(ITEM *OrderedCells, int count);

void qs(ITEM *OrderedCells, int left, int right);

void ReadChannelState(char *Path, DATE *Current, Channel *Head);

void ReadMetRecord(OPTIONSTRUCT *Options, DATE *Current, int NSoilLayers,
		   FILES *InFile, MET *MetRecord);

void ReadPRISMMap(DATE *Current, int Dt, char *HDFFileName);

void ResetAggregate(LAYER *Soil, LAYER *Veg, AGGREGATED *Total,
                    OPTIONSTRUCT *Options);

void ResetValues(MAPSIZE *Map, SOILPIX **SoilMap);

int Round(double x);

void RouteSubSurface(int Dt, MAPSIZE *Map, TOPOPIX **TopoMap,
		     VEGTABLE *VType, VEGPIX **VegMap,
		     NETSTRUCT **Network, SOILTABLE *SType,
		     SOILPIX **SoilMap, CHANNEL *ChannelData, 
		     TIMESTRUCT *Time, OPTIONSTRUCT *Options, 
		     char *DumpPath);

void RouteSubSurfaceSpinup(int Dt, MAPSIZE *Map, TOPOPIX **TopoMap,
                           VEGTABLE *VType, VEGPIX **VegMap,
                           NETSTRUCT **Network, SOILTABLE *SType,
                           SOILPIX **SoilMap, OPTIONSTRUCT *Options,
                           float **SubFlowGrad, unsigned char ***SubDir, unsigned int **SubTotalDir);

void RouteSurface(MAPSIZE * Map, TIMESTRUCT * Time, TOPOPIX ** TopoMap,
  SOILPIX ** SoilMap, OPTIONSTRUCT *Options,
  DUMPSTRUCT *Dump, VEGPIX ** VegMap, VEGTABLE *VType, LAKETABLE *LType, SOILTABLE *SType, CHANNEL *ChannelData,
  float Tair, float Rh);

float SatVaporPressure(float Temperature);

int ScanInts(FILE *FilePtr, int *X, int N);

int ScanDoubles(FILE *FilePtr, double *X, int N);

int ScanFloats(FILE *FilePtr, float *X, int N);

uchar ScanUChars(FILE *FilePtr, uchar *X, int N);

void SkipHeader(FILES *InFile, int NLines);

void SkipLines(FILES *InFile, int NLines);

void StoreChannelState(char *Path, DATE *Current, Channel *Head);

void StoreChannelStateExtra(char *Path, DATE *Current, Channel *Head);

void StoreModelState(char *Path, DATE *Current, MAPSIZE *Map,
		     OPTIONSTRUCT *Options, TOPOPIX **TopoMap, PRECIPPIX **PrecipMap, 
             SNOWPIX **SnowMap, VEGPIX **VegMap, 
             LAYER *Veg, SOILPIX **SoilMap, LAYER *Soil, NETSTRUCT **Network, 
		     CHANNEL *ChannelData);

void SnowStats(DATE *Now, MAPSIZE *Map, OPTIONSTRUCT *Options, 
        TOPOPIX **TopoMap, SNOWPIX **Snow, int Dt);

void UpdateVegMap(DATE *Current, OPTIONSTRUCT * Options, MAPSIZE * Map,
                LAYER *Veg, VEGPIX *** VegMap, VEGTABLE *VType, DYNAVEG *DVeg);

float viscosity(float Tair, float Rh);

/* functions used in canopy gapping */
void AggregateCanopyGap(CanopyGapStruct **Gap, VEGPIX *LocalVeg,
  SOILPIX *LocalSoil, SNOWPIX *LocalSnow, EVAPPIX *LocalEvap,
  PRECIPPIX *LocalPrecip, PIXRAD *LocalRad, double weight, int NSoil, int NVeg, int NVegLayer);

float AreaIntegral(float Extn, float Lmax, float SolarAltitude, float R,
  float xmax, float xmin, float Rsb, float Rsd, float Albedo);

void CalcCanopyGapAerodynamic(CanopyGapStruct **Gap, int NVegLayers,
  float *Height);

void CalcCanopyGapET(CanopyGapStruct **Gap, int MaxSoilLayer, VEGTABLE *VType,
  VEGPIX *LocalVeg, SOILTABLE *SType, SOILPIX *LocalSoil, PIXMET *LocalMet,
  EVAPPIX *LocalEvap, NETSTRUCT *LocalNetwork, int Dt, float UpperRa,
  float LowerRa, float DX, float DY, int x, int y, CHANNEL *ChannelData);

void CalcGapSurroudingET(int Dt, CanopyGapStruct **Gap,
  SOILTABLE *SType, VEGTABLE *VType, PIXRAD *LocalRad, PIXMET *LocalMet,
  SOILPIX *LocalSoil, NETSTRUCT *LocalNetwork, float UpperRa, float LowerRa,
  VEGPIX *LocalVeg, float DX, float DY, int x, int y, CHANNEL *ChannelData);

void CanopyGapInterception(OPTIONSTRUCT *Options, CanopyGapStruct **Gap,
  int HeatFluxOption, int y, int x, int Dt, int NVegLActual,
  float DX, float DY, float UpperRa, float UpperWind, VEGTABLE *VType,
  SOILPIX *LocalSoil, VEGPIX *LocalVeg, SNOWPIX *LocalSnow,
  PRECIPPIX *LocalPrecip, PIXRAD *LocalRad, PIXMET *LocalMet);

void CanopyGapInterceptionStorage(int NAct, float *MaxInt, float *Fract,
  float *Int, float *Precip);

void CanopyGapRadiation(CanopyGapStruct **Gap, float SunAngle, float Rs,
  float Rsb, float Rsd, float Ld, float TSurf, float Tcanopy, float SoilAlbedo,
  VEGTABLE *VType, SNOWPIX *LocalSnow, PIXRAD *LocalRad, float Gapping, VEGPIX *LocalVeg);

float CanopyGapShortRadiation(int Understory, float GapView, float h, float dm,
  float SunAngle, float Rsb, float Rsd, float Extn, float SoilAlbedo, VEGTABLE *VType,
  SNOWPIX *LocalSnow, float Vf);

void CanopyGapLongRadiation(CanopyGapStruct *Gap, float h, float dm, float Ld,
  float Tsurf, float Vf);

void CanopyGapSnowMelt(OPTIONSTRUCT *Options, int y, int x, int Dt,
  CanopyGapStruct **Gap, float DX, float DY, VEGTABLE *VType, VEGPIX *LocalVeg,
  SNOWPIX *LocalSnow, PRECIPPIX *LocalPrecip, PIXRAD *LocalRad, PIXMET *LocalMet);

void GapSurroundingLongRadiation(CanopyGapStruct *Forest, float Ld, float Vf, float F,
  float Tcanopy, float Tsurf);

void GapSurroundingShortRadiation(CanopyGapStruct *Forest, VEGTABLE *VType,
  SNOWPIX *LocalSnow, float SoilAlbedo, float SineSolarAltitude, float Rs, VEGPIX *LocalVeg);

void InitCanopyGapMap(OPTIONSTRUCT *Options, LISTPTR Input, MAPSIZE *Map,
  LAYER *Soil, LAYER *Veg, VEGTABLE *VType, VEGPIX ***VegMap, SOILTABLE *SType,
  SOILPIX ***SoilMap);

float NonGapShortRadiation(float Rs, float SunAngle, float SoilAlbedo,
  CanopyGapStruct *Forest, VEGTABLE *VType, SNOWPIX *LocalSnow);

void CalcGapSurroudingIntercept(OPTIONSTRUCT *Options, int HeatFluxOption,
  int y, int x, int Dt, int NVegLActual, CanopyGapStruct **Gap, VEGTABLE *VType,
  PIXRAD *LocalRad, PIXMET *LocalMet, float UpperRa, float UpperWind, VEGPIX *LocalVeg);

float CalcGapView(float R, float H, float Vf);

#endif
