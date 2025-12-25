
#ifndef SOILMOISTURE_H
#define SOILMOISTURE_H

#define NO_CUT -10

void AdjustStorage(int NSoilLayers, float TotalDepth, float *RootDepth,
		   float Area, float DX, float DY, float BankHeight, 
		   float *PercArea, float *Adjust, int *CutBankZone);

float CalcAvailableWater(int NRootLayers, float TotalDepth, float *RootDepth,
			 float *Moisture, float *FCap, float *Moist, float TableDepth,
			 float *Adjust);

float CalcTotalWater(int NSoilLayers, float TotalDepth, float *RootDepth,
		     float *Moist, float *Adjust);

void CutBankGeometry(int i, float RootDepth, float TopZone, float BankHeight,
		     float Area, float DX, float DY, float *PercArea, 
		     float *Adjust, int *CutBankZone);

void DistributeSatflow(int Dt, float DX, float DY, float SatFlow,
                       int NSoilLayers, float TotalDepth, float *RootDepth,
                       float *Porosity, float *FCap, float *Adjust,
                       float *TableDepth, float *Runoff, float *Moist);

void UnsaturatedFlow(int Dt, float DX, float DY, float Infiltration, 
		     float SatFlow, int NSoilLayers, 
		     float TotalDepth, float Area, float *RootDepth, float *Ks, 
		     float *PoreDist, float *Porosity, float *FCap, float *Perc, 
		     float *PercArea, float *Adjust, int CutBankZone, float BankHeight,
			   float *TableDepth, float *Runoff, float *Moist, int InfiltOption);

float WaterTableDepth(int NRootLayers, float TotalDepth, float *RootDepth,
		      float *Porosity, float *FCap, float *Adjust,
		      float *Moist);

#endif
