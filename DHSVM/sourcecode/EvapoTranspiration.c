
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "massenergy.h"
#include "constants.h"
#include "functions.h"

/*****************************************************************************
EvapoTranspiration()
*****************************************************************************/
void EvapoTranspiration(int Layer, int ImpvRad, int Dt, PIXMET *Met,
  float NetRad, float Rp, VEGTABLE *VType, SOILTABLE *SType,
  float MoistureFlux, float *Moist, float *SoilTemp, float *Int,
  float *EPot, float *EInt, float **ESoil, float *EAct, float *ETot,
  float *Adjust, float Ra, VEGPIX *LocalVeg, SOILPIX *LocalSoil,
  OPTIONSTRUCT *Options)
{
  float *Rc;			/* canopy resistance associated with
                        conditions in each soil layer (s/m) */
  float *RootWt;      /* fraction of transpiration drawn from each soil layer */
  float RcPlant;      /* plant-level canopy resistance, Options only (s/m) */
  float DryEvapTime;	/* amount of time remaining during a timestep
                        after the interception storage is depleted (s) */
  float F;			    /* Fractional coverage by vegetation layer */
  float SoilMoisture, PlantAvailableMoist;	/* Amount of water in each soil layer (m) */
  float WetArea;		/* relative leaf area wetted by interception storage */
  float WetEvapRate;	/* evaporation rate from wetted fraction per unit ground area (m/s) */
  float WetEvapTime;	/* amount of time needed to evaporate the amount of water
                        in interception storage (s) */
  CANOPYHYD Hyd;
  int i;


  F = LocalVeg->Fract[Layer];
  NetRad /= F;
  Rp /= F;

  /* Convert the water amounts related to partial canopy cover to a pixel depth
  as if the entire pixel is covered. These depths will be converted back later on. */
  *Int /= F;
  MoistureFlux /= F;
  LocalVeg->MaxInt[Layer] /= F;

  /* allocate memory for the canopy resistance array */
  if (!(Rc = (float *)calloc(VType->NSoilLayers, sizeof(float))))
    ReportError("EvapoTranspiration()", 1);
  if (!(RootWt = (float *)calloc(VType->NSoilLayers, sizeof(float))))
    ReportError("EvapoTranspiration()", 1);

  /* Calculate the evaporation rate in m/s */
  EPot[Layer] = (Met->Slope * NetRad + Met->AirDens * CP * Met->Vpd / Ra) /
    (WATER_DENSITY * Met->Lv * (Met->Slope + Met->Gamma));

  /* The potential evaporation rate accounts for the amount of moisture that
  the atmosphere can absorb. If we do not account for the amount of
  evaporation from overlying evaporation, we can end up with the situation
  that all vegetation layers and the soil layer transpire/evaporate at the
  potential rate, resulting in an overprediction of the actual evaporation
  rate. Thus we subtract the amount of evaporation that has already
  been calculated for overlying layers from the potential evaporation.
  Another mechanism that could be used to account for this would be to
  decrease the vapor pressure deficit while going down through the canopy
  (not implemented here) */
  EPot[Layer] -= MoistureFlux / Dt;

  if (EPot[Layer] < 0)
    EPot[Layer] = 0;

  /* WetArea = pow(*Int/VType->MaxInt[Layer], (double) 2.0/3.0); */

  WetArea = cbrt(*Int / LocalVeg->MaxInt[Layer]);
  WetArea = MIN(WetArea, 1);
  WetArea = WetArea * WetArea;
  
  /* calculate the amount of water that can evaporate from the interception
  storage.  Given this evaporation rate, calculate the amount of time it
  will take to evaporate the entire amount of intercepted water.  If this
  time period is shorter than the length of the time interval, the
  previously wetted leaves can transpire during the remaining part of the
  time interval */

  /* WORK IN PROGRESS:  the amount of interception storage can be replenished
  during (low intensity) rainstorms */

  WetEvapRate = WetArea * EPot[Layer];
  if (WetEvapRate > 0) {
    WetEvapTime = *Int / WetEvapRate;
    if (WetEvapTime > Dt) {
      WetEvapTime = Dt;
    }
  }
  else if (*Int > 0)
    WetEvapTime = Dt;
  else
    WetEvapTime = 0;

  if (WetEvapRate > 0) {
    if (WetEvapTime < Dt) {
      EInt[Layer] = *Int;
      *Int = 0.0;
      DryEvapTime = Dt - WetEvapTime;
    }
    else {
      EInt[Layer] = Dt * WetEvapRate;
      *Int -= EInt[Layer];
      DryEvapTime = 0.0;
    }
  }
  else {
    EInt[Layer] = 0.0;
    if (*Int > 0)
      DryEvapTime = 0.0;
    else
      DryEvapTime = Dt;
  }

  /* Correct the evaporation from interception and the interception storage for
  the fractional overstory coverage */
  EInt[Layer] *= F;
  *ETot += EInt[Layer];
  *Int *= F;
  LocalVeg->MaxInt[Layer] *= F;

  /* Canopy resistance and the distribution of root uptake among soil layers.

     JARVIS: resistance is evaluated separately for each layer's moisture and
     uptake follows the prescribed root fractions (Wigmosta et al. 1994).

     Everything else goes through the single seam in CanopyResistance.c: one
     plant-level resistance, and per-layer uptake weights that come either
     from the empirical beta (HYD_NONE) or from the Krs/SUF macroscopic root
     model (HYD_KRSSUF; Couvreur 2012, Vanderborght 2021). */
  /* Root-zone soil water potential, computed for EVERY scheme.  It is a
     property of the soil state, so it is the one variable on which Jarvis,
     the empirical-beta paths and the hydraulic schemes can be compared
     directly.  Costs one Brooks-Corey evaluation per soil layer. */
  LocalVeg->PsiSoil[Layer] =
    CanopyRootZonePotential(Layer, VType, SType, LocalSoil, Moist);

  if (Options->StomScheme == JARVIS) {
    for (i = 0; i < VType->NSoilLayers; i++) {
      Rc[i] = CanopyResistance(LocalVeg->LAI[Layer], VType->RsMin[Layer],
                               VType->RsMax[Layer], VType->Rpc[Layer],
                               VType->VpdThres[Layer], VType->MoistThres[Layer],
                               SType->WP[i], SoilTemp[i], Moist[i],
                               Met->Vpd, Rp);
      RootWt[i] = VType->RootFract[Layer][i];
    }
  }
  else {
    RcPlant = CanopyResistanceScheme(Options->StomScheme, Options->Hydraulics,
                LocalVeg->LAI[Layer], Layer, Rp, NetRad, Met, VType, SType,
                LocalSoil, Moist, LocalVeg->PhotoDormancy, &Hyd);

    /* ProfitMax, ProfitMax2 and SOX solve transpiration themselves, with a
       leaf energy balance and the leaf-to-air VPD; the number they return is
       also the one their leaf water potentials, PLC and root uptake were
       computed from.  Penman-Monteith below must therefore return that
       transpiration, not re-derive one from Rc and the AIR VPD -- which is
       what produced near-potential rates at dawn and dusk whenever the leaf
       VPD was small.  Invert PM for the resistance that reproduces Escheme
       (docs/provenance.md departure #1, now implemented).  The Medlyn paths
       keep conductance as their state variable and leave Rc as computed. */
    if (Hyd.EConsistent) {
      float RcMin = (LocalVeg->LAI[Layer] > 0.0)
                    ? VType->RsMin[Layer] / LocalVeg->LAI[Layer]
                    : VType->RsMax[Layer];
      RcPlant = CanopyResistanceFromFlux(Hyd.Escheme, EPot[Layer], Met->Slope,
                                         Met->Gamma, Ra, RcMin,
                                         VType->RsMax[Layer]);
      Hyd.Rc = RcPlant;
    }

    for (i = 0; i < VType->NSoilLayers; i++) {
      Rc[i] = RcPlant;
      RootWt[i] = Hyd.Efrac[i];
    }

    /* Per canopy layer, so the overstory no longer overwrites the
       understory (or the other way round). */
    LocalVeg->PsiLeaf[Layer]       = Hyd.PsiLeaf;
    LocalVeg->PsiRoot[Layer]       = Hyd.PsiRoot;
    LocalVeg->PsiSoil[Layer]       = Hyd.Heff;
    LocalVeg->PLC[Layer]           = Hyd.PLC;
    LocalVeg->SafetyMargin[Layer]  = Hyd.SafetyMargin;
    LocalVeg->HydStress[Layer]     = Hyd.Beta;
    LocalVeg->Tleaf[Layer]         = Hyd.Tleaf;
    LocalVeg->VpdLeaf[Layer]       = Hyd.VpdLeaf;
    LocalVeg->AnCanopy[Layer]      = Hyd.An;
    LocalVeg->Escheme[Layer]       = Hyd.Escheme;
    LocalVeg->Ecrit[Layer]         = Hyd.Ecrit;
    LocalVeg->Tsupply[Layer]       = Hyd.Tsupply;
    LocalVeg->SupplyLimited[Layer] = (Hyd.Tsupply >= 0.0f);
  }

  LocalVeg->Rc[Layer] = Rc[0]; /* Plant-level under Medlyn/Sperry, layer 0 under Jarvis */

  /* Calculate the transpiration rate for the current vegetation layer,
  and adjust the soil moisture content in each of the soil layers */
  for (i = 0; i < VType->NSoilLayers; i++) {
    ESoil[Layer][i] = ((Met->Slope + Met->Gamma) /
                      (Met->Slope + Met->Gamma * (1 + Rc[i] / Ra))) *
                        RootWt[i] * EPot[Layer] * Adjust[i];

    /* Calculate the amounts of water transpired during each timestep based
    on the evaporation and transpiration rates.  While there is still water
    in interception storage only the area that is not covered by intercep-
    ted water will transpire.  When all of the interception storage has
    disappeared all leaves will contribute to the transpiration */
    ESoil[Layer][i] *= WetEvapTime * (1 - WetArea) + DryEvapTime;

    SoilMoisture = Moist[i] * VType->RootDepth[i] * Adjust[i];
    PlantAvailableMoist = (Moist[i] - SType->WP[i]) * VType->RootDepth[i] * Adjust[i];
    if (PlantAvailableMoist < 0.0)
      PlantAvailableMoist = 0.0;
    if (PlantAvailableMoist < ESoil[Layer][i])
      ESoil[Layer][i] = PlantAvailableMoist;

    /* Correct the evaporation for the fractional overstory coverage and update
    the soil moisture */
    ESoil[Layer][i] *= F;
    SoilMoisture -= ESoil[Layer][i];

    Moist[i] = SoilMoisture / (VType->RootDepth[i] * Adjust[i]);

  }

  /* Supply limit.  When the root Dirichlet switch fires (Leitner et al. 2025
     sec. 2.3) the soil cannot deliver what the atmosphere is demanding, and
     Penman-Monteith has no way to know that.  Scale the layer fluxes down to
     what the root system can actually supply.  Without this the model
     transpires water the hydraulics have already declared unavailable.
     Tsupply < 0 means unlimited.  Runs AFTER the plant-available-water clamp
     so the tighter of the two limits wins. */
  if (Options->StomScheme != JARVIS && Hyd.Tsupply >= 0.0f) {
    float EDemand = 0.0f, ESupply, Scale;
    for (i = 0; i < VType->NSoilLayers; i++) EDemand += ESoil[Layer][i];
    ESupply = Hyd.Tsupply * (float)Dt * F;
    if (EDemand > ESupply && EDemand > 0.0f) {
      Scale = ESupply / EDemand;
      for (i = 0; i < VType->NSoilLayers; i++) {
        /* Give the water back to the soil that the clamp already removed. */
        Moist[i] += (ESoil[Layer][i] * (1.0f - Scale)) /
                    (VType->RootDepth[i] * Adjust[i]);
        ESoil[Layer][i] *= Scale;
      }
    }
  }

  for (i = 0, EAct[Layer] = 0; i < VType->NSoilLayers; i++) {
    *ETot += ESoil[Layer][i];
    EAct[Layer] += ESoil[Layer][i];
  }
  
  free(Rc);
  free(RootWt);   /* Phase V: was leaking since the Medlyn path went in */
}

