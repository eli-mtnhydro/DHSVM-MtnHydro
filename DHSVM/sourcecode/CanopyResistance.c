
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "massenergy.h"
#include "constants.h"
#include "photosynthesis.h"

/*****************************************************************************
  CanopyResistance()
*****************************************************************************/
float CanopyResistance(float LAI, float RsMin, float RsMax, float Rpc,
  float VpdThres, float MoistThres, float WP,
  float TSoil, float SoilMoisture, float Vpd, float Rp)
{
  float MoistFactor;	/* multiplier for resistance due to soil moisture feed-back */
  float Resistance;		/* Canopy resistance (s/m) */
  float RpFactor;		/* multiplier for resistance due to light level feed-back */
  float TFactor;		/* multiplier for resistance due to soil temperaure feed-back */
  float VpdFactor;		/* multiplier for resistance due to vapor pressure deficit feed-back */

  if (TSoil <= 0) {
    Resistance = DHSVM_HUGE;
    return Resistance;
  }

  /* for OBS */
  TFactor = 1.0 / (0.176 + 0.0770 * TSoil - 0.0018 * TSoil * TSoil);

  if (TFactor <= 0) {
    Resistance = DHSVM_HUGE;
    return Resistance;
  }

  /* equation 14, Wigmosta et al [1994] */

  if (Vpd >= VpdThres) {
    Resistance = DHSVM_HUGE;
    return Resistance;
  }
  else
    VpdFactor = 1.0 / (1 - Vpd / VpdThres);

  /* equation 15, Wigmosta et al [1994 */

  RpFactor = 1.0 / ((RsMin / RsMax + Rp / Rpc) / (1 + Rp / Rpc));

  /* equation 16, Wigmosta et al [1994] */

  if (SoilMoisture <= WP) {
    Resistance = DHSVM_HUGE;
    return Resistance;
  }
  else if (SoilMoisture < MoistThres)
    MoistFactor = (MoistThres - WP) / (SoilMoisture - WP);
  else
    MoistFactor = 1.0;

  Resistance = TFactor * VpdFactor * RpFactor * MoistFactor * RsMin / LAI;

  return Resistance;

}

/*****************************************************************************
 CanopyResistancePhoto()
 *****************************************************************************/
float CanopyResistancePhoto(float Lai, float Vcmax25, float G1, float G0,
  float Dormancy, float Beta, int Layer, float RsMax, float Rp, PIXMET *Met)
{
  float FBeam, ParBeam, ParDiff, GsCanopy;

  if (Lai <= 0.0 || Vcmax25 <= 0.0)
    return DHSVM_HUGE;

  /* Layer 0 is always the topmost canopy layer and sees the direct beam.
     Layer 1 sits under an overstory, where transmitted light is nearly all
     diffuse -- using the above-canopy beam fraction there would badly
     overestimate the sunlit fraction. */
  if (Layer == 0 && (Met->SinBeam + Met->SinDiffuse) > 0.0)
    FBeam = Met->SinBeam / (Met->SinBeam + Met->SinDiffuse);
  else
    FBeam = 0.0;

  ParBeam = Rp * FBeam;
  ParDiff = Rp * (1.0 - FBeam);

  GsCanopy = PhotoCanopyConductance(Vcmax25, G1, G0, Lai,
    Met->SineSolarAltitude, ParBeam, ParDiff, Met->Tair,
    Met->Vpd * PHOTO_PA_TO_KPA, PHOTO_CA, Met->Press, PHOTO_GB, Beta, Dormancy, NULL);

  if (GsCanopy <= 1.0 / RsMax)
    return RsMax;

  return 1.0 / GsCanopy;
}
