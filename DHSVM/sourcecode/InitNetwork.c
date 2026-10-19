
#include <stdio.h>
#include <stdlib.h>
#include "constants.h"
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "settings.h"
#include "soilmoisture.h"
#include "DHSVMChannel.h"

 /*****************************************************************************
   Function name: InitNetwork()

   Purpose      : Initialize channel.  Memory is allocated, and the
                  necessary adjustments for the soil profile are calculated

 *****************************************************************************/
void InitNetwork(int NY, int NX, float DX, float DY, TOPOPIX **TopoMap,
  SOILPIX **SoilMap, VEGPIX **VegMap, int SoilLayers,
  NETSTRUCT ***Network, CHANNEL *ChannelData)
{
  const char *Routine = "InitNetwork";
  int x;			/* column counter */
  int y;			/* row counter */
  /* Allocate memory for network structure */

  if (!(*Network = (NETSTRUCT **)calloc(NY, sizeof(NETSTRUCT *))))
    ReportError((char *)Routine, 1);

  for (y = 0; y < NY; y++) {
    if (!((*Network)[y] = (NETSTRUCT *)calloc(NX, sizeof(NETSTRUCT))))
      ReportError((char *)Routine, 1);
  }

  for (y = 0; y < NY; y++) {
    for (x = 0; x < NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
        if (!((*Network)[y][x].Adjust =
          (float *)calloc((SoilLayers + 1), sizeof(float))))
          ReportError((char *)Routine, 1);

        if (!((*Network)[y][x].PercArea =
          (float *)calloc((SoilLayers + 1), sizeof(float))))
          ReportError((char *)Routine, 1);
      }
    }
  }
  
  /* Read the channel Network information, and calculate the storage
     adjustment factors */
  for (y = 0; y < NY; y++) {
    for (x = 0; x < NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
        ChannelCut(y, x, ChannelData, &((*Network)[y][x]));

        AdjustStorage(SoilLayers,
          SoilMap[y][x].Depth,
          SoilMap[y][x].LayerThickness,
          (*Network)[y][x].Area, DX, DY,
          (*Network)[y][x].BankHeight,
          (*Network)[y][x].PercArea,
          (*Network)[y][x].Adjust,
          &((*Network)[y][x].CutBankZone));
      }
    }
  }
}
