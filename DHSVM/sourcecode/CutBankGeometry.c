
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "settings.h"
#include "soilmoisture.h"

/*****************************************************************************
  Function name: CutBankGeometry()

  Purpose      : This subroutine calculates corrections to adjust for the
                 effects of channels in grid cells.  It
                 also add precip and updates the upper zone soil
                 moisture.  If the water table is below the
                 channel bed precip is added to the coresponding zone.
	
  Required     :
    int i            - Number of the soil layer being processed.
    float RootDepth  - Depth of the soil layer being processed (m)
    float TopZone    - Distance from the ground surface to top of zone i (m) 
    float BankHeight - Distance from ground surface to channel bed (m) 
    float Area       - Area of channel surface (m)
    float DX         - Grid cell width (m)
    float DY         - Grid cell width (m)

  Modifies     :
    float *PercArea  - Area of the bottom of zone i, for perc (m/m)
    float *Adjust    - Corrects for loss of soil storage due to channel
                       Multiplied with RootDepth to give the zone
                       thickness for use in calculating soil moisture
    int *CutBankZone - Number of the soil layer containing the bottom of the
                       cut-bank.  Used in UnsaturatedFlow to check for
                       surface runoff.  If BankHeight = 0.; CutBankZone =
                       NO_CUT 

    Schematic: 
    
    <-----------------------------------DX---------------------------------->
    |====================|     -                     |======================|
    |      ^             |     |                     |   |                  |
    |      TopZone[0]    |<----|------Area---------->|   |- RootDepth[0]    |
    |                    |     |                     |   |                  |
    |<---PercArea*DX---->|     |                     |<--|--PercArea*DX---->|
    |                    |     |-BankHeight          |   V                  |
    |====================|     |                     |======================|
    |       ^            |     |                     |   |                  |
    |       TopZone[1]   |     |                     |   |-RootDepth[1]     |
    |                    |     V                     |   |                  |
    |                    |---------------------------|   |                  |
    |                                                    |                  |
    |                             CutBankZone            V                  |
    |=======================================================================|
    |                                                                       |
    |                                                                       |
    |=======================================================================|    

*****************************************************************************/
void CutBankGeometry(int i, float RootDepth, float TopZone, float BankHeight,
		     float Area, float DX, float DY, float *PercArea, 
		     float *Adjust, int *CutBankZone)
{
  *PercArea = 1.0;
  *Adjust = 1.0;

  if (BankHeight > 0.0) {
    if (BankHeight <= TopZone) {
      /* below cut depth - full area */
      *PercArea = 1.0;
      *Adjust = 1.0;
    }
    else {
      if (BankHeight <= (TopZone + RootDepth)) {
	/* cut depth in this zone partial area */
	*PercArea = 1.0;
	*Adjust = 1.0 - (Area * (BankHeight - TopZone) / (RootDepth * DX * DY));
	*CutBankZone = i;
      }
      else {
	/* above cut depth - less than full area  */
	if (DX*DY < Area)
	  Area = DX*DY;
	*PercArea = 1 - Area / (DX * DY);
	*Adjust = *PercArea;
      }
    }
  }
}
