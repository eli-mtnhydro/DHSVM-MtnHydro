/*
 * SUMMARY:      UpdateVegMaps() - Update Vegetation Maps at selected date
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Zhuoran Duan
 * ORG:          Pacific Northwest National Laboratory, Hydrology Group
 * E-MAIL:       zhuoran.duan@pnnl.gov
 * ORIG-DATE:    Mar-2020
 * DESCRIPTION:  Update vegetation map for type, fractional cover, LAI, and height
 * DESCRIP-END.
 * FUNCTIONS:    UpdateVegMaps()
 *               InitTopoMap()
 *               InitSoilMap()
 *               InitVegMap()
 * COMMENTS:
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "fileio.h"
#include "functions.h"
#include "constants.h"
#include "getinit.h"
#include "sizeofnt.h"
#include "varid.h"
#include "Calendar.h"


/*******************************************************************************
  Function name: InitVegUpdate()

  Purpose      : Initialize the vegetation updates.  This information is in the
         [VEGETATION] section of the input file

  Required     :
    LISTPTR Input         - Linked list with input strings
    int NUpdate           - Number of vegetation layers to update
    DATE **DUpdate        - Array with update dates

  Returns      : void

  Modifies     : DUpdates and its members

  Comments     :
*****************************************************************************/
void InitVegUpdate(LISTPTR Input, int NUpdate, DATE ** DUpdate)
{
  char *Routine = "InitVegUpdate";
  int i;			/* counter */
  char KeyName[BUFSIZE + 1];
  char *KeyStr = "UPDATE DATE";
  char *SectionName = "VEGETATION";
  char VarStr[BUFSIZE + 1];

  if (!(*DUpdate = (DATE *)calloc(NUpdate, sizeof(DATE))))
    ReportError(Routine, 1);

  for (i = 0; i < NUpdate; i++) {  
    sprintf(KeyName, "%s %d", KeyStr, i + 1);
    GetInitString(SectionName, KeyName, "", VarStr,
      (unsigned long)BUFSIZE, Input);
    if (!SScanDate(VarStr, &((*DUpdate)[i])))
      ReportError(KeyName, 51);
  }
}

/*****************************************************************************
  IsVegDate()
*****************************************************************************/
uchar IsVegDate(DATE *Current, DYNAVEG *DVeg){
  //char *Routine = "IsVegDate";
  int i;			/* counter */

  for (i = 0; i < DVeg->NUpdate; i++) {
    if (IsEqualTime(Current, &(DVeg->DUpdate[i])))
    {
      return TRUE;
      break;
    }
  } 
  return FALSE;
}

/*****************************************************************************
  Function name: UpdateVegMap()

  Purpose      : Update Vegetation Maps at user defined date, this updates:
                  VEGPIX *** VegMap
  Required     :
    VEGPIX *** VegMap       - Updates vegetation information
    DYNAVEG *DVeg           - Dynamic Vegetaion, with input path and dates 

  Returns      : void

  Modifies     : (see list of required above)

  Comments     :
*****************************************************************************/
void UpdateVegMap(DATE *Current, OPTIONSTRUCT * Options, LISTPTR Input, MAPSIZE * Map,
                LAYER *Veg, VEGPIX *** VegMap, VEGTABLE *VType, DYNAVEG *DVeg)
{  
  const char *Routine = "UpdateVegMap";
  char VarName[BUFSIZE + 1];
  char Path[BUFSIZE + 1];
  char FileName[NAMESIZE + 1];
  char Str[NAMESIZE + 1];
  int i;			/* counter */
  int x;			/* counter */
  int y;			/* counter */
  int flag;
  int NumberType;		/* number type */
  unsigned char *Type;		/* Vegetation type */
  float *FC = NULL;		/* Vegetation Fractional Coverage */
  float *LAIMonthly= NULL; /* Vegetation Leaf Area Index, monthly */
  float *Height = NULL; /*Vegetaion Height*/
  float *Gap;		/* gap diameter */

  int NSet; /*Counter for LAI map month*/

  
  strcpy(Path, DVeg->DynaVegPath);

  printf("Updating vegetation maps...\n");

  /*Vegetation Maps must be named as Veg.XXX.MM.DD.YYYY.HH.MM.SS.ext*/
  sprintf(Str, "%02d.%02d.%02d.%02d.%02d.%02d", Current->Month, Current->Day,
  Current->Year, Current->Hour, Current->Min, Current->Sec);
    
  /*Update Vegetation Type Map*/
  
  sprintf(FileName, "%sVegetation.Type.%s%s", Path, Str, fileext);

  if (access(FileName, F_OK) == 0){
    /* Read the vegetation type */
    GetVarName(005, 0, VarName);
    GetVarNumberType(005, &NumberType);
    if (!(Type = (unsigned char *)calloc(Map->NX * Map->NY,
      SizeOfNumberType(NumberType))))
      ReportError((char *)Routine, 1);
    flag = Read2DMatrix(FileName, Type, NumberType, Map, 0, VarName, 0);
    
    if ((Options->FileFormat == NETCDF && flag == 0)
      || (Options->FileFormat == BIN))
    {
      for (y = 0, i = 0; y < Map->NY; y++) {
        for (x = 0; x < Map->NX; x++, i++) {
          (*VegMap)[y][x].Veg = Type[i];
        }
      }
    }
    else if (Options->FileFormat == NETCDF && flag == 1) {
      for (y = Map->NY - 1, i = 0; y >= 0; y--) {
        for (x = 0; x < Map->NX; x++, i++) {
          (*VegMap)[y][x].Veg = Type[i];
        }
      }
    }
    free(Type);
  }
  else printf("Spatial TYPE map NOT provided, vegetation type NOT updated\n");


  /*Update: after update of vegetation type, associated spatial vegetation properties 
  also needs to be updated, for example: FC, Height and LAI. Note on 03/30/2020*/

  /*Update Vegetation Fractional Cover Map*/
  sprintf(FileName, "%sVegetation.FC.%s%s", Path, Str, fileext);
  /* Read the vegetation fractional coverage map */
  GetVarName(010, 0, VarName);
  GetVarNumberType(010, &NumberType);

  if (access(FileName, F_OK) == 0) {
    printf("Spatial FC map provided, updating...\n");
    if (!(FC = (float *)calloc(Map->NX * Map->NY,
      SizeOfNumberType(NumberType))))
      ReportError((char *)Routine, 1);
    flag = Read2DMatrix(FileName, FC, NumberType, Map, 0, VarName, 0);

    if ((Options->FileFormat == NETCDF && flag == 0)
      || (Options->FileFormat == BIN))
    {
      for (y = 0, i = 0; y < Map->NY; y++) {
        for (x = 0; x < Map->NX; x++, i++) {
          if ( VType[(*VegMap)[y][x].Veg - 1].OverStory == TRUE) {
            if (FC[i] > 0.0)
              (*VegMap)[y][x].Fract[0] = FC[i];
            else
              (*VegMap)[y][x].Fract[0] = VType[(*VegMap)[y][x].Veg - 1].Fract[0];
           
            /*If understory exists, set default understory FC=1.0*/
            if (VType[(*VegMap)[y][x].Veg - 1].UnderStory == TRUE)
              (*VegMap)[y][x].Fract[1] = 1.0;
          }
          else{
            if (VType[(*VegMap)[y][x].Veg - 1].UnderStory == TRUE)
              (*VegMap)[y][x].Fract[0] = 1.0;
          }
        }
      }
    }
    else if (Options->FileFormat == NETCDF && flag == 1) {
      for (y = Map->NY - 1, i = 0; y >= 0; y--) {
        for (x = 0; x < Map->NX; x++, i++) {
          /*Allocate memory*/
          if ( VType[(*VegMap)[y][x].Veg - 1].OverStory == TRUE) {  
            if (FC[i] > 0.0)
              (*VegMap)[y][x].Fract[0] = FC[i];
            else
              (*VegMap)[y][x].Fract[0] = VType[(*VegMap)[y][x].Veg - 1].Fract[0];
           
            /*If understory exists, set default understory FC=1.0*/
            if (VType[(*VegMap)[y][x].Veg - 1].UnderStory == TRUE)
              (*VegMap)[y][x].Fract[1] = 1.0;
          }
          else{
            if ( VType[(*VegMap)[y][x].Veg - 1].UnderStory == TRUE)
              (*VegMap)[y][x].Fract[0] = 1.0;	   
          }
        }
      }
    }
    else ReportError((char *)Routine, 57);
    free(FC);
  }
  else printf("Vegetation FC was not updated!\n");
    
  /*Calculate Vf */
  if (Options->ImprovRadiation == TRUE) {
    for (y = 0, i = 0; y < Map->NY; y++) {
        for (x = 0; x < Map->NX; x++, i++) {
          if ( VType[(*VegMap)[y][x].Veg - 1].NVegLayers >0) 
            (*VegMap)[y][x].Vf = (*VegMap)[y][x].Fract[0] * VType[(*VegMap)[y][x].Veg - 1].VfAdjust;
        }
    }
  }

  /*Update Vegetation LAI Map*/
  sprintf(FileName, "%sVegetation.LAI.%s%s", Path, Str, fileext);
 
  if (access(FileName, F_OK) == 0) {
    printf("Spatial LAI map provided, updating...\n");
    GetVarName(011, 0, VarName);
    GetVarNumberType(011, &NumberType);

  /*Read data monthy by month*/
    for (NSet = 0; NSet < 12; NSet++) {
      if (!(LAIMonthly = (float *)calloc(Map->NX * Map->NY,
        SizeOfNumberType(NumberType))))
        ReportError((char *)Routine, 1);
      flag = Read2DMatrix(FileName, LAIMonthly, NumberType, Map, NSet, VarName, 0);
          
      if ((Options->FileFormat == NETCDF && flag == 0)
        || (Options->FileFormat == BIN))
      {
        for (y = 0, i = 0; y < Map->NY; y++) {
          for (x = 0; x < Map->NX; x++, i++) {
            if ( VType[(*VegMap)[y][x].Veg - 1].OverStory == TRUE) {
            if (LAIMonthly[i] > 0.0)
                (*VegMap)[y][x].LAIMonthly[0][NSet] = LAIMonthly[i];
            else
              (*VegMap)[y][x].LAIMonthly[0][NSet] = VType[(*VegMap)[y][x].Veg - 1].LAIMonthly[0][NSet];

              if ( VType[(*VegMap)[y][x].Veg - 1].UnderStory  == TRUE )
                (*VegMap)[y][x].LAIMonthly[1][NSet] = VType[(*VegMap)[y][x].Veg - 1].LAIMonthly[1][NSet];
            }
            else{
              if ( VType[(*VegMap)[y][x].Veg - 1].UnderStory == TRUE)
                (*VegMap)[y][x].LAIMonthly[0][NSet] = VType[(*VegMap)[y][x].Veg - 1].LAIMonthly[0][NSet];
            }
          }
        }
      }
      else if (Options->FileFormat == NETCDF && flag == 1) {
        for (y = Map->NY - 1, i = 0; y >= 0; y--) {
          for (x = 0; x < Map->NX; x++, i++) {
          
            if ( VType[(*VegMap)[y][x].Veg - 1].OverStory == TRUE) {
              if (LAIMonthly[i] > 0.0)
                (*VegMap)[y][x].LAIMonthly[0][NSet] = LAIMonthly[i];
              else
              (*VegMap)[y][x].LAIMonthly[0][NSet] = VType[(*VegMap)[y][x].Veg - 1].LAIMonthly[0][NSet];
        
              if (VType[(*VegMap)[y][x].Veg - 1].UnderStory == TRUE)
                (*VegMap)[y][x].LAIMonthly[1][NSet] = VType[(*VegMap)[y][x].Veg - 1].LAIMonthly[1][NSet];
            }
            else{
              if (VType[(*VegMap)[y][x].Veg - 1].UnderStory == TRUE)
                (*VegMap)[y][x].LAIMonthly[0][NSet] = VType[(*VegMap)[y][x].Veg - 1].LAIMonthly[0][NSet];
            }
          }
        }
      }
      else ReportError((char *)Routine, 57);  

      free(LAIMonthly);     
    }   
  }
  else printf("Spatial LAI file not provided, LAI values not updated!\n");


  /*Update Vegetation Height Map*/
  sprintf(FileName, "%sVegetation.Height.%s%s", Path, Str, fileext);
  
   if (access(FileName, F_OK) == 0) {
    printf("Spatial tree height map provided, updating...\n");
    GetVarName(014, 0, VarName);
    GetVarNumberType(014, &NumberType);

    for (NSet = 0; NSet < 2; NSet++) {
      if (!(Height = (float *)calloc(Map->NX * Map->NY, SizeOfNumberType(NumberType))))
        ReportError((char *)Routine, 1);
      
      flag = Read2DMatrix(FileName, Height, NumberType, Map, NSet, VarName, 0);

      if ((Options->FileFormat == NETCDF && flag == 0)
        || (Options->FileFormat == BIN))
      {
        for (y = 0, i = 0; y < Map->NY; y++) {
          for (x = 0; x < Map->NX; x++, i++) {
              if (Height[i] > 0.0)
                (*VegMap)[y][x].Height[NSet] = Height[i];
              else
                (*VegMap)[y][x].Height[NSet] = VType[(*VegMap)[y][x].Veg - 1].Height[NSet];
          }
        }
      }
      else if (Options->FileFormat == NETCDF && flag == 1) {
        for (y = Map->NY - 1, i = 0; y >= 0; y--) {
          for (x = 0; x < Map->NX; x++, i++) {
              if (Height[i] > 0.0)
                (*VegMap)[y][x].Height[NSet] = Height[i];
              else
                (*VegMap)[y][x].Height[NSet] = VType[(*VegMap)[y][x].Veg - 1].Height[NSet];
          }
        }
      }
      else ReportError((char *)Routine, 57);
      free(Height);
    } 
  }
  else printf("Spatial HEIGHT file not provided, tree height values not updated!\n");
   
  if (Options->CanopyGapping) {
    /*Update Vegetation Gap Map*/
    sprintf(FileName, "%sVegetation.Gap.%s%s", Path, Str, fileext);
   
    if (access(FileName, F_OK) == 0) {
      printf("Spatial GAP map provided, updating...\n");
     
      GetVarName(007, 0, VarName);
      GetVarNumberType(007, &NumberType);
    
      if (!(Gap = (float *)calloc(Map->NX * Map->NY,
        SizeOfNumberType(NumberType))))
        ReportError((char *)Routine, 1);
      flag = Read2DMatrix(FileName, Gap, NumberType, Map, 0, VarName, 0);
      /* if NetCDF, may need to reverse the matrix */
      if ((Options->FileFormat == NETCDF && flag == 0)
        || (Options->FileFormat == BIN))
      {
        for (y = 0, i = 0; y < Map->NY; y++) {
          for (x = 0; x < Map->NX; x++, i++) {
            (*VegMap)[y][x].Gapping = Gap[i];
            /* set gapping to false for cells with no overstory */
            if (VType[(*VegMap)[y][x].Veg - 1].OverStory == FALSE)
              (*VegMap)[y][x].Gapping = 0.0;
          }
        }
      }
      else if (Options->FileFormat == NETCDF && flag == 1) {
        for (y = Map->NY - 1, i = 0; y >= 0; y--) {
          for (x = 0; x < Map->NX; x++, i++) {
            (*VegMap)[y][x].Gapping = Gap[i];
            /* set gapping to false for cells with no overstory */
            if (VType[(*VegMap)[y][x].Veg - 1].OverStory == FALSE)
              (*VegMap)[y][x].Gapping = 0.0;
            /* set gapping to false given glacier cell */
            if (VType[(*VegMap)[y][x].Veg - 1].Index == GLACIER)
              (*VegMap)[y][x].Gapping = 0.0;
          }
        }
      }
      else ReportError((char *)Routine, 57);
      free(Gap);
    }
    else printf("Spatial gap file not provided, gap values not updated!\n");

  }
}