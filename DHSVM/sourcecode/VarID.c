
/*
 * DESCRIPTION:  Maintains a structure that acts as a database with info on each
 *               variable, and provides functions to query this database
 
 * COMMENTS:     If the number of IDs increases it might be worthwhile to use a
 *               better, faster search.  This is not done here, because in the 
 *               overall scheme of DHSVM it is not worth the programming effort
 *               right now.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "sizeofnt.h"
#include "varid.h"

extern char fileext[];

struct {
  int ID;
  char Name[BUFSIZE + 1];
  char LongName[BUFSIZE + 1];
  char Format[BUFSIZE + 1];
  char Units[BUFSIZE + 1];
  char FileLabel[BUFSIZE + 1];
  int NumberType;
  int IsMultiLayer;
  int IsVegLayer;
  int IsSoilLayer;
  int AddLayer;
} varinfo[] = {
  {
  001, "Basin.DEM",
      "DEM", "%.3f",
      "m", "Digital Elevation Model", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  002, "Basin.Mask",
      "Basin mask", "%d", "", "Basin mask", NC_BYTE, FALSE, FALSE, FALSE, 0}, {
  003, "Soil.Type",
      "Soil type", "%d", "", "Soil type", NC_BYTE, FALSE, FALSE, FALSE, 0}, {
  004, "Soil.Depth",
      "Soil depth", "%.3f",
      "m", "Total soil depth", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  005, "Veg.Type",
      "Vegetation type", "%d",
      "", "Vegetation type", NC_BYTE, FALSE, FALSE, FALSE, 0}, {
  006, "Travel.Time",
      "Travel time", "%d",
      "hours", "Travel time", NC_SHORT, FALSE, FALSE, FALSE, 0}, {
  007, "Veg.CanopyGap",
       "Canopy Gap", "%.2f", 
       "", "Canopy Gap", NC_FLOAT, FALSE, FALSE, FALSE, 0 },{
  010, "Veg.Fract",
       "Overstory Fractional Coverage", "%.2f", 
       "", "Overstory Fractional Coverage", NC_FLOAT, FALSE, FALSE, FALSE, 0 },{
  011, "Veg.LAI",
       "Overstory Leaf Area Index", "%.2f", 
       "", "Overstory Leaf Area Index", NC_FLOAT, FALSE, FALSE, FALSE, 0 },{
  012, "Soil.KsLat",
       "Soil Lateral Conductivity", "%.6f", 
       "", "Soil Lateral Conductivity", NC_FLOAT, FALSE, FALSE, FALSE, 0 },{
  013, "Soil.Porosity",
       "Soil Porosity", "%.3f", 
       "", "Soil Porosity", NC_FLOAT, TRUE, FALSE, FALSE, 0 },{
  014, "Veg.Height",
       "Tree Height", "%.2f", 
       "", "Tree Height", NC_FLOAT, FALSE, FALSE, FALSE, 0 },{
  015, "Soil.FCap",
       "Soil Field Capacity", "%.3f", 
       "", "Soil Field Capacity", NC_FLOAT, TRUE, FALSE, FALSE, 0 },{  
  016, "Soil.KsLatExp",
       "Soil Exponential Decrease", "%.6f", 
       "", "Soil Exponential Decrease", NC_FLOAT, FALSE, FALSE, FALSE, 0 },{  
  017, "Lake.ID",
      "Lake Map ID", "%d", 
      "", "Lake Map ID", NC_BYTE, FALSE, FALSE, FALSE, 0 }, {
  100, "Met.PrecipMultiplier",
	   "PptMultiplier", "%.8f",
	   "", "Precipitation Multiplier", NC_FLOAT, FALSE, FALSE, FALSE, 0 },{
  101, "Evap.ETot",
      "Evapotranspiration (Total)", "%.4g",
      "m/timestep", "Total amount of evapotranspiration",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  102, "Evap.EPot",
      "Potential Evapotranspiration", "%.4g",
      "m/timestep", "Potential evaporation/transpiration",
      NC_FLOAT, TRUE, TRUE, FALSE, 1}, {
  103, "Evap.EInt",
      "Interception Evaporation", "%.4g",
      "m/timestep", "Evaporation from interception",
      NC_FLOAT, TRUE, TRUE, FALSE, 1}, {
  104, "Evap.ESoil",
      "Not implemented yet", "%.4g",
      "", "Not implemented yet", NC_FLOAT, TRUE, TRUE, FALSE, 0}, {
  105, "Evap.EAct",
      "Evaporation", "%.4g",
      "m/timestep", "Actual evaporation/transpiration",
      NC_FLOAT, TRUE, TRUE, FALSE, 1}, {
  201, "Precip",
      "Precipitation", "%.4g",
      "m/timestep", "Precipitation", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  202, "Precip.IntRain",
      "Interception Storage (liquid)", "%.4g",
      "m", "Interception storage (liquid)", NC_FLOAT, TRUE, TRUE, FALSE, 0}, {
  203, "Precip.IntSnow",
      "Interception Storage (frozen)", "%.4g",
      "m", "Interception storage (frozen)", NC_FLOAT, TRUE, TRUE, FALSE, 0}, {
  204, "Temp.Instor",
      "Temporary interception storage for top vegetation layer", "%.4g",
      "m", "Temporary interception storage for top vegetation layer",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  205, "PRISM.Precip",
      "PRISM Precipitation", "%.4g",
      "mm/month", "PRISM precipitation", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  206, "SumPrecip",
      "SumPrecipitation", "%.4g",
      "m", "Accumulated Precipitation", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  207, "Snow.Pattern",
      "SnowPattern", "%.4g",
      "", "Relative Snow Pattern", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  301, "Rad.ISW",
      "Incoming ShortWave Radiation", "%.4g",
      "W/m2", "Incoming ShortWave Radiation",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  302, "Rad.NSW",
      "Net Shortwave Radiation", "%.4g",
      "W/m2", "Net Shortwave solar radiation", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  303, "Rad.Beam",
      "Net Beam Radiation", "%.4g",
      "W/m2", "Net Beam Radiation", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  304, "Shade.Factor",
      "Shade Factor", "%d",
      "", "Shade Factor", NC_BYTE, FALSE, FALSE, FALSE, 0}, {
  305, "SkyView.Factor",
      "SkyView Factor", "%.4g",
      "-", "Skyview Factor", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  401, "Snow.HasSnow", "Snow Presence/Absence", "%1d", "", "Snow cover flag",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  402, "Snow.SnowCoverOver",
      "Overstory Snow Flag", "%1d", "", "Flag overstory can be covered",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  403, "Snow.LastSnow",
      "Last Snowfall", "%4d", "days", "Days since last snowfall",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  404, "Snow.Swq",
      "Snow Water Equivalent", "%.4g",
      "m", "Snow water equivalent", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  405, "Snow.Melt",
      "Snow Melt", "%.4g",
      "m/timestep", "Snow Melt", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  406, "Snow.PackWater",
      "Liquid Water Content (Deep Layer)", "%.4g",
      "m", "Liquid water content of snow pack",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  407, "Snow.TPack",
      "Snow Temperature (Deep Layer)", "%.4g",
      "C", "Temperature of snow pack", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  408, "Snow.SurfWater",
      "Liquid Water Content (Surface Layer)", "%.4g",
      "m", "Liquid water content of surface layer",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  409, "Snow.TSurf",
      "Snow Temperature (Surface Layer)", "%.4g",
      "C", "Temperature of snow pack surface layer",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  410, "Snow.ColdContent",
      "Snow Cold Content", "%.4g",
      "J", "Cold content of snow pack", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  411, "Snow.Albedo",
      "Snow Albedo", "%.4g",
      " ", "Albedo of snow pack surface", NC_FLOAT, FALSE, FALSE, FALSE, 0 }, {
  412, "Snow.MaxSwe",
      "Peak SWE", "%.4g",
      " ", "Peak SWE of current water year", NC_FLOAT, FALSE, FALSE, FALSE, 0 }, {
  413, "Snow.MaxSweDate",
      "Peak SWE Date", "%d",
      " ", "Peak SWE Date of current water year", NC_INT, FALSE, FALSE, FALSE, 0 }, {
  414, "Snow.MeltOutDate",
      "Melt out date", "%d",
      " ", "Snow disappearance date of current water year", NC_INT, FALSE, FALSE, FALSE, 0 }, {
  415, "Snow.SumAccum",
    "Snow Accumulation", "%.4g",
    " ", "Total snow accumulation for current water year", NC_FLOAT, FALSE, FALSE, FALSE, 0 }, {
  501, "Soil.Moist",
      "Soil Moisture Content", "%.4g",
      "", "Soil moisture for layer %d", NC_FLOAT, TRUE, FALSE, TRUE, 0}, {
  502, "Soil.Perc",
      "Percolation", "%.4g",
      "m/timestep", "Percolation", NC_FLOAT, TRUE, FALSE, TRUE, 0}, {
  503, "Soil.TableDepth",
      "Water Table Depth", "%.4g",
      "m below surface", "Depth of water table",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  504, "Soil.NetFlux",
      "Net Water Flux", "%.4g",
      "m/timestep", "Net flux of water", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  505, "Soil.TSurf",
      "Surface Temperature", "%.4g",
      "C", "Soil surface temperature", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  506, "Soil.Qnet",
      "Net Radiation", "%.4g",
      "W/m2", "Net radiation exchange at surface",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  507, "Soil.Qs",
      "Sensible Heat Flux", "%.4g",
      "W/m2", "Sensible heat exchange", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  508, "Soil.Qe",
      "Latent Heat Flux", "%.4g",
      "W/m2", "Latent heat exchange", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  509, "Soil.Qg",
      "Ground Heat Flux", "%.4g",
      "W/m2", "Ground heat exchange", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  510, "Soil.Qst",
      "Ground Heat Storage", "%.4g",
      "W/m2", "Ground heat storage", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  511, "Soil.Temp",
      "Soil Temperature", "%.4g",
      "C", "Soil Temperature", NC_FLOAT, TRUE, FALSE, TRUE, 0}, {
  512, "Soil.Runoff",
      "Surface Ponding", "%.4g",
      "m", "Surface Ponding", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  513, "SoilMap.IExcess",
      "Surface runoff from HOF and Return Flow", "%.4g",
      "m", "Surface runoff from HOF and Return Flow",
       NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  514, "SoilMap.InfiltAcc",
      "Infiltration Accumulation", "%.4g",
      "m", "Accumulated water in top layer",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  602, "Precip.Lapse",
      "Precipitation Lapse Rate", "%.5f",
      "", "Precipitation Lapse Rate", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  800, "Ts",
      "Snow Temperature Threshold", "%.4f", "", "Snow Temperature Threshold", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  801, "Tr",
      "Rain Temperature Threshold", "%.4f", "", "Rain Temperature Threshold", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  802, "Snow.amax",
      "Fresh Snow Albedo", "%.4f", "", "Fresh Snow Albedo", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {        
  803, "Snow.LamdaAcc",
      "Albedo lambda during accumulation", "%.4f", "", 
      "Albedo decay lambda during accumulation", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  804, "Snow.LamdaMelt",
      "Albedo lambda during melt", "%.4f", "", 
      "Albedo decay lambda during melt", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  805, "Snow.MinAlbedoAcc",
      "Min Albedo during accumulation", "%.4f", "", 
      "Min Albedo during accumulation", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {        
  806, "Snow.MinAlbedoMelt",
      "Min Albedo during melt", "%.4f", "", 
      "Min Albedo during melt", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {   
  ENDOFLIST, "", "", "", "", "",
      ENDOFLIST, ENDOFLIST, ENDOFLIST, ENDOFLIST, ENDOFLIST}
};

/*****************************************************************************
  GetVarAttr()
*****************************************************************************/
void GetVarAttr(MAPDUMP * DMap)
{
  GetVarName(DMap->ID, DMap->Layer, DMap->Name);
  GetVarLongName(DMap->ID, DMap->Layer, DMap->LongName);
  GetVarUnits(DMap->ID, DMap->Units);
  GetVarFileName(DMap->ID, DMap->Layer, DMap->Resolution, DMap->FileName);
  GetVarFileLabel(DMap->ID, DMap->FileLabel);
  GetVarNumberType(DMap->ID, &(DMap->NumberType));
}

/******************************************************************************/
/*				  GetVarName()                                */
/******************************************************************************/
void GetVarName(int ID, int Layer, char *Name)
{
  char *Routine = "GetVarName";
  int i;

  i = 0;
  while (varinfo[i].ID != ENDOFLIST) {
    if (varinfo[i].ID == ID) {
      if (varinfo[i].IsMultiLayer == TRUE)
	sprintf(Name, "%d.%s", Layer, varinfo[i].Name);
      else
	strcpy(Name, varinfo[i].Name);
      return;
    }
    i++;
  }
  ReportError((char *) Routine, 26);
}

/******************************************************************************/
/*				GetVarLongName()                              */
/******************************************************************************/
void GetVarLongName(int ID, int Layer, char *LongName)
{
  char *Routine = "GetVarLongName";
  int i;

  i = 0;
  while (varinfo[i].ID != ENDOFLIST) {
    if (varinfo[i].ID == ID) {
      if (varinfo[i].IsMultiLayer == TRUE)
	sprintf(LongName, "%s (Layer %d)", varinfo[i].LongName, Layer);
      else
	strcpy(LongName, varinfo[i].LongName);
      return;
    }
    i++;
  }
  ReportError((char *) Routine, 26);
}

/******************************************************************************/
/*				  GetVarUnits()                               */
/******************************************************************************/
void GetVarUnits(int ID, char *Units)
{
  char *Routine = "GetVarUnits";
  int i;

  i = 0;
  while (varinfo[i].ID != ENDOFLIST) {
    if (varinfo[i].ID == ID) {
      strcpy(Units, varinfo[i].Units);
      return;
    }
    i++;
  }
  ReportError((char *) Routine, 26);
}

/******************************************************************************/
/*				 GetVarFileName()                             */
/******************************************************************************/
void GetVarFileName(int ID, int Layer, unsigned char Resolution, char *FileName)
{
  char *Routine = "GetVarFileName";
  char Name[BUFSIZE + 1];
  char Str[BUFSIZE + 7];
  int i;

  i = 0;
  while (varinfo[i].ID != ENDOFLIST) {
    if (varinfo[i].ID == ID) {
      GetVarName(ID, Layer, Name);
      if (Resolution == MAP_OUTPUT) {
	sprintf(Str, "%sMap.%s%s", FileName, Name, fileext);
      }
      else
	ReportError((char *) Routine, 21);
      strncpy(FileName, Str, BUFSIZE + 7);
      return;
    }
    i++;
  }
  ReportError((char *) Routine, 26);
}

/******************************************************************************/
/*				GetVarFileLabel()                             */
/******************************************************************************/
void GetVarFileLabel(int ID, char *FileLabel)
{
  char *Routine = "GetVarFileLabel";
  int i;

  i = 0;
  while (varinfo[i].ID != ENDOFLIST) {
    if (varinfo[i].ID == ID) {
      strcpy(FileLabel, varinfo[i].FileLabel);
      return;
    }
    i++;
  }
  ReportError((char *) Routine, 26);
}

/******************************************************************************/
/*			       GetVarNumberType()                             */
/******************************************************************************/
void GetVarNumberType(int ID, int *NumberType)
{
  char *Routine = "GetVarNumberType";
  int i;

  i = 0;
  while (varinfo[i].ID != ENDOFLIST) {
    if (varinfo[i].ID == ID) {
      *NumberType = varinfo[i].NumberType;
      return;
    }
    i++;
  }
  ReportError((char *) Routine, 26);
}

/******************************************************************************/
/*				    IsValidID()                               */
/******************************************************************************/
unsigned char IsValidID(int ID)
{
  int i;

  i = 0;
  while (varinfo[i].ID != ENDOFLIST) {
    if (varinfo[i].ID == ID)
      return TRUE;
    i++;
  }
  return FALSE;
}

/******************************************************************************/
/*				 IsMultiLayer()                               */
/******************************************************************************/
unsigned char IsMultiLayer(int ID)
{
  char *Routine = "IsMultiLayer";
  int i;

  i = 0;
  while (varinfo[i].ID != ENDOFLIST) {
    if (varinfo[i].ID == ID)
      return varinfo[i].IsMultiLayer;
    i++;
  }
  ReportError((char *) Routine, 26);
  return FALSE;
}

/******************************************************************************/
/*			 	 GetVarNLayers()                              */
/******************************************************************************/
int GetVarNLayers(int ID, int MaxSoilLayers, int MaxVegLayers)
{
  char *Routine = "GetVarNLayers";
  int NLayers = -1;
  int i;

  i = 0;
  while (varinfo[i].ID != ENDOFLIST) {
    if (varinfo[i].ID == ID) {
      if (varinfo[i].IsVegLayer == TRUE)
	NLayers = MaxVegLayers + varinfo[i].AddLayer;
      else if (varinfo[i].IsSoilLayer == TRUE)
	NLayers = MaxSoilLayers + varinfo[i].AddLayer;
      else
	NLayers = 1;
      return NLayers;
    }
    i++;
  }
  ReportError((char *) Routine, 26);
  return NLayers;
}
