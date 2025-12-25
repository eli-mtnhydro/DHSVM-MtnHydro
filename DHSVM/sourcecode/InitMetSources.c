
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include "settings.h"
#include "data.h"
#include "Calendar.h"
#include "fileio.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "getinit.h"
#include "constants.h"
#include "rad.h"

 /*******************************************************************************
   Function name: InitMetSources()

   Purpose      : Initialize and configure the model to process meteorological
                  data from various different sources
                  Processes the following section in the input file:
                  [METEOROLOGY]

   Required     :
     LISTPTR Input           - Linked list with input strings
     OPTIONSTRUCT *Options   - Structure with different program options
     MAPSIZE *Map            - Coverage and resolution of model area
     int NSoilLayers         - Number of soil layers
     TIMESTRUCT *Time        - Begin and end times, model timestep
     INPUTFILES *InFiles     - Various input filenames
     int *NStats             - Number of meteorological stations
     METLOCATION **Stat      - Station information

 *******************************************************************************/
void InitMetSources(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map,
  TOPOPIX **TopoMap, int NSoilLayers, TIMESTRUCT *Time, INPUTFILES *InFiles,
  int *NStats, METLOCATION **Stat)
{
  
  if (Options->Outside == TRUE) {
    printf("\nAll met stations in list will be included \n");
    if (Options->Prism == TRUE) {
      printf("WARNING: PRISM Option is also on\n");
      printf("Make sure file .prism files exist\n\n");
    }
  }
  
  InitStations(Input, Map, Time->NDaySteps, Options, NStats, Stat);
}

/*******************************************************************************
  Function name: InitStations()

  Purpose      : Read the station information from the options file.  This
                 information is in the [METEOROLOGY] section

  Required     :
    LISTPTR Input       - Linked list with input strings
    MAPSIZE *Map        - Information about the basin area
    int NDaysSteps      - Number of time steps in a day
    int *NStats         - Number of met stations
    METLOCATION **Stat  - Information about each met station

*****************************************************************************/
void InitStations(LISTPTR Input, MAPSIZE *Map, int NDaySteps,
  OPTIONSTRUCT *Options, int *NStats, METLOCATION **Stat)
{
  char *Routine = "InitStations";
  int i;
  int j;
  int k;
  char tempfilename[BUFSIZE * 2 + 7];
  char KeyName[station_file + 1][BUFSIZE + 1];
  char *KeyStr[] = {
    "STATION NAME",
    "NORTH COORDINATE",
    "EAST COORDINATE",
    "ELEVATION",
    "STATION FILE"
  };
  char *SectionName = "METEOROLOGY";
  char VarStr[station_file + 1][BUFSIZE + 1];
  float East;
  float North;
  FILE *PrismStatFile;
  FILE *SnowPatternStatFile;

  /* Get the number of different stations */
  GetInitString(SectionName, "NUMBER OF STATIONS", "", VarStr[0],
    (unsigned long)BUFSIZE, Input);
  if (!CopyInt(NStats, VarStr[0], 1))
    ReportError("NUMBER OF STATIONS", 51);

  if (*NStats <= 0)
    ReportError("Input Options File", 6);

  printf("\nEvaluating %d Met stations for inclusion\n", *NStats);

  /* Allocate memory for the stations */
  if (!(*Stat = (METLOCATION *)calloc(*NStats, sizeof(METLOCATION))))
    ReportError(Routine, 1);

  /* Read key-entry pairs for each station from the input file */
  /* for each potential station, up to NStats, read in the data and */
  /* determine if it is in the current model bounding box */
  /* If it is then put it into memory, otherwise, forget about it */
  /* unless Outside option is TRUE, then include it anyway */
  /* use temp counter k to track number of valid stations */
  k = 0;
  for (i = 0; i < *NStats; i++) {
    for (j = 0; j <= station_file; j++) {
      sprintf(KeyName[j], "%s %d", KeyStr[j], i + 1);
      GetInitString(SectionName, KeyName[j], "", VarStr[j],
        (unsigned long)BUFSIZE, Input);
    }

    /* Assign the entries to the variables */
    if (IsEmptyStr(VarStr[station_name]))
      ReportError(KeyName[station_name], 51);
    strcpy((*Stat)[k].Name, VarStr[station_name]);

    if (!CopyFloat(&North, VarStr[station_north], 1))
      ReportError(KeyName[station_north], 51);

    if (!CopyFloat(&East, VarStr[station_east], 1))
      ReportError(KeyName[station_east], 51);

    (*Stat)[k].Loc.N = Round(((Map->Yorig - 0.5 * Map->DY) - North) / Map->DY);
    (*Stat)[k].Loc.E = Round((East - (Map->Xorig + 0.5 * Map->DX)) / Map->DX);

    if (!CopyFloat(&((*Stat)[k].Elev), VarStr[station_elev], 1))
      ReportError(KeyName[station_elev], 51);

    if (IsEmptyStr(VarStr[station_file]))
      ReportError(KeyName[station_file], 51);
    strcpy((*Stat)[k].MetFile.FileName, VarStr[station_file]);

    OpenFile(&((*Stat)[k].MetFile.FilePtr), (*Stat)[k].MetFile.FileName, "r", FALSE);

    /* check to see if the stations are inside the bounding box */
    if (((*Stat)[k].Loc.N >= Map->NY || (*Stat)[k].Loc.N < 0 ||
      (*Stat)[k].Loc.E >= Map->NX || (*Stat)[k].Loc.E < 0)
      && Options->Outside == FALSE){
      k = k;
	  //printf("Station %d outside bounding box: %s ignored\n", i + 1, (*Stat)[k].Name);
	}
    else
      k = k + 1;
  }
  if (Options->Outside == FALSE)
    printf("Final number of stations in bounding box is %d \n\n", k);
  else
    printf("Forced to include all %d stations \n", k);
  *NStats = k;

  if (Options->Outside == TRUE && Options->Prism == TRUE) {

    for (i = 0; i < *NStats; i++) {
      sprintf(tempfilename, "%s.prism", (*Stat)[i].MetFile.FileName);
      OpenFile(&PrismStatFile, tempfilename, "rt", FALSE);
      for (k = 0; k < 12; k++) {
        if (fscanf(PrismStatFile, "%f ", &(*Stat)[i].PrismPrecip[k]) == EOF)
          ReportError(tempfilename, 2);
      }
      fclose(PrismStatFile);
    }
  }
  
  if (Options->SnowPattern == TRUE) {
    for (i = 0; i < *NStats; i++) {
      sprintf(tempfilename, "%s.snowpattern", (*Stat)[i].MetFile.FileName);
      OpenFile(&SnowPatternStatFile, tempfilename, "rt", FALSE);
      if (fscanf(SnowPatternStatFile, "%f ", &(*Stat)[i].SnowPatternBase) == EOF)
        ReportError(tempfilename, 2);
      fclose(SnowPatternStatFile);
    }
  }
}
