
/*
 * DESCRIPTION:  Store the state of the channel.  The channel state file
                 has one line per segment: the unique channel ID, the
                 storage in the segment in m3, and optionally the restart
                 memory (written by every dump): the routing storage
                 constant K (1/s, a moving average over time steps) and the
                 water depth at the top of the segment (m).
		 This is the content of the fields  
		   _channel_rec_
		     SegmentID id
		     float storage
		     float K
		     float top_water_depth
		 A file with only ID and storage (e.g. a cold start) starts K
		 from its initial value (channel_routing_parameters(), hydraulic
		 radius 3/4 of the bank height) and the water depth uniform.
		 Storage, K and depth are written with 9 significant digits so
		 that they are read back exactly.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h> 
#include <math.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "fileio.h"
#include "functions.h"
#include "constants.h"
#include "sizeofnt.h"
#include "channel.h"

typedef struct _RECORDSTRUCT {
  SegmentID id;
  float storage;
  float K;
  float top_water_depth;
  int memory;			/* TRUE if K and top_water_depth were read */
} RECORDSTRUCT;

int CompareRecord(const void *record1, const void *record2);
int CompareRecordID(const void *key, const void *record);

/*****************************************************************************
  ReadChannelState()

  Read the state of the channel from a previous run.  Currently just read an
  ASCII file, with the unique channel IDs in the first column and the amount
  of storage in the second column (m3), optionally followed by the routing
  storage constant K and the water depth at the top of the segment
*****************************************************************************/
void ReadChannelState(char *Path, DATE *Now, int deltat, Channel *Head)
{
  char InFileName[BUFSIZ + 15] = "";
  char Str[BUFSIZ + 1] = "";
  Channel *Current = NULL;
  FILE *InFile = NULL;
  int i = 0;
  int NLines = 0;
  int max_seg = 0;
  RECORDSTRUCT *Match = NULL;
  RECORDSTRUCT *Record = NULL;
  char Line[BUFSIZ + 1];
  int NFields;

  /* Re-create the storage file name and open it */
  sprintf(Str, "%02d.%02d.%04d.%02d.%02d.%02d", Now->Month, Now->Day,
	  Now->Year, Now->Hour, Now->Min, Now->Sec);
  snprintf(InFileName, sizeof(InFileName), "%sChannel.State.%s", Path, Str);
  OpenFile(&InFile, InFileName, "r", TRUE);
  NLines = CountLines(InFile);
  rewind(InFile);

  /* Allocate memory and read the file */
  Record = (RECORDSTRUCT *) calloc(NLines, sizeof(RECORDSTRUCT));
  if (Record == NULL)
    ReportError("ReadChannelState", 1);
  for (i = 0; i < NLines; i++) {
    if (fgets(Line, BUFSIZ, InFile) == NULL)
      ReportError(InFileName, 2);
    NFields = sscanf(Line, "%hu %f %f %f", &(Record[i].id), &(Record[i].storage),
                     &(Record[i].K), &(Record[i].top_water_depth));
    if (NFields != 2 && NFields != 4)
      ReportError(InFileName, 2);
    Record[i].memory = (NFields == 4);
  }
  qsort(Record, NLines, sizeof(RECORDSTRUCT), CompareRecord);

  /* Assign the storages to the correct IDs */
  Current = Head;
  while (Current) {
    Match = bsearch(&(Current->id), Record, NLines, sizeof(RECORDSTRUCT),
		    CompareRecordID);
    if (Current->id > max_seg)
      max_seg = Current->id;
    if (Match == NULL)
      ReportError("ReadChannelState", 55);
    Current->storage = Match->storage;
    
    /* Initialize depth uniformly in each segment */
    Current->top_water_depth = Current->storage / (Current->class2->width * Current->length);
    Current->bottom_water_depth = Current->top_water_depth;

    /* Restore the routing memory if present */
    if (Match->memory) {
      Current->K = Match->K;
      Current->X = exp(-Current->K * deltat);
      Current->top_water_depth = Match->top_water_depth;
    }
    
    Current = Current->next;
  }

  /* Clean up */
  if (Record)
    free(Record);
  fclose(InFile);
}

/*****************************************************************************
  StoreChannelState()

  Store the current state of the channel, i.e. the storage in each channel 
  segment.

*****************************************************************************/
void StoreChannelState(char *Path, DATE * Now, Channel * Head)
{
  char OutFileName[BUFSIZ + 15] = "";
  char Str[BUFSIZ + 1] = "";
  Channel *Current = NULL;
  FILE *OutFile = NULL;

  printf("Storing channel state\n");
  /* Create storage file */
  sprintf(Str, "%02d.%02d.%04d.%02d.%02d.%02d", Now->Month, Now->Day,
	  Now->Year, Now->Hour, Now->Min, Now->Sec);
  snprintf(OutFileName, sizeof(OutFileName), "%sChannel.State.%s", Path, Str);
  OpenFile(&OutFile, OutFileName, "w", TRUE);

  /* Store data */
  Current = Head;
  while (Current) {
    fprintf(OutFile, "%12hu ", Current->id);
    fprintf(OutFile, "%16.9g %16.9g %16.9g\n", Current->storage, Current->K,
            Current->top_water_depth);
    Current = Current->next;
  }

  /* Close file */
  fclose(OutFile);
}

/*****************************************************************************
  StoreChannelStateExtra()
  
  Store extra information about the channel state,
  i.e. storage, inflow, lateral inflow, infiltration, evaporation, and outflow
  
*****************************************************************************/
void StoreChannelStateExtra(char *Path, DATE * Now, Channel * Head)
{
  char OutFileName[BUFSIZ + 21] = "";
  char Str[BUFSIZ + 1] = "";
  Channel *Current = NULL;
  FILE *OutFile = NULL;
  
  printf("Storing channel state with extra information\n");
  /* Create storage file */
  sprintf(Str, "%02d.%02d.%04d.%02d.%02d.%02d", Now->Month, Now->Day,
          Now->Year, Now->Hour, Now->Min, Now->Sec);
  snprintf(OutFileName, sizeof(OutFileName), "%sChannel.State.Extra.%s", Path, Str);
  OpenFile(&OutFile, OutFileName, "w", TRUE);
  
  fprintf(OutFile, "%12s %12s %12s %12s %12s %12s %12s\n",
          "ID", "Storage", "Inflow", "LateralFlow",
          "Infiltration", "Evaporation", "Outflow");
  
  /* Store data */
  Current = Head;
  while (Current) {
    fprintf(OutFile, "%12hu ", Current->id);
    fprintf(OutFile, "%12g ", Current->storage);
    fprintf(OutFile, "%12g ", Current->inflow);
    fprintf(OutFile, "%12g ", Current->lateral_inflow);
    fprintf(OutFile, "%12g ", Current->infiltration);
    fprintf(OutFile, "%12g ", Current->evaporation);
    fprintf(OutFile, "%12g\n", Current->outflow);
    Current = Current->next;
  }
  
  /* Close file */
  fclose(OutFile);
}

/*****************************************************************************
  CompareRecord()

  Compare two RECORDSTRUCT elements for qsort
*****************************************************************************/
int CompareRecord(const void *record1, const void *record2)
{
  RECORDSTRUCT *x = NULL;
  RECORDSTRUCT *y = NULL;

  x = (RECORDSTRUCT *) record1;
  y = (RECORDSTRUCT *) record2;

  return (int) x->id - y->id;
}

/*****************************************************************************
  CompareRecordID()

  Compare RECORDSTRUCT element with an ID to see if the RECORDSTRUCT has the
  right ID for bsearch
*****************************************************************************/
int CompareRecordID(const void *key, const void *record)
{
  SegmentID *x = NULL;
  RECORDSTRUCT *y = NULL;

  x = (SegmentID *) key;
  y = (RECORDSTRUCT *) record;

  return (int) (*x - y->id);
}
