
/*
 * DESCRIPTION:  State of the channel segments in a snapshot: the file
 *               State.<stamp>.Segments.csv has one row per segment, in the
 *               order of the segment table:
 *
 *               Segment,Channel.Storage,Channel.RoutingK,Channel.TopWaterDepth
 *               -,m3,1/s,m
 *               <id>,<storage>,<K>,<top_water_depth>
 *
 *               Storage (m3), the routing storage constant K (1/s, a moving
 *               average over time steps) and the water depth at the top of
 *               the segment (m) are written with 9 significant digits, so
 *               that they are read back exactly.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "fileio.h"
#include "functions.h"
#include "constants.h"
#include "channel.h"

#define SEGMENT_COLUMNS "Segment,Channel.Storage,Channel.RoutingK,Channel.TopWaterDepth"
#define SEGMENT_UNITS "-,m3,1/s,m"

/*****************************************************************************
  ColdStartChannelState()

  Every segment holds water of depth Depth over the channel width and the
  segment length; the routing constants keep their initial values
  (channel_routing_parameters())
*****************************************************************************/
void ColdStartChannelState(float Depth, Channel *Head)
{
  Channel *Current;

  for (Current = Head; Current != NULL; Current = Current->next) {
    Current->storage = Depth * Current->class2->width * Current->length;
    Current->top_water_depth = Current->storage / (Current->class2->width * Current->length);
    Current->bottom_water_depth = Current->top_water_depth;
  }
}

/*****************************************************************************
  ReadChannelState()

  Reads the segment state of a snapshot. The rows must list the segments in
  the order of the segment table.
*****************************************************************************/
void ReadChannelState(char *FileName, int deltat, Channel *Head)
{
  FILE *InFile = NULL;
  Channel *Current = NULL;
  char Line[BUFSIZE * 2 + 1];
  char *Fields[5];
  float Value[3];
  int LineNo = 2;
  int i;
  long Id;

  if (!(InFile = fopen(FileName, "rb")))
    ReportErrorCode(112, "%s: not found", FileName);
  if (!ReadTextLine(InFile, Line, sizeof(Line), FileName, 1, 112) ||
      strcmp(Line, SEGMENT_COLUMNS) != 0)
    ReportErrorCode(112, "%s: line 1: the columns must be %s", FileName, SEGMENT_COLUMNS);
  if (!ReadTextLine(InFile, Line, sizeof(Line), FileName, 2, 112) ||
      strcmp(Line, SEGMENT_UNITS) != 0)
    ReportErrorCode(112, "%s: line 2: the units must be %s", FileName, SEGMENT_UNITS);

  for (Current = Head; Current != NULL; Current = Current->next) {
    LineNo++;
    if (!ReadTextLine(InFile, Line, sizeof(Line), FileName, LineNo, 112) ||
        IsBlankLine(Line))
      ReportErrorCode(112, "%s: line %d: no row, expected Segment %d", FileName, LineNo,
                      (int) Current->id);
    if (SplitFields(Line, Fields, 5) != 4)
      ReportErrorCode(112, "%s: line %d: expected 4 fields", FileName, LineNo);
    if (!TextToInteger(Fields[0], &Id) || Id != (long) Current->id)
      ReportErrorCode(112, "%s: line %d: Segment %s, expected %d", FileName, LineNo,
                      Fields[0], (int) Current->id);
    for (i = 0; i < 3; i++) {
      if (!TextToFloat(Fields[i + 1], &(Value[i])))
        ReportErrorCode(112, "%s: line %d: column %d: \"%s\" is not a number", FileName,
                        LineNo, i + 2, Fields[i + 1]);
    }

    Current->storage = Value[0];
    Current->top_water_depth = Current->storage / (Current->class2->width * Current->length);
    Current->bottom_water_depth = Current->top_water_depth;
    Current->K = Value[1];
    Current->top_water_depth = Value[2];
  }

  while (ReadTextLine(InFile, Line, sizeof(Line), FileName, ++LineNo, 112)) {
    if (!IsBlankLine(Line))
      ReportErrorCode(112, "%s: line %d: a row after the last segment", FileName, LineNo);
  }
  fclose(InFile);
}

/*****************************************************************************
  StoreChannelState()

  Writes the segment state of a snapshot
*****************************************************************************/
void StoreChannelState(char *FileName, Channel *Head)
{
  Channel *Current = NULL;
  FILE *OutFile = NULL;

  if (!(OutFile = fopen(FileName, "wb")))
    ReportErrorCode(113, "%s: cannot write", FileName);
  fprintf(OutFile, "%s\n%s\n", SEGMENT_COLUMNS, SEGMENT_UNITS);
  for (Current = Head; Current != NULL; Current = Current->next)
    fprintf(OutFile, "%d,%.9g,%.9g,%.9g\n", (int) Current->id, Current->storage,
            Current->K, Current->top_water_depth);
  if (ferror(OutFile) || fclose(OutFile) != 0)
    ReportErrorCode(113, "%s: cannot write", FileName);
}
