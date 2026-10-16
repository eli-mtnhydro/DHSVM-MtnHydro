
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "fileio.h"
#include "functions.h"
#include "constants.h"
#include "tables.h"

/*****************************************************************************
  ReadMetRecord()

  Reads the record of the time step that starts at Current from a station
  data file: its time must be the step start. The lapse rate is the
  station's lapse rate of the month of the step.
*****************************************************************************/
void ReadMetRecord(DATE *Current, METLOCATION *Station)
{
  char *Path = Station->MetFile.FileName;
  char Found[BUFSIZE + 1];
  char Expected[BUFSIZE + 1];
  float Array[NMETCOLUMNS];	/* Values of the record, in the order of the columns */
  DATE MetDate;			/* Date of meteorological record */
  MET *MetRecord = &(Station->Data);
  int k;

  SPrintIsoDate(Current, Expected);
  if (Station->Pending)
    Station->Pending = FALSE;
  else {
    if (!ReadTextLine(Station->MetFile.FilePtr, Station->Record, TABLE_LINESIZE, Path,
                      ++(Station->LineNo), 107))
      ReportErrorCode(107, "%s: line %d: no record, expected %s", Path, Station->LineNo,
                      Expected);
    if (IsBlankLine(Station->Record))
      ReportErrorCode(107, "%s: line %d is blank, expected %s", Path, Station->LineNo,
                      Expected);
    if (SplitFields(Station->Record, Station->Fields, Station->NFields) != Station->NFields)
      ReportErrorCode(107, "%s: line %d has a different number of fields than line 1", Path,
                      Station->LineNo);
  }

  if (!SScanIsoDate(Station->Fields[Station->Column[0]], &MetDate))
    ReportErrorCode(107, "%s: line %d: time \"%s\" is not YYYY-MM-DD HH:MM", Path,
                    Station->LineNo, Station->Fields[Station->Column[0]]);
  if (MinutesOf(&MetDate) != MinutesOf(Current)) {
    SPrintIsoDate(&MetDate, Found);
    ReportErrorCode(107, "%s: line %d has time %s, expected %s", Path, Station->LineNo,
                    Found, Expected);
  }

  for (k = 1; k < NMETCOLUMNS; k++) {
    if (!TextToFloat(Station->Fields[Station->Column[k]], &Array[k]))
      ReportErrorCode(107, "%s: line %d: field %d: \"%s\" is not a number", Path,
                      Station->LineNo, Station->Column[k] + 1,
                      Station->Fields[Station->Column[k]]);
  }

  MetRecord->Tair = Array[1];
  MetRecord->Wind = Array[2];
  MetRecord->Rh = Array[3];
  MetRecord->Sin = Array[4];
  MetRecord->Lin = Array[5];
  MetRecord->Precip = Array[6];

  MetRecord->TempLapse = Station->LapseRate[Current->Month - 1];
}
