/*
 * DESCRIPTION:  Comma-separated tables (the class, channel, lake and station
 *               tables). Line 1 holds the column names, line 2 the units of
 *               each column, data start on line 3. Columns are found by
 *               name; the units of every column used must be the expected
 *               ones. No quoting: no field holds a comma or a double quote.
 *               A leading UTF-8 byte order mark is skipped, LF and CR LF
 *               line ends are read, blank lines at the end are ignored.
 *               Errors stop the run with exit code 106 naming the file (and
 *               the line and column); columns never looked up are warned.
 */
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "DHSVMerror.h"
#include "fileio.h"
#include "mapio.h"
#include "tables.h"

static const char *MonthTokens[12] = {
  "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
};

static char *CopyText(const char *Text)
{
  char *Copy = malloc(strlen(Text) + 1);

  if (Copy == NULL)
    ReportError("Tables", 1);
  strcpy(Copy, Text);
  return Copy;
}

/*****************************************************************************
  ReadTable()

  Reads the whole table into memory and checks its form: names and units
  lines, the same number of fields on every line, no blank line between data
  lines, unique column names, no double quote.
*****************************************************************************/
void ReadTable(const char *Path, TABLE *Table)
{
  FILE *File;
  char *Line;
  char **Fields;
  int MaxFields = 1;
  int LineNo = 0;
  int NBlank = 0;
  int Capacity = 0;
  int n, i, j;

  memset(Table, 0, sizeof(TABLE));
  if (strlen(Path) >= sizeof(Table->Path))
    ReportErrorCode(104, "%.40s...: path longer than %d characters", Path,
                    (int) sizeof(Table->Path) - 1);
  strcpy(Table->Path, Path);
  if (!(File = fopen(Path, "rb")))
    ReportErrorCode(104, "%s: cannot open", Path);
  if (!(Line = malloc(TABLE_LINESIZE)))
    ReportError("ReadTable", 1);

  /* Column names */
  if (!ReadTextLine(File, Line, TABLE_LINESIZE, Path, ++LineNo, 106) || IsBlankLine(Line))
    ReportErrorCode(106, "%s: line 1 (column names) is missing", Path);
  if (strchr(Line, '"'))
    ReportErrorCode(106, "%s: line 1 holds a double quote", Path);
  for (n = 1, i = 0; Line[i] != '\0'; i++)
    if (Line[i] == ',')
      n++;
  Table->NColumns = n;
  if (!(Table->Names = calloc(n, sizeof(char *))) ||
      !(Table->Units = calloc(n, sizeof(char *))) ||
      !(Table->Known = calloc(n, sizeof(char))) ||
      !(Fields = calloc(n, sizeof(char *))))
    ReportError("ReadTable", 1);
  SplitFields(Line, Fields, n);
  for (j = 0; j < n; j++) {
    if (Fields[j][0] == '\0')
      ReportErrorCode(106, "%s: column %d has no name", Path, j + 1);
    for (i = 0; i < j; i++) {
      if (strcmp(Fields[i], Fields[j]) == 0)
        ReportErrorCode(106, "%s: column %s repeated", Path, Fields[j]);
    }
    Table->Names[j] = CopyText(Fields[j]);
  }

  /* Units */
  if (!ReadTextLine(File, Line, TABLE_LINESIZE, Path, ++LineNo, 106) || IsBlankLine(Line))
    ReportErrorCode(106, "%s: line 2 (units) is missing", Path);
  if (strchr(Line, '"'))
    ReportErrorCode(106, "%s: line 2 holds a double quote", Path);
  if (SplitFields(Line, Fields, n) != n)
    ReportErrorCode(106, "%s: line 2 has a different number of fields than line 1", Path);
  for (j = 0; j < n; j++)
    Table->Units[j] = CopyText(Fields[j]);
  free(Fields);
  MaxFields = n;

  /* Data */
  while (ReadTextLine(File, Line, TABLE_LINESIZE, Path, ++LineNo, 106)) {
    if (IsBlankLine(Line)) {
      NBlank++;
      continue;
    }
    if (NBlank > 0)
      ReportErrorCode(106, "%s: line %d is blank", Path, LineNo - NBlank);
    if (strchr(Line, '"'))
      ReportErrorCode(106, "%s: line %d holds a double quote", Path, LineNo);
    if (Table->NRows == Capacity) {
      Capacity = (Capacity == 0) ? 64 : 2 * Capacity;
      if (!(Table->Lines = realloc(Table->Lines, Capacity * sizeof(char *))) ||
          !(Table->Fields = realloc(Table->Fields, Capacity * sizeof(char **))))
        ReportError("ReadTable", 1);
    }
    Table->Lines[Table->NRows] = CopyText(Line);
    if (!(Table->Fields[Table->NRows] = calloc(MaxFields, sizeof(char *))))
      ReportError("ReadTable", 1);
    if (SplitFields(Table->Lines[Table->NRows], Table->Fields[Table->NRows], MaxFields) != n)
      ReportErrorCode(106, "%s: line %d has a different number of fields than line 1",
                      Path, LineNo);
    Table->NRows++;
  }
  free(Line);
  fclose(File);
}

/*****************************************************************************
  FreeTable()
*****************************************************************************/
void FreeTable(TABLE *Table)
{
  int i;

  for (i = 0; i < Table->NColumns; i++) {
    free(Table->Names[i]);
    free(Table->Units[i]);
  }
  for (i = 0; i < Table->NRows; i++) {
    free(Table->Lines[i]);
    free(Table->Fields[i]);
  }
  free(Table->Names);
  free(Table->Units);
  free(Table->Known);
  free(Table->Lines);
  free(Table->Fields);
  memset(Table, 0, sizeof(TABLE));
}

/*****************************************************************************
  TableColumn()

  Index of the column Name, its units checked; -1 if it is absent and not
  Required
*****************************************************************************/
int TableColumn(TABLE *Table, const char *Name, const char *Units, int Required)
{
  int j;

  for (j = 0; j < Table->NColumns; j++) {
    if (strcmp(Table->Names[j], Name) == 0) {
      Table->Known[j] = TRUE;
      if (strcmp(Table->Units[j], Units) != 0)
        ReportErrorCode(106, "%s: column %s has units \"%s\", expected \"%s\"",
                        Table->Path, Name, Table->Units[j], Units);
      return j;
    }
  }
  if (Required)
    ReportErrorCode(106, "%s: column %s missing", Table->Path, Name);
  return -1;
}

/* Name of a layer column: <Base>_1 .. <Base>_L, <Base>_Deep, or
   <Base>_Jan .. <Base>_Dec */
static void LayerColumnName(char *Buffer, size_t Size, const char *Base, int LayerSet,
                            int SoilLayers, int k)
{
  if (LayerSet == MAP_MONTH)
    snprintf(Buffer, Size, "%s_%s", Base, MonthTokens[k]);
  else if (k < SoilLayers)
    snprintf(Buffer, Size, "%s_%d", Base, k + 1);
  else
    snprintf(Buffer, Size, "%s_Deep", Base);
}

/*****************************************************************************
  TableLayerColumns()

  Indexes of the layer columns of Base: 12 months (MAP_MONTH), the soil
  layers 1..L and Deep (MAP_SOIL), or the soil layers 1..L only (MAP_ONE,
  for the root zone columns); -1 where absent and not Required
*****************************************************************************/
void TableLayerColumns(TABLE *Table, const char *Base, int LayerSet, int SoilLayers,
                       const char *Units, int Required, int *Columns)
{
  char Name[BUFSIZE + 1];
  int n = (LayerSet == MAP_MONTH) ? 12 : (LayerSet == MAP_SOIL ? SoilLayers + 1 : SoilLayers);
  int k;

  for (k = 0; k < n; k++) {
    LayerColumnName(Name, sizeof(Name), Base, LayerSet, SoilLayers, k);
    Columns[k] = TableColumn(Table, Name, Units, Required);
  }
}

/*****************************************************************************
  TableKnown(), TableKnownLayers()

  Mark columns the model knows but does not use under the run's options
  (read silently and ignored)
*****************************************************************************/
void TableKnown(TABLE *Table, const char *Name)
{
  int j;

  for (j = 0; j < Table->NColumns; j++) {
    if (strcmp(Table->Names[j], Name) == 0)
      Table->Known[j] = TRUE;
  }
}

void TableKnownLayers(TABLE *Table, const char *Base, int LayerSet, int SoilLayers)
{
  char Name[BUFSIZE + 1];
  int n = (LayerSet == MAP_MONTH) ? 12 : (LayerSet == MAP_SOIL ? SoilLayers + 1 : SoilLayers);
  int k;

  for (k = 0; k < n; k++) {
    LayerColumnName(Name, sizeof(Name), Base, LayerSet, SoilLayers, k);
    TableKnown(Table, Name);
  }
}

/*****************************************************************************
  TableIsEmpty(), TableText()
*****************************************************************************/
int TableIsEmpty(TABLE *Table, int Row, int Column)
{
  return (Column < 0 || Table->Fields[Row][Column][0] == '\0');
}

const char *TableText(TABLE *Table, int Row, int Column)
{
  if (TableIsEmpty(Table, Row, Column))
    ReportErrorCode(106, "%s: line %d: column %s is empty", Table->Path, Row + 3,
                    Table->Names[Column]);
  return Table->Fields[Row][Column];
}

/*****************************************************************************
  TableFloat(), TableDouble(), TableInt()

  A field as a number: (float) strtod, strtod, base-10 strtol; the whole
  field must be the number
*****************************************************************************/
float TableFloat(TABLE *Table, int Row, int Column)
{
  return (float) TableDouble(Table, Row, Column);
}

double TableDouble(TABLE *Table, int Row, int Column)
{
  const char *Text = TableText(Table, Row, Column);
  char *End;
  double x;

  x = strtod(Text, &End);
  if (End == Text || *End != '\0')
    ReportErrorCode(106, "%s: line %d: column %s: \"%s\" is not a number", Table->Path,
                    Row + 3, Table->Names[Column], Text);
  return x;
}

int TableInt(TABLE *Table, int Row, int Column)
{
  const char *Text = TableText(Table, Row, Column);
  char *End;
  long x;

  errno = 0;
  x = strtol(Text, &End, 10);
  if (End == Text || *End != '\0' || errno == ERANGE || x < INT_MIN || x > INT_MAX)
    ReportErrorCode(106, "%s: line %d: column %s: \"%s\" is not an integer", Table->Path,
                    Row + 3, Table->Names[Column], Text);
  return (int) x;
}

/*****************************************************************************
  TableCheckRows()

  A class table: at least one row, at most MaxRows, and the IdColumn of row
  i equal to i
*****************************************************************************/
void TableCheckRows(TABLE *Table, const char *IdColumn, int MaxRows)
{
  int Column = TableColumn(Table, IdColumn, "-", TRUE);
  int i;

  if (Table->NRows < 1)
    ReportErrorCode(106, "%s: no data rows", Table->Path);
  if (Table->NRows > MaxRows)
    ReportErrorCode(106, "%s: %d rows, at most %d allowed", Table->Path, Table->NRows, MaxRows);
  for (i = 0; i < Table->NRows; i++) {
    if (TableInt(Table, i, Column) != i + 1)
      ReportErrorCode(106, "%s: line %d: %s is %s, expected %d", Table->Path, i + 3,
                      IdColumn, TableText(Table, i, Column), i + 1);
  }
}

/*****************************************************************************
  TableWarnUnknown()

  Warns about every column never looked up
*****************************************************************************/
void TableWarnUnknown(TABLE *Table)
{
  int j;

  for (j = 0; j < Table->NColumns; j++) {
    if (!Table->Known[j])
      fprintf(stderr, "DHSVM WARNING: %s: column %s not used\n", Table->Path, Table->Names[j]);
  }
}
