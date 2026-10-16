/*
 * DESCRIPTION:  What to write when: the [OUTPUT] keys (variable lists,
 *               classes, streamflow, network dynamics interval, map times),
 *               the snapshot instants of [RUN], the checks of the map times
 *               and the network dynamics interval against the time steps,
 *               and the point table.
 */
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "fileio.h"
#include "functions.h"
#include "constants.h"
#include "getinit.h"
#include "output.h"
#include "registry.h"
#include "tables.h"

static const char *ListBases[NOUTLISTS] = { "Basin", "Point", "Class", "Channels", "Maps" };

static void ReadListKeys(LISTPTR Input, const char *Base, OUTLIST *List);
static void ReadRequestTable(const char *Path, DUMPSTRUCT *Dump);
static void AddListKey(OUTLIST *List, const char *Where, const char *Value);

/*****************************************************************************
  InitOutputKeys()

  [OUTPUT]: each key present is read and its form checked (exit 103). The
  variable lists Basin, Point, Class, Channels and Maps, each continued by
  <Key> 2, <Key> 3, ... up to the first missing number; Classes (vegetation
  class numbers); Streamflow (TRUE or FALSE); Network Dynamics Interval
  (hours); Map Time 1, 2, ... (instants); Request Table (a table of
  variables whose requests follow the list keys of each kind; empty = none).
  Map Time 1 is required when a list of Maps names a variable. Every key is
  optional.
*****************************************************************************/
void InitOutputKeys(LISTPTR Input, DUMPSTRUCT *Dump)
{
  char Value[BUFSIZE + 1];
  char Text[BUFSIZE + 1];
  char KeyName[BUFSIZE + 1];
  char *Token;
  char *End;
  long Class;
  int i, n;

  for (i = 0; i < NOUTLISTS; i++) {
    memset(&(Dump->List[i]), 0, sizeof(OUTLIST));
    ReadListKeys(Input, ListBases[i], &(Dump->List[i]));
  }

  /* Vegetation classes with a Class_<N>.csv */
  Dump->ClassesText[0] = '\0';
  if (GetConfigValue(Input, "OUTPUT", "Classes", Value)) {
    strcpy(Text, Value);
    for (Token = strtok(Text, " \t"); Token != NULL; Token = strtok(NULL, " \t")) {
      errno = 0;
      Class = strtol(Token, &End, 10);
      if (End == Token || *End != '\0' || errno == ERANGE || Class < -1000000 || Class > 1000000)
        ReportErrorCode(103, "[OUTPUT] Classes = %s: expected vegetation class numbers", Value);
    }
    strcpy(Dump->ClassesText, Value);
  }

  /* Streamflow.csv */
  Dump->Streamflow = FALSE;
  if (GetConfigValue(Input, "OUTPUT", "Streamflow", Value)) {
    RequireConfigValue(Input, "OUTPUT", "Streamflow", Value);
    Dump->Streamflow = ConfigBool("OUTPUT", "Streamflow", Value);
  }

  /* NetworkDynamics.csv: the interval in hours (checked against the time
     step by CheckOutputTimes()) */
  Dump->NetworkText[0] = '\0';
  Dump->NetworkSteps = 0;
  if (GetConfigValue(Input, "OUTPUT", "Network Dynamics Interval", Value)) {
    RequireConfigValue(Input, "OUTPUT", "Network Dynamics Interval", Value);
    ConfigDouble("OUTPUT", "Network Dynamics Interval", Value);
    strcpy(Dump->NetworkText, Value);
  }

  /* Layer times of the map stacks */
  for (n = 0; ; n++) {
    snprintf(KeyName, sizeof(KeyName), "Map Time %d", n + 1);
    if (!GetConfigValue(Input, "OUTPUT", KeyName, Value))
      break;
  }
  Dump->NMapTimes = n;
  Dump->MapTime = NULL;
  if (n > 0 && !(Dump->MapTime = (DATE *) calloc(n, sizeof(DATE))))
    ReportError("InitOutputKeys", 1);
  for (i = 0; i < n; i++) {
    snprintf(KeyName, sizeof(KeyName), "Map Time %d", i + 1);
    RequireConfigValue(Input, "OUTPUT", KeyName, Value);
    if (!SScanIsoDate(Value, &(Dump->MapTime[i])))
      ReportErrorCode(103, "[OUTPUT] %s = %s: expected YYYY-MM-DD HH:MM", KeyName, Value);
  }

  /* Requests of the request table, after the list keys of each kind */
  if (GetConfigValue(Input, "OUTPUT", "Request Table", Value) && !IsEmptyStr(Value))
    ReadRequestTable(Value, Dump);

  for (i = 0; i < Dump->List[OUT_MAPS].NKeys; i++) {
    if (!IsEmptyStr(Dump->List[OUT_MAPS].Value[i]) && Dump->NMapTimes == 0)
      ReportErrorCode(103, "[OUTPUT] Map Time 1 is missing");
  }
}

/* The keys Base, Base 2, Base 3, ... of a variable list */
static void ReadListKeys(LISTPTR Input, const char *Base, OUTLIST *List)
{
  char Value[BUFSIZE + 1];
  char KeyName[BUFSIZE + 1];
  char Where[BUFSIZE + 1];
  int n;

  if (GetConfigValue(Input, "OUTPUT", Base, Value)) {
    snprintf(Where, sizeof(Where), "[OUTPUT] %s", Base);
    AddListKey(List, Where, Value);
  }
  for (n = 2; ; n++) {
    snprintf(KeyName, sizeof(KeyName), "%s %d", Base, n);
    if (!GetConfigValue(Input, "OUTPUT", KeyName, Value))
      break;
    snprintf(Where, sizeof(Where), "[OUTPUT] %s %d", Base, n);
    AddListKey(List, Where, Value);
  }
}

/*****************************************************************************
  ReadRequestTable()

  The request table: column Variable (one name of a variable list) and the
  flag columns Basin, Point, Class, Channels and Maps (1 = requested, 0 or
  an empty field = not; a missing column is 0 in every row). Each row
  flagged for a kind adds its variable to that kind's list, in file order,
  named in messages by its line.
*****************************************************************************/
static void ReadRequestTable(const char *Path, DUMPSTRUCT *Dump)
{
  TABLE Table;
  char Where[2 * BUFSIZE + 1];
  int Column[NOUTLISTS];
  int Variable;
  int Flag;
  int i, k;

  ReadTable(Path, &Table);
  Variable = TableColumn(&Table, "Variable", "-", TRUE);
  for (k = 0; k < NOUTLISTS; k++)
    Column[k] = TableColumn(&Table, ListBases[k], "-", FALSE);
  TableWarnUnknown(&Table);

  for (i = 0; i < Table.NRows; i++) {
    for (k = 0; k < NOUTLISTS; k++) {
      if (TableIsEmpty(&Table, i, Column[k]))
        continue;
      Flag = TableInt(&Table, i, Column[k]);
      if (Flag != 0 && Flag != 1)
        ReportErrorCode(106, "%s: line %d: column %s is %d, expected 0 or 1", Path, i + 3,
                        ListBases[k], Flag);
      if (Flag == 1) {
        snprintf(Where, sizeof(Where), "%s: line %d", Path, i + 3);
        AddListKey(&(Dump->List[k]), Where, TableText(&Table, i, Variable));
      }
    }
  }
  FreeTable(&Table);
}

/* Adds a value to a variable list; Where names its origin in messages
   ([OUTPUT] <Key>, or <request table>: line <l>) */
static void AddListKey(OUTLIST *List, const char *Where, const char *Value)
{
  if (!(List->Key = realloc(List->Key, (List->NKeys + 1) * sizeof(*(List->Key)))) ||
      !(List->Value = realloc(List->Value, (List->NKeys + 1) * sizeof(*(List->Value)))))
    ReportError("AddListKey", 1);
  snprintf(List->Key[List->NKeys], sizeof(*(List->Key)), "%s", Where);
  snprintf(List->Value[List->NKeys], BUFSIZE + 1, "%s", Value);
  List->NKeys++;
}

/*****************************************************************************
  Function name: InitStateDump()

  Purpose      : Read the snapshot instants [RUN] State Time 1, 2, ... (the
                 list ends at the first missing number). Each lies on the
                 step grid between Start and End, after the one before.

  Required     :
    LISTPTR Input         - Linked list with input strings
    TIMESTRUCT *Time      - Run period and time step
    DUMPSTRUCT *Dump      - Snapshot instants and their step numbers

*****************************************************************************/
void InitStateDump(LISTPTR Input, TIMESTRUCT *Time, DUMPSTRUCT *Dump)
{
  char *Routine = "InitStateDump";
  int i;			/* counter */
  int n;
  char KeyName[BUFSIZE + 1];
  char VarStr[BUFSIZE + 1];
  long long Seconds;
  DATE Instant;

  /* Count the keys */
  for (n = 0; ; n++) {
    snprintf(KeyName, sizeof(KeyName), "State Time %d", n + 1);
    if (!GetConfigValue(Input, "RUN", KeyName, VarStr))
      break;
  }
  Dump->NStates = n;
  Dump->DState = NULL;
  Dump->StateStep = NULL;
  if (n == 0)
    return;

  if (!(Dump->DState = (DATE *)calloc(n, sizeof(DATE))) ||
      !(Dump->StateStep = (int *)calloc(n, sizeof(int))))
    ReportError(Routine, 1);

  for (i = 0; i < n; i++) {
    snprintf(KeyName, sizeof(KeyName), "State Time %d", i + 1);
    RequireConfigValue(Input, "RUN", KeyName, VarStr);
    if (!SScanIsoDate(VarStr, &Instant))
      ReportErrorCode(103, "[RUN] %s = %s: expected YYYY-MM-DD HH:MM", KeyName, VarStr);
    Seconds = (MinutesOf(&Instant) - MinutesOf(&(Time->Start))) * SECPMIN;
    if (Seconds < 0 || Seconds > (long long) Time->NSteps * Time->Dt)
      ReportErrorCode(110, "[RUN] %s = %s: not between Start and End", KeyName, VarStr);
    if (Seconds % Time->Dt != 0)
      ReportErrorCode(110, "[RUN] %s = %s: not on the time step grid", KeyName, VarStr);
    CopyDate(&((Dump->DState)[i]), &Instant);
    Dump->StateStep[i] = (int) (Seconds / Time->Dt);
    if (i > 0 && Dump->StateStep[i] <= Dump->StateStep[i - 1])
      ReportErrorCode(110, "[RUN] %s = %s: not after State Time %d", KeyName, VarStr, i);
  }
}

/*****************************************************************************
  CheckOutputTimes()

  Every Map Time on the step grid, after the one before, between Start
  and End, and before End when a requested map variable is a value over a
  time step (its layer stamped T holds the step starting at T). The
  network dynamics interval, converted to seconds as the time step, a
  positive whole number of time steps (0 = none); an interval longer than
  the run writes the first step only.
*****************************************************************************/
void CheckOutputTimes(OPTIONSTRUCT *Options, MAPSIZE *Map, TIMESTRUCT *Time,
                      DUMPSTRUCT *Dump)
{
  OUTLIST *Maps = &(Dump->List[OUT_MAPS]);
  char Text[NAMESIZE + 1];
  const char *Interval = NULL;
  long long Seconds;
  long long Last = -1;
  double Hours;
  double x;
  long long n;
  int i;

  /* The map variables, as far as their names are valid (the names are
     checked with the other output requests) */
  ExpandList(Maps, KIND_MAP, "a map", Options, Map->SoilLayers, 0);
  for (i = 0; i < Maps->N && Interval == NULL; i++) {
    if (RegistryIsInterval(Maps->Vars[i].Var))
      Interval = Maps->Vars[i].Name;
  }

  for (i = 0; i < Dump->NMapTimes; i++) {
    SPrintIsoDate(&(Dump->MapTime[i]), Text);
    Seconds = (MinutesOf(&(Dump->MapTime[i])) - MinutesOf(&(Time->Start))) * SECPMIN;
    if (Seconds % Time->Dt != 0)
      ReportErrorCode(110, "[OUTPUT] Map Time %d = %s: not on the time step grid", i + 1, Text);
    if (Seconds < 0 || Seconds > (long long) Time->NSteps * Time->Dt)
      ReportErrorCode(110, "[OUTPUT] Map Time %d = %s: not between Start and End", i + 1, Text);
    if (Seconds <= Last)
      ReportErrorCode(110, "[OUTPUT] Map Time %d = %s: not after Map Time %d", i + 1, Text, i);
    if (Interval != NULL && Seconds == (long long) Time->NSteps * Time->Dt)
      ReportErrorCode(110, "[OUTPUT] Map Time %d = %s: %s is a value over a time step, its "
                      "last layer can start at End - Time Step", i + 1, Text, Interval);
    Last = Seconds;
  }

  if (!IsEmptyStr(Dump->NetworkText)) {
    Hours = strtod(Dump->NetworkText, NULL);
    if (Hours != 0.0) {
      x = Hours * 3600.0;
      if (!(x > 0.0 && x < 1e15))
        ReportErrorCode(110, "[OUTPUT] Network Dynamics Interval = %s: must be a multiple of "
                        "the time step", Dump->NetworkText);
      n = llround(x);
      if (!(fabs(x - (double) n) < 0.01) || n <= 0 || n % Time->Dt != 0)
        ReportErrorCode(110, "[OUTPUT] Network Dynamics Interval = %s: must be a multiple of "
                        "the time step", Dump->NetworkText);
      if (n / Time->Dt < Time->NSteps)
        Dump->NetworkSteps = (int) (n / Time->Dt);
      else
        Dump->NetworkSteps = Time->NSteps;
    }
  }
}

/*****************************************************************************
  ReadPointTable()

  The point table ([INPUTS] Point Table), read when [OUTPUT] Point names a
  variable: columns Name, North, East (m). The names and cells are checked
  with the other output requests.
*****************************************************************************/
void ReadPointTable(LISTPTR Input, DUMPSTRUCT *Dump)
{
  TABLE Table;
  char Path[BUFSIZE + 1];
  int Name, North, East;
  int i;

  Dump->NPoints = 0;
  Dump->Points = NULL;
  if (!OutputListNamesAny(&(Dump->List[OUT_POINT])) ||
      !GetConfigValue(Input, "INPUTS", "Point Table", Path) || IsEmptyStr(Path))
    return;

  ReadTable(Path, &Table);
  Name = TableColumn(&Table, "Name", "-", TRUE);
  North = TableColumn(&Table, "North", "m", TRUE);
  East = TableColumn(&Table, "East", "m", TRUE);
  TableWarnUnknown(&Table);
  Dump->NPoints = Table.NRows;
  if (Table.NRows > 0 && !(Dump->Points = (POINTOUT *) calloc(Table.NRows, sizeof(POINTOUT))))
    ReportError("ReadPointTable", 1);
  for (i = 0; i < Table.NRows; i++) {
    snprintf(Dump->Points[i].Name, sizeof(Dump->Points[i].Name), "%s",
             TableText(&Table, i, Name));
    Dump->Points[i].Line = i + 3;
    Dump->Points[i].North = TableDouble(&Table, i, North);
    Dump->Points[i].East = TableDouble(&Table, i, East);
  }
  FreeTable(&Table);
}

/* TRUE if a variable list names anything */
int OutputListNamesAny(OUTLIST *List)
{
  int i;

  for (i = 0; i < List->NKeys; i++) {
    if (!IsEmptyStr(List->Value[i]))
      return TRUE;
  }
  return FALSE;
}
