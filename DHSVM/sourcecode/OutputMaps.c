/*
 * DESCRIPTION:  Map stacks: one map file per requested variable layer,
 *               <Output Directory>/<layer name>.map, one layer per
 *               [OUTPUT] Map Time, named by the time. A layer stamped T
 *               holds a state's value at T (written after the time step
 *               ending at T, or at Start before the first step) or a value
 *               over the time step starting at T (written after that step).
 *               -9999 outside the mask and where a cell lacks the layer.
 */
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "fileio.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"
#include "mapio.h"
#include "output.h"
#include "registry.h"

/*****************************************************************************
  CreateMapStacks()

  Writes the header of every map stack (all layer names) and the time step
  after which each layer is written
*****************************************************************************/
void CreateMapStacks(MAPSIZE *Map, TIMESTRUCT *Time, DUMPSTRUCT *Dump)
{
  OUTLIST *List = &(Dump->List[OUT_MAPS]);
  char (*Names)[MAP_NAMESIZE] = NULL;
  long long Seconds;
  FILE *File;
  int i, j;

  Dump->Maps = NULL;
  if (List->N == 0)
    return;
  if (!(Dump->Maps = (MAPOUT *) calloc(List->N, sizeof(MAPOUT))) ||
      !(Names = calloc(Dump->NMapTimes, MAP_NAMESIZE)))
    ReportError("CreateMapStacks", 1);
  for (j = 0; j < Dump->NMapTimes; j++)
    SPrintIsoDate(&(Dump->MapTime[j]), Names[j]);

  for (i = 0; i < List->N; i++) {
    MAPOUT *Stack = &(Dump->Maps[i]);
    MakePath(Stack->FileName, sizeof(Stack->FileName), "%s%s.map", Dump->Path,
             List->Vars[i].Name);
    if (!(Stack->Step = (int *) calloc(Dump->NMapTimes, sizeof(int))))
      ReportError("CreateMapStacks", 1);
    for (j = 0; j < Dump->NMapTimes; j++) {
      Seconds = (MinutesOf(&(Dump->MapTime[j])) - MinutesOf(&(Time->Start))) * SECPMIN;
      Stack->Step[j] = (int) (Seconds / Time->Dt) + (RegistryIsInterval(List->Vars[i].Var) ? 1 : 0);
    }
    File = CreateMap(Stack->FileName, List->Vars[i].Name, RegistryRow(List->Vars[i].Var)->Units,
                     Dump->NMapTimes, Names, NULL, Map);
    CloseMapFile(File, Stack->FileName);
  }
  free(Names);
}

/*****************************************************************************
  WriteMapLayers()

  Appends the layers due after StepsDone time steps
*****************************************************************************/
void WriteMapLayers(int StepsDone, MODELVIEW *Model, DUMPSTRUCT *Dump)
{
  OUTLIST *List = &(Dump->List[OUT_MAPS]);
  MAPSIZE *Map = Model->Map;
  OUTVAR *Var;
  FILE *File;
  float *Array = NULL;
  float Value;
  int Has;
  int i, j, x, y;

  for (i = 0; i < List->N; i++) {
    for (j = 0; j < Dump->NMapTimes; j++) {
      if (Dump->Maps[i].Step[j] != StepsDone)
        continue;
      if (Array == NULL && !(Array = (float *) calloc(Map->NY * Map->NX, sizeof(float))))
        ReportError("WriteMapLayers", 1);
      Var = &(List->Vars[i]);
      for (y = 0; y < Map->NY; y++) {
        for (x = 0; x < Map->NX; x++) {
          Value = NA;
          if (INBASIN(Model->Topo[y][x].Mask)) {
            Value = CellValue(Var->Code, Var->Layer, Model, y, x, &Has);
            if (!Has)
              Value = NA;
          }
          Array[y * Map->NX + x] = Value;
        }
      }
      if (!(File = fopen(Dump->Maps[i].FileName, "ab")))
        ReportErrorCode(113, "%s: cannot write", Dump->Maps[i].FileName);
      WriteMapLayer(File, Array, Map, Dump->Maps[i].FileName);
      CloseMapFile(File, Dump->Maps[i].FileName);
    }
  }
  free(Array);
}
