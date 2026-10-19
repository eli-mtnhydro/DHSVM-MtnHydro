#ifndef TABLES_H
#define TABLES_H

#include "settings.h"

#define TABLE_LINESIZE 65536	/* Longest line of a table file */
#define TABLE_MAXROWS  254	/* Most rows of a class table (the class maps are BYTE, 255 = no-data) */

/* A comma-separated table: line 1 the column names, line 2 their units,
   data from line 3. Columns are found by name. */
typedef struct {
  char Path[2 * BUFSIZE + 1];
  int NColumns;
  char **Names;			/* Column names */
  char **Units;			/* Units of each column */
  char *Known;			/* TRUE once a column has been looked up */
  int NRows;			/* Data rows */
  char ***Fields;		/* Fields[row][column] */
  char **Lines;			/* Storage of the data lines */
} TABLE;

void ReadTable(const char *Path, TABLE *Table);
void FreeTable(TABLE *Table);
int TableColumn(TABLE *Table, const char *Name, const char *Units, int Required);
void TableLayerColumns(TABLE *Table, const char *Base, int LayerSet, int SoilLayers,
                       const char *Units, int Required, int *Columns);
void TableKnown(TABLE *Table, const char *Name);
void TableKnownLayers(TABLE *Table, const char *Base, int LayerSet, int SoilLayers);
int TableIsEmpty(TABLE *Table, int Row, int Column);
const char *TableText(TABLE *Table, int Row, int Column);
float TableFloat(TABLE *Table, int Row, int Column);
double TableDouble(TABLE *Table, int Row, int Column);
int TableInt(TABLE *Table, int Row, int Column);
void TableCheckRows(TABLE *Table, const char *IdColumn, int MinRows, int MaxRows);
void TableWarnUnknown(TABLE *Table);

#endif
