
#ifndef FILEIO_H
#define FILEIO_H

#include "data.h"

void InitFileIO(void);

/* global file extension string */
extern char fileext[];

/* function pointers for 2D file IO */

void CreateMapFile(char *FileName, char *FileLabel, MAPSIZE *Map);

int Read2DMatrix(char *FileName, void *Matrix, int NumberType, 
                 MAPSIZE *Map, int NDataSet, char *VarName, int index);

int Write2DMatrix(char *FileName, void *Matrix, int NumberType, 
                  MAPSIZE *Map, MAPDUMP *DMap, int index);


/* generic file functions */
void OpenFile(FILE **FilePtr, char *FileName, char *Mode,
	      unsigned char OverWrite);

#endif
