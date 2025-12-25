
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fileio.h"
#include "fifobin.h"
#include "DHSVMerror.h"

void (*CreateMapFileFmt) (char *FileName, ...);
int (*Read2DMatrixFmt) (char *FileName, void *Matrix, int NumberType, int NY, int NX, int NDataSet, ...);
int (*Write2DMatrixFmt) (char *FileName, void *Matrix, int NumberType, int NY, int NX, ...);

/*******************************************************************************
  Function name: InitFileIO()

  Purpose      : Initialize function pointers for file IO
 
  Comments     :

  This function sets the file pointers for file I/O functions to the
  functions that implement the desired file format.  By using file pointers,
  the main routines do not need to be changed if a new file format is to be
  supported in the future (at least that is what I am trying to accomplish).
  The only thing that will need to be done is write the necessary I/O
  functions for the new file format, and add the additional options to this 
  initialization routine.

  ONLY SUPPORTS BINARY GOING FORWARD

  Information is stored in all files in the following way:
  fastest varying dimension: X (West to East)
  next fastest dimension   : Y (North to South)
          .                : Variable (if more than one)
  slowest varying dimension: Time (if more than one timestep)
  All the information is written out or read in for one timestep at a time.
            
  The following terminology is used (this is not based on anything
  "official", and may (will) not correspond to either dictionary or
  scientific definitions, but there is only so much one can do with variable
  and function names).  Some of the distinctions being made may not always
  seem to make sense for certain file formats, and partly result from the
  fact that the first version of the model used the HDF file format:
   - 2DMatrix :
       a map layer with X and Y dimension.  
   - 2DImage :
       a map layer with X and Y dimension in which the data are stored as 
       unsigned char, i.e. values in the interval [a, b] are mapped to numbers 
       in the range [0, 255].

*******************************************************************************/
void InitFileIO(void)
{
  printf("Initializing file IO\n");

  strcpy(fileext, ".bin");
  CreateMapFileFmt = CreateMapFileBin;
  Read2DMatrixFmt = Read2DMatrixBin;
  Write2DMatrixFmt = Write2DMatrixBin;
}

/******************************************************************************/
/*                            CreateMapFile                                   */
/******************************************************************************/
void CreateMapFile(char *FileName, char *FileLabel, MAPSIZE *Map)
{
  CreateMapFileFmt(FileName, FileLabel, Map);
}

/******************************************************************************/
/*                              Read2DMatrix                                  */
/******************************************************************************/

int Read2DMatrix(char *FileName, void *Matrix, int NumberType, MAPSIZE *Map,
                 int NDataSet, char *VarName, int index)
{
  int result;
  result = Read2DMatrixFmt(FileName, Matrix, NumberType,
                           Map->NY, Map->NX, NDataSet, VarName, index);
  return result;
}

/******************************************************************************/
/*                              Write2DMatrix                                  */
/******************************************************************************/
int Write2DMatrix(char *FileName, void *Matrix, int NumberType, MAPSIZE *Map,
                  MAPDUMP *DMap, int index)
{
  int result;
  result = Write2DMatrixFmt(FileName, Matrix, NumberType, 
                            Map->NY, Map->NX, DMap, index);
  return result;
}
