
#ifndef FIFOBIN_H
#define FIFOBIN_H

void CreateMapFileBin(char *FileName, ...);
int Read2DMatrixBin(char *FileName, void *Matrix, int NumberType, int NY,
		    int NX, int NDataSet, ...); 
int Write2DMatrixBin(char *FileName, void *Matrix, int NumberType, int NY,
		     int NX, ...);

#endif
