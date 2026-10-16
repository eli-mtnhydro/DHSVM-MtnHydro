
#ifndef FILEIO_H
#define FILEIO_H

#include "data.h"

/* generic file functions */
void OpenFile(FILE **FilePtr, char *FileName, char *Mode,
	      unsigned char OverWrite);

int IsBlankLine(const char *Line);

int IsDirectory(const char *Path);

void MakePath(char *Buffer, size_t Size, const char *Format, ...);

int ReadTextLine(FILE *File, char *Buffer, int Size, const char *Path, int LineNo,
                 int Code);

int SplitFields(char *Line, char **Fields, int MaxFields);

int TextToFloat(const char *Text, float *Value);

int TextToInteger(const char *Text, long *Value);

#endif
