
#include <ctype.h>
#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"
#include "fileio.h"

/*****************************************************************************
  OpenFile()
*****************************************************************************/
void OpenFile(FILE **FilePtr, char *FileName, char *Mode, uchar OverWrite)
{
  struct stat FileInfo;

  /* if OverWrite is FALSE and the Mode is not read, check whether the file 
     already exists */

  if (!OverWrite && strstr(Mode, "w")) {
    if (!(stat(FileName, &FileInfo)))
      ReportError(FileName, 4);
  }
  if (!(*FilePtr = fopen(FileName, Mode)))
    ReportError(FileName, 3);

}

/*****************************************************************************
  MakePath()

  Composes a path into Buffer (Size >= 2 * BUFSIZE + 1); a path that does
  not fit is an error, never cut short
*****************************************************************************/
void MakePath(char *Buffer, size_t Size, const char *Format, ...)
{
  va_list Args;
  int n;

  va_start(Args, Format);
  n = vsnprintf(Buffer, Size, Format, Args);
  va_end(Args);
  if (n < 0 || (size_t) n >= Size) {
    if (Size > 41)
      Buffer[40] = '\0';
    ReportErrorCode(104, "%s...: path longer than %d characters", Buffer, (int) Size - 1);
  }
}

/*****************************************************************************
  IsDirectory()

  TRUE if Path names an existing directory (trailing separators allowed)
*****************************************************************************/
int IsDirectory(const char *Path)
{
  char Dir[2 * BUFSIZE + 1];
  struct stat Info;
  size_t n;

  if (strlen(Path) >= sizeof(Dir))
    return FALSE;
  strcpy(Dir, Path);
  n = strlen(Dir);
  while (n > 1 && (Dir[n - 1] == '/' || Dir[n - 1] == '\\') &&
         !(n == 3 && Dir[1] == ':'))
    Dir[--n] = '\0';
  if (stat(Dir, &Info) != 0)
    return FALSE;
  return S_ISDIR(Info.st_mode) ? TRUE : FALSE;
}

/*****************************************************************************
  ReadTextLine()

  Reads the next line of a text file into Buffer, without its line end (LF
  or CR LF; the last line may lack it) and, on line 1, without a UTF-8 byte
  order mark. Returns FALSE at the end of the file. A line that does not fit
  into Buffer stops the run with exit code Code.
*****************************************************************************/
int ReadTextLine(FILE *File, char *Buffer, int Size, const char *Path, int LineNo,
                 int Code)
{
  int c;
  int n = 0;

  c = fgetc(File);
  if (c == EOF)
    return FALSE;
  while (c != EOF && c != '\n') {
    if (n >= Size - 1)
      ReportErrorCode(Code, "%s: line %d is longer than %d characters", Path, LineNo,
                      Size - 1);
    Buffer[n++] = (char) c;
    c = fgetc(File);
  }
  if (n > 0 && Buffer[n - 1] == '\r')
    n--;
  Buffer[n] = '\0';
  if (LineNo == 1 && n >= 3 && (unsigned char) Buffer[0] == 0xEF &&
      (unsigned char) Buffer[1] == 0xBB && (unsigned char) Buffer[2] == 0xBF)
    memmove(Buffer, Buffer + 3, (size_t) n - 2);
  return TRUE;
}

/*****************************************************************************
  IsBlankLine()
*****************************************************************************/
int IsBlankLine(const char *Line)
{
  for (; *Line != '\0'; Line++) {
    if (!isspace((unsigned char) *Line))
      return FALSE;
  }
  return TRUE;
}

/*****************************************************************************
  SplitFields()

  Splits a comma-separated line in place; returns the number of fields
  (Fields holds the first MaxFields)
*****************************************************************************/
int SplitFields(char *Line, char **Fields, int MaxFields)
{
  int n = 0;
  char *p = Line;

  for (;;) {
    if (n < MaxFields)
      Fields[n] = p;
    n++;
    p = strchr(p, ',');
    if (p == NULL)
      break;
    *p++ = '\0';
  }
  return n;
}

/*****************************************************************************
  TextToFloat(), TextToInteger()

  Convert a whole field ((float) strtod; base-10 strtol); FALSE if the field
  is empty or holds anything else, or the integer is out of range
*****************************************************************************/
int TextToFloat(const char *Text, float *Value)
{
  char *End;
  double x;

  x = strtod(Text, &End);
  if (End == Text || *End != '\0')
    return FALSE;
  *Value = (float) x;
  return TRUE;
}

int TextToInteger(const char *Text, long *Value)
{
  char *End;

  errno = 0;
  *Value = strtol(Text, &End, 10);
  return (End != Text && *End == '\0' && errno != ERANGE);
}
