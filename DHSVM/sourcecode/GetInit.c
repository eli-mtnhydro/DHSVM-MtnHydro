
/*
 * DESCRIPTION:  Functions are designed to behave more or less the same as the
 *               GetPrivateProfileString and GetPrivateProfileInt functions 
 *               that are part of MFC (the functions here are a little less 
 *               general).  The functions search in a specified file for the 
 *               entry that is associated with a certain key.
 *               The file is assumed to be organized in the same way as
 *               windows .ini files, i.e. the file consists of sections and
 *               key-value pairs.  In addition there can be comments.  Thus
 *               input files have the format:
 *
 *               # comment
 *               [section]
 *               key=value          # comment
 */
#define _CRT_SECURE_NO_DEPRECATE
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "DHSVMerror.h"
#include "fileio.h"
#include "getinit.h"

/*#####################################################################################
 This function is used to find the matching key word in the input file for the "key" 
 specified in the fucntion: InitVegTable( )
 #####################################################################################*/
unsigned char LocateKey(const char *Key, char *Entry, LISTPTR Input)
{
  unsigned char Found = FALSE;
  char Buffer[BUFSIZE + 1];
  char KeyBuffer[BUFSIZE + 1];
  char EntryBuffer[BUFSIZE + 1];
  char *StrPtr = NULL;

  /* Find a key in the current section */

  if (Input) {
    strncpy(Buffer, Input->Str, BUFSIZE + 1);
    while (!IsSection(Buffer)) {

      /* Check whether the current line contains a key-entry pair */

      if (IsKeyEntryPair(Buffer)) {
	StrPtr = strchr(Buffer, SEPARATOR);
	*StrPtr = '\0';
	/* strcpy(KeyBuffer, Buffer); */
	memmove(KeyBuffer, Buffer, strlen(Buffer)+1);
	++StrPtr;
    /* strcpy(EntryBuffer, StrPtr); */
	memmove(EntryBuffer, StrPtr, strlen(StrPtr)+1);
	Strip(KeyBuffer);
	MakeKeyString(KeyBuffer);
	if (strcmp(Key, KeyBuffer) == 0) {
	  Found = TRUE;
	  Input->Read = TRUE;
	  Strip(EntryBuffer);
	  /* strcpy(Entry, EntryBuffer); */
	  memmove(Entry, EntryBuffer, strlen(EntryBuffer)+1);
	  break;
	}
      }
      /* Get the next line */
      Input = Input->Next;
      if (Input)
	strncpy(Buffer, Input->Str, BUFSIZE + 1);
      else
	break;
    }
  }

  return Found;
}
/*#####################################################################################*/
LISTPTR LocateSection(const char *Section, LISTPTR Input)
{
  char Buffer[BUFSIZE + 1];
  char *StartPtr = NULL;
  char *EndPtr = NULL;

  while (Input) {
    strncpy(Buffer, Input->Str, BUFSIZE + 1);
    if (IsSection(Buffer)) {
      if (Buffer[0] == OPENSECTION) {
	StartPtr = &Buffer[1];
	EndPtr = strchr(Buffer, CLOSESECTION);
	*EndPtr = '\0';
	/* strcpy(Buffer, StartPtr);*/
	memmove(Buffer, StartPtr, strlen(StartPtr)+1);
	Strip(Buffer);
	MakeKeyString(Buffer);
	if (strcmp(Section, Buffer) == 0) {
	  Input = Input->Next;
	  break;
	}
      }
    }
    Input = Input->Next;
  }

  return Input;
}
/*#####################################################################################*/
unsigned char IsKeyEntryPair(char *Buffer)
{
  char *StrSeparator = NULL;

  StrSeparator = strchr(Buffer, SEPARATOR);
  if (StrSeparator == NULL)
    return FALSE;

  return TRUE;
}
/*#####################################################################################*/
unsigned char IsSection(char *Buffer)
{
  char *StrEndSection = NULL;
  char *StrStartComment = NULL;

  if (Buffer[0] != OPENSECTION)
    return FALSE;

  StrStartComment = strchr(Buffer, OPENCOMMENT);

  StrEndSection = strchr(Buffer, CLOSESECTION);
  if (StrEndSection == NULL)
    return FALSE;

  if (StrStartComment != NULL && StrStartComment < StrEndSection)
    return FALSE;

  return TRUE;
}
/*#####################################################################################*/
void Strip(char *Buffer)
{
  char *StrEnd = NULL;
  char *StrStart = Buffer;

  /* remove leading whitespace */
  while (*StrStart != '\0' && isspace((int) *StrStart))
    ++StrStart;


  /* remove comment */
  StrEnd = strchr(Buffer, OPENCOMMENT);
  if (StrEnd != NULL)
    *StrEnd = '\0';

  /* remove trailing whitespace */

  StrEnd = &Buffer[strlen(Buffer)];
  --StrEnd;
  while (StrEnd >= StrStart && isspace((int) *StrEnd)) {
    *StrEnd = '\0';
    --StrEnd;
  }
  /* strcpy(Buffer, StrStart); */
  /* commented out by Ning b/c source and destination overlap */
  memmove(Buffer, StrStart, strlen(StrStart)+1);
}
/*#####################################################################################
 Uppercases a key or section name and turns every run of white space (spaces
 and tabs) into one space
 #####################################################################################*/
void MakeKeyString(char *Buffer)
{
  char Str[BUFSIZE + 1];
  char *PtrStr = Str;
  char *PtrBuffer = Buffer;

  while (*PtrBuffer != '\0') {
    if (isspace((int) *PtrBuffer)) {
      *PtrStr++ = ' ';
      while (isspace((int) *PtrBuffer))
        PtrBuffer++;
    }
    else
      *PtrStr++ = (char) toupper(*PtrBuffer++);
  }
  *PtrStr = '\0';
  strcpy(Buffer, Str);
}
/*#####################################################################################*/
int CopyDouble(double *Value, char *Str, const int NValues)
{
  char *EndPtr = NULL;
  int i;

  for (i = 0; i < NValues; i++) {
    Value[i] = strtod(Str, &EndPtr);
    if (EndPtr == Str)
      return FALSE;
    Str = EndPtr;
  }

  if (EndPtr && *EndPtr != '\0')
    return FALSE;

  return TRUE;
}
/*#####################################################################################*/
int CopyFloat(float *Value, char *Str, const int NValues)
{
  char *EndPtr = NULL;
  int i;

  for (i = 0; i < NValues; i++) {
    Value[i] = (float) strtod(Str, &EndPtr);
    if (EndPtr == Str)
      return FALSE;
    Str = EndPtr;
  }

  if (EndPtr && *EndPtr != '\0')
    return FALSE;

  return TRUE;
}
/*#####################################################################################*/
int CopyInt(int *Value, char *Str, const int NValues)
{
  char *EndPtr = NULL;
  long x;
  int i;

  for (i = 0; i < NValues; i++) {
    errno = 0;
    x = strtol(Str, &EndPtr, 10);
    if (EndPtr == Str || errno == ERANGE || x < INT_MIN || x > INT_MAX)
      return FALSE;
    Value[i] = (int) x;
    Str = EndPtr;
  }

  if (EndPtr && *EndPtr != '\0')
    return FALSE;

  return TRUE;
}
/*#####################################################################################*/
int CopyUChar(unsigned char *Value, char *Str, const int NValues)
{
  char *EndPtr = NULL;
  long x;
  int i;

  for (i = 0; i < NValues; i++) {
    errno = 0;
    x = strtol(Str, &EndPtr, 10);
    if (EndPtr == Str || errno == ERANGE || x < 0 || x > UCHAR_MAX)
      return FALSE;
    Value[i] = (unsigned char) x;
    Str = EndPtr;
  }

  if (EndPtr && *EndPtr != '\0')
    return FALSE;

  return TRUE;
}
/*#####################################################################################*/
int IsEmptyStr(char *Str)
{
  if (Str == NULL)
    return TRUE;
  if (Str[0] == '\0')
    return TRUE;
  return FALSE;
}
/*#####################################################################################
 Reads the config into a list of its section and key = value lines. Each
 line is read whole (LF or CRLF; the last line also without a line end). A
 line longer than BUFSIZE characters (comments included), or a line that is
 not blank, a comment, a section or a key = value line, stops the run; so
 does a section other than the six of a config, a repeated section, or a key
 line before the first section.
 #####################################################################################*/
void ReadInitFile(char *TemplateFileName, LISTPTR * Input)
{
  const char *Sections[] = { "RUN", "GRID", "OPTIONS", "INPUTS", "PARAMETERS", "OUTPUT" };
  int Seen[6] = { FALSE, FALSE, FALSE, FALSE, FALSE, FALSE };
  int InSection = FALSE;		/* TRUE after the first section line */
  FILE *InFile = NULL;			/* File with input information */
  char Line[BUFSIZE + 3];		/* Line as read */
  char Buffer[BUFSIZE + 1];		/* Line without comment and outer blanks */
  char Name[BUFSIZE + 1];		/* Section name */
  int c;						/* Character read */
  int Last;						/* Last character of the line */
  int Length;					/* Length of the line */
  int LineNo = 0;				/* Line number */
  int i;
  LISTPTR Current = NULL;		/* pointer to current node in list */
  LISTPTR Head = NULL;			/* pointer to the start of the list */

  if (!(InFile = fopen(TemplateFileName, "rb")))
    ReportErrorCode(101, "%s: cannot open", TemplateFileName);

  c = fgetc(InFile);
  while (c != EOF) {
    LineNo++;
    Length = 0;
    Last = 0;
    while (c != EOF && c != '\n') {
      if (Length < BUFSIZE + 2)
        Line[Length] = (char) c;
      Length++;
      Last = c;
      c = fgetc(InFile);
    }
    if (c == '\n')
      c = fgetc(InFile);
    if (Last == '\r')
      Length--;
    if (Length > BUFSIZE)
      ReportErrorCode(102, "%s: line %d is longer than %d characters",
                      TemplateFileName, LineNo, BUFSIZE);
    Line[Length] = '\0';

    strcpy(Buffer, Line);
    Strip(Buffer);
    if (Buffer[0] == '\0')
      continue;
    if (!IsSection(Buffer) && !IsKeyEntryPair(Buffer))
      ReportErrorCode(102, "%s: line %d is not a section, a key = value line or a comment",
                      TemplateFileName, LineNo);
    if (IsSection(Buffer)) {
      strcpy(Name, Buffer + 1);
      *strchr(Name, CLOSESECTION) = '\0';
      Strip(Name);
      MakeKeyString(Name);
      for (i = 0; i < 6 && strcmp(Name, Sections[i]) != 0; i++)
        ;
      if (i == 6)
        ReportErrorCode(102, "%s: line %d: unknown section [%s]", TemplateFileName, LineNo, Name);
      if (Seen[i])
        ReportErrorCode(102, "%s: line %d: section [%s] repeated", TemplateFileName, LineNo, Name);
      Seen[i] = TRUE;
      InSection = TRUE;
    }
    else if (!InSection)
      ReportErrorCode(102, "%s: line %d: key before the first section", TemplateFileName, LineNo);
    if (Head == NULL) {
      Head = CreateNode();
      Current = Head;
      *Input = Head;
    }
    else {
      Current->Next = CreateNode();
      Current = Current->Next;
    }
    strncpy(Current->Str, Buffer, BUFSIZE + 1);
  }

  fclose(InFile);

  return;
}
/*#####################################################################################*/
LISTPTR CreateNode(void)
{
  LISTPTR NewNode = NULL;

  NewNode = calloc(1, sizeof(INPUTSTRUCT));
  if (NewNode == NULL)
    ReportError("CreateNode", 1);
  NewNode->Next = NULL;

  return NewNode;
}
/*#####################################################################################*/
void DeleteList(LISTPTR Head)
{
  LISTPTR Current = NULL;

  Current = Head;
  while (Current != NULL) {
    Head = Head->Next;
    free(Current);
    Current = Head;
  }

  return;
}
/*#####################################################################################
 Warns about every key in the input file that was never read.
 Call after initialization, before DeleteList().
 #####################################################################################*/
void WarnUnreadKeys(LISTPTR Input)
{
  char Section[BUFSIZE + 1] = "";
  char Key[BUFSIZE + 1];
  char *StrPtr = NULL;
  int NUnread = 0;

  for (; Input != NULL; Input = Input->Next) {
    if (IsSection(Input->Str)) {
      strncpy(Section, Input->Str, BUFSIZE + 1);
      continue;
    }
    if (!Input->Read && IsKeyEntryPair(Input->Str)) {
      strncpy(Key, Input->Str, BUFSIZE + 1);
      StrPtr = strchr(Key, SEPARATOR);
      *StrPtr = '\0';
      Strip(Key);
      Warning("config key not used: %s %s", Section, Key);
      NUnread++;
    }
  }
  if (NUnread > 0)
    Warning("%d config keys not used (misspelled, duplicated, in the wrong section, or for an "
            "option that is off)", NUnread);
}
/*#####################################################################################*/
/*#####################################################################################
 GetConfigValue()
 Looks up a key in a section, both named as written (e.g. "RUN", "Time Step");
 returns TRUE and the value if the key is present (its first occurrence)
 #####################################################################################*/
int GetConfigValue(LISTPTR Input, const char *Section, const char *Key, char *Value)
{
  char SectionKey[BUFSIZE + 1];
  char KeyKey[BUFSIZE + 1];
  LISTPTR SectionHead = NULL;

  strncpy(SectionKey, Section, BUFSIZE);
  SectionKey[BUFSIZE] = '\0';
  MakeKeyString(SectionKey);
  strncpy(KeyKey, Key, BUFSIZE);
  KeyKey[BUFSIZE] = '\0';
  MakeKeyString(KeyKey);

  Value[0] = '\0';
  if ((SectionHead = LocateSection(SectionKey, Input)) == NULL)
    return FALSE;
  return LocateKey(KeyKey, Value, SectionHead);
}
/*#####################################################################################
 RequireConfigValue()
 As GetConfigValue(), for a key that must be present with a value
 #####################################################################################*/
void RequireConfigValue(LISTPTR Input, const char *Section, const char *Key, char *Value)
{
  if (!GetConfigValue(Input, Section, Key, Value))
    ReportErrorCode(103, "[%s] %s is missing", Section, Key);
  if (IsEmptyStr(Value))
    ReportErrorCode(103, "[%s] %s is empty", Section, Key);
}
/*#####################################################################################
 ConfigFloat(), ConfigDouble(), ConfigInt(), ConfigBool()
 Convert a whole config value ((float) strtod or strtod to a number,
 base-10 strtol in the int range, TRUE or FALSE in any case); anything else
 stops the run naming the key
 #####################################################################################*/
float ConfigFloat(const char *Section, const char *Key, char *Value)
{
  float x;

  if (!CopyFloat(&x, Value, 1))
    ReportErrorCode(103, "[%s] %s = %s: expected a number", Section, Key, Value);
  return x;
}

double ConfigDouble(const char *Section, const char *Key, char *Value)
{
  double x;

  if (!CopyDouble(&x, Value, 1))
    ReportErrorCode(103, "[%s] %s = %s: expected a number", Section, Key, Value);
  return x;
}

int ConfigInt(const char *Section, const char *Key, char *Value)
{
  char *EndPtr = NULL;
  long x;

  errno = 0;
  x = strtol(Value, &EndPtr, 10);
  if (EndPtr == Value || *EndPtr != '\0' || errno == ERANGE || x < INT_MIN || x > INT_MAX)
    ReportErrorCode(103, "[%s] %s = %s: expected an integer", Section, Key, Value);
  return (int) x;
}

int ConfigBool(const char *Section, const char *Key, char *Value)
{
  if (IsChoice(Value, "TRUE"))
    return TRUE;
  if (IsChoice(Value, "FALSE"))
    return FALSE;
  ReportErrorCode(103, "[%s] %s = %s: expected TRUE or FALSE", Section, Key, Value);
  return FALSE;
}
/*#####################################################################################
 IsChoice()
 TRUE if Value equals Choice ignoring case
 #####################################################################################*/
int IsChoice(const char *Value, const char *Choice)
{
  while (*Value != '\0' && *Choice != '\0') {
    if (toupper((unsigned char) *Value) != toupper((unsigned char) *Choice))
      return FALSE;
    Value++;
    Choice++;
  }
  return (*Value == '\0' && *Choice == '\0');
}
