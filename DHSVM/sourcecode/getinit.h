#ifndef GETINIT_H
#define GETINIT_H

#ifndef BUFSIZE
#define BUFSIZE 255
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#define SEPARATOR    '='
#define OPENCOMMENT  '#'
#define OPENSECTION  '['
#define CLOSESECTION ']'

typedef struct _INPUTSTRUCT *LISTPTR;

typedef struct _INPUTSTRUCT {
  char Str[BUFSIZE + 1];
  unsigned char Read;		/* TRUE once the model has read this key */
  LISTPTR Next;
} INPUTSTRUCT;

int ConfigBool(const char *Section, const char *Key, char *Value);
double ConfigDouble(const char *Section, const char *Key, char *Value);
float ConfigFloat(const char *Section, const char *Key, char *Value);
int ConfigInt(const char *Section, const char *Key, char *Value);

int CopyDouble(double *Value, char *Str, const int NValues);
int CopyFloat(float *Value, char *Str, const int NValues);
int CopyInt(int *Value, char *Str, const int NValues);
int CopyUChar(unsigned char *Value, char *Str, const int NValues);

LISTPTR CreateNode(void);

void DeleteList(LISTPTR StartNode);

int GetConfigValue(LISTPTR Input, const char *Section, const char *Key, char *Value);


int IsChoice(const char *Value, const char *Choice);

int IsEmptyStr(char *Str);

unsigned char IsKeyEntryPair(char *Buffer);

unsigned char IsSection(char *Buffer);

unsigned char LocateKey(const char *Key, char *Entry, LISTPTR Input);

LISTPTR LocateSection(const char *Section, LISTPTR Input);

void MakeKeyString(char *Buffer);

void ReadInitFile(char *TemplateFileName, LISTPTR * Input);

void RequireConfigValue(LISTPTR Input, const char *Section, const char *Key, char *Value);

void Strip(char *Buffer);

void WarnUnreadKeys(LISTPTR Input);

#endif
