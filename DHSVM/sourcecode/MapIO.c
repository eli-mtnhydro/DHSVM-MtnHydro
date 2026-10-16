/*
 * DESCRIPTION:  Map files (.map). A map file is a text header of HeaderBytes
 *               bytes (a multiple of 4096) followed by the values:
 *
 *               DHSVM-MAP 1 <HeaderBytes>
 *               Key = Value         (Name, Rows, Columns, Layers, NumberType,
 *               ...                  ByteOrder, CellSize, North, West, NoData,
 *               Layer <k> = <name>   Units, Time, Writer, Layer 1..Layers)
 *               End
 *               <spaces up to HeaderBytes - 1, then a newline>
 *
 *               Values are little-endian, layer after layer, each layer row
 *               after row from north to south, each row from west to east.
 *               NumberType is BYTE (unsigned char, no-data 255) or FLOAT
 *               (float, no-data -9999).
 */
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "getinit.h"
#include "sizeofnt.h"
#include "mapio.h"

#define MAP_HEADER_BLOCK 4096

static const char *MonthTokens[12] = {
  "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
};

/* Header keys other than "Layer <k>"; the first NREQUIREDKEYS are required */
static const char *HeaderKeys[] = {
  "Name", "Rows", "Columns", "Layers", "NumberType", "ByteOrder", "CellSize",
  "North", "West", "NoData", "Units", "Time", "Writer"
};
#define NHEADERKEYS 13
#define NREQUIREDKEYS 11

/*****************************************************************************
  Local helpers
*****************************************************************************/
static int HostIsLittleEndian(void)
{
  unsigned short One = 1;
  return *((unsigned char *) &One) == 1;
}

static int SeekTo(FILE *File, long long Offset)
{
#ifdef _WIN32
  return _fseeki64(File, Offset, SEEK_SET);
#else
  return fseeko(File, (off_t) Offset, SEEK_SET);
#endif
}

static long long SizeOfFile(FILE *File)
{
#ifdef _WIN32
  if (_fseeki64(File, 0, SEEK_END) != 0)
    return -1;
  return _ftelli64(File);
#else
  if (fseeko(File, 0, SEEK_END) != 0)
    return -1;
  return (long long) ftello(File);
#endif
}

static const char *TypeName(int NumberType)
{
  return (NumberType == NC_BYTE) ? "BYTE" : "FLOAT";
}

/* A positive decimal integer without sign or leading zero, else -1 */
static long long PositiveInteger(const char *Str)
{
  long long x = 0;

  if (*Str < '1' || *Str > '9')
    return -1;
  for (; *Str != '\0'; Str++) {
    if (!isdigit((unsigned char) *Str) || x > 100000000000LL)
      return -1;
    x = 10 * x + (*Str - '0');
  }
  return x;
}

/* Index of Key among the header lines, -1 if absent */
static int FindKey(char **Keys, int NLines, const char *Key)
{
  int i;

  for (i = 0; i < NLines; i++)
    if (strcmp(Keys[i], Key) == 0)
      return i;
  return -1;
}

/* Number written with the fewest digits that read back as the same double */
static void HeaderNumber(char *Buffer, double x)
{
  int Digits;

  if (x == floor(x) && fabs(x) < 1e15) {
    sprintf(Buffer, "%.0f", x);
    return;
  }
  for (Digits = 1; Digits <= 17; Digits++) {
    sprintf(Buffer, "%.*g", Digits, x);
    if (strtod(Buffer, NULL) == x)
      return;
  }
}

/*****************************************************************************
  MapLayerCount(), MapLayerNames()

  Number and names of the layers of a map with a layer set: one layer named
  like the map; soil layers <Name>.1 .. <Name>.<L>, <Name>.Deep; months
  <Name>.Jan .. <Name>.Dec; azimuths <Name>.Az000 .. <Name>.Az350 (every
  10 degrees clockwise from north)
*****************************************************************************/
int MapLayerCount(int LayerSet, int SoilLayers)
{
  switch (LayerSet) {
  case MAP_SOIL:
    return SoilLayers + 1;
  case MAP_MONTH:
    return 12;
  case MAP_AZIMUTH:
    return MAP_NAZIMUTHS;
  default:
    return 1;
  }
}

void MapLayerNames(const char *Name, int LayerSet, int SoilLayers,
                   char (*LayerNames)[MAP_NAMESIZE])
{
  int k;

  switch (LayerSet) {
  case MAP_SOIL:
    for (k = 0; k < SoilLayers; k++)
      snprintf(LayerNames[k], MAP_NAMESIZE, "%s.%d", Name, k + 1);
    snprintf(LayerNames[SoilLayers], MAP_NAMESIZE, "%s.Deep", Name);
    break;
  case MAP_MONTH:
    for (k = 0; k < 12; k++)
      snprintf(LayerNames[k], MAP_NAMESIZE, "%s.%s", Name, MonthTokens[k]);
    break;
  case MAP_AZIMUTH:
    for (k = 0; k < MAP_NAZIMUTHS; k++)
      snprintf(LayerNames[k], MAP_NAMESIZE, "%s.Az%03d", Name, k * 360 / MAP_NAZIMUTHS);
    break;
  default:
    snprintf(LayerNames[0], MAP_NAMESIZE, "%s", Name);
    break;
  }
}

/*****************************************************************************
  CheckMap()

  Opens a map file and checks, in this order: the tag, the version and the
  header length; the header keys (all required keys, no unknown or repeated
  key, an End line); the grid against [GRID]; Name, NumberType, ByteOrder,
  Units, NoData and Layers against the expected ones; the layer names; Time
  when Time is not NULL; the file size. Stops at the first failure with exit
  code Code (105 for input maps, 112 for snapshots); a file that cannot be
  opened is error 104. The file stays open for ReadMapLayer().
*****************************************************************************/
void CheckMap(MAPFILE *File, const char *Path, const char *Name, int NumberType,
              const char *Units, int NLayers, char (*LayerNames)[MAP_NAMESIZE],
              const char *Time, MAPSIZE *Map, int Code)
{
  FILE *InFile;
  char Line1[MAP_NAMESIZE + 1];
  char KeyName[MAP_NAMESIZE];
  char *Header;
  char *Line;
  char *Next;
  char *Eq;
  char *EndPtr;
  char **Keys;
  char **Values;
  int NLines = 0;
  int EndFound = FALSE;
  int i, j, k;
  size_t n;
  long long HeaderBytes = 0;
  long long FileLayers;
  long long Size;
  long long Expected;
  double x;
  double NoData;

  if (!HostIsLittleEndian())
    ReportErrorCode(Code, "%s: maps are little-endian; this computer is not", Path);
  if (!(InFile = fopen(Path, "rb")))
    ReportErrorCode(104, "%s: cannot open", Path);

  /* Step 1: tag, version and header length */
  n = fread(Line1, 1, MAP_NAMESIZE, InFile);
  Line1[n] = '\0';
  Next = (char *) memchr(Line1, '\n', n);
  if (Next == NULL || memchr(Line1, '\0', (size_t) (Next - Line1)) != NULL)
    ReportErrorCode(Code, "%s: not a DHSVM-MAP 1 file", Path);
  *Next = '\0';
  if (strncmp(Line1, "DHSVM-MAP 1 ", 12) != 0 ||
      (HeaderBytes = PositiveInteger(Line1 + 12)) < 0 ||
      HeaderBytes % MAP_HEADER_BLOCK != 0)
    ReportErrorCode(Code, "%s: not a DHSVM-MAP 1 file", Path);
  Size = SizeOfFile(InFile);
  if (Size < HeaderBytes)
    ReportErrorCode(Code, "%s: %.0f bytes, expected at least the %.0f header bytes (truncated)",
                    Path, (double) Size, (double) HeaderBytes);

  if (!(Header = (char *) malloc((size_t) HeaderBytes + 1)))
    ReportError("CheckMap", 1);
  if (SeekTo(InFile, 0) != 0 ||
      fread(Header, 1, (size_t) HeaderBytes, InFile) != (size_t) HeaderBytes)
    ReportErrorCode(Code, "%s: cannot read the header", Path);
  Header[HeaderBytes] = '\0';
  if (strlen(Header) != (size_t) HeaderBytes)
    ReportErrorCode(Code, "%s: the header holds a NUL byte", Path);

  /* Step 2: header keys */
  if (!(Keys = (char **) calloc((size_t) HeaderBytes / 4 + 1, sizeof(char *))) ||
      !(Values = (char **) calloc((size_t) HeaderBytes / 4 + 1, sizeof(char *))))
    ReportError("CheckMap", 1);
  Line = strchr(Header, '\n') + 1;
  while ((Next = strchr(Line, '\n')) != NULL) {
    *Next = '\0';
    if (strcmp(Line, "End") == 0) {
      EndFound = TRUE;
      break;
    }
    if ((Eq = strstr(Line, " = ")) == NULL)
      ReportErrorCode(Code, "%s: header line \"%s\" is not Key = Value", Path, Line);
    *Eq = '\0';
    Keys[NLines] = Line;
    Values[NLines] = Eq + 3;
    NLines++;
    Line = Next + 1;
  }
  if (!EndFound)
    ReportErrorCode(Code, "%s: the header has no End line", Path);
  for (i = 0; i < NLines; i++) {
    for (j = 0; j < NHEADERKEYS; j++)
      if (strcmp(Keys[i], HeaderKeys[j]) == 0)
        break;
    if (j == NHEADERKEYS &&
        (strncmp(Keys[i], "Layer ", 6) != 0 || PositiveInteger(Keys[i] + 6) < 0))
      ReportErrorCode(Code, "%s: unknown header key %s", Path, Keys[i]);
    if (FindKey(Keys, i, Keys[i]) >= 0)
      ReportErrorCode(Code, "%s: header key %s repeated", Path, Keys[i]);
  }
  for (j = 0; j < NREQUIREDKEYS; j++)
    if (FindKey(Keys, NLines, HeaderKeys[j]) < 0)
      ReportErrorCode(Code, "%s: header key %s missing", Path, HeaderKeys[j]);
  if ((FileLayers = PositiveInteger(Values[FindKey(Keys, NLines, "Layers")])) < 0 ||
      FileLayers > 100000)
    ReportErrorCode(Code, "%s: Layers = %s is not a positive integer", Path,
                    Values[FindKey(Keys, NLines, "Layers")]);
  for (i = 0; i < NLines; i++)
    if (strncmp(Keys[i], "Layer ", 6) == 0 && PositiveInteger(Keys[i] + 6) > FileLayers)
      ReportErrorCode(Code, "%s: unknown header key %s", Path, Keys[i]);
  for (k = 1; k <= FileLayers; k++) {
    snprintf(KeyName, sizeof(KeyName), "Layer %d", k);
    if (FindKey(Keys, NLines, KeyName) < 0)
      ReportErrorCode(Code, "%s: header key %s missing", Path, KeyName);
  }

  /* Step 3: grid */
  i = FindKey(Keys, NLines, "Rows");
  if (PositiveInteger(Values[i]) != Map->NY)
    ReportErrorCode(Code, "%s: Rows = %s, [GRID] has %d", Path, Values[i], Map->NY);
  i = FindKey(Keys, NLines, "Columns");
  if (PositiveInteger(Values[i]) != Map->NX)
    ReportErrorCode(Code, "%s: Columns = %s, [GRID] has %d", Path, Values[i], Map->NX);
  i = FindKey(Keys, NLines, "CellSize");
  x = strtod(Values[i], &EndPtr);
  if (EndPtr == Values[i] || *EndPtr != '\0' || !(fabs(x - Map->DX) <= 0.001))
    ReportErrorCode(Code, "%s: CellSize = %s, [GRID] has %g", Path, Values[i], Map->DX);
  i = FindKey(Keys, NLines, "North");
  x = strtod(Values[i], &EndPtr);
  if (EndPtr == Values[i] || *EndPtr != '\0' || !(fabs(x - Map->Yorig) <= 0.001))
    ReportErrorCode(Code, "%s: North = %s, [GRID] has %.17g", Path, Values[i], Map->Yorig);
  i = FindKey(Keys, NLines, "West");
  x = strtod(Values[i], &EndPtr);
  if (EndPtr == Values[i] || *EndPtr != '\0' || !(fabs(x - Map->Xorig) <= 0.001))
    ReportErrorCode(Code, "%s: West = %s, [GRID] has %.17g", Path, Values[i], Map->Xorig);

  /* Step 4: name, number type, units, no-data value, layer count */
  i = FindKey(Keys, NLines, "Name");
  if (strcmp(Values[i], Name) != 0)
    ReportErrorCode(Code, "%s: Name is %s, expected %s", Path, Values[i], Name);
  i = FindKey(Keys, NLines, "NumberType");
  if (strcmp(Values[i], TypeName(NumberType)) != 0)
    ReportErrorCode(Code, "%s: NumberType is %s, expected %s", Path, Values[i],
                    TypeName(NumberType));
  i = FindKey(Keys, NLines, "ByteOrder");
  if (strcmp(Values[i], "LITTLE") != 0)
    ReportErrorCode(Code, "%s: ByteOrder is %s, expected LITTLE", Path, Values[i]);
  i = FindKey(Keys, NLines, "Units");
  if (strcmp(Values[i], Units) != 0)
    ReportErrorCode(Code, "%s: Units is %s, expected %s", Path, Values[i], Units);
  i = FindKey(Keys, NLines, "NoData");
  NoData = (NumberType == NC_BYTE) ? 255.0 : -9999.0;
  x = strtod(Values[i], &EndPtr);
  if (EndPtr == Values[i] || *EndPtr != '\0' || x != NoData)
    ReportErrorCode(Code, "%s: NoData is %s, expected %g", Path, Values[i], NoData);
  if (FileLayers != NLayers)
    ReportErrorCode(Code, "%s: Layers is %d, expected %d", Path, (int) FileLayers, NLayers);
  if (Time != NULL) {
    if ((i = FindKey(Keys, NLines, "Time")) < 0)
      ReportErrorCode(Code, "%s: header key Time missing", Path);
    if (strcmp(Values[i], Time) != 0)
      ReportErrorCode(Code, "%s: Time is %s, expected %s", Path, Values[i], Time);
  }

  /* Step 5: layer names */
  for (k = 1; k <= NLayers; k++) {
    snprintf(KeyName, sizeof(KeyName), "Layer %d", k);
    i = FindKey(Keys, NLines, KeyName);
    if (strcmp(Values[i], LayerNames[k - 1]) != 0)
      ReportErrorCode(Code, "%s: layer %d is named %s, expected %s", Path, k,
                      Values[i], LayerNames[k - 1]);
  }

  /* Step 6: file size */
  Expected = HeaderBytes + (long long) NLayers * Map->NY * Map->NX *
    (long long) SizeOfNumberType(NumberType);
  if (Size != Expected)
    ReportErrorCode(Code, "%s: %.0f bytes, expected %.0f (truncated or extra data)",
                    Path, (double) Size, (double) Expected);

  free(Keys);
  free(Values);
  free(Header);

  snprintf(File->Path, sizeof(File->Path), "%s", Path);
  File->Ptr = InFile;
  File->HeaderBytes = HeaderBytes;
  File->NumberType = NumberType;
  File->Layers = NLayers;
  File->ErrorCode = Code;
}

/*****************************************************************************
  OpenInputMap()

  Opens and checks an input map whose layers follow a layer set
*****************************************************************************/
void OpenInputMap(MAPFILE *File, const char *Path, const char *Name,
                  int NumberType, const char *Units, int LayerSet, MAPSIZE *Map)
{
  int NLayers = MapLayerCount(LayerSet, Map->SoilLayers);
  char (*LayerNames)[MAP_NAMESIZE];

  if (!(LayerNames = calloc(NLayers, MAP_NAMESIZE)))
    ReportError("OpenInputMap", 1);
  MapLayerNames(Name, LayerSet, Map->SoilLayers, LayerNames);
  CheckMap(File, Path, Name, NumberType, Units, NLayers, LayerNames, NULL, Map, 105);
  free(LayerNames);
}

/*****************************************************************************
  ReadMapLayer()

  Reads layer Layer (from 0) of an open map into Values (NY x NX values of
  the map's number type)
*****************************************************************************/
void ReadMapLayer(MAPFILE *File, int Layer, void *Values, MAPSIZE *Map)
{
  size_t NCells = (size_t) Map->NY * Map->NX;
  size_t Bytes = SizeOfNumberType(File->NumberType);

  if (Layer < 0 || Layer >= File->Layers ||
      SeekTo(File->Ptr, File->HeaderBytes + (long long) Layer * NCells * Bytes) != 0 ||
      fread(Values, Bytes, NCells, File->Ptr) != NCells)
    ReportErrorCode(File->ErrorCode, "%s: cannot read layer %d", File->Path, Layer + 1);
}

/*****************************************************************************
  CloseMap()
*****************************************************************************/
void CloseMap(MAPFILE *File)
{
  if (File->Ptr != NULL)
    fclose(File->Ptr);
  File->Ptr = NULL;
}

/*****************************************************************************
  ReadInputMap()

  Reads every layer of the input map named by the [INPUTS] key Name (the
  map's name) into Values, one NY x NX block per layer
*****************************************************************************/
void ReadInputMap(LISTPTR Input, const char *Name, int NumberType,
                  const char *Units, int LayerSet, MAPSIZE *Map, void *Values)
{
  char Path[BUFSIZE + 1];
  MAPFILE File;
  size_t Block = (size_t) Map->NY * Map->NX * SizeOfNumberType(NumberType);
  int k;

  RequireConfigValue(Input, "INPUTS", Name, Path);
  OpenInputMap(&File, Path, Name, NumberType, Units, LayerSet, Map);
  for (k = 0; k < File.Layers; k++)
    ReadMapLayer(&File, k, (char *) Values + k * Block, Map);
  CloseMap(&File);
}

/*****************************************************************************
  CreateMap()

  Creates a FLOAT map file and writes its header (Time is left out when
  NULL); the NLayers layers follow with WriteMapLayer()
*****************************************************************************/
FILE *CreateMap(const char *Path, const char *Name, const char *Units,
                int NLayers, char (*LayerNames)[MAP_NAMESIZE],
                const char *Time, MAPSIZE *Map)
{
  FILE *OutFile;
  char *Body;
  char *p;
  char Line1[MAP_NAMESIZE];
  char Number[3][40];
  long long Content;
  long long HeaderBytes = MAP_HEADER_BLOCK;
  long long New;
  long long k;

  if (!HostIsLittleEndian())
    ReportErrorCode(113, "%s: maps are little-endian; this computer is not", Path);
  if (!(Body = (char *) malloc(1024 + (size_t) NLayers * (MAP_NAMESIZE + 24))))
    ReportError("CreateMap", 1);
  HeaderNumber(Number[0], (double) Map->DX);
  HeaderNumber(Number[1], Map->Yorig);
  HeaderNumber(Number[2], Map->Xorig);
  p = Body;
  p += sprintf(p, "Name = %s\nRows = %d\nColumns = %d\nLayers = %d\n", Name,
               Map->NY, Map->NX, NLayers);
  p += sprintf(p, "NumberType = FLOAT\nByteOrder = LITTLE\n");
  p += sprintf(p, "CellSize = %s\nNorth = %s\nWest = %s\n", Number[0], Number[1], Number[2]);
  p += sprintf(p, "NoData = -9999\nUnits = %s\n", Units);
  if (Time != NULL)
    p += sprintf(p, "Time = %s\n", Time);
  p += sprintf(p, "Writer = DHSVM %s\n", DHSVM_VERSION);
  for (k = 0; k < NLayers; k++)
    p += sprintf(p, "Layer %d = %s\n", (int) k + 1, LayerNames[k]);
  p += sprintf(p, "End\n");

  /* HeaderBytes = 4096 x ceiling((Content + 1) / 4096); line 1 holds it */
  for (;;) {
    sprintf(Line1, "DHSVM-MAP 1 %.0f\n", (double) HeaderBytes);
    Content = (long long) strlen(Line1) + (long long) (p - Body);
    New = MAP_HEADER_BLOCK * ((Content + MAP_HEADER_BLOCK) / MAP_HEADER_BLOCK);
    if (New == HeaderBytes)
      break;
    HeaderBytes = New;
  }

  if (!(OutFile = fopen(Path, "wb")))
    ReportErrorCode(113, "%s: cannot write", Path);
  if (fputs(Line1, OutFile) == EOF || fputs(Body, OutFile) == EOF)
    ReportErrorCode(113, "%s: cannot write", Path);
  for (k = Content; k < HeaderBytes - 1; k++)
    if (fputc(' ', OutFile) == EOF)
      ReportErrorCode(113, "%s: cannot write", Path);
  if (fputc('\n', OutFile) == EOF)
    ReportErrorCode(113, "%s: cannot write", Path);

  free(Body);
  return OutFile;
}

/*****************************************************************************
  WriteMapLayer(), CloseMapFile()
*****************************************************************************/
void WriteMapLayer(FILE *File, float *Values, MAPSIZE *Map, const char *Path)
{
  size_t NCells = (size_t) Map->NY * Map->NX;

  if (fwrite(Values, sizeof(float), NCells, File) != NCells)
    ReportErrorCode(113, "%s: cannot write", Path);
}

void CloseMapFile(FILE *File, const char *Path)
{
  if (fclose(File) != 0)
    ReportErrorCode(113, "%s: cannot write", Path);
}
