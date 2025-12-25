
#ifndef VARID_H
#define VARID_H

#ifndef ENDOFLIST
#define ENDOFLIST -1
#endif

void GetVarAttr(MAPDUMP * DMap);
void GetVarFileLabel(int ID, char *FileLabel);
void GetVarFileName(int ID, int Layer, unsigned char Resolution,
		    char *FileName);
void GetVarLongName(int ID, int Layer, char *LongName);
int GetVarNLayers(int ID, int MaxSoilLayers, int MaxVegLayers);
void GetVarName(int ID, int Layer, char *Name);
void GetVarNumberType(int ID, int *NumberType);
void GetVarUnits(int ID, char *Units);
unsigned char IsMultiLayer(int ID);
unsigned char IsValidID(int ID);

#endif
