
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "fileio.h"
#include "functions.h"
#include "constants.h"
#include "getinit.h"
#include "sizeofnt.h"
#include "slopeaspect.h"
#include "varid.h"

/*****************************************************************************
InitParameterMaps()
*****************************************************************************/
void InitParameterMaps(OPTIONSTRUCT *Options, MAPSIZE *Map, int Id, 
  char *FileName, SNOWPIX ***SnowMap, int ParamType, float temp) {


  const char *Routine = "InitParameterMaps";
  char VarName[BUFSIZE + 1];	/* Variable name */
  int i;			/* Counter */
  int x;			/* Counter */
  int y;			/* Counter */
  int NumberType;		/* Number type of data set */
  float *Array;

  /* Read the map */

  if (ParamType == MAP) {
	GetVarName(Id, 0, VarName);
	GetVarNumberType(Id, &NumberType);
	if (!(Array = (float *)calloc(Map->NX * Map->NY, SizeOfNumberType(NumberType))))
	  ReportError((char *)Routine, 1);

	Read2DMatrix(FileName, Array, NumberType, Map, 0, VarName, 0);
	
	for (y = 0, i = 0; y < Map->NY; y++) {
	  for (x = 0; x < Map->NX; x++, i++) {
	    switch (Id) {
	    case 800:
	      (*SnowMap)[y][x].Ts = Array[i];
	      break;
	    case 801:
	      (*SnowMap)[y][x].Tr = Array[i];
	      break;
	    case 802:
	      (*SnowMap)[y][x].amax = Array[i];
	      break;
	    case 803:
	      (*SnowMap)[y][x].LamdaAcc = Array[i];
	      break;
	    case 804:
	      (*SnowMap)[y][x].LamdaMelt = Array[i];
	      break;
	    case 805:
	      (*SnowMap)[y][x].AccMin = Array[i];
	      break;
	    case 806:
	      (*SnowMap)[y][x].MeltMin = Array[i];
	      break;
	    default:
	      printf("%s: Map ID %d not found", Routine, Id);
	    exit(74);
	    } /* end switch (n) */
	  }
	}

	free(Array);
  }
  /* assign a constant parameter to all model grids*/
  else if (ParamType == CONSTANT) {
	for (y = 0; y < Map->NY; y++) {
	  for (x = 0; x < Map->NX; x++) {
		switch (Id) {
		case 800:
		  (*SnowMap)[y][x].Ts = temp;
		  break;
		case 801:
		  (*SnowMap)[y][x].Tr = temp;
		  break;
		case 802:
		  (*SnowMap)[y][x].amax = temp;
		  break;
		case 803:
		  (*SnowMap)[y][x].LamdaAcc = temp;
		  break;
		case 804:
		  (*SnowMap)[y][x].LamdaMelt = temp;
		  break;
		case 805:
		  (*SnowMap)[y][x].AccMin = temp;
		  break;
		case 806:
		  (*SnowMap)[y][x].MeltMin = temp;
		  break;
		default:
		  printf("%s: Map ID %d not found", Routine, Id);
		  exit(0);
		}
	  }
	}
  }
  else {
	printf("%s: Parameter type %d not found", Routine, ParamType);
	exit(75);
  }
}