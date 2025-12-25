
#include <stdlib.h>
#include <stdio.h>
#include "functions.h"

void InitCharArray(char *Array, int Size)
{
  int i;			/* counter */

  for (i = 0; i < Size; i++)
    Array[i] = '\0';
}
