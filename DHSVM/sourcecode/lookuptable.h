
#ifndef LOOKUP_TABLE_H
#define LOOKUP_TABLE_H

typedef struct {
  unsigned long Size;		/* Number of elements in lookup table */
  float Offset;			/* Value of key of first entry in the table */
  float Delta;			/* Interval between keys */
  float *Data;			/* Pointer to array with entries */
} FLOATTABLE;

float FloatLookup(float x, FLOATTABLE * Table);
void InitFloatTable(unsigned long Size, float Offset, float Delta,
		    float (*Function) (float), FLOATTABLE * Table);

#endif
