
#ifndef SETTINGS_H
#define SETTINGS_H

#ifndef _AIX			/* AIX 3.2.5 already defines this */
typedef unsigned char uchar;
#endif
typedef unsigned short unshort;
typedef unsigned int unint;

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define INBASIN(x) ((x) != OUTSIDEBASIN)
#ifndef ABSVAL
#define ABSVAL(x)  ( (x) < 0 ? -(x) : (x) )
#endif

#ifndef TRUE
#define TRUE           1
#endif
#ifndef FALSE
#define FALSE          0
#endif

#define DHSVM_HUGE     (1e20)

#define DHSVM_VERSION  "X.3.1"	/* the [RUN] Model Version this executable reads */

/* When compiling DHSVM using Microsoft C++ define MSC++ at compile time.
   Microsoft C++ treats all numeric constants as doubles, and issues a
   warning if this constant is assigned to a float without an explicit type
   cast.   This warning can becme somewhat annoying, and the following pragma
   directive disables the warning */
#ifdef MS_C_PLUSPLUS
#pragma warning(disable : 4244)
#endif

/* Default value for not applicable */
#define NOT_APPLICABLE -9999

/* Options for flow gradient calculation */
#define TOPOGRAPHY     1
#define WATERTABLE     2

/* Options for infiltration */
#define STATIC 1
#define DYNAMIC 2

#define TINY       1e-20
#define DEBUG      FALSE

#define HEADERLINES    5
#define BUFSIZE      255
#define MAXUCHAR     255	/* Maximum value of a 1-byte uchar */
#define MAXSTRING    255
#define NAMESIZE     127

#define MAXDIRS        8	/* Maximum number of directions in which water can flow, must equal 8 */
#define NNEIGHBORS     8  /* Number of directions in which water can flow based on fine grid, must equal 8 */
#define MAXVEGLAYERS   2  /* Vegetation layers a cell can hold: overstory and understory */


#define NA          -9999	/* Not applicable */

#define MIN_SWE 0.005

#endif
