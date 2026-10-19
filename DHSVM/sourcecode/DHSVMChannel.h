
#ifndef _DHSVMChannel_h_
#define _DHSVMChannel_h_

#include "settings.h"		/* for data.h */
#include "data.h"
#include "getinit.h"
#include "channel.h"
#include "channel_grid.h"

/* -------------------------------------------------------------
   struct CHANNEL
   ------------------------------------------------------------- */
typedef struct {
  ChannelClass *stream_class;
  Channel *streams;
  ChannelMapPtr **stream_map;
} CHANNEL;

/* -------------------------------------------------------------
   struct CHANNELTABLES
   The channel class, segment and cell tables as read (rows in table
   order; cell rows and columns from 1)
   ------------------------------------------------------------- */
typedef struct {
  char ClassPath[BUFSIZE + 1];
  char SegmentPath[BUFSIZE + 1];
  char CellPath[BUFSIZE + 1];
  int NClasses;
  int *ClassID;
  float *Width;
  float *BankHeight;
  float *Friction;
  int NSegments;
  int *SegmentID;
  int *Order;
  float *Slope;
  float *Length;
  int *Class;
  int *Outlet;
  char **SaveName;
  int NCells;
  int *Row;
  int *Column;
  int *Segment;
  float *CellLength;
  float *CutHeight;
  float *CutWidth;
} CHANNELTABLES;

/* -------------------------------------------------------------
   available functions
   ------------------------------------------------------------- */
void ReadChannelTables(LISTPTR Input, CHANNELTABLES *Tables);
void ValidateChannels(OPTIONSTRUCT *Options, MAPSIZE *Map, TOPOPIX **TopoMap,
                      SOILPIX **SoilMap, CHANNELTABLES *Tables, LAKETABLE *LType);
void InitChannel(MAPSIZE *Map, int deltat, CHANNEL *channel, CHANNELTABLES *Tables,
	    SOILTABLE *SType, SOILPIX ** SoilMap,
	    LAKETABLE *LType, TOPOPIX **TopoMap,
	    int *MaxStreamID, OPTIONSTRUCT *Options);
double ChannelCulvertFlow(int y, int x, CHANNEL *ChannelData);
void RouteChannel(CHANNEL *ChannelData, TIMESTRUCT *Time, MAPSIZE *Map,
		  TOPOPIX **TopoMap, SOILPIX **SoilMap,
		  OPTIONSTRUCT *Options, NETSTRUCT **Network, SOILTABLE *SType,
		  VEGPIX **VegMap, EVAPPIX **Evap, LAKETABLE *LType);
void ChannelCut(int y, int x, CHANNEL *ChannelData, NETSTRUCT *Network);
void ChannelLimitVegFC(int y, int x, float DXDY, CHANNEL * ChannelData,
                       VEGPIX *LocalVeg);
uchar ChannelFraction(TOPOPIX *topo, ChannelMapRec *rds);

#endif
