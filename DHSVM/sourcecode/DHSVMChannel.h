
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
  FILE *streamout;
  FILE *streamflowout;
} CHANNEL;

/* -------------------------------------------------------------
   available functions
   ------------------------------------------------------------- */
void InitChannel(LISTPTR Input, MAPSIZE *Map, int deltat, CHANNEL *channel,
	    SOILTABLE *SType, SOILPIX ** SoilMap, VEGTABLE *VType, VEGPIX **VegMap,
	    LAKETABLE *LType, TOPOPIX **TopoMap,
	    int *MaxStreamID, OPTIONSTRUCT *Options);
void InitChannelDump(OPTIONSTRUCT *Options, CHANNEL *channel, char *DumpPath);
double ChannelCulvertFlow(int y, int x, CHANNEL *ChannelData);
void RouteChannel(CHANNEL *ChannelData, TIMESTRUCT *Time, MAPSIZE *Map,
		  TOPOPIX **TopoMap, SOILPIX **SoilMap, AGGREGATED *Total, 
		  OPTIONSTRUCT *Options, NETSTRUCT **Network, SOILTABLE *SType,
		  VEGTABLE *VType, VEGPIX **VegMap, EVAPPIX **Evap, LAKETABLE *LType);
void ChannelCut(int y, int x, CHANNEL *ChannelData, NETSTRUCT *Network);
void ChannelLimitVegFC(int y, int x, float DXDY, CHANNEL * ChannelData,
                       VEGTABLE *VType, VEGPIX *LocalVeg);
uchar ChannelFraction(TOPOPIX *topo, ChannelMapRec *rds);

#endif
