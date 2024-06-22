/* -------------------------------------------------------------
   file: channel_grid.h

   This module provides the necessary interface between the watershed
   model and the channel routing module.
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Battelle Memorial Institute
   Pacific Northwest Laboratory
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Created January  4, 1996 by  William A Perkins
   $Id: channel_grid.h,v 1.7 2004/05/03 03:28:48 colleen Exp $
   ------------------------------------------------------------- */

#ifndef _channel_grid_h_
#define _channel_grid_h_

#include "channel.h"
#include "settings.h"
#include "data.h"

/* -------------------------------------------------------------
   struct ChannelMapRec
   This is used to locate the channel segment located within a grid
   cell.  And to determine if the channel network has a sink in any of
   all of the segments which pass thru the cell
   ------------------------------------------------------------- */

struct _channel_map_rec_ {
  float length;			/* channel length within cell (m) */
  float aspect;			/* channel aspect within cell (radians) */
  float cut_height;		/* channel cut depth (m) */
  float cut_width;		/* "effective" cut width (m) */
  char sink;			/* is this cell a channel sink? */
  float azimuth;        /* channel azimuth */
  float infiltration_rate; /* infiltration rate out of the bottom of the channel (m/s) */
  float infiltration; /* amount of water (m^3) that could infiltrate if available */
  float avail_water; /* amount of water (m^3) in segment that is derived from uphill */
  float satflow; /* amount of water (m^3) flowing laterally into channel from soil */
  Channel *channel;		/* pointer to segment record */

  struct _channel_map_rec_ *next;
  struct _channel_map_rec_ *next_seg; /* Subset of cells in a particular segment only */
};
typedef struct _channel_map_rec_ ChannelMapRec;
typedef struct _channel_map_rec_ *ChannelMapPtr;

/* -------------------------------------------------------------
   externally available routines
   ------------------------------------------------------------- */

				/* Module Functions */

void channel_grid_init(int cols, int rows);
void channel_grid_done(void);

				/* Input Functions */

ChannelMapPtr **channel_grid_read_map(Channel *net, const char *file,
				      SOILTABLE *SType, SOILPIX **SoilMap, VEGTABLE *VType, VEGPIX **VegMap);

				/* Query Functions */

int channel_grid_has_channel(ChannelMapPtr **map, int col, int row);
int channel_grid_has_sink(ChannelMapPtr **map, int col, int row);
double channel_grid_cell_length(ChannelMapPtr **map, int col, int row);
double channel_grid_cell_width(ChannelMapPtr **map, int col, int row);
double channel_grid_cell_bankht(ChannelMapPtr **map, int col, int row);
float channel_grid_cell_maxbankht(ChannelMapPtr **map, int col, int row);

float channel_grid_calc_satflow(ChannelMapPtr ** map, int col, int row,
                                float TableDepth,
                                float Transmissivity, float AvailableWater,
                                float DX, float DY, float Dt);
void channel_grid_satflow(ChannelMapPtr ** map, int col, int row);

void channel_grid_inc_inflow(ChannelMapPtr **map, int col, int row, float mass);
void channel_grid_inc_melt(ChannelMapPtr **map, int col, int row, float mass);                                                                              
double channel_grid_outflow(ChannelMapPtr **map, int col, int row);
double channel_grid_flowlength(ChannelMapPtr **map, int col, int row, 
			       float floslope);
double channel_grid_flowslope(ChannelMapPtr **map, int col, int row);
ChannelClass* channel_grid_class(ChannelMapPtr **map, int col, int row);

void channel_grid_calc_infiltration(ChannelMapPtr ** map, int col, int row, int deltat,
                                    float TableDepth, float MaxInfiltrationCap);
float channel_grid_infiltration(ChannelMapPtr ** map, int col, int row);

float channel_grid_evaporation(ChannelMapPtr ** map, int col, int row,
                               float EPot, float MaxEvapCap);

void channel_grid_free_map(ChannelMapPtr **map);

/* new functions for RBM model */
void channel_grid_inc_other(ChannelMapPtr **map, int col, int row, PIXRAD *LocalRad , 
							PIXMET *LocalMet, float skyview);
void Init_segment_ncell(TOPOPIX **TopoMap, ChannelMapPtr **map, int NY, int NX, Channel *net);
void channel_grid_avg (Channel *Channel);
#endif
