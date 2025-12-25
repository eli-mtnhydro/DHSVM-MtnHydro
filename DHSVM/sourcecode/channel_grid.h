
#ifndef _channel_grid_h_
#define _channel_grid_h_

#include "channel.h"
#include "settings.h"
#include "data.h"

/* -------------------------------------------------------------
   struct ChannelMapRec
   This is used to locate the channel segment located within a gridcell
   ------------------------------------------------------------- */

struct _channel_map_rec_ {
  float length;			/* channel length within cell (m) */
  float cut_height;		/* channel cut depth (m) */
  float cut_width;		/* "effective" cut width (m) */
  float table_depth; /* local water table depth in portion of grid cell below channel */
  float infiltration_rate; /* infiltration rate out of the bottom of the channel (m/s) */
  float infiltration; /* amount of water (m^3) that could infiltrate if available */
  float evaporation; /* amount of water (m^3) that could evaporate if available */
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

				/* Input Functions */
				
void channel_grid_init(int cols, int rows);

ChannelMapPtr **channel_grid_read_map(Channel *net, const char *file,
				      SOILTABLE *SType, SOILPIX **SoilMap, VEGTABLE *VType, VEGPIX **VegMap);

void channel_combine_map_network(Channel * net, ChannelMapPtr ** map, MAPSIZE * Map);

				/* Query Functions */

int channel_grid_has_channel(ChannelMapPtr **map, int col, int row);
double channel_grid_cell_length(ChannelMapPtr **map, int col, int row);
double channel_grid_cell_width(ChannelMapPtr **map, int col, int row);
double channel_grid_cell_bankht(ChannelMapPtr **map, int col, int row);
float channel_grid_cell_maxbankht(ChannelMapPtr **map, int col, int row);
float channel_grid_cell_water_depth(ChannelMapPtr ** map, int col, int row);

float channel_grid_lateral_outflow(ChannelMapPtr ** map, int col, int row,
                                   float TableDepth,
                                   float Transmissivity,  float SoilDeficit,
                                   float DX, float DY, float Dt);

float channel_grid_calc_satflow(ChannelMapPtr ** map, int col, int row,
                                float TableDepth,
                                float Transmissivity, float AvailableWater,
                                float DX, float DY, float Dt);
void channel_grid_satflow(ChannelMapPtr ** map, int col, int row);

void channel_grid_inc_inflow(ChannelMapPtr **map, int col, int row, float mass);

void channel_grid_init_table(ChannelMapPtr ** map, int col, int row,
                             float GridTableDepth);
float channel_grid_table_depth(ChannelMapPtr ** map, int col, int row, int deltat,
                               float GridTableDepth, float Transmissivity,
                               float SoilDeficit, float DX);

void channel_grid_calc_infiltration(ChannelMapPtr ** map, int col, int row, int deltat,
                                    float TableDepth, float MaxInfiltrationCap, float DX);
float channel_grid_infiltration(ChannelMapPtr ** map, int col, int row);

void channel_grid_calc_evaporation(ChannelMapPtr ** map, int col, int row,
                                   float EPot, float MaxEvapCap);
float channel_grid_evaporation(ChannelMapPtr ** map, int col, int row);
float channel_grid_dry_evaporation(ChannelMapPtr ** map, int col, int row,
                                   float EPot, float MaxEvapCap, float DXDY,
                                   float Dt, float *Porosity, float *FCap, float *Ks,
                                   float *Press, float *m, float LayerThickness,
                                   float *MoistContent, float *Adjust, int CutBankZone);

void channel_grid_free_map(ChannelMapPtr **map);

#endif
