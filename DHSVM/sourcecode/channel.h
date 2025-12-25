
#ifndef _channel_h_
#define _channel_h_

#include "settings.h"

typedef struct LAKETABLE LAKETABLE;

typedef unsigned short int SegmentID, ClassID;

/* -------------------------------------------------------------
   struct ChannelClass
   ------------------------------------------------------------- */

typedef struct _channel_class_rec_ {
  ClassID id;			/* unique identifier */
  float width;			/* channel width */
  float bank_height;		/* bank height for streams */
  float friction;		        /* Manning's n for the channel*/
  struct _channel_class_rec_ *next;

} ChannelClass;

/* -------------------------------------------------------------
   struct Channel
   This is the basic unit of channel information.
   ------------------------------------------------------------- */
struct _channel_rec_ {
  SegmentID id;
  SegmentID outid;
  unsigned order;		/* determines computation order */
  char *record_name;	/* The name this segment is to have in the output, if output is recorded */
  char record;			/* TRUE if outflow values are to be saved by channel_save_outflow */
  float length;			/* Parameters */
  float slope;              /* Effective water surface slope, m/m (recomputed each timestep) */
  float ground_slope;       /* Ground slope, m/m (equal to slope if water depth is uniform) */
  float top_water_depth;    /* Depth of water storage (m) at top of channel segment */
  float bottom_water_depth; /* Depth of water storage (m) at bottom of channel segment */
  float K;              /* Travel time constant, a function of slope */
  float X;              /* Weighting factor (0~1), exponential function of K */
  ChannelClass *class2;	/* ChannelClass identifier */
  uchar IntersectsLake; /* Logical flag for whether the channel enters a lake */
  LAKETABLE *lake;        /* Pointer to lake that is intersected by this channel */

  /* necessary routing terms */
  float lateral_inflow;	/* cubic meters */
  float last_inflow;	/* cubic meters */
  float last_outflow;	/* cubic meters */
  float last_storage;	/* cubic meters */
  float inflow;			/* cubic meters */
  float lake_inflow; /* cubic meters */
  float outflow;		/* cubic meters */
  float storage;		/* cubic meters */
  float infiltration;		/* cubic meters */
  float remaining_infil; /* cubic meters */
  float evaporation; /* cubic meters */
  float remaining_evap; /* cubic meters */
  float last_lateral_inflow;

  struct _channel_rec_ *outlet;	/* NULL if does not drain to another segment */
  struct _channel_rec_ *next;
  struct _channel_map_rec_ *grid; /* Pointer to first map cell in current segment */
};
typedef struct _channel_rec_ Channel, *ChannelPtr;

/* -------------------------------------------------------------
   externally available routines
   ------------------------------------------------------------- */

/* ChannelClass */
ChannelClass *channel_read_classes(const char *file, int ChanType);
void channel_free_classes(ChannelClass *head);

/* Channel */
Channel *channel_read_network(const char *file, ChannelClass * class_list, int *MaxID);
void channel_routing_parameters(Channel *net, int deltat);
void channel_update_routing_parameters(Channel *network, int deltat, int max_order);
Channel *channel_find_segment(Channel *net, SegmentID id);
int channel_step_initialize_network(Channel *net);
int channel_incr_lat_inflow(Channel *segment, float linflow);
void channel_segment_infil_evap(Channel * segment);
int channel_route_network(Channel *net, int deltat);
int channel_save_outflow(double time, Channel * net, FILE *file, FILE *file2);
int channel_save_outflow_text(char *tstring, Channel *net, FILE *out,
			      FILE *out2, int flag);
void channel_free_network(Channel *net);

#endif
