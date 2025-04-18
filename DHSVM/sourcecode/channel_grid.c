/* -------------------------------------------------------------
   file: channel_grid.c
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Battelle Memorial Institute
   Pacific Northwest Laboratory
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Created January  5, 1996 by  William A Perkins
   $Id: channel_grid.c,v 3.1.2 2014/1/2 Ning Exp $
   ------------------------------------------------------------- */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifndef PI
#define PI 3.14159265358979323846
#endif
#include <errno.h>
#include <string.h>

#include "channel_grid.h"
#include "tableio.h"
#include "errorhandler.h"
#include "settings.h"
#include "data.h"
#include "DHSVMChannel.h"
#include "constants.h"
#include "functions.h"
#include "massenergy.h"

/* -------------------------------------------------------------
   local function prototype
   ------------------------------------------------------------- */
static ChannelMapRec *alloc_channel_map_record(void);
static ChannelMapPtr **channel_grid_create_map(int cols, int rows);
Channel *Find_First_Segment(ChannelMapPtr **map, int col, int row, float SlopeAspect, 
			    char *Continue);
char channel_grid_has_intersection(ChannelMapPtr **map, int Currid, int Nextid, int row, 
				   int col, int Flag);

/* -------------------------------------------------------------
   local module variables
   ------------------------------------------------------------- */
static int channel_grid_cols = 0;
static int channel_grid_rows = 0;
static char channel_grid_initialized = FALSE;


/* -------------------------------------------------------------
   Find_First_Segment
   ------------------------------------------------------------- */

Channel *Find_First_Segment(ChannelMapPtr ** map, int col, int row, float SlopeAspect, char *Continue)
{
  ChannelMapPtr cell = map[col][row];
  float test;
  float DeltaAspect;
  Channel *Ptr = NULL;
  
  DeltaAspect = 2.*PI;

  while (cell != NULL) {
    test = fabs(SlopeAspect - cell->aspect);

    if(test > PI) {
      if(SlopeAspect < cell->aspect)
	test = fabs(SlopeAspect - (cell->aspect - 2.*PI));
      else
	test = fabs(cell->aspect - (SlopeAspect - 2.*PI));
    }
    if(test < 0. || test > PI) {
      printf("Problem in Find_First_Segment\n");
      exit(0);
    }
    
    if( test < DeltaAspect) {
      Ptr = cell->channel;
      DeltaAspect = test;
    }
    cell = cell->next;
  }
  if(DeltaAspect <= 70.*PI/180.)
    *Continue = TRUE;
  else
    *Continue = FALSE;

  return (Ptr);
}

/* -------------------------------------------------------------
   alloc_channel_map_record
   ------------------------------------------------------------- */
static ChannelMapRec *alloc_channel_map_record(void)
{
  ChannelMapRec *p;
  if ((p = (ChannelMapRec *) malloc(sizeof(ChannelMapRec))) == NULL) {
    error_handler(ERRHDL_FATAL,
		  "alloc_channel_map_record: %s", strerror(errno));
  }
  p->length = 0.0;
  p->aspect = 0.0;
  p->infiltration_rate = 0.0;
  p->sink = FALSE;
  p->channel = NULL;
  p->next = NULL;
  p->next_seg = NULL;
  return (p);
}

/* -------------------------------------------------------------
   channel_grid_create_map
   ------------------------------------------------------------- */
static ChannelMapPtr **channel_grid_create_map(int cols, int rows)
{
  ChannelMapPtr **map;
  int row, col;
  ChannelMapPtr *junk;

  if ((map = (ChannelMapPtr **) malloc(cols * sizeof(ChannelMapPtr *))) == NULL) {
    error_handler(ERRHDL_FATAL, "channel_grid_create_map: malloc failed: %s",
		  strerror(errno));
  }
  if ((junk =
       (ChannelMapPtr *) malloc(rows * cols * sizeof(ChannelMapPtr))) == NULL) {
    free(map);
    error_handler(ERRHDL_FATAL,
		  "channel_grid_create_map: malloc failed: %s",
		  strerror(errno));
  }

  for (col = 0; col < cols; col++) {
    if (col == 0) {
      map[col] = junk;
    }
    else {
      map[col] = &(map[0][col * rows]);
    }
    for (row = 0; row < rows; row++) {
      map[col][row] = NULL;
    }
  }
  return (map);
}

/* -------------------------------------------------------------
   free_channel_map_record
   ------------------------------------------------------------- */
static void free_channel_map_record(ChannelMapRec * cell)
{
  if (cell->next != NULL) {
    free_channel_map_record(cell->next);
  }
  free(cell);
}

/* -------------------------------------------------------------
   channel_grid_free_map
   ------------------------------------------------------------- */
void channel_grid_free_map(ChannelMapPtr ** map)
{
  int c, r;
  for (c = 0; c < channel_grid_cols; c++) {
    for (r = 0; r < channel_grid_rows; r++) {
      if (map[c][r] != NULL) {
	free_channel_map_record(map[c][r]);
      }
    }
  }
  free(map[0]);
  free(map);
}

/* -------------------------------------------------------------
   ------------------- Input Functions -------------------------
   ------------------------------------------------------------- */

/* -------------------------------------------------------------
   channel_grid_read_map
   ------------------------------------------------------------- */
ChannelMapPtr **channel_grid_read_map(Channel *net, const char *file,
				      SOILTABLE *SType, SOILPIX ** SoilMap, VEGTABLE *VType, VEGPIX **VegMap)
{
  ChannelMapPtr **map;
  static const int fields = 8;
  static char *sink_words[2] = {
    "SINK", "\n"
  };
  static TableField map_fields[8] = {
    {"Column", TABLE_INTEGER, TRUE, FALSE, {0.0}, "", NULL},
    {"Row", TABLE_INTEGER, TRUE, FALSE, {0.0}, "", NULL},
    {"Segment ID", TABLE_INTEGER, TRUE, FALSE, {0.0}, "", NULL},
    {"Segment Length", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
    {"Cut Height", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
    {"Cut Width", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
    {"Segment Azimuth", TABLE_REAL, FALSE, FALSE, {0.0}, "", NULL},
    {"Sink?", TABLE_WORD, FALSE, FALSE, {0.0}, "", sink_words}
  };
  int done, err = 0;
  int iSoil;
  float Depth;

  if (!channel_grid_initialized) {
    error_handler(ERRHDL_ERROR,
		  "channel_grid_read_map: channel_grid module not initialized");
    return NULL;
  }

  error_handler(ERRHDL_STATUS,
		"channel_grid_read_map: reading file \"%s\"", file);

  if (table_open(file) != 0) {
    error_handler(ERRHDL_ERROR,
		  "channel.grid_read_map: unable to read file \"%s\"", file);
    return NULL;
  }

  map = channel_grid_create_map(channel_grid_cols, channel_grid_rows);

  done = FALSE;
  while (!done) {
    int i;
    int row = 0, col = 0;
    int rec_err = 0;
    ChannelMapPtr cell;

    done = (table_get_fields(fields, map_fields) < 0);
    if (done) {
      for (i = 0; i < fields; i++) {
	if (map_fields[i].read)
	  break;
      }
      if (i >= fields)
	continue;
    }

    if (map_fields[0].read) {
      if (map_fields[0].value.integer < 0 ||
	  map_fields[0].value.integer >= channel_grid_cols) {
	rec_err++;
      }
      else {
	col = map_fields[0].value.integer;
      }
    }
    else {
      rec_err++;
    }
    if (map_fields[1].read) {
      if (map_fields[1].value.integer < 0 ||
	  map_fields[1].value.integer >= channel_grid_rows) {
	rec_err++;
      }
      else {
	row = map_fields[1].value.integer;
      }
    }
    else {
      rec_err++;
    }

    if (rec_err) {
      error_handler(ERRHDL_ERROR,
		    "%s: line %d: bad coordinates", file, table_lineno());
      err++;
      continue;
    }

    if (map[col][row] != NULL) {
      cell = map[col][row];
      while (cell->next != NULL)
        cell = cell->next;
      cell->next = alloc_channel_map_record();
      cell = cell->next;
    }
    else {
      map[col][row] = alloc_channel_map_record();
      cell = map[col][row];
    }
    
    for (i = 2; i < fields; i++) {
      if (map_fields[i].read) {
	switch (i) {
	case 2:
	  if ((cell->channel =
	       channel_find_segment(net,
				    map_fields[i].value.integer)) == NULL) {
	    error_handler(ERRHDL_ERROR,
			  "%s, line %d: unable to locate segment %d", file,
			  table_lineno(), map_fields[i].value.integer);
	    err++;
	  }
	  break;
	case 3:
	  cell->length = map_fields[i].value.real;
	  if (cell->length < 0.0) {
	    error_handler(ERRHDL_ERROR,
			  "%s, line %d: bad length", file, table_lineno());
	    err++;
	  }
	  break;
	case 4:
	  cell->cut_height = map_fields[i].value.real;
	  if (cell->cut_height > SoilMap[row][col].Depth) {
	    /*printf("warning overriding cut depths with 0.95 soil depth \n");*/
	    cell->cut_height = SoilMap[row][col].Depth*0.95;
	  }
	  if (cell->cut_height < 0.0
	      || cell->cut_height > SoilMap[row][col].Depth) {
	    error_handler(ERRHDL_ERROR, "%s, line %d: bad cut_depth", file,
			  table_lineno());
	    err++;
	  }
	  break;
	case 5:
	  cell->cut_width = map_fields[i].value.real;
	  if (cell->cut_width < 0.0) {
	    error_handler(ERRHDL_ERROR,
			  "%s, line %d: bad cut_width", file, table_lineno());
	    err++;
	  }
	  break;
	case 6:
	  /* road aspect is read in degrees and
	     stored in radians */
	  cell->azimuth = (float)(map_fields[i].value.real);
	  cell->aspect = map_fields[i].value.real * PI / 180.0;
	  break;
	case 7:
	  cell->sink = TRUE;
	  break;
	default:
	  error_handler(ERRHDL_FATAL,
			"channel_grid_read_map: this should not happen");
	  break;
	  }
    }
    }
    
    /* Set infiltration rate to vertical hydraulic conductivity
       of soil layer containing bottom of channel cut */
    Depth = 0.0;
    for (iSoil = 0; iSoil < SType[SoilMap[row][col].Soil - 1].NLayers && Depth < cell->cut_height; iSoil++) {
      if (VType[VegMap[row][col].Veg - 1].RootDepth[iSoil] < (SoilMap[row][col].Depth - Depth))
        Depth += VType[VegMap[row][col].Veg - 1].RootDepth[iSoil];
      else
        Depth = SoilMap[row][col].Depth;
    }
    if (Depth > cell->cut_height)
      cell->infiltration_rate = SoilMap[row][col].KsVert[iSoil - 1];
    else /* Channel bottom is below root zone */
      cell->infiltration_rate = SoilMap[row][col].KsVert[iSoil];
    
    /* Initialize available storage for re-infiltration */
    cell->avail_water = 0.0;
  }

  table_errors += err;
  error_handler(ERRHDL_STATUS,
		"channel_grid_read_map: %s: %d errors, %d warnings",
		file, table_errors, table_warnings);

  table_close();

  error_handler(ERRHDL_STATUS,
		"channel_grid_read_map: done reading file \"%s\"", file);

  if (table_errors) {
    error_handler(ERRHDL_ERROR,
		  "channel_grid_read_map: %s: too many errors", file);
    channel_grid_free_map(map);
    map = NULL;
  }

  return (map);
}

/* -------------------------------------------------------------
 channel_combine_map_network
 ------------------------------------------------------------- */
void channel_combine_map_network(Channel * net, ChannelMapPtr ** map, MAPSIZE * Map)
{
  Channel *segment;
  ChannelMapPtr cell, cell2;
  int j, k, col, row;
  
  for (segment = net; segment != NULL; segment = segment->next) {
    /* Check all grid cells to see if they contain the current segment */
    for (k = (Map->NumCells - 1); k > -1;  k--) {
      row = Map->OrderedCells[k].y;
      col = Map->OrderedCells[k].x;
      cell = map[col][row];
      while (cell != NULL) {
        if (cell->channel->id == segment->id){
          
          /* Entry point from network to map */
          if (segment->grid == NULL)
            segment->grid = cell;
          
          /* Find next downstream cell which also matches channel id */
          for (j = (k - 1); (cell->next_seg == NULL && j > -1);  j--) {
            row = Map->OrderedCells[j].y;
            col = Map->OrderedCells[j].x;
            cell2 = map[col][row];
            while (cell2 != NULL && cell->next_seg == NULL) {
              if (cell2->channel->id == segment->id)
                cell->next_seg = cell2;
              cell2 = cell2->next;
            }
          } /* End of finding downstream cell */
        }
        cell = cell->next;
      }
    }
  }
}

/* -------------------------------------------------------------
   ---------------------- Query Functions ---------------------
   ------------------------------------------------------------- */

/* -------------------------------------------------------------
   channel_grid_has_channel
   ------------------------------------------------------------- */
int channel_grid_has_channel(ChannelMapPtr ** map, int col, int row)
{
  if (map != NULL)
    return (map[col][row] != NULL);
  else
    return FALSE;
}

/* -------------------------------------------------------------
   channel_grid_has_sink
   ------------------------------------------------------------- */
int channel_grid_has_sink(ChannelMapPtr ** map, int col, int row)
{
  ChannelMapPtr cell = map[col][row];
  char test = FALSE;

  while (cell != NULL) {
    test = (test || cell->sink);
    cell = cell->next;
  }
  return (test);
}

/* -------------------------------------------------------------
   channel_grid_has_intersection
   ------------------------------------------------------------- */
char channel_grid_has_intersection(ChannelMapPtr ** map, int Currid, int Nextid,  int row, int col, int Flag)
{
  ChannelMapPtr cell = map[col][row];
  char Next = FALSE;
  char Current = FALSE;
  char Intersection = FALSE;
  
  if( Flag < 2) { /* search for intersecting cells given Currid and Nextid
		      occur in the cell of intersection */
    while (cell != NULL) {
      if(cell->channel->id == Currid)
	Current = TRUE;
      if(cell->channel->id == Nextid)
	Next = TRUE;
      cell = cell->next;
    }
    
    if(Current && Next)
      Intersection = TRUE;
  }

  else if(Flag == 2) { /* if Flag == 2 , only search for the nearest cell that has id of next
	    segment */
    while (cell != NULL) {
      if(cell->channel->id == Nextid)
	Next = TRUE;
      cell = cell->next;
    }
    
    if(Next)
      Intersection = TRUE;
  }
  else { /* if Flag == 3 , only search for the nearest cell that has id of current
	    segment */
    while (cell != NULL) {
      if(cell->channel->id == Currid)
	Current = TRUE;
      cell = cell->next;
    }
    
    if(Current)
      Intersection = TRUE;
  }

  return (Intersection);
}

/* -------------------------------------------------------------
   channel_grid_cell_length
   returns the total length of channel(s) in the cell.
   ------------------------------------------------------------- */
double channel_grid_cell_length(ChannelMapPtr ** map, int col, int row)
{
  ChannelMapPtr cell = map[col][row];
  double len = 0.0;

  while (cell != NULL) {
    len += cell->length;
    cell = cell->next;
  }
  return len;
}

/* -------------------------------------------------------------
   channel_grid_cell_width
   returns a length-weighted average of the channel widths in the cell
   ------------------------------------------------------------- */
double channel_grid_cell_width(ChannelMapPtr ** map, int col, int row)
{
  ChannelMapPtr cell = map[col][row];
  double len = channel_grid_cell_length(map, col, row);
  double width = 0.0;

  if (len > 0.0) {
    while (cell != NULL) {
      width += cell->cut_width * cell->length;
      cell = cell->next;
    }
    width /= len;
  }

  return width;
}

/* -------------------------------------------------------------
   channel_grid_cell_bankheight
   ------------------------------------------------------------- */
double channel_grid_cell_bankht(ChannelMapPtr ** map, int col, int row)
{
  ChannelMapPtr cell = map[col][row];
  double len = channel_grid_cell_length(map, col, row);
  double height = 0.0;

  if (len > 0.0) {
    while (cell != NULL) {
      height += cell->cut_height * cell->length;
      cell = cell->next;
    }
    height /= len;
  }
  return (height);
}

/* -------------------------------------------------------------
 channel_grid_cell_maxbankheight
 ------------------------------------------------------------- */
float channel_grid_cell_maxbankht(ChannelMapPtr ** map, int col, int row)
{
  ChannelMapPtr cell = map[col][row];
  float height = 0.0;
  
  while (cell != NULL) {
    if (cell->cut_height > height)
      height = cell->cut_height;
    cell = cell->next;
  }
  return (height);
}

/* -------------------------------------------------------------
   channel_grid_cell_water_depth
   ------------------------------------------------------------- */
float channel_grid_cell_water_depth(ChannelMapPtr ** map, int col, int row)
{
  ChannelMapPtr cell = map[col][row];
  float len = channel_grid_cell_length(map, col, row);
  float water_depth;
  float water_depth_avg = 0.0;

  if (len > 0.0) {
    while (cell != NULL) {
      water_depth = ((cell->channel->storage + cell->channel->last_storage) / 2.0) /
                    (cell->channel->class2->width * cell->channel->length);
      water_depth_avg += water_depth * cell->length;
      cell = cell->next;
    }
    water_depth_avg /= len;
  }
  return (water_depth_avg);
}

/* -------------------------------------------------------------
   channel_grid_lateral_outflow
   New method for calculating channel outflow to subsurface
   taking into account the different cut depths so that only
   channels with cut depth < water table lose water
   ------------------------------------------------------------- */
float channel_grid_lateral_outflow(ChannelMapPtr ** map, int col, int row,
                                   float TableDepth,
                                   float Transmissivity,  float SoilDeficit,
                                   float DX, float DY, float Dt)
{
  ChannelMapPtr cell = map[col][row];
  float outflow; /* From one channel segment */
  float cell_outflow = 0.0; /* From all channel segments in cell */
  float max_outflow;
  float eff_dist;
  float drop;
  float grad;
  float water_depth;
  
  while (cell != NULL) {
    
    water_depth = ((cell->channel->storage + cell->channel->last_storage) / 2.0) /
                  (cell->channel->class2->width * cell->channel->length);
    
    if ((cell->cut_height - water_depth) < TableDepth) {
      
      /* Compute gradient from halfway between edge of cell and edge of channel */
      /* But enforce upper limit of 1 m for steepest gradient calculation */
      /* Assumes DX == DY */
      eff_dist = (DX - cell->cut_width) / 4;
      if (eff_dist < 1.0)
        eff_dist = 1.0;
      
      /* Water can flow out from both sides of each channel segment */
      if (water_depth > cell->cut_height)
        water_depth = cell->cut_height;
      drop = (TableDepth - (cell->cut_height - water_depth)) / eff_dist;
      grad = drop * (cell->length * 2);
      
      outflow = Transmissivity * grad * Dt;
      
      /* Limit outflow to amount that would equilibrate channel water level and TableDepth */
      /* Let TableDepthNew = cut_height - WaterDepthNew */
      /* TableDepthNew = TableDepth - outflow / (SoilDeficit * DX * DY) */
      /* WaterDepthNew = water_depth - outflow / (channel width * channel length) */
      if (SoilDeficit < 0.00001)
        max_outflow = 0.0;
      else
        max_outflow = (TableDepth - (cell->cut_height - water_depth)) /
          (1.0 / (SoilDeficit * DX * DY) + 1.0 / (cell->channel->class2->width * cell->channel->length));
      if (max_outflow > cell->channel->storage)
        max_outflow = cell->channel->storage;
      if (outflow > max_outflow)
        outflow = max_outflow;
      if (outflow < 0.0)
        outflow = 0.0;
      
      cell->satflow -= outflow;
      cell_outflow += outflow;
      TableDepth -= outflow / (SoilDeficit * DX * DY);
    }
    cell = cell->next;
  }
  return (cell_outflow);
}
/* -------------------------------------------------------------
   channel_grid_calc_satflow
   New method for calculating channel inflow from subsurface
   taking into account the different cut depths so that only
   channels with cut depth > water table get water
   ------------------------------------------------------------- */
float channel_grid_calc_satflow(ChannelMapPtr ** map, int col, int row,
                                float TableDepth,
                                float Transmissivity, float AvailableWater,
                                float DX, float DY, float Dt)
{
  ChannelMapPtr cell = map[col][row];
  float inflow; /* Into one channel segment */
  float cell_inflow = 0.0; /* Into all channel segments in cell */
  float max_inflow;
  float eff_dist;
  float drop;
  float grad;
  float water_depth;
  
  max_inflow = AvailableWater * DX * DY;
  
  while (cell != NULL) {
    
    water_depth = ((cell->channel->storage + cell->channel->last_storage) / 2.0) /
                  (cell->channel->class2->width * cell->channel->length);
    
    if ((cell->cut_height - water_depth) > TableDepth) {
      
      /* Compute gradient from halfway between edge of cell and edge of channel */
      /* But enforce upper limit of 1 m for steepest gradient calculation */
      /* Assumes DX == DY */
      eff_dist = (DX - cell->cut_width) / 4;
      if (eff_dist < 1.0)
        eff_dist = 1.0;
      
      /* Water can flow in from both sides of each channel segment */
      drop = (cell->cut_height - TableDepth - water_depth) / eff_dist;
      grad = drop * (cell->length * 2);
      
      inflow = Transmissivity * grad * Dt;
      
      if (inflow > max_inflow)
        inflow = max_inflow;
      if (inflow < 0.0)
        inflow = 0.0;
      
      cell->satflow += inflow;
      cell_inflow += inflow;
      max_inflow -= inflow;
    }
    cell = cell->next;
  }
  return (cell_inflow);
}

/* -------------------------------------------------------------
   channel_grid_satflow
   Transfer satflow to lateral inflow
   ------------------------------------------------------------- */
void channel_grid_satflow(ChannelMapPtr ** map, int col, int row)
{
  ChannelMapPtr cell = map[col][row];
  
  while (cell != NULL) {
    
    cell->channel->lateral_inflow += cell->satflow;
    cell->satflow = 0.0;
    
    /* Updates the amount of water available for infiltration;
       only water originating from inflow above the current cell
       is eligible for infiltration at the current location. */
    if (cell->avail_water < 0.0)
      cell->avail_water = 0.0;
    cell->avail_water += cell->channel->lateral_inflow;
    
    cell = cell->next;
  }
}

/* -------------------------------------------------------------
   channel_grid_inc_inflow
   Given a flow (or actually mass), this function increases the inflow
   any channel(s) in the cell in proportion to their length within the
   cell.
   ------------------------------------------------------------- */
void channel_grid_inc_inflow(ChannelMapPtr ** map, int col, int row, float mass)
{
  ChannelMapPtr cell = map[col][row];
  float len = channel_grid_cell_length(map, col, row);

  /* 
     if (mass > 0 && len <= 0.0) {
     error_handler(ERRHDL_ERROR,
     "channel_grid_inc_inflow: attempt to add flow in cell with no channels! (col=%d, row=%d)", 
     col, row);
     return;
     }
   */

  while (cell != NULL) {
    cell->channel->lateral_inflow += mass * cell->length / len;
    cell = cell->next;
  }
}
/* -------------------------------------------------------------
channel_grid_inc_melt
Given a flow (or actually mass), this function increases the inflow contributed
from melt to any channel(s) in the cell in proportion to their length within the
cell.
------------------------------------------------------------- */
void channel_grid_inc_melt(ChannelMapPtr ** map, int col, int row, float mass)
{
  ChannelMapPtr cell = map[col][row];
  float len = channel_grid_cell_length(map, col, row);

  /*
  if (mass > 0 && len <= 0.0) {
  error_handler(ERRHDL_ERROR,
  "channel_grid_inc_inflow: attempt to add flow in cell with no channels! (col=%d, row=%d)",
  col, row);
  return;
  }
  */

  while (cell != NULL) {
	cell->channel->melt += mass * cell->length / len;
	cell = cell->next;
  }
}
/* -------------------------------------------------------------
   channel_grid_outflow
   If the channel(s) within the cell are marked as ``sinks'', this
   function totals the mass from the channels(s) and returns the total
   mass.
   ------------------------------------------------------------- */
double channel_grid_outflow(ChannelMapPtr ** map, int col, int row)
{
  ChannelMapPtr cell = map[col][row];
  double mass = 0.0;

  while (cell != NULL) {
    if (cell->sink) {
      mass += cell->channel->outflow;
    }
    cell = cell->next;
  }
  return mass;
}

/* -------------------------------------------------------------
 channel_grid_init_table
 ------------------------------------------------------------- */
void channel_grid_init_table(ChannelMapPtr ** map, int col, int row,
                             float GridTableDepth)
{
  ChannelMapPtr cell = map[col][row];
  
  while (cell != NULL) {
    cell->table_depth = GridTableDepth;
    cell = cell->next;
  }
}

/* -------------------------------------------------------------
 channel_grid_table_depth
 Update local water table below channel based on infiltration rate
 and return shallowest local table depth
 ------------------------------------------------------------- */
float channel_grid_table_depth(ChannelMapPtr ** map, int col, int row, int deltat,
                               float GridTableDepth, float Transmissivity,
                               float SoilDeficit, float DX)
{
  ChannelMapPtr cell = map[col][row];
  float table_depth_min = GridTableDepth;
  float eff_dist, drop, grad, flow;
  
  while (cell != NULL) {
    if (GridTableDepth > cell->cut_height) {
      
      /* Compute gradient from halfway between edge of cell and edge of channel */
      /* But enforce upper limit of 1 m for steepest gradient calculation */
      /* Assumes DX == DY */
      eff_dist = (DX - cell->cut_width) / 4;
      if (eff_dist < 1.0)
        eff_dist = 1.0;
      
      /* Water flows away from both sides of infiltration zone */
      drop = (GridTableDepth - cell->table_depth) / eff_dist;
      if (drop < 0.0)
        drop = 0.0;
      grad = drop * (cell->length * 2);
      flow = Transmissivity * grad * deltat;
      
      flow /= cell->length * cell->cut_width; /* m^3 to m water depth */
      cell->table_depth += flow / SoilDeficit; /* water depth to table depth */
    }
    else
      cell->table_depth = cell->cut_height;
    
    /* Must remain between bottom of channel cut and grid cell water table */
    if (cell->table_depth > GridTableDepth)
      cell->table_depth = GridTableDepth;
    if (cell->table_depth < cell->cut_height)
      cell->table_depth = cell->cut_height;
    
    if (cell->table_depth < table_depth_min)
      table_depth_min = cell->table_depth;
    cell = cell->next;
  }
  return (table_depth_min);
}

/* -------------------------------------------------------------
 channel_grid_calc_infiltration
 If the channel(s) within the cell are above the water table,
 this function calculates the POTENTIAL infiltration from the channel(s)
 so that it can be subtracted during stream network routing.
 ------------------------------------------------------------- */
void channel_grid_calc_infiltration(ChannelMapPtr ** map, int col, int row, int deltat,
                                    float TableDepth, float MaxInfiltrationCap, float DX)
{
  ChannelMapPtr cell = map[col][row];
  float infiltration; /* From one channel segment */
  float water_depth, gradient, max_infiltration;
  float eff_dist_x, eff_dist_z, eff_dist;
  
  while (cell != NULL) {
    
    if (MaxInfiltrationCap > 0.0) {
      
      /* Infiltration rate is assumed equal to vertical hydraulic conductivity */
      /* Infiltration volume = conductivity * gradient * area * time */
      /* Channel head = water surface elev in channel - water table elev */
      water_depth = ((cell->channel->storage + cell->channel->last_storage) / 2.0) /
                    (cell->channel->class2->width * cell->channel->length);
      if (water_depth > cell->cut_height)
        water_depth = cell->cut_height;
      
      /* Compute gradient from halfway between edge of cell and edge of channel */
      /* But enforce upper limit of 1 m for steepest gradient calculation */
      /* Assumes DX == DY */
      eff_dist_x = (DX - cell->cut_width) / 4;
      eff_dist_z = TableDepth - (cell->cut_height - water_depth);
      eff_dist = sqrt(eff_dist_x*eff_dist_x + eff_dist_z*eff_dist_z);
      if (eff_dist < 1.0)
        eff_dist = 1.0;
      
      /* Unlike analogous gradients used in other functions,
         the infiltration gradient is along a diagonal direction
         to account for variations in the relative importance
         of vertical or lateral flow paths from the channel bottom
         to the rest of the grid cell */
      gradient = eff_dist_z / eff_dist;
      
      infiltration = cell->infiltration_rate * gradient * cell->length * cell->cut_width * deltat;
      
      /* Only allow enough infiltration to raise water table to channel bottom
         for the area directly below channel so that infiltration rates
         are not dependent on the choice of grid resolution */
      max_infiltration = MaxInfiltrationCap * cell->length * cell->cut_width;
      
      if (infiltration > max_infiltration)
        infiltration = max_infiltration;
      if (infiltration < 0.0)
        infiltration = 0.0;
      
      cell->infiltration = infiltration;
      
      /* Water table gets shallower below channel in proportion to how much of the capacity was filled */
      cell->table_depth -= (infiltration / max_infiltration) * (cell->table_depth - cell->cut_height);
      
    } else
      cell->infiltration = 0.0;
    
    cell = cell->next;
  }
}

/* -------------------------------------------------------------
 channel_grid_infiltration
 This function returns the sum of all ACTUAL infiltration
 after the updated routing method enforces constraints on
 the maximum amount of water available from inflow and storage.
 ------------------------------------------------------------- */
float channel_grid_infiltration(ChannelMapPtr ** map, int col, int row)
{
  ChannelMapPtr cell = map[col][row];
  float tot_infiltration = 0.0; /* From all channel segments in cell */
  
  while (cell != NULL) {
    
    if (cell->channel->remaining_infil > 0.0){
      
      if (cell->infiltration > cell->channel->remaining_infil)
        cell->infiltration = cell->channel->remaining_infil;
      
      cell->channel->remaining_infil -= cell->infiltration;
      cell->avail_water -= cell->infiltration;
      tot_infiltration += cell->infiltration;
    } else
      cell->infiltration = 0.0;
    
    cell = cell->next;
  }
  return tot_infiltration;
}

/* -------------------------------------------------------------
 channel_grid_calc_evaporation
 If the channel(s) within the cell have nonzero water storage,
 this function totals the mass of evaporation from the channel(s),
 assuming that sufficient water is available (handled later)
 ------------------------------------------------------------- */
void channel_grid_calc_evaporation(ChannelMapPtr ** map, int col, int row,
                                    float EPot, float MaxEvapCap)
{
  ChannelMapPtr cell = map[col][row];
  
  while (cell != NULL) {
    
    cell->evaporation = EPot * cell->cut_width * cell->length;
    
    if (cell->evaporation > MaxEvapCap)
      cell->evaporation = MaxEvapCap;
    if (cell->evaporation < 0.0)
      cell->evaporation = 0.0;
    
    /* Keep track of remaining total grid cell EPot */
    MaxEvapCap -= cell->evaporation;
    if (MaxEvapCap < 0.0)
      MaxEvapCap = 0.0;
    
    cell = cell->next;
  }
}


/* -------------------------------------------------------------
 channel_grid_evaporation
 This function returns the sum of all ACTUAL evaporation
 after the updated routing method enforces constraints on
 the maximum amount of water available from inflow and storage.
 ------------------------------------------------------------- */
float channel_grid_evaporation(ChannelMapPtr ** map, int col, int row)
{
  ChannelMapPtr cell = map[col][row];
  float tot_evap = 0.0; /* From all channel segments in cell */
  
  while (cell != NULL) {
    
    if (cell->channel->remaining_evap > 0.0){
      
      if (cell->evaporation > cell->channel->remaining_evap)
        cell->evaporation = cell->channel->remaining_evap;
      
      cell->channel->remaining_evap -= cell->evaporation;
      cell->avail_water -= cell->evaporation;
      tot_evap += cell->evaporation;
    } else
      cell->evaporation = 0.0;
    
    cell = cell->next;
  }
  return tot_evap;
}

/* -------------------------------------------------------------
 channel_grid_dry_evaporation
 If any channel(s) within the cell have ZERO water storage,
 this function totals the soil evaporation from each segment.
 ------------------------------------------------------------- */
float channel_grid_dry_evaporation(ChannelMapPtr ** map, int col, int row,
                                   float EPot, float MaxEvapCap, float DXDY,
                                   float Dt, float *Porosity, float *FCap, float *Ks,
                                   float *Press, float *m, float LayerThickness,
                                   float *MoistContent, float *Adjust, int CutBankZone)
{
  ChannelMapPtr cell = map[col][row];
  /* Note: unlike the similar channel_grid_evaporation function,
     SoilEvap and TotalSoilEvap are in depth units (m), not volumes */ 
  float SoilEvap; /* From one channel segment */
  float TotalSoilEvap = 0.0; /*  From all channel segments in cell */
  float DesorptionVolume;
  float tmp;
  float SoilMoisture;
  
  while (cell != NULL) {
    
    if (fequal(cell->channel->storage, 0.0) && 
        MoistContent[CutBankZone] > FCap[CutBankZone]) {
      
      DesorptionVolume = Desorption(Dt, MoistContent[CutBankZone], Porosity[CutBankZone],
                                    Ks[CutBankZone], Press[CutBankZone], m[CutBankZone]);
      
      SoilEvap = MIN(EPot, DesorptionVolume);
      
      if (SoilEvap > MaxEvapCap)
        SoilEvap = MaxEvapCap;
      
      /* Re-scale from channel to grid cell depth */
      SoilEvap *= (cell->cut_width * cell->length) / DXDY;
      
      SoilMoisture = MoistContent[CutBankZone] * LayerThickness * Adjust[CutBankZone];
      tmp = FCap[CutBankZone] * LayerThickness * Adjust[CutBankZone];
      
      if (SoilEvap > SoilMoisture - tmp) {
        SoilEvap = SoilMoisture - tmp;
        MoistContent[CutBankZone] = FCap[CutBankZone];
      }
      else {
        SoilMoisture -= SoilEvap;
        MoistContent[CutBankZone] = SoilMoisture / (LayerThickness * Adjust[CutBankZone]);
      }
      
      TotalSoilEvap += SoilEvap;
      
      /* Keep track of remaining total grid cell EPot */
      MaxEvapCap -= SoilEvap;
      if (MaxEvapCap < 0.0)
        MaxEvapCap = 0.0;
    }
    cell = cell->next;
  }
  return TotalSoilEvap;
}

/* -------------------------------------------------------------
   channel_grid_flowlength
   returns the flowlength along a road surface in a channel
   if there is more than one road in a grid cell, the road
   with the greatest surface area is used to calculate the 
   flowlength.
   This can result in a flolen that is greater than the 
   horizontal length of the road in the cell. 
  ------------------------------------------------------------- */
double channel_grid_flowlength(ChannelMapPtr ** map, int col, int row, float floslope)
{
  ChannelMapPtr cell = map[col][row];
  double flolen = 0.0;
  double area;
  double maxarea = 0.0;

  while (cell != NULL) {
    area = cell->length * cell->cut_width;
    if(area > maxarea){
      flolen = cell->cut_width * (floslope/ROADCROWN)*sqrt(1+pow(ROADCROWN,2));
      maxarea = area;
    }
    if(flolen < cell->cut_width)
      flolen = cell->cut_width;
    /* If crowned, only one half goes to ditch. */ 
    if (cell->channel->class2->crown == CHAN_CROWNED) 
      flolen *= 0.5;
    
    cell = cell->next;
  }
  return flolen;
}

/* -------------------------------------------------------------
   channel_grid_flowslope
   returns the flowlength along a road surface in a grid cell
   if there is more than one road in a grid cell, the road
   with the greatest surface area is used to calculate the 
   flowslope
   ------------------------------------------------------------- */

double channel_grid_flowslope(ChannelMapPtr ** map, int col, int row)
{
  ChannelMapPtr cell = map[col][row];
  float floslope = 0.0;
  double area;
  double maxarea = 0.0;

  while (cell != NULL) {
    area = cell->length * cell->cut_width;
    if(area > maxarea){
      floslope = sqrt(pow(ROADCROWN, 2) + pow((cell->channel->slope),2));
      maxarea = area;
    }
    cell = cell->next;
  }
  return floslope;
}

/* -------------------------------------------------------------
   channel_grid_class
   returns the erodibility coeffienct of the road surface in a
   gird cell. if there is more than one road in a grid cell, the road
   with the greatest surface area is used
   ------------------------------------------------------------- */

ChannelClass* channel_grid_class(ChannelMapPtr ** map, int col, int row)
{
  ChannelMapPtr cell = map[col][row];
  ChannelClass *pntr = NULL;
  double area;
  double maxarea = 0.0;

  while (cell != NULL) {
    area = cell->length * cell->cut_width;
    if(area > maxarea){
      pntr = cell->channel->class2;
      maxarea = area;
    }
    cell = cell->next;
  }
  return pntr;
}

/* -------------------------------------------------------------
   ---------------------- Module Functions ---------------------
   ------------------------------------------------------------- */

/* -------------------------------------------------------------
   channel_grid_init
   ------------------------------------------------------------- */
void channel_grid_init(int cols, int rows)
{
  channel_grid_cols = cols;
  channel_grid_rows = rows;
  channel_grid_initialized = 1;
}

/* -------------------------------------------------------------
   channel_grid_done
   ------------------------------------------------------------- */
void channel_grid_done(void)
{
  /* ? */
}

#ifdef TEST_MAIN
/* -------------------------------------------------------------
   interpolate
   ------------------------------------------------------------- */
static float interpolate(int n, float *x, float *y, float x0)
{
  int i;
  if (x0 <= x[0]) {
    return ((x0 - x[0]) / (x[1] - x[0]) * (y[1] - y[0]) + y[0]);
  }
  for (i = 0; i < n - 1; i++) {
    if (x0 < x[i + 1]) {
      return ((x0 - x[i]) / (x[i + 1] - x[i]) * (y[i + 1] - y[i]) + y[i]);
    }
  }
  return ((x0 - x[i - 1]) / (x[i] - x[i - 1]) * (y[i] - y[i - 1]) + y[i]);
}

/* -------------------------------------------------------------
   Main Program
   ------------------------------------------------------------- */
int main(int argc, char **argv)
{
  static const int columns = 5;
  static const int rows = 10;

  int r, c;

  ChannelClass *class;
  Channel *simple = NULL, *current;
  ChannelMapPtr **map = NULL;

  static int interval = 3600;	/* seconds */
  static float timestep = 1.0;	/* hour */
  static float endtime = 144.0;
#define TIMES 6
  static float bndflow[TIMES] = { 0.0, 0.0, 300.0, 300.0, 0.0, 0.0 };
  static float bndtime[TIMES] = { 0.0, 12.0, 36.0, 48.0, 60.0, 1000.0 };
  float time;

  /* module initialization */

  error_handler_init(argv[0], NULL, ERRHDL_DEBUG);
  channel_init();
  channel_grid_init(columns, rows);

  /* read channel classes */

  if ((class = channel_read_classes("example_classes.dat")) == NULL) {
    error_handler(ERRHDL_FATAL, "example_classes.dat: trouble reading file");
  }

  /* read a network */

  if ((simple = channel_read_network("example_network.dat", class)) == NULL) {
    error_handler(ERRHDL_FATAL, "example_network.dat: trouble reading file");
  }

  /* read channel map */

  if ((map = channel_grid_read_map(simple, "example_map.dat")) == NULL) {
    error_handler(ERRHDL_FATAL, "example_map.dat: trouble reading file");
  }

  /* check channel_grid_read_map */

  printf("channel_grid_read_map check:\n");
  for (r = rows - 1; r >= 0; r--) {
    for (c = 0; c < columns; c++) {
      ChannelMapPtr cell = map[c][r];
      int count;
      for (count = 0; cell != NULL; cell = cell->next) {
	count++;
      }
      printf("%3d", count);
    }
    printf("\n");
  }
  printf("\n");

  /* check channel_grid_cell_length */

  printf("channel_grid_cell_length check:\n");
  for (r = rows - 1; r >= 0; r--) {
    for (c = 0; c < columns; c++) {
      printf(" %8.2g", channel_grid_cell_length(map, c, r));
    }
    printf("\n");
  }
  printf("\n");

  /* use routing example to test channel_grid_inc_inflow and
     channel_grid_outflow */

  /* initialize flows */

  for (current = simple; current != NULL; current = current->next) {
    current->inflow = bndflow[0] * timestep;
    current->outflow = bndflow[0] * timestep;
  }

  /* time loop */

  for (time = 0.0; time <= endtime; time += timestep) {
    float inflow = interpolate(TIMES, bndtime, bndflow, time) * interval;
    float outflow;

    channel_step_initialize_network(simple);
    channel_grid_inc_inflow(map, 2, 0, inflow);
    (void) channel_route_network(simple, interval);
    outflow = channel_grid_outflow(map, 2, 6);
    channel_save_outflow(time * interval, simple, stdout);
    printf("outflow: %8.3g\n", outflow);
  }

  /* deallocate memory */

  channel_grid_free_map(map);
  channel_free_network(simple);
  channel_free_classes(class);

  /* module shutdown */

  channel_grid_done();
  channel_done();
  error_handler_done();
  exit(0);

}
#endif

/* -----------------------------------------------------------------------------------
   Function name : channel_grid_inc_other ()
   Author: Nathalie Voisin Aug 2010
   Function usage:
   This function generates outputs required by John's RBM model

   Given a mass/flux, this function increases the mass/flux in any
   channel(s) is in proportion to their length within the cell 
   BECAUSE THOSE ARE NRG fluxes, similar to each segment
   --------------------------------------------------------------------------------- */

void channel_grid_inc_other(ChannelMapPtr ** map, int col, int row, PIXRAD * LocalRad, 
							PIXMET * LocalMet, float skyview)
{
  ChannelMapPtr cell = map[col][row];
  
  while (cell != NULL ) {
	/* ISW is the total incoming shortwave radiation (VIC outputs) */
	cell->channel->ISW += LocalRad->ObsShortIn;
	
	cell->channel->NSW += LocalRad->RBMNetShort;
	cell->channel->Beam += LocalRad->PixelBeam;
	cell->channel->Diffuse += LocalRad->PixelDiffuse;

    cell->channel->ILW += LocalRad->PixelLongIn;
	cell->channel->NLW += LocalRad->RBMNetLong;

    cell->channel->VP += LocalMet->Eact;
    cell->channel->WND += LocalMet->Wind;
    cell->channel->ATP += LocalMet->Tair;

	cell->channel->azimuth += 
		cell->azimuth*cell->length /cell->channel->length;

	cell->channel->skyview += skyview;

    cell = cell->next;
  }
}
/*************************************************************************************
void channel_grid_avg( ): average the heat budget variables by the total cell numbers
*************************************************************************************/
void channel_grid_avg(Channel *Channel)
{  
  while (Channel) {
    if (Channel->Ncells > 0 ) {
	  Channel->ISW /= Channel->Ncells ;
	  
	  Channel->NSW /= Channel->Ncells;
	  Channel->Beam /= Channel->Ncells;
	  Channel->Diffuse /= Channel->Ncells;

	  Channel->ILW /= Channel->Ncells ;
	  Channel->NLW /= Channel->Ncells ;
      
      Channel->VP  /= Channel->Ncells;
      Channel->WND /= Channel->Ncells;
	  Channel->ATP /= Channel->Ncells;

	  Channel->skyview /= Channel->Ncells;
    }
	Channel = Channel->next; 
  }
}
/*********************************************************************************
Init_segment_ncell : computes the number of grid cell contributing to one segment
**********************************************************************************/
void Init_segment_ncell(TOPOPIX **TopoMap, ChannelMapPtr ** map, int NY, 
						int NX, Channel* net)
{
  int y,x;
  ChannelMapPtr cell; 

  for (y = 0; y < NY; y++) {
    for (x = 0; x < NX; x++) {      
	  if (INBASIN(TopoMap[y][x].Mask)) {
		if (channel_grid_has_channel(map, x, y)){	
           cell = map[x][y];
		   while (cell != NULL) {
			 cell->channel->Ncells++;
             cell = cell->next;
		   }   
        }
      }
    }
  }

  // then check all segments
  for (; net != NULL; net = net->next) {
    if (net->Ncells == 0 ) {
      error_handler(ERRHDL_ERROR,"Init_segment_ncells: write error:%s", strerror(errno));
    }
  }
} 



