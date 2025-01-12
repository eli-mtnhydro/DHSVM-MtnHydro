/* -------------------------------------------------------------
file: channel.c
------------------------------------------------------------- */
/* -------------------------------------------------------------
Battelle Memorial Institute
Pacific Northwest Laboratory
------------------------------------------------------------- */
/* -------------------------------------------------------------
Created October 24, 1995 by  William A Perkins
$Id: channel.c,v3.1.2 2014/1/2 Ning Exp $
------------------------------------------------------------- */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include "errorhandler.h"
#include "DHSVMerror.h"
#include "channel.h"
#include "channel_grid.h"
#include "constants.h"
#include "tableio.h"
#include "settings.h"

/* for test msw */
#define TEST_MAIN 0
/* end test */

/* -------------------------------------------------------------
-------------- ChannelClass Functions -----------------------
------------------------------------------------------------- */

/* -------------------------------------------------------------
alloc_channel_class
------------------------------------------------------------- */
static ChannelClass *alloc_channel_class(void)
{
  ChannelClass *p;

  if ((p = (ChannelClass *) malloc(sizeof(ChannelClass))) == NULL) {
    error_handler(ERRHDL_ERROR, "alloc_channel_class: malloc failed: %s",
      strerror(errno));
    return NULL;
  }

  p->id = 0;
  p->width = 0.0;
  p->bank_height = 0.0;
  p->friction = 0.0;
  p->infiltration = 0.0;
  p->crown = CHAN_OUTSLOPED;
  p->next = (ChannelClass *) NULL;

  return p;
}

/* -------------------------------------------------------------
channel_free_classes
------------------------------------------------------------- */
void channel_free_classes(ChannelClass * head)
{
  if (head->next != NULL) {
    channel_free_classes(head->next);
  }
  free(head);
}

/* -------------------------------------------------------------
find_channel_class
------------------------------------------------------------- */
static ChannelClass *find_channel_class(ChannelClass * list, ClassID id)
{
  while (list != (ChannelClass *) NULL) {
    if (list->id == id)
      break;
    list = list->next;
  }
  return list;
}

/* -------------------------------------------------------------
channel_read_classes
This function opens and reads the specified file and returns a
linked list of ChannelClass structs.  If anything goes wrong, NULL
is returned and any ChannelClass structs are destroyed.
------------------------------------------------------------- */
ChannelClass *channel_read_classes(const char *file, int ChanType)
{
  ChannelClass *head = NULL, *current = NULL;
  static const int fields = 6;
  int done;
  int err = 0;
  static char *crown_words[4] = {
    "OUTSLOPED", "CROWNED", "INSLOPED", NULL
  };

  static TableField class_fields[6] = {
    {"ID", TABLE_INTEGER, TRUE, FALSE, {0}, "", NULL},
    {"Channel Width", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
    {"Bank (stream) or Cut Height (road)", TABLE_REAL, TRUE, FALSE, {0.0}, "",
    NULL},
    {"Friction Coefficient (Manning's n)", TABLE_REAL, TRUE, FALSE, {0.0}, "",
    NULL},
    {"Maximum Road Infiltration Rate (m/s)", TABLE_REAL, FALSE, FALSE, {0.0}, 
    "", NULL},
    {"Road Crown Type", TABLE_WORD, FALSE, FALSE, {0}, "", crown_words}
  };

  // Extra fields are required if we're dealing with a road
  if (ChanType == road_class) {
    // max infiltration rate and crown type are required
    class_fields[4].required = TRUE;
    class_fields[5].required = TRUE;
  }

  error_handler(ERRHDL_STATUS,
    "channel_read_classes: reading file \"%s\"", file);

  if (table_open(file) != 0) {
    error_handler(ERRHDL_ERROR,
      "channel_read_classes: unable to open file \"%s\": %s",
      file, strerror(errno));
    return NULL;
  }

  done = FALSE;
  while (!done) {
    int i;
    done = (table_get_fields(fields, class_fields) < 0);
    if (done) {
      for (i = 0; i < fields; i++) {
        if (class_fields[i].read)
          break;
      }
      if (i >= fields)
        continue;
    }

    if (head == NULL) {
      head = alloc_channel_class();
      current = head;
    }
    else {
      current->next = alloc_channel_class();
      current = current->next;
    }

    for (i = 0; i < fields; i++) {
      if (class_fields[i].read) {
        switch (i) {
        case 0:
          current->id = class_fields[i].value.integer;
          if (current->id <= 0) {
            error_handler(ERRHDL_ERROR,
              "%s: class %d: class id invalid", file, current->id);
            err++;
          }
          break;
          // We do not want these 3 values to be 0.0 EVER
          // Exit if this is the case.
        case 1:
          if (class_fields[i].value.real>0.0)
            current->width = class_fields[i].value.real;
          else error_handler(ERRHDL_FATAL,
            "channel_read_classes: %s: width cannot be 0.0", file);
          break;
        case 2:
          if (class_fields[i].value.real>0.0)
            current->bank_height = class_fields[i].value.real;
          else error_handler(ERRHDL_FATAL,
            "channel_read_classes: %s: bank cannot be 0.0", file);
          break;
        case 3:
          if (class_fields[i].value.real>0.0)
            current->friction = class_fields[i].value.real;
          else error_handler(ERRHDL_FATAL,
            "channel_read_classes: %s: friction cannot be 0.0", file);
          break;
        case 4:
          current->infiltration = class_fields[i].value.real;
          break;
        case 5:
          switch (class_fields[i].value.integer) {
          case -1:
            error_handler(ERRHDL_ERROR,
              "channel_read_classes: %s: unknown road crown type: %s",
              file, class_fields[i].field);
            err++;
            break;
          case 0:
            current->crown = CHAN_OUTSLOPED;
            break;
          case 1:
            current->crown = CHAN_CROWNED;
            break;
          case 2:
            current->crown = CHAN_INSLOPED;
            break;
          default:
            error_handler(ERRHDL_FATAL,
              "channel_read_classes: this should not happen");
          }
          break;
        default:
          error_handler(ERRHDL_FATAL,
            "channel_read_classes: this should not happen either");
        }
      }
    }
  }
  error_handler(ERRHDL_STATUS,
    "channel_read_classes: %s: %d errors, %d warnings",
    file, table_errors, table_warnings);

  table_close();

  error_handler(ERRHDL_STATUS,
    "channel_read_classes: done reading file \"%s\"", file);

  if (table_errors) {
    error_handler(ERRHDL_ERROR,
      "channel_read_classes: %s: too many errors", file);
    channel_free_classes(head);
    head = NULL;
  }

  return (head);
}

/* -------------------------------------------------------------
---------------------- Channel Functions --------------------
------------------------------------------------------------- */

/* -------------------------------------------------------------
alloc_channel_segment
------------------------------------------------------------- */
static Channel *alloc_channel_segment(void)
{
  Channel *seg;
  if ((seg = (Channel *) malloc(sizeof(Channel))) == NULL) {
    error_handler(ERRHDL_ERROR, "alloc_channel_segment: malloc failed: %s",
      strerror(errno));
    return NULL;
  }
  seg->id = 0;
  seg->outid = 0;
  seg->order = 0;
  seg->record_name = NULL;
  seg->record = FALSE;
  seg->length = 0.0;
  seg->slope = 0.0;
  seg->class2 = NULL;
  seg->lateral_inflow = 0.0;
  seg->melt = 0.0;                  
  seg->last_inflow = 0.0;
  seg->last_outflow = 0.0;
  seg->inflow = 0.0;
  seg->outflow = 0.0;
  seg->storage = 0.0;
  seg->infiltration = 0.0;
  seg->remaining_infil = 0.0;
  seg->evaporation = 0.0;
  seg->outlet = NULL;
  seg->next = NULL;
  seg->grid = NULL;

  /* Initialize the variables required by John's RBM model */
  seg->ISW = 0.;   /* incident shortwave radiation */
  seg->Beam = 0.;
  seg->Diffuse = 0.;
  seg->ILW = 0.;   /* incident longwave radiation */
  seg->NSW = 0.;   /* net shortwave radiation */
  seg->NLW = 0.;   /* net longwave radiation */
  seg->VP = 0.;     /* actual vapor pressure */
  seg->WND = 0.;
  seg->ATP = 0.;
  seg->Ncells = 0; /* not used for now */
  seg->azimuth = 0;
  seg->skyview = 0;

  return seg;
}

/* -------------------------------------------------------------
channel_find_segment
A simple linear search of the channel network to find a segment
with the given id
------------------------------------------------------------- */
Channel *channel_find_segment(Channel * head, SegmentID id)
{
  for (; head != NULL; head = head->next) {
    if (head->id == id)
      break;
  }
  if (head == NULL) {
    error_handler(ERRHDL_WARNING,
      "channel_find_segment: unable to find segment %d", id);
  }
  else {
    error_handler(ERRHDL_DEBUG, "channel_find_segment: found segment %d", id);
  }

  return head;
}

/* -------------------------------------------------------------
channel_routing_parameters
------------------------------------------------------------- */
void channel_routing_parameters(Channel *network, int deltat)
{
  float Rh;
  Channel *segment;

  for (segment = network; segment != NULL; segment = segment->next) {
    
    /* Initialized assuming uniform depth */
    segment->slope = segment->ground_slope;
    
    Rh = segment->class2->bank_height * 0.75;
    segment->K = sqrt(segment->slope) * pow((double) Rh, 2.0 / 3.0) /
                 (segment->class2->friction * segment->length);
    segment->X = exp(-segment->K * deltat);
  }
}

/* -------------------------------------------------------------
 channel_update_routing_parameters
 ------------------------------------------------------------- */
void channel_update_routing_parameters(Channel *network, int deltat, int max_order)
{
  float WaterDepth; /* Length-average depth in channel (m) */
  float Kold, Rh;
  Channel *segment;
  int order;
  
  for (order = max_order; order >= 1; order--) {
    segment = network;
    while (segment != NULL) {
      if ((int) segment->order == order) {
        
        Kold = segment->K;
        
        if (segment->storage > 0.0) {
          
          WaterDepth = ((segment->storage + segment->last_storage) / 2.0) / (segment->class2->width * segment->length);
          
          /* Assume uniform water depth in outlet segment, and back-propagate uphill */
          if (segment->outlet != NULL)
            segment->bottom_water_depth = segment->outlet->top_water_depth;
          else
            segment->bottom_water_depth = WaterDepth;
          
          /* Assume bottom depth is no greater than cut depth
           (perhaps extra flow expands the channel width into a floodplain) */
          if (segment->bottom_water_depth > segment->class2->bank_height)
            segment->bottom_water_depth = segment->class2->bank_height;
          
          /* Enforce topographic constraint so the network doesn't "flatten" upstream */
          if (segment->bottom_water_depth < 0.0)
            segment->bottom_water_depth = 0.0;
          
          segment->top_water_depth = 2 * WaterDepth - segment->bottom_water_depth;
          
          segment->slope = segment->ground_slope;
          segment->slope += (segment->top_water_depth - segment->bottom_water_depth) / segment->length;
          
          /* Don't allow slope to exceed 45 degrees unless ground slope is steeper */
          if (segment->slope > 1.0 && segment->slope > segment->ground_slope)
            segment->slope = segment->ground_slope;
          
          Rh = (WaterDepth * segment->class2->width) / (WaterDepth * 2 + segment->class2->width);
          
          if (segment->slope > 0.0)
            segment->K = sqrt(segment->slope) * pow((double) Rh, 2.0 / 3.0) /
              (segment->class2->friction * segment->length);
          else
            segment->K = 0.0;
        } else
          segment->K = 0.0;
        
        /* Use moving average to prevent numeric oscillation across timesteps */
        segment->K = (segment->K + Kold) / 2.0;
        
        if (segment->K < MINSTORAGEK)
          segment->K = MINSTORAGEK;
        
        segment->X = exp(-segment->K * deltat);
      }
      segment = segment->next;
    }
  }
}

/* -------------------------------------------------------------
channel_read_network
------------------------------------------------------------- */
Channel *channel_read_network(const char *file, ChannelClass *class_list, int *MaxID)
{
  Channel *head = NULL, *current = NULL;
  int err = 0;
  int done;
  static const int fields = 8;
  static char *save_words[2] = {
    "SAVE", "\0"
  };
  static TableField chan_fields[8] = {
    {"ID", TABLE_INTEGER, TRUE, FALSE, {0}, "", NULL},
    {"Order", TABLE_INTEGER, TRUE, FALSE, {0}, "", NULL},
    {"Slope", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
    {"Length", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
    {"Class", TABLE_INTEGER, TRUE, FALSE, {0}, "", NULL},
    {"Outlet ID", TABLE_INTEGER, FALSE, FALSE, {0}, "", NULL},
    {"Save Flag", TABLE_WORD, FALSE, FALSE, {0}, "", save_words},
    {"Save Name", TABLE_STRING, FALSE, FALSE, {0.0}, "", NULL}
  };

  error_handler(ERRHDL_STATUS,
    "channel_read_network: reading file \"%s\"", file);

  if (table_open(file) != 0) {
    error_handler(ERRHDL_ERROR,
      "channel_read_network: unable to open file \"%s\": %s",
      file, strerror(errno));
    return NULL;
  }

  *MaxID = 0;
  done = FALSE;
  while (!done) {
    int i;

    done = (table_get_fields(fields, chan_fields) < 0);

    if (done) {
      for (i = 0; i < fields; i++) {
        if (chan_fields[i].read)
          break;
      }
      if (i >= fields)
        continue;
    }

    if (head == NULL) {
      head = alloc_channel_segment();
      current = head;
    }
    else {
      current->next = alloc_channel_segment();
      current = current->next;
    }

    for (i = 0; i < fields; i++) {
      if (chan_fields[i].read) {
        switch (i) {
        case 0:
          current->id = chan_fields[i].value.integer;
          if(current->id > *MaxID) *MaxID = current->id;
          if (current->id <= 0) {
            error_handler(ERRHDL_ERROR,
              "%s: segment %d: channel id invalid",
              file, current->id);
            err++;
          }
          break;
        case 1:
          if (chan_fields[i].value.integer > 0) {
            current->order = chan_fields[i].value.integer;
          }
          else {
            error_handler(ERRHDL_ERROR,
              "%s: segment %d: channel order (%d) invalid",
              file, current->id, chan_fields[i].value.integer);
            err++;
          }
          break;
        case 2:
          if (chan_fields[i].value.real > 0) {
            current->ground_slope = chan_fields[i].value.real;
          }
          else {
            error_handler(ERRHDL_ERROR,
              "%s: segment %d: channel slope (%f) invalid",
              file, current->id, chan_fields[i].value.real);
            err++;
          }
          break;
        case 3:
          if (chan_fields[i].value.real > 0) {
            current->length = chan_fields[i].value.real;
          }
          else {
            error_handler(ERRHDL_ERROR,
              "%s: segment %d: channel length (%f) invalid",
              file, current->id, chan_fields[i].value.real);
            err++;
          }
          break;
        case 4:
          if ((current->class2 =
            find_channel_class(class_list, chan_fields[i].value.integer)
            ) == NULL) {
              error_handler(ERRHDL_ERROR,
                "%s: segment %d: channel class %d not found",
                file, current->id, chan_fields[i].value.integer);
              err++;
          }
          break;
        case 5:
          current->outid = chan_fields[i].value.integer;
          break;
        case 6:
          current->record = TRUE;
          break;
        case 7:
          current->record_name = (char *) strdup(chan_fields[i].field);
          break;
        default:
          error_handler(ERRHDL_FATAL,
            "channel_read_network: what is this field %d?", i);
          break;
        }
      }
    }
  }

  table_close();

  /* find segment outlet segments, if
  specified */

  for (current = head; current != NULL; current = current->next) {
    if (current->outid != 0) {
      current->outlet = channel_find_segment(head, current->outid);
      if (current->outlet == NULL) {
        error_handler(ERRHDL_ERROR,
          "%s: cannot find outlet (%d) for segment %d",
          file, current->outid, current->id);
        err++;
      }
    }
  }

  table_errors += err;

  error_handler(ERRHDL_STATUS,
    "channel_read_network: %s: %d errors, %d warnings",
    file, table_errors, table_warnings);

  if (table_errors) {
    error_handler(ERRHDL_ERROR,
      "channel_read_network: %s: too many errors", file);
    channel_free_network(head);
    head = NULL;
  }

  return (head);
}

/* -------------------------------------------------------------
 channel_segment_infil_evap
------------------------------------------------------------- */
void channel_segment_infil_evap(Channel * segment)
{
  ChannelMapPtr cell = segment->grid;
  float max_avail;
  
  /* Prevent issues with floating point arithmetic */
  if (segment->storage < 1e-6 * (segment->inflow + segment->lateral_inflow) &&
      segment->storage < 1.0)
    segment->storage = 0.0;
  
  segment->infiltration = 0.0;
  segment->evaporation = 0.0;
  max_avail = segment->storage + segment->inflow + segment->lateral_inflow;
  
  while (cell != NULL) {
    
    /* Update available storage with water entering segment from upstream */
    cell->avail_water += segment->inflow;
    
    /* Limit available water to a maximum of the latest inflow and storage,
       only considering amount originating uphill of current grid cell */
    if (cell->avail_water > max_avail)
      cell->avail_water = max_avail;
    
    /* Cell at bottom of channel must have all water available 
       (otherwise rounding error can cause water to get "stuck") */
    if (cell->next_seg == NULL)
      cell->avail_water = max_avail;
    
    /* If insufficient water is available to completely satisfy the pre-calculated
       infiltration capacity and potential evaporation, equally divide the
       available water between infiltration and evaporation,
       while ensuring that neither exceeds the maximum capacity/potential */
    if ((cell->infiltration + cell->evaporation) > cell->avail_water) {
      if (cell->infiltration > cell->evaporation) {
        cell->evaporation = MIN(cell->avail_water / 2.0, cell->evaporation);
        cell->infiltration = cell->avail_water - cell->evaporation;
      } else {
        cell->infiltration = MIN(cell->avail_water / 2.0, cell->infiltration);
        cell->evaporation = cell->avail_water - cell->infiltration;
      }
    }
    
    segment->infiltration += cell->infiltration;
    segment->evaporation += cell->evaporation;
    cell = cell->next_seg;
  }
  
  if (segment->infiltration < 0.0)
    segment->infiltration = 0.0;
  if (segment->evaporation < 0.0)
    segment->evaporation = 0.0;
  
  if ((segment->infiltration + segment->evaporation) > max_avail) {
    if (segment->infiltration > segment->evaporation) {
      segment->evaporation = MIN(max_avail / 2.0, segment->evaporation);
      segment->infiltration = max_avail - segment->evaporation;
    } else {
      segment->infiltration = MIN(max_avail / 2.0, segment->infiltration);
      segment->evaporation = max_avail - segment->infiltration;
    }
  }
  
  segment->remaining_infil = segment->infiltration;
  segment->remaining_evap = segment->evaporation;
}

/* -------------------------------------------------------------
channel_route_segment
------------------------------------------------------------- */
static int channel_route_segment(Channel * segment, int deltat)
{
  float net_balance, outflow, storage, storage_loss;
  int err = 0;
  
  if (segment->inflow < 0.0)
    segment->inflow = 0.0;
  if (segment->storage < 0.0)
    segment->storage = 0.0;
  
  /* Account for infiltration and evaporation prior to routing each segment */
  channel_segment_infil_evap(segment);
  
  /* Subtract any necessary water from storage instead of inflow */
  storage_loss = (segment->infiltration + segment->evaporation) -
                 (segment->inflow + segment->lateral_inflow);
  if (storage_loss < 0.0)
    storage_loss = 0.0;
  if (storage_loss > segment->storage)
    storage_loss = segment->storage;
  segment->storage -= storage_loss;
  
  /* Add up mass balance components and change masses to rates */
  net_balance = segment->inflow + segment->lateral_inflow -
                (segment->infiltration + segment->evaporation - storage_loss);
  net_balance /= deltat;
  
  if (segment->K > 1e-10)
    storage = (net_balance / segment->K) +
              (segment->storage - net_balance / segment->K) * segment->X;
  else
    storage = net_balance;
  
  if (storage < 0.0)
    storage = 0.0;
  
  outflow = net_balance - (storage - segment->storage) / deltat;
  
  if (outflow < 0.0)
    outflow = 0.0;
  
  segment->outflow = outflow * deltat;
  segment->storage = storage;
  
  if (segment->outlet != NULL)
    segment->outlet->inflow += segment->outflow;
  
  return (err);
}

/* -------------------------------------------------------------
channel_route_network
------------------------------------------------------------- */
int channel_route_network(Channel * net, int deltat)
{
  int order;
  int order_count;
  int err = 0;
  Channel *current;
  
  for (order = 1;; order++) {
    order_count = 0;
    current = net;
    while (current != NULL) {
      if ((int) current->order == order) {
        
        err += channel_route_segment(current, deltat);
        order_count++;
      }
      current = current->next;
    }
    if (order_count == 0)
      break;
  }
  
  channel_update_routing_parameters(net, deltat, (order - 1));
  
  return (err);
}

/* -------------------------------------------------------------
channel_step_initialize_network
------------------------------------------------------------- */
int channel_step_initialize_network(Channel *net)
{
  if (net != NULL) {
    
    net->last_inflow = net->inflow;
    net->inflow = 0.0;
    net->lateral_inflow = 0.0;
    net->melt = 0.0;                 
    net->last_outflow = net->outflow;
    net->last_storage = net->storage;
    net->infiltration = 0.0;
    net->remaining_infil = 0.0;
    net->evaporation = 0.0;

    /* Initialzie variables for John's RBM model */ 
    net->ILW = 0.; /* incident longwave radiation */
    net->NLW = 0.; /* net longwave radiation */
    net->ISW = 0.; /* incident shortwave radiation */
    net->Beam = 0.;
    net->Diffuse = 0.;
    net->NSW = 0.0;            /* net shortwave radiation with both topo and canopy shading */

    net->VP = 0.;              /* actual vapor pressure */
    net->WND = 0.;
    net->ATP = 0.;
    net->azimuth = 0;
    net->skyview = 0;
    //net->Ncells = 0; /* not used for now */

    channel_step_initialize_network(net->next);
  }
  return (0);
}

/* -------------------------------------------------------------
channel_save_outflow
This routine saves the channel output
------------------------------------------------------------- */
int channel_save_outflow(double time, Channel * net, FILE * out, FILE * out2)
{
  char buffer[16];
  sprintf(buffer, "%12.5g", time);
  return (channel_save_outflow_text(buffer, net, out, out2, 0));
}

/* -------------------------------------------------------------
channel_save_outflow_wtext
Saves the channel outflow using a text string as the time field
------------------------------------------------------------- */
int
  channel_save_outflow_text(char *tstring, Channel * net, FILE * out,
  FILE * out2, int flag)
{
  int err = 0;
  float total_outflow = 0.0;
  float total_lateral_inflow = 0.0;
  float total_storage = 0.0;
  float total_storage_change = 0.0;
  float total_error = 0.0;
  float total_infiltration = 0.0;
  float total_evaporation = 0.0;

  if (flag == 1) {
    fprintf(out2, "DATE ");
    for (; net != NULL; net = net->next) {
      total_lateral_inflow += net->lateral_inflow;
      if (net->outlet == NULL) {
        total_outflow += net->outflow;
      }
      if (net->record)
        fprintf(out2, "%s ", net->record_name);
    }
    fprintf(out2, "\n");
  }

  //tsstring = date in the form of 01.01.1915-00:00:00 
  if (fprintf(out2, "%15s ", tstring) == EOF) {
    error_handler(ERRHDL_ERROR,
      "channel_save_outflow: write error:%s", strerror(errno));
    err++;
  }

  for (; net != NULL; net = net->next) {
    total_lateral_inflow += net->lateral_inflow;
    if (net->outlet == NULL) {
      total_outflow += net->outflow;
    }
    total_storage += net->storage;
    total_storage_change += net->storage - net->last_storage;
    total_infiltration += net->infiltration;
    total_evaporation += net->evaporation;

    if (net->record) {
      if (fprintf(out, "%15s %10d %12.5g %12.5g %12.5g %12.5g %12.5g %12.5g",
        tstring, net->id, net->inflow, net->lateral_inflow,
        net->outflow, net->storage - net->last_storage, net->infiltration, net->evaporation) == EOF) {
          error_handler(ERRHDL_ERROR,
            "channel_save_outflow: write error:%s", strerror(errno));
          err++;
      }
      if (fprintf(out2, "%12.5g ", net->outflow) == EOF) {
        error_handler(ERRHDL_ERROR,
          "channel_save_outflow: write error:%s", strerror(errno));
        err++;
      }
      if (net->record_name != NULL) {
        if (fprintf(out, "   \"%s\"\n", net->record_name) == EOF) {
          error_handler(ERRHDL_ERROR,
            "channel_save_outflow: write error:%s",
            strerror(errno));
          err++;
        }

      }
      else {
        if (fprintf(out, "\n") == EOF) {
          error_handler(ERRHDL_ERROR,
            "channel_save_outflow: write error:%s",
            strerror(errno));
          err++;
        }
      }
    }
  }
  total_error = total_storage_change - total_lateral_inflow + total_outflow + total_infiltration + total_evaporation;
  if (fprintf(out, "%15s %10d %12.5g %12.5g %12.5g %12.5g %12.5g %12.5g %12.5g \"Totals\"\n",
    tstring, 0, total_lateral_inflow,
    total_outflow, total_storage,
    total_storage_change, total_infiltration, total_evaporation, total_error) == EOF) {
      error_handler(ERRHDL_ERROR,
        "channel_save_outflow: write error:%s", strerror(errno));
      err++;
  }
  fprintf(out2, "\n");

  return (err);
}

/* -------------------------------------------------------------
channel_free_network
------------------------------------------------------------- */
void channel_free_network(Channel * net)
{
  if (net->next != NULL) {
    channel_free_network(net->next);
  }
  free(net);
}

/* -------------------------------------------------------------
channel_init
------------------------------------------------------------- */
void channel_init(void)
{
  /* do nothing */
  return;
}

/* -------------------------------------------------------------
channel_done
------------------------------------------------------------- */
void channel_done(void)
{
  /* do nothing */
  return;
}

#if TEST_MAIN

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
  static int interval = 3600;	/* timestep in seconds */
  static int timestep;
  static int endtime = 144;
#define TIMES 6
  static float bndflow[TIMES] = { 0.0, 0.0, 300.0, 300.0, 0.0, 0.0 };
  static float bndtime[TIMES] = { 0.0, 12.0, 36.0, 48.0, 60.0, 1000.0 };

  float time;
  ChannelClass *class;
  Channel *simple = NULL, *current, *tail;

  error_handler_init(argv[0], NULL, ERRHDL_ERROR);
  channel_init();

  /* read classes */

  if ((class = channel_read_classes("example_classes.dat")) == NULL) {
    error_handler(ERRHDL_FATAL, "example_classes.dat: trouble reading file");
  }

  /* read a network */

  if ((simple = channel_read_network("example_network.dat", class)) == NULL) {
    error_handler(ERRHDL_FATAL, "example_network.dat: trouble reading file");
  }

  /* initialize flows */

  for (current = simple; current != NULL; current = current->next) {
    current->inflow = bndflow[0];
    current->outflow = bndflow[0];
    current->outlet = current->next;
    tail = current;
  }

  /* time loop */

  for (timestep = 0; timestep <= endtime; timestep++) {
    float inflow = interpolate(TIMES, bndtime, bndflow, timestep) * interval;
    float outflow;

    channel_step_initialize_network(simple);
    simple->inflow = inflow;
    (void) channel_route_network(simple, interval);
    outflow = tail->outflow / interval;
    channel_save_outflow(timestep * interval, simple, stdout);
  }

  channel_free_network(simple);
  channel_free_classes(class);
  channel_done();
  error_handler_done();
  exit(0);
}
#endif

/* -------------------------------------------------------------
channel_read_rveg_param
------------------------------------------------------------- */
int channel_read_rveg_param(Channel *head, const char *file, int *MaxID)
{
  Channel *current = NULL;
  int err = 0;
  int done;
  static const int fields = 18;
  static TableField rveg_fields[18] = {
    {"ID", TABLE_INTEGER, TRUE, FALSE, {0}, "", NULL},
    {"Height", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
    {"BufferWidth", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
    {"ExtnCoeff1", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
    {"ExtnCoeff2", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
    {"ExtnCoeff3", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
    {"ExtnCoeff4", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
    {"ExtnCoeff5", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
    {"ExtnCoeff6", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
    {"ExtnCoeff7", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
    {"ExtnCoeff8", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
    {"ExtnCoeff9", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
    {"ExtnCoeff10", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
    {"ExtnCoeff11", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
    {"ExtnCoeff12", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
    {"Dist", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
    {"Overhang", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
    {"StreamWidth", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
  };

  error_handler(ERRHDL_STATUS,
    "channel_read_rveg_param: reading file \"%s\"", file);

  if (table_open(file) != 0) {
    error_handler(ERRHDL_ERROR,
      "channel_read_rveg_param: unable to open file \"%s\": %s",
      file, strerror(errno));
    exit(3);
  }

  *MaxID = 0;
  done = FALSE;
  while (!done) {
    int i;
    done = (table_get_fields(fields, rveg_fields) < 0);
    if (done) {
      for (i = 0; i < fields; i++) {
        if (rveg_fields[i].read)
          break;
      }
      if (i >= fields)
        continue;
    }

    if (current == NULL) {
      current = head;
    }
    else {
      current = current->next;
    }


    for (i = 0; i < fields; i++) {
      if (rveg_fields[i].read) {
        switch (i) {
        case 0:
          current->id = rveg_fields[i].value.integer;
          if (current->id > *MaxID) 
            *MaxID = current->id;
          if (current->id <= 0) {
            error_handler(ERRHDL_ERROR, "%s: segment %d: channel id invalid",
              file, current->id);
            err++;
          }
          break;
        case 1:
          if (rveg_fields[i].value.real >= 0) {
            current->rveg.TREEHEIGHT = rveg_fields[i].value.real;
          }
          else {
            error_handler(ERRHDL_ERROR, "%s: segment %d: tree height (%f) invalid",
              file, current->id, rveg_fields[i].value.real);
            err++;
          }
          break;
        case 2:
          if (rveg_fields[i].value.real >= 0) {
            current->rveg.BUFFERWIDTH = rveg_fields[i].value.real;
          }
          else {
            error_handler(ERRHDL_ERROR, "%s: segment %d: buffer width (%f) invalid",
              file, current->id, rveg_fields[i].value.real);
            err++;
          }
          break;
        case 3:
          if (rveg_fields[i].value.real >= 0) {
            current->rveg.ExtnCoeff[0] = rveg_fields[i].value.real;
          }
          else {
            error_handler(ERRHDL_ERROR, "%s: segment %d: extinction coeff in month[1] (%f) invalid",
              file, current->id, rveg_fields[i].value.real);
            err++;
          }
          break;
        case 4:
          if (rveg_fields[i].value.real >= 0) {
            current->rveg.ExtnCoeff[1] = rveg_fields[i].value.real;
          }
          else {
            error_handler(ERRHDL_ERROR, "%s: segment %d: extinction coeff in month[2] (%f) invalid",
              file, current->id, rveg_fields[i].value.real);
            err++;
          }
          break;
        case 5:
          if (rveg_fields[i].value.real >= 0) {
            current->rveg.ExtnCoeff[2] = rveg_fields[i].value.real;
          }
          else {
            error_handler(ERRHDL_ERROR, "%s: segment %d: extinction coeff in month[3] (%f) invalid",
              file, current->id, rveg_fields[i].value.real);
            err++;
          }
          break;
        case 6:
          if (rveg_fields[i].value.real >= 0) {
            current->rveg.ExtnCoeff[3] = rveg_fields[i].value.real;
          }
          else {
            error_handler(ERRHDL_ERROR, "%s: segment %d: extinction coeff in month[4] (%f) invalid",
              file, current->id, rveg_fields[i].value.real);
            err++;
          }
          break;
        case 7:
          if (rveg_fields[i].value.real >= 0) {
            current->rveg.ExtnCoeff[4] = rveg_fields[i].value.real;
          }
          else {
            error_handler(ERRHDL_ERROR, "%s: segment %d: extinction coeff in month[5] (%f) invalid",
              file, current->id, rveg_fields[i].value.real);
            err++;
          }
          break;
        case 8:
          if (rveg_fields[i].value.real >= 0) {
            current->rveg.ExtnCoeff[5] = rveg_fields[i].value.real;
          }
          else {
            error_handler(ERRHDL_ERROR, "%s: segment %d: extinction coeff in month[6] (%f) invalid",
              file, current->id, rveg_fields[i].value.real);
            err++;
          }
          break;
        case 9:
          if (rveg_fields[i].value.real >= 0) {
            current->rveg.ExtnCoeff[6] = rveg_fields[i].value.real;
          }
          else {
            error_handler(ERRHDL_ERROR, "%s: segment %d: extinction coeff in month[7] (%f) invalid",
              file, current->id, rveg_fields[i].value.real);
            err++;
          }
          break;
        case 10:
          if (rveg_fields[i].value.real >= 0) {
            current->rveg.ExtnCoeff[7] = rveg_fields[i].value.real;
          }
          else {
            error_handler(ERRHDL_ERROR, "%s: segment %d: extinction coeff in month[8] (%f) invalid",
              file, current->id, rveg_fields[i].value.real);
            err++;
          }
          break;
        case 11:
          if (rveg_fields[i].value.real >= 0) {
            current->rveg.ExtnCoeff[8] = rveg_fields[i].value.real;
          }
          else {
            error_handler(ERRHDL_ERROR, "%s: segment %d: extinction coeff in month[9] (%f) invalid",
              file, current->id, rveg_fields[i].value.real);
            err++;
          }
          break;
        case 12:
          if (rveg_fields[i].value.real >= 0) {
            current->rveg.ExtnCoeff[9] = rveg_fields[i].value.real;
          }
          else {
            error_handler(ERRHDL_ERROR, "%s: segment %d: extinction coeff in month[10] (%f) invalid",
              file, current->id, rveg_fields[i].value.real);
            err++;
          }
          break;
        case 13:
          if (rveg_fields[i].value.real >= 0) {
            current->rveg.ExtnCoeff[10] = rveg_fields[i].value.real;
          }
          else {
            error_handler(ERRHDL_ERROR, "%s: segment %d: extinction coeff in month[11] (%f) invalid",
              file, current->id, rveg_fields[i].value.real);
            err++;
          }
          break;
        case 14:
          if (rveg_fields[i].value.real >= 0) {
            current->rveg.ExtnCoeff[11] = rveg_fields[i].value.real;
          }
          else {
            error_handler(ERRHDL_ERROR, "%s: segment %d: extinction coeff in month[12] (%f) invalid",
              file, current->id, rveg_fields[i].value.real);
            err++;
          }
          break;
        case 15:
          if (rveg_fields[i].value.real >= 0) {
            current->rveg.CanopyBankDist = rveg_fields[i].value.real;
          }
          else {
            error_handler(ERRHDL_ERROR, "%s: segment %d: distance to bank (%f) invalid",
              file, current->id, rveg_fields[i].value.real);
            err++;
          }
          break;
        case 16:
          if (rveg_fields[i].value.real >= 0) {
            current->rveg.OvhCoeff = rveg_fields[i].value.real;
          }
          else {
            error_handler(ERRHDL_ERROR, "%s: segment %d: overhanging coeff (%f) invalid",
              file, current->id, rveg_fields[i].value.real);
            err++;
          }
          break;
        case 17:
          if (rveg_fields[i].value.real >= 0) {
            current->rveg.StreamWidth = rveg_fields[i].value.real;
          }
          else {
            error_handler(ERRHDL_ERROR, "%s: segment %d: segment width (%f) invalid",
              file, current->id, rveg_fields[i].value.real);
            err++;
          }
          break;
        default:
          error_handler(ERRHDL_FATAL, "channel_read_rveg_param: what is this field %d?", i);
          break;
        }
      }
    }
  }

  table_close();
  table_errors += err;

  error_handler(ERRHDL_STATUS,
    "channel_read_rveg_param: %s: %d errors, %d warnings",
    file, table_errors, table_warnings);

  if (table_errors) {
    error_handler(ERRHDL_ERROR,
      "channel_read_rveg_param: %s: too many errors", file);
    channel_free_network(current);
    current = NULL;
  }

  return (err);
}
