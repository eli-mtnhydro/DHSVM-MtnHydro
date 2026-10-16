
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
#include "settings.h"

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
channel_create_classes
Builds the linked list of ChannelClass structs from the rows of the
channel class table (in table order)
------------------------------------------------------------- */
ChannelClass *channel_create_classes(int N, int *ID, float *Width, float *BankHeight,
                                     float *Friction)
{
  ChannelClass *head = NULL, *current = NULL;
  int i;

  for (i = 0; i < N; i++) {
    if (head == NULL) {
      head = alloc_channel_class();
      current = head;
    }
    else {
      current->next = alloc_channel_class();
      current = current->next;
    }
    current->id = (ClassID) ID[i];
    current->width = Width[i];
    current->bank_height = BankHeight[i];
    current->friction = Friction[i];
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
  seg->ground_slope = 0.0;
  seg->top_water_depth = 0.0;
  seg->bottom_water_depth = 0.0;
  seg->K = 0.0;
  seg->class2 = NULL;
  seg->IntersectsLake = FALSE;
  seg->lake = NULL;
  seg->lateral_inflow = 0.0;
  seg->last_inflow = 0.0;
  seg->last_outflow = 0.0;
  seg->last_storage = 0.0;
  seg->inflow = 0.0;
  seg->lake_inflow = 0.0;
  seg->outflow = 0.0;
  seg->storage = 0.0;
  seg->infiltration = 0.0;
  seg->remaining_infil = 0.0;
  seg->evaporation = 0.0;
  seg->remaining_evap = 0.0;
  seg->last_lateral_inflow = 0.0;
  seg->outlet = NULL;
  seg->next = NULL;
  seg->grid = NULL;

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
      }
      segment = segment->next;
    }
  }
}

/* -------------------------------------------------------------
channel_create_network
Builds the linked list of channel segments from the rows of the channel
segment table (in table order); a segment with a save name is recorded
in the outputs under that name
------------------------------------------------------------- */
Channel *channel_create_network(int N, int *ID, int *Order, float *Slope, float *Length,
                                int *Class, int *Outlet, char **SaveName,
                                ChannelClass *class_list, int *MaxID)
{
  Channel *head = NULL, *current = NULL;
  int i;

  *MaxID = 0;
  for (i = 0; i < N; i++) {
    if (head == NULL) {
      head = alloc_channel_segment();
      current = head;
    }
    else {
      current->next = alloc_channel_segment();
      current = current->next;
    }

    current->IntersectsLake = FALSE;
    current->id = (SegmentID) ID[i];
    if (current->id > *MaxID)
      *MaxID = current->id;
    current->order = Order[i];
    current->ground_slope = Slope[i];
    current->length = Length[i];
    current->class2 = find_channel_class(class_list, (ClassID) Class[i]);
    current->outid = (SegmentID) Outlet[i];
    if (SaveName[i][0] != '\0') {
      current->record = TRUE;
      current->record_name = (char *) strdup(SaveName[i]);
    }
  }

  /* find segment outlet segments, if specified */
  for (current = head; current != NULL; current = current->next) {
    if (current->outid != 0)
      current->outlet = channel_find_segment(head, current->outid);
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
  double k, s0;
  int err = 0;
  
  segment->inflow += segment->lake_inflow;
  
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
  
  /* Linear reservoir over the step, S0 exp(-K dt) + net (1 - exp(-K dt)) / K,
     in a form without cancellation for small K; K = 0 keeps the net inflow */
  k = segment->K;
  s0 = segment->storage;
  if (k > 0.0)
    storage = (float) (s0 + ((double) net_balance - k * s0) * (-expm1(-k * deltat) / k));
  else
    storage = (float) (s0 + (double) net_balance * deltat);
  
  if (storage < 0.0)
    storage = 0.0;
  
  outflow = net_balance - (storage - segment->storage) / deltat;
  
  if (outflow < 0.0)
    outflow = 0.0;
  
  segment->outflow = outflow * deltat;
  segment->storage = storage;
  
  if (segment->IntersectsLake) {
    segment->outflow += segment->storage;
    segment->storage = 0.0;
    segment->lake_inflow -= segment->outflow;
    segment->lake->Inflow += segment->outflow / segment->lake->Area;
  } else {
    if (segment->outlet != NULL)
      segment->outlet->inflow += segment->outflow;
  }
  
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
    net->lake_inflow = 0.0;
    net->lateral_inflow = 0.0;
    net->last_outflow = net->outflow;
    net->last_storage = net->storage;
    net->infiltration = 0.0;
    net->remaining_infil = 0.0;
    net->evaporation = 0.0;
    
    channel_step_initialize_network(net->next);
  }
  return (0);
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

