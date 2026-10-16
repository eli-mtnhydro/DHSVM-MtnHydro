/*
 * DESCRIPTION:  Where the surface runoff of the impervious fraction of a
 *               cell goes: the cell its D8 flow path reaches first among the
 *               channel cells (cells with a row in the channel cell table).
 *               The D8 direction of every mask cell is the neighbour with
 *               the steepest drop on the whole elevation grid (neighbours
 *               outside the mask count, neighbours outside the grid do
 *               not), visited in the order NE, E, SE, S, SW, W, NW, N; a
 *               neighbour replaces the choice only with a strictly larger
 *               drop, so a tie keeps the first, and a cell with no lower
 *               neighbour points east. From each mask cell the path follows
 *               these directions; a move that would leave the mask turns
 *               clockwise to the next direction until it stays in the mask.
 *               The path ends on a channel cell, or after 100 moves
 *               wherever it is; a cell with no neighbour in the mask is the
 *               end of its own path.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"
#include "channel_grid.h"

#define NROUTEDIRS 8			/* Neighbour directions */
#define MAX_ROUTE_MOVES 100		/* Longest flow path */

/* Directions clockwise from east: E, SE, S, SW, W, NW, N, NE */
static const int RouteRow[NROUTEDIRS] = { 0, 1, 1, 1, 0, -1, -1, -1 };
static const int RouteCol[NROUTEDIRS] = { 1, 1, 0, -1, -1, -1, 0, 1 };

/* The order in which the D8 direction visits them: NE, E, SE, S, SW, W, NW, N */
static const int D8Order[NROUTEDIRS] = { 7, 0, 1, 2, 3, 4, 5, 6 };

static int InMask(MAPSIZE *Map, TOPOPIX **TopoMap, int y, int x)
{
  return (y >= 0 && y < Map->NY && x >= 0 && x < Map->NX && INBASIN(TopoMap[y][x].Mask));
}

/*****************************************************************************
  InitSurfaceRouting()

  Sets drains_y, drains_x of every mask cell: the end of its flow path
*****************************************************************************/
void InitSurfaceRouting(MAPSIZE *Map, TOPOPIX **TopoMap, CHANNEL *ChannelData)
{
  const char *Routine = "InitSurfaceRouting";
  unsigned char *Dir;
  double Diagonal = (double) Map->DX * sqrt(2.0);
  double Distance, Drop, Best;
  int x, y, yn, xn, n, k, D, Move, Turn;

  Message("Deriving the impervious surface routing\n");

  if (!(Dir = (unsigned char *) calloc((size_t) Map->NY * Map->NX, sizeof(unsigned char))))
    ReportError((char *) Routine, 1);

  /* D8 direction of every mask cell */
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (!INBASIN(TopoMap[y][x].Mask))
        continue;
      Best = 0.0;
      D = 0;
      for (n = 0; n < NROUTEDIRS; n++) {
        k = D8Order[n];
        yn = y + RouteRow[k];
        xn = x + RouteCol[k];
        if (yn < 0 || yn >= Map->NY || xn < 0 || xn >= Map->NX)
          continue;
        Distance = (RouteRow[k] != 0 && RouteCol[k] != 0) ? Diagonal : (double) Map->DX;
        Drop = ((double) TopoMap[y][x].Dem - (double) TopoMap[yn][xn].Dem) / Distance;
        if (Drop > Best) {
          Best = Drop;
          D = k;
        }
      }
      Dir[y * Map->NX + x] = (unsigned char) D;
    }
  }

  /* Flow path of every mask cell to a channel cell */
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (!INBASIN(TopoMap[y][x].Mask))
        continue;
      yn = y;
      xn = x;
      for (Move = 0; Move < MAX_ROUTE_MOVES &&
           !channel_grid_has_channel(ChannelData->stream_map, xn, yn); Move++) {
        D = Dir[yn * Map->NX + xn];
        for (Turn = 0; Turn < NROUTEDIRS; Turn++) {
          k = (D + Turn) % NROUTEDIRS;
          if (InMask(Map, TopoMap, yn + RouteRow[k], xn + RouteCol[k]))
            break;
        }
        /* No neighbour in the mask: only possible for the first cell of a
           path (every later cell has the one it came from) */
        if (Turn == NROUTEDIRS)
          break;
        yn += RouteRow[k];
        xn += RouteCol[k];
      }
      TopoMap[y][x].drains_y = yn;
      TopoMap[y][x].drains_x = xn;
    }
  }
  free(Dir);
}
