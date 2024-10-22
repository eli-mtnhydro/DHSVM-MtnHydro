/*
 * SUMMARY:      CalcSnowAlbedo.c - Calculate snow albedo
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Calculate the snow albedo as a function of snow age
 * DESCRIP-END.
 * FUNCTIONS:    CalcSnowAlbedo()
 * COMMENTS:
 * $Id: CalcSnowAlbedo.c,v 1.4 2003/07/01 21:26:10 olivier Exp $     
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "constants.h"
#include "data.h"
#include "Calendar.h"
#include "functions.h"

/*****************************************************************************
  CalcSnowAlbedo()

  Source:
  Laramie, R. L., and J. C. Schaake, Jr., Simulation of the continuous
  snowmelt process, Ralph M. Parsons Laboratory, Mass. Inst. of Technol.,
  1972

  Snow albedo is calculated as a function of the number of days since the
  last observed snow fall. There are separete albedo curves for the freeze
  and thaw conditions.
*****************************************************************************/
float CalcSnowAlbedo(SNOWPIX *LocalSnow, int StepsPerDay)
{
  /* Laramie and Schaake (1972) */
  /* Updated based on Storck (2000) */
  /* Updated by Eli Boardman 2024 to fix various issues */
  float Last, Albedo;
  
  Last = LocalSnow->LastSnow / (float) StepsPerDay;
  
  if (Last > (float) DAYPYEAR)
    Last = (float) DAYPYEAR;
  
  if (LocalSnow->AccumSeason == TRUE) {
    /* Accumulation season */
    Albedo = LocalSnow->amax * pow(LocalSnow->LamdaAcc, pow(Last, 0.58));
    if (Albedo < LocalSnow->AccMin)
      Albedo = LocalSnow->AccMin;
  } else {
    /* Melt season */
    Albedo = LocalSnow->amax * pow(LocalSnow->LamdaMelt, pow(Last, 0.46));
    if (Albedo < LocalSnow->MeltMin)
      Albedo = LocalSnow->MeltMin;
  }
  
  return Albedo;
}
