
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "functions.h"
#include "Calendar.h"
#include "DHSVMerror.h"


/*****************************************************************************
  DayOfYear()
*****************************************************************************/
int DayOfYear(int Year, int Month, int Day)
{
  int i;
  int Jday;
  int DaysPerMonth[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

  if (IsLeapYear(Year))
    DaysPerMonth[1] = 29;
  else
    DaysPerMonth[1] = 28;

  for (i = 0, Jday = 0; i < (Month - 1); i++)
    Jday += DaysPerMonth[i];

  Jday += Day;

  return Jday;
}

/*****************************************************************************
  IsLeapYear()
*****************************************************************************/
uchar IsLeapYear(int Year)
{
  if ((Year % 4 == 0 && Year % 100 != 0) || Year % 400 == 0)
    return TRUE;
  return FALSE;
}

/*****************************************************************************
  IsEqualTime()
*****************************************************************************/
uchar IsEqualTime(DATE * Day1, DATE * Day2)
{
  return dequal(Day1->Julian, Day2->Julian);
}

/*****************************************************************************
  SScanIsoDate()

  Reads an instant written YYYY-MM-DD HH:MM (exactly 16 characters, a real
  calendar date and time of day). Returns FALSE for anything else.
*****************************************************************************/
int SScanIsoDate(const char *Str, DATE * Day)
{
  const char *Form = "dddd-dd-dd dd:dd";
  int DaysPerMonth[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
  int i;

  if (Str == NULL || strlen(Str) != 16)
    return FALSE;
  for (i = 0; i < 16; i++) {
    if (Form[i] == 'd') {
      if (!isdigit((unsigned char) Str[i]))
        return FALSE;
    }
    else if (Str[i] != Form[i])
      return FALSE;
  }

  Day->Year = (Str[0] - '0') * 1000 + (Str[1] - '0') * 100 + (Str[2] - '0') * 10 + (Str[3] - '0');
  Day->Month = (Str[5] - '0') * 10 + (Str[6] - '0');
  Day->Day = (Str[8] - '0') * 10 + (Str[9] - '0');
  Day->Hour = (Str[11] - '0') * 10 + (Str[12] - '0');
  Day->Min = (Str[14] - '0') * 10 + (Str[15] - '0');
  Day->Sec = 0;

  if (IsLeapYear(Day->Year))
    DaysPerMonth[1] = 29;
  if (Day->Year < 1 || Day->Month < 1 || Day->Month > MONTHPYEAR)
    return FALSE;
  if (Day->Day < 1 || Day->Day > DaysPerMonth[Day->Month - 1])
    return FALSE;
  if (Day->Hour > 23 || Day->Min > 59)
    return FALSE;

  Day->JDay = DayOfYear(Day->Year, Day->Month, Day->Day);
  Day->Julian =
    GregorianToJulianDay(Day->Year, Day->Month, Day->Day, Day->Hour, Day->Min,
			 Day->Sec);

  return TRUE;
}

/*****************************************************************************
  SPrintIsoDate(), SPrintIsoStamp()

  An instant as text (YYYY-MM-DD HH:MM) and in file names (YYYY-MM-DD_HHMM).
  Buffer holds at least 17 characters.
*****************************************************************************/
void SPrintIsoDate(DATE * Day, char *Buffer)
{
  sprintf(Buffer, "%04d-%02d-%02d %02d:%02d", Day->Year, Day->Month, Day->Day,
          Day->Hour, Day->Min);
}

void SPrintIsoStamp(DATE * Day, char *Buffer)
{
  sprintf(Buffer, "%04d-%02d-%02d_%02d%02d", Day->Year, Day->Month, Day->Day,
          Day->Hour, Day->Min);
}

/*****************************************************************************
  DaysFromCivil()

  Days from 1970-01-01 to a Gregorian date
*****************************************************************************/
static long long DaysFromCivil(int Year, int Month, int Day)
{
  long long y = Year - (Month <= 2);
  long long Era = (y >= 0 ? y : y - 399) / 400;
  long long YearOfEra = y - Era * 400;
  long long DayOfYr = (153 * (Month + (Month > 2 ? -3 : 9)) + 2) / 5 + Day - 1;
  long long DayOfEra = YearOfEra * 365 + YearOfEra / 4 - YearOfEra / 100 + DayOfYr;

  return Era * 146097 + DayOfEra - 719468;
}

/*****************************************************************************
  MinutesOf()

  Minutes from 1970-01-01 00:00 to an instant (instants are whole minutes)
*****************************************************************************/
long long MinutesOf(DATE * Day)
{
  return DaysFromCivil(Day->Year, Day->Month, Day->Day) * MINPDAY +
    Day->Hour * MINPHOUR + Day->Min;
}

/*****************************************************************************
  DayOfWaterYear()

  Day of the water year of a date, October 1 = 1
*****************************************************************************/
int DayOfWaterYear(DATE * Day)
{
  int StartYear = (Day->Month >= 10) ? Day->Year : Day->Year - 1;

  return (int) (DaysFromCivil(Day->Year, Day->Month, Day->Day) -
                DaysFromCivil(StartYear, 10, 1)) + 1;
}


/*****************************************************************************
  CopyDate()
*****************************************************************************/
void CopyDate(DATE * Copy, DATE * Original)
{
  Copy->Year = Original->Year;
  Copy->Month = Original->Month;
  Copy->Day = Original->Day;
  Copy->Hour = Original->Hour;
  Copy->Min = Original->Min;
  Copy->Sec = Original->Sec;
  Copy->JDay = Original->JDay;
  Copy->Julian = Original->Julian;
}

/*****************************************************************************
  IsNewMonth()
*****************************************************************************/
uchar IsNewMonth(DATE *Now, int Interval)
{
  int Year;
  int Month;
  int Day;
  int Hour;
  int Min;
  double Sec;
  double Julian;

  Julian = Now->Julian - ((double) Interval) / SECPDAY;
  JulianDayToGregorian(Julian, &Year, &Month, &Day, &Hour, &Min, &Sec);
  if (Month != Now->Month)
    return TRUE;
  else
    return FALSE;
}

/*****************************************************************************
  IsNewWaterYear()
*****************************************************************************/
uchar IsNewWaterYear(DATE *Now)
{
  int Year;
  int Month;
  int Day;
  int Hour;
  int Min;
  double Sec;
  double Julian;

  Julian = Now->Julian;
  JulianDayToGregorian(Julian, &Year, &Month, &Day, &Hour, &Min, &Sec);
  if (Month == 10 && Day == 1 && Hour == 0)
    return TRUE;
  else
    return FALSE;
}

/*****************************************************************************
  IsNewDay()
*****************************************************************************/
uchar IsNewDay(int DayStep)
{
  if (DayStep == 0)
    return TRUE;
  else
    return FALSE;
}

/*****************************************************************************
  Before()
*****************************************************************************/
uchar Before(DATE * Day1, DATE * Day2)
{
  if (Day1->Julian < Day2->Julian)
    return TRUE;
  else
    return FALSE;
}

/*****************************************************************************
  After()
*****************************************************************************/
uchar After(DATE *Day1, DATE *Day2)
{
  if (Day1->Julian > Day2->Julian)
    return TRUE;
  else
    return FALSE;
}

/*****************************************************************************
  IncreaseTime()
*****************************************************************************/
void IncreaseTime(TIMESTRUCT * Time)
{
  double Sec;

  (Time->Step)++;
  Time->DayStep = (Time->DayStep + 1) % Time->NDaySteps;
  Time->Current.Julian = Time->Start.Julian +
    ((double) Time->Step) * (((double) Time->Dt) / SECPDAY);
  JulianDayToGregorian(Time->Current.Julian, &(Time->Current.Year),
		       &(Time->Current.Month), &(Time->Current.Day),
		       &(Time->Current.Hour), &(Time->Current.Min), &Sec);
  Time->Current.Sec = (int) Sec;
  Time->Current.JDay = DayOfYear(Time->Current.Year, Time->Current.Month,
				 Time->Current.Day);
}

/*****************************************************************************
 IncreaseVariableTime()
 *****************************************************************************/
void IncreaseVariableTime(TIMESTRUCT *Time, float VariableDT, TIMESTRUCT *NextTime)
{
  double Sec;
  
  Time->Current.Julian = Time->Current.Julian + VariableDT/SECPDAY;
  
  if(After(&(Time->Current), &(NextTime->Current)) || 
     IsEqualTime(&(Time->Current), &(NextTime->Current))) {
    
    (Time->Step)++;
    Time->DayStep = (Time->DayStep + 1) % Time->NDaySteps;
  }
  
  JulianDayToGregorian(Time->Current.Julian, &(Time->Current.Year),
                       &(Time->Current.Month), &(Time->Current.Day),
                       &(Time->Current.Hour), &(Time->Current.Min), 
                       &(Sec));
                       Time->Current.Sec = (int) Sec;
                       Time->Current.JDay = DayOfYear(Time->Current.Year, Time->Current.Month,
                                                      Time->Current.Day);
                       
}

/******************************************************************************/
/*				   InitTime()                                 */
/******************************************************************************/
/* Start and End are instants on the step grid with End after Start by a
   whole number of steps (checked by the caller); Dt divides a day */
void InitTime(TIMESTRUCT * Time, DATE * Start, DATE * End, int Dt)
{
  Time->Dt = Dt;		/* timestep in seconds */

  CopyDate(&(Time->Start), Start);
  CopyDate(&(Time->End), End);
  CopyDate(&(Time->Current), &(Time->Start));

  Time->Step = 0;
  Time->NDaySteps = SECPDAY / Time->Dt;
  Time->DayStep = ((Time->Start.Hour * MINPHOUR + Time->Start.Min) * SECPMIN) / Time->Dt;
  Time->NSteps = (int) ((MinutesOf(&(Time->End)) - MinutesOf(&(Time->Start))) *
                        SECPMIN / Time->Dt);
}

/******************************************************************************/
/*			     GregorianToJulianDay()                           */
/******************************************************************************/
/*
** Takes a date, and returns a Julian day. A Julian day is the number of
** days since some base date  (in the very distant past).
** Handy for getting date of x number of days after a given Julian date
** (use jdate to get that from the Gregorian date).
** Author: Robert G. Tantzen, translator: Nat Howard
** Translated from the algol original in Collected Algorithms of CACM
** (This and jdate are algorithm 199).
** Copied from the xmgr auxiliary functions Mon Feb  1 14:16:40 1999
*/
double GregorianToJulianDay(int year, int mon, int day, int h, int mi,
			    double se)
{
  long m = mon;
  long d = day;
  long y = year;
  long c;
  long ya;
  long j;
  double seconds = h * 3600.0 + mi * 60 + se;

  if (m > 2)
    m -= 3;
  else {
    m += 9;
    --y;
  }
  c = y / 100L;
  ya = y - (100L * c);
  j = (146097L * c) / 4L + (1461L * ya) / 4L +
    (153L * m + 2L) / 5L + d + 1721119L;
  if (seconds < 12 * 3600.0) {
    j--;
    seconds += 12.0 * 3600.0;
  }
  else {
    seconds = seconds - 12.0 * 3600.0;
  }
  return (j + (seconds / 3600.0) / 24.0);
}

/******************************************************************************/
/*			     JulianDayToGregorian()                           */
/******************************************************************************/
/* Julian date converter. Takes a julian date (the number of days since
** some distant epoch or other), and returns an int pointer to static space.
** ip[0] = month;
** ip[1] = day of month;
** ip[2] = year (actual year, like 1977, not 77 unless it was  77 a.d.);
** ip[3] = day of week (0->Sunday to 6->Saturday)
** These are Gregorian.
** Copied from Algorithm 199 in Collected algorithms of the CACM
** Author: Robert G. Tantzen, Translator: Nat Howard
** Copied from the xmgr auxiliary functions Mon Feb  1 14:16:40 1999
*/
void JulianDayToGregorian(double jd, int *y, int *m, int *d, int *h, int *mi,
			  double *sec)
{
  long j = jd;
  double tmp, frac = jd - j;

  /* The following four lines are added so that we round to the whole seconds, 
     Bart Nijssen Wed Feb  3 08:31:50 1999 */
  if (rint(frac * SECPDAY) != floor(frac * SECPDAY)) {
    jd += 1. / SECPDAY;
    j = jd;
    frac = jd - j;
  }
  /* end of added code */

  if (frac >= 0.5) {
    frac = frac - 0.5;
    j++;
  }
  else {
    frac = frac + 0.5;
  }

  j -= 1721119L;
  *y = (4L * j - 1L) / 146097L;
  j = 4L * j - 1L - 146097L * *y;
  *d = j / 4L;
  j = (4L * *d + 3L) / 1461L;
  *d = 4L * *d + 3L - 1461L * j;
  *d = (*d + 4L) / 4L;
  *m = (5L * *d - 3L) / 153L;
  *d = 5L * *d - 3 - 153L * *m;
  *d = (*d + 5L) / 5L;
  *y = 100L * *y + j;
  if (*m < 10)
    *m += 3;
  else {
    *m -= 9;
    *y += 1;
  }
  tmp = 3600.0 * (frac * 24.0);
  *h = (int) (tmp / 3600.0);
  tmp = tmp - *h * 3600.0;
  *mi = (int) (tmp / 60.0);
  *sec = tmp - *mi * 60.0;
}
