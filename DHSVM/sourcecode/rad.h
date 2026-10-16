
#ifndef RAD_H
#define RAD_H

#define ALBEDO  0.15		/* WORK IN PROGRESS, See InitNewStep() */

void SeparateRadiation(float TotalSolar, float ClearIndex,
		       float *Beam, float *Diffuse);

void SolarDay(int DayOfYear, float Longitude, float Latitude,
	      float StandardMeridian, float *NoonHour, float *Declination,
	      float *HalfDayLength, float *Sunrise, float *Sunset,
	      float *TimeAdjustment, float *SunEarthDist);

void SolarHour(float Latitude, float LocalHour, float Dt, float NoonHour,
	       float Declination, float Sunrise, float Sunset,
	       float TimeAdjustment, float SunEarthDist,
	       float *SineSolarAltitude, int *DayLight, float *SolarTimeStep,
	       float *SunMax, float *SolarAzimuth);

#endif
