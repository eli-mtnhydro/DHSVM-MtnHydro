
#ifndef CONSTANTS_H
#define CONSTANTS_H

#define CELL_PARTITION 2        /* Number of veg type in a grid cell */
#define CH_ICE     (2100.0e3)	/* Volumetric heat capacity (J/(m3*C) of ice (0C) */
#define CH_WATER   (4186.8e3)	/* Volumetric heat capacity (J/(m3*C) of water */
#define CP         1013.0		/* Specific heat of moist air at constant pressure (J/(kg*C)) */
#define DELTAT     50.		    /* Used in SensibleHeatFlux to bracket the effective surface temperature (C) */
#define DEGPRAD     57.29578	/* degree per radian */
#define D0_MULTIPLIER  0.63		/* Multiplier for vegetation height to get displacement height (m) */
#define DZ_TOP      0.1			/* Thickness of soil surface layer for which heat stoarge change is calculated (m) */
#define EPS         0.622		/* ratio of molecular weight of water vapor to that for dry air */
#define G           9.81		/* gravitational accelleration (m/(s^2)) */
#define GRAMSPKG    1000.	    /* grams per kilogram */
#define JOULESPCAL  4.1868	    /* Joules per calorie */
#define KhH2O       0.58		/* Thermal conductivity of water (W/(mk)) */
#define LF            (333.7e3) /* latent heat of fusion (J/kg) */
#define MINPDEG        4.		/* minutes per degree longitude */
#undef PI
#define PI   3.14159265358979323846
#define RADPHOUR    0.2617994   /* radians per hour: Earth's Rotation (2 PI rad/day) * (1 day/24 h) */
#define RADPDEG     (PI/180.0)	/* radians per degree */
#define SOLARCON    1360.	    /* Solar constant (W/m^2) */
#define STEFAN    (5.6696e-8)	/* Stefan-Boltzmann constant (W/(M^2*C^4) */
#define VISFRACT    0.5		    /* part of shortwave that is in the visible range of the spectrum */
#define VON_KARMAN  0.4		    /* Von Karman's constant */
#define WATER_DENSITY 1000.		/* Density of water in kg/m3 */
#define Z0_MULTIPLIER 0.13		/* Multiplier for vegetation height to get roughness length (m) */
#define MINSTORAGEK (1e-10)   /* Minimum allowed value of the linear storage parameter for numerical stability */
#define MTHRESH        0.9   /* Saturation extent is based on the number of pixels with a water table 
                                that is at least MTHRESH of soil depth. */

/**************** extern constants - see globals.c ****************/

extern int NDIRS;               /* How many neighbors are used in surface/subsurface routing */
extern int xdirection4[];
extern int ydirection4[];
extern int xdirection8[];
extern int ydirection8[];
extern int *xdirection, *ydirection;

extern float LAI_SNOW_MULTIPLIER;		/* Multiplier to calculate the amount of available snow interception as a function of LAI */
extern float LAI_WATER_MULTIPLIER;		/* Multiplier to determine maximum interception storage as a function of LAI  */
extern float LIQUID_WATER_CAPACITY;		/* Water holding capacity of snow as a fraction of snow-water-equivalent */
extern float MAX_SNOW_TEMP;				/* Maximum temperature at which snow can occur (C) */
extern float MIN_INTERCEPTION_STORAGE;	/* The amount of snow on the canopy that can only be melted off (m) */
extern float MIN_SNOW_RESET_ALBEDO; /* Minimum snow (m) required to completely reset snow albedo */
extern float LAMBDA_FOREST_OFFSET; /* Additional correction applied to albedo decay rate in forest */
extern float MIN_RAIN_TEMP;				/* Minimum temperature at which rain can occur (C) */
extern unsigned char OUTSIDEBASIN;		/* Mask value indicating outside the basin */
extern float MINELEV; /* Smallest elevation of all grid cells (m) */
extern float SNOWPAT_WEIGHT; /* Fractional weight between snow pattern and precip pattern for snowfall */
extern float TEMPLAPSE;					/* Temperature lapse rate in C/m */
extern float TEMPERATURE_OFFSET; /* Uniform offset added to input air temperature (C) */
extern float LAPSE_RATE_BIAS; /* Additional lapse rate relative to LAPSE_BIAS_ELEV (C/m) */
extern float LAPSE_BIAS_ELEV; /* Elevation where temperature bias is zero (m) */
extern float SOIL_DEPTH_ADJ; /* Calibration param - additive */
extern float SOIL_KSAT_ADJ; /* Calibration param - multiplicative */
extern float SOIL_EXPDEC_ADJ; /* Calibration param - multiplicative */
extern float SOIL_POROSITY_ADJ; /* Calibration param - multiplicative */
extern float SOIL_FIELDCAP_ADJ; /* Calibration param - multiplicative */
extern float VEG_LAI_ADJ; /* Calibration param - multiplicative */
extern float VEG_FC_ADJ; /* Calibration param - multiplicative */
extern float Z0_GROUND;					/* Roughness length for bare soil (m) */
extern float Z0_SNOW;					/* Roughness length for snow (m) */
extern float Zref;						/* Reference height (m) */

/* Snow albedo decay curve */
extern float ALB_MAX;                   /* Fresh snow albedo */                                                               
extern float ALB_ACC_LAMBDA;            /* Snow freeze albedo curve control parameters */
extern float ALB_MELT_LAMBDA;           /* Snow thaw albedo curve control parameters */
extern float ALB_ACC_MIN;
extern float ALB_MELT_MIN;
extern float PRECIP_MULTIPLIER;        /* Precipitation multiplier */
extern float MAX_SURFACE_SWE; 	       /* Maximum depth of the surface layer in snow water equivalent (m) */

extern float GAPWIND_FACTOR;
extern int TotNumGap;                  /* Total number of grid cells with a gap structure */

extern float SNOWSLIDE1;               /* First Parameter in Snowslide equation */
extern float SNOWSLIDE2;               /* Second Parameter in Snowslide equation */
#endif
