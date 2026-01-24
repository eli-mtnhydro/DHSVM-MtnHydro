
int NDIRS;  /* How many neighbors are used in surface/subsurface routing */
/* These indices are so neighbors can be looked up quickly */
int xdirection4[] = {0,  1,  0, -1};
int ydirection4[] = {-1,  0,  1,  0};
int xdirection8[] = {-1,  0,  1,  1,  1,  0, -1, -1};
int ydirection8[] = {1,  1,  1,  0, -1, -1, -1,  0};
int *xdirection;
int *ydirection;

float LAI_SNOW_MULTIPLIER;		/* Multiplier to calculate the amount of available snow interception as a function of LAI */
float LAI_WATER_MULTIPLIER;		/* Multiplier to determine maximum interception storage as a function of LAI  */
float LIQUID_WATER_CAPACITY;		/* Water holding capacity of snow as a fraction of snow-water-equivalent */
float MAX_SNOW_TEMP;				/* Maximum temperature at which snow can occur (C) */
float MIN_INTERCEPTION_STORAGE;	/* The amount of snow on the canopy that can only be melted off (m) */
float MIN_SNOW_RESET_ALBEDO; /* Minimum snow (m) required to completely reset snow albedo */
float LAMBDA_FOREST_OFFSET; /* Additional correction applied to albedo decay rate in forest */
float MIN_RAIN_TEMP;				/* Minimum temperature at which rain can occur (C) */
unsigned char OUTSIDEBASIN;		/* Mask value indicating outside the basin */
float MINELEV; /* Smallest elevation of all grid cells (m) */
float SNOWPAT_WEIGHT; /* Fractional weight between snow pattern and precip pattern for snowfall */
float TEMPLAPSE;					/* Temperature lapse rate in C/m */
float TEMPERATURE_OFFSET; /* Uniform offset added to input air temperature (C) */
float LAPSE_RATE_BIAS; /* Additional lapse rate relative to LAPSE_BIAS_ELEV (C/m) */
float LAPSE_BIAS_ELEV; /* Elevation where temperature bias is zero (m) */
float SOIL_DEPTH_ADJ; /* Calibration param - additive */
float SOIL_KSAT_ADJ; /* Calibration param - multiplicative */
float SOIL_EXPDEC_ADJ; /* Calibration param - multiplicative */
float SOIL_POROSITY_ADJ; /* Calibration param - multiplicative */
float SOIL_FIELDCAP_ADJ; /* Calibration param - multiplicative */
float VEG_LAI_ADJ; /* Calibration param - multiplicative */
float VEG_FC_ADJ; /* Calibration param - multiplicative */
float Z0_GROUND;					/* Roughness length for bare soil (m) */
float Z0_SNOW;					/* Roughness length for snow (m) */
float Zref;						/* Reference height (m) */

/* Snow albedo decay curve */
float ALB_MAX;                   /* Fresh snow albedo */                                                               
float ALB_ACC_LAMBDA;            /* Snow freeze albedo curve control parameters */
float ALB_MELT_LAMBDA;           /* Snow thaw albedo curve control parameters */
float ALB_ACC_MIN;
float ALB_MELT_MIN;
float PRECIP_MULTIPLIER;        /* Precipitation multiplier */
float MAX_SURFACE_SWE; 	       /* Maximum depth of the surface layer in snow water equivalent (m) */

float GAPWIND_FACTOR;
int TotNumGap;                  /* Total number of grid cells with a gap structure */

float SNOWSLIDE1;               /* First Parameter in Snowslide equation */
float SNOWSLIDE2;               /* Second Parameter in Snowslide equation */
