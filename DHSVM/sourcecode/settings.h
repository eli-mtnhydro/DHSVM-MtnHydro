
#ifndef SETTINGS_H
#define SETTINGS_H

#ifndef _AIX			/* AIX 3.2.5 already defines this */
typedef unsigned char uchar;
#endif
typedef unsigned short unshort;
typedef unsigned int unint;

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define INBASIN(x) ((x) != OUTSIDEBASIN)
#ifndef ABSVAL
#define ABSVAL(x)  ( (x) < 0 ? -(x) : (x) )
#endif

#ifndef TRUE
#define TRUE           1
#endif
#ifndef FALSE
#define FALSE          0
#endif

#define DHSVM_HUGE     (1e20)

/* When compiling DHSVM using Microsoft C++ define MSC++ at compile time.
   Microsoft C++ treats all numeric constants as doubles, and issues a
   warning if this constant is assigned to a float without an explicit type
   cast.   This warning can becme somewhat annoying, and the following pragma
   directive disables the warning */
#ifdef MS_C_PLUSPLUS
#pragma warning(disable : 4244)
#endif

/* Default value for not applicable */
#define NOT_APPLICABLE -9999

/* Options for flow gradient calculation */
#define TOPOGRAPHY     1
#define WATERTABLE     2

/* Options for meterological interpolation */
#define INVDIST        1
#define NEAREST        2
#define VARCRESS       3
#define UNIFORM        4

/* Options for model extent */
#define POINT 1
#define BASIN 2

/* Options for lapse rate and spatial parameters */
#define CONSTANT 1
#define VARIABLE 2
#define MAP 3

/* Options for infiltration */
#define STATIC 1
#define DYNAMIC 2

/* Options for canopy radiation attenuation */
#define FIXED    1
#define VARIABLE 2

#define TINY       1e-20
#define DEBUG      FALSE

#define HEADERLINES    5
#define BUFSIZE      255
#define MAXUCHAR     255	/* Maximum value of a 1-byte uchar */
#define MAXSTRING    255
#define NAMESIZE     127

#define MAXDIRS        8	/* Maximum number of directions in which water can flow, must equal 8 */
#define NNEIGHBORS     8  /* Number of directions in which water can flow based on fine grid, must equal 8 */


#define NA          -9999	/* Not applicable */

#define MAP_OUTPUT 1

#define MIN_SWE 0.005
#define MELTOUT_SWE 0.05

// Canopy type used in canopy gapping option
enum CanopyType {
  Opening,
  Forest
};

enum StomatalSchemeType {
  JARVIS = 0,      /* Wigmosta et al. 1994, multiplicative                 */
  MEDLYN,          /* Medlyn et al. 2011                                   */
  PROFITMAX,       /* Sperry et al. 2017                                   */
  PROFITMAX2,      /* Wang et al. 2020                                     */
  SOX              /* Eller et al. 2018 / 2020                             */
};

/* Plant hydraulics.  Orthogonal to the scheme above; valid combinations are
   enforced in InitConstants.c.  PROFITMAX2 must be parsed BEFORE PROFITMAX:
   strncmp on "PROFITMAX" with length 9 matches both. */
enum PlantHydraulicsType {
  HYD_NONE = 0,    /* empirical beta from soil moisture                    */
  HYD_KRSSUF       /* Couvreur 2012 / Vanderborght 2021 macroscopic roots  */
};

enum KEYS {
/* Options *//* list order must match order in InitConstants.c */
  extent = 0, gradient, routing_neighbors, routing_mfd,
  sensible_heat_flux, routing, lakedyna, interflow, vertksatsource, infiltration,
  interpolation, max_interp_dist, prism, snowpattern,
  canopy_radatt, shading, outside, rhoverride, 
  temp_lapse, cressman_radius, cressman_stations,
  prism_data_path, prism_data_ext, snowpattern_data_path,
  shading_data_path, shading_data_ext, skyview_data_path, 
  improv_radiation, stomatal_scheme, plant_hydraulics, gapping, snowslide, sepr, 
  snowstats, dynaveg, streamdata, streamtime, gw_spinup, gw_spinup_yrs, gw_spinup_recharge,
  /* Area */
  coordinate_system, extreme_north, extreme_west, center_latitude,
  center_longitude, time_zone_meridian, number_of_rows,
  number_of_columns, grid_spacing, point_north, point_east,
  /* Time */
  time_step, model_start, model_end, 
  /* Constants */
  ground_roughness, snow_roughness, 
  snow_water_capacity, reference_height, rain_lai_multiplier,
  snow_lai_multiplier, min_intercepted_snow, min_snow_reset_albedo, albedo_lambda_forest,
  outside_basin,
  temp_lapse_rate, max_swe, snowslide_parameter1, snowslide_parameter2,
  gapwind_adj, snowpattern_weight, temperature_offset, lapse_bias, lapse_elev,
  soil_depth_adj, soil_ksat_adj, soil_expdec_adj, soil_porosity_adj, soil_fieldcap_adj, veg_lai_adj, veg_fc_adj,
  atmos_co2, atmos_co2_file,
  /* Constants that can vary spatially */
  rain_threshold = 0,
  snow_threshold,
  fresh_alb,
  alb_acc_lambda,
  alb_melt_lambda,
  alb_acc_min,
  alb_melt_min,
  multiplier,
  multiplier_melt,
  /* Station information */
  station_name = 0, station_north, station_east, station_elev, station_file,
  /* Soil information */
  soil_description = 0, lateral_ks, exponent, depth_thresh, anisotropy, max_infiltration, capillary_drive, deepflux,
  soil_albedo, manning, number_of_layers, porosity, pore_size, bubbling_pressure, field_capacity,
  wilting_point, bulk_density, vertical_ks, solids_thermal, thermal_capacity,
  rhizosphere_ks, soil_last_key = rhizosphere_ks,
  /* Vegetation information */
  veg_description = 0, overstory, understory, fraction, hemifraction, trunk_space,
  aerodynamic_att, beam_attn, diff_attn, clumping_factor, leaf_angle_a, leaf_angle_b,
  scat, snow_int_cap, mass_drip_ratio, snow_int_eff, imperv_frac, detention_frac, 
  detention_decay, height, max_resistance, min_resistance,
  max_carboxylation, stomatal_slope, stomatal_intercept, xylem_p50,
  moisture_threshold, vpd, rpc,
  number_of_root_zones, root_zone_depth, overstory_fraction, understory_fraction, 
  monextn, vf_adj, overstory_monlai, understory_monlai, overstory_monalb, understory_monalb, 
  /* optional physiology / hydraulics keys (defaults when absent) */
  xylem_shape, xylem_curve, xylem_conductance, root_conductance, root_compensation,
  collar_pressure_min, max_stomatal_conductance, leaf_width, jmax_vcmax_ratio,
  cica_target, rd_vcmax_ratio, root_length_index, fine_root_radius,
  veg_last_key = fine_root_radius,
  /* terrain information */
  demfile = 0, maskfile, lakefile, 
  lake_name = 0, lake_outlet, lake_scale, lake_exponent,
  soiltype_file = 0, soildepth_file, kslat_file, expdec_file, porosity_file, fc_file, 
  vegtype_file = 0, vegfc_file, veglai_file, vegheight_file, 
  dynaveg_path, dynaveg_num, 
  /* DHSVM channel keys */
  stream_network = 0, stream_map, stream_class,
  /* number of each type of output */
  output_path =
    0, initial_state_path, npixels, nstates, nmapvars,
  /* pixel information */
  north = 0, east, name,
  /* state information */
  state_date = 0,
  /* map information */
  map_variable = 0, map_layer, nmaps, map_date
};

#endif
