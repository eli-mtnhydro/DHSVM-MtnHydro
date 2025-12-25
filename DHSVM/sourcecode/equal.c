
/******************************************************************************/
/*				    INCLUDES                                  */
/******************************************************************************/

#include <float.h>
/* #define NO_FLOAT_H */
#include <math.h>
#include "functions.h"

/******************************************************************************/
/*				GLOBAL VARIABLES                              */
/******************************************************************************/
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif
#ifndef FLT_DIGITS
#define FLT_DIGITS 7		/* default sig. digits for float data */
#endif
#ifndef DBL_DIGITS
#define DBL_DIGITS 15		/* default sig. digits for double data */
#endif
static double double_eps;
static float float_eps;
static int initeps = 0;

/******************************************************************************/
/******************************************************************************/
/*			           FUNCTIONS                                  */
/******************************************************************************/
/******************************************************************************/

static double double_epsilon(void);
static float float_epsilon(void);
static void init_epsilons(void);

/******************************************************************************/
/*				     dequal                                   */
/******************************************************************************/
unsigned char dequal(double a, double b)
{
  if (!initeps) {		/* make sure epsilons get initialized */
    init_epsilons();
    initeps = 1;
  }

  /* Two double values only need to be equal to within machine precision */
  if ((a > 0) == (b > 0) &&	/* prevents potential overflow */
      (ABSVAL(a - b) <= ABSVAL(double_eps * b)))
    return TRUE;
  else
    return FALSE;
}

/******************************************************************************/
/*				     fequal                                   */
/******************************************************************************/
unsigned char fequal(float a, float b)
{
  if (!initeps) {		/* make sure epsilons get initialized */
    init_epsilons();
    initeps = 1;
  }

  /* Two float values anly need to be equal to within machine precision */
  if ((a > 0) == (b > 0) &&	/* prevents potential overflow */
      (ABSVAL(a - b) <= ABSVAL(float_eps * b)))
    return TRUE;
  else
    return FALSE;
}

/******************************************************************************/
/*				 double_epsilon                               */
/******************************************************************************/
static double double_epsilon(void)
{
  double double_eps;
#ifndef NO_FLOAT_H
  double_eps = DBL_EPSILON;
#else /* NO_FLOAT_H */
  {
    double etop, ebot, eps;
    double one = 1.0;
    double two = 2.0;
    etop = 1.0;
    ebot = 0.0;
    eps = ebot + (etop - ebot) / two;
    while (eps != ebot && eps != etop) {
      double epsp1;

      epsp1 = one + eps;
      if (epsp1 > one)
	etop = eps;
      else
	ebot = eps;
      eps = ebot + (etop - ebot) / two;
    }
    double_eps = two * etop;
  }
#endif /* NO_FLOAT_H */
  return double_eps;
}

/******************************************************************************/
/*				 float_epsilon                                */
/******************************************************************************/
static float float_epsilon(void)
{
  float float_eps;
#ifndef NO_FLOAT_H
  float_eps = FLT_EPSILON;
#else /* NO_FLOAT_H */
  {
    float etop, ebot, eps;
    float one = 1.0;
    float two = 2.0;
    etop = 1.0;
    ebot = 0.0;
    eps = ebot + (etop - ebot) / two;
    while (eps != ebot && eps != etop) {
      float epsp1;

      epsp1 = one + eps;
      if (epsp1 > one)
	etop = eps;
      else
	ebot = eps;
      eps = ebot + (etop - ebot) / two;
    }
    float_eps = two * etop;
  }
#endif /* NO_FLOAT_H */
  return float_eps;
}

/******************************************************************************/
/*				 init_epsilons                                */
/******************************************************************************/
static void init_epsilons(void)
{
  float_eps = float_epsilon();
  double_eps = double_epsilon();
}
