
#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#include "errorhandler.h"

static FILE *LOG = NULL;
static const char *Program = "unknown program";
static ErrorLevel Level;

/* -------------------------------------------------------------
   error_handler_init
   ------------------------------------------------------------- */
int
error_handler_init(const char *program, const char *logfile,
		   ErrorLevel debug_level)
{
  if (program != NULL)
    Program = program;

  if (debug_level < ERRHDL_ERROR) {
    error_handler(ERRHDL_WARNING,
		  "error_handler_init: specified debug level (%d) too low changeing to ERRHDL_ERROR",
		  debug_level);
    Level = ERRHDL_ERROR;
  }
  else {
    Level = debug_level;
  }

  if (logfile != NULL) {
    if ((LOG = fopen(logfile, "w")) == NULL) {
      error_handler(ERRHDL_ERROR,
		    "error_handler_init: unable to open log file \"%s\": %s",
		    logfile, strerror(errno));
      LOG = NULL;
    }
  }
  return (0);
}

/* -------------------------------------------------------------
   error_handler
   ------------------------------------------------------------- */
void error_handler(ErrorLevel debug_level, const char *fmt, ...)
{
  va_list ap;
  FILE *out = (LOG == NULL) ? stderr : LOG;
  char buffer[256];

  if (Level < debug_level)
    return;

  /* Errors stop the run with exit code 1 */
  if (debug_level <= ERRHDL_ERROR)
    snprintf(buffer, sizeof(buffer), "DHSVM ERROR 1: %s\n", fmt);
  else if (debug_level == ERRHDL_WARNING)
    snprintf(buffer, sizeof(buffer), "DHSVM WARNING: %s\n", fmt);
  else
    snprintf(buffer, sizeof(buffer), "%s: %s\n", Program, fmt);

  va_start(ap, fmt);
  if (vfprintf(out, buffer, ap) <= 0) {
    if (out != stderr)
      fprintf(stderr,
	      "%s: print_error_msg: error writing to log file: %s",
	      Program, (char *) strerror(errno));
    else
      abort();
  }
  va_end(ap);
  if (debug_level <= ERRHDL_ERROR) {
    fflush(out);
    error_handler_done();
    exit(1);
  }
  fflush(out);

  return;
}

/* -------------------------------------------------------------
   error_handler_done
   ------------------------------------------------------------- */
int error_handler_done(void)
{
  if (LOG != NULL) {
    if (fclose(LOG) != 0) {
      fprintf(stderr,
	      "%s: error_handler_done: unable to close log file: %s",
	      Program, (char *) strerror(errno));
      return (-1);
    }
  }
  LOG = NULL;
  return (0);
}
