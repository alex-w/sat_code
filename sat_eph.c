#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>
#include "watdefs.h"
#include "afuncs.h"
#include "date.h"
#include "norad.h"
#include "mpc_func.h"
#include "observe.h"

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923


typedef struct
{
   double lat, lon, alt, rho_sin_phi, rho_cos_phi;
   double jd_start, jd_end, step_size;
   int n_steps;
   const char *desig;
} ephem_t;

static int verbose = 0;

static char *fgets_trimmed( char *buff, const int buffsize, FILE *ifile)
{
   char *rval = fgets( buff, buffsize, ifile);

   if( rval)
      {
      size_t i = 0;

      while( rval[i] != 10 && rval[i] != 13 && rval[i])
         i++;
      while( i && rval[i - 1] == ' ')
         i--;        /* drop trailing spaces */
      rval[i] = '\0';
      }
   return( rval);
}

static void show_base_60( char *buff, const unsigned n_millisec)
{
   snprintf( buff, 15, "%03u %02u %02u.%03u",
               n_millisec / 3600000u, (n_millisec / 60000u) % 60u,
               (n_millisec / 1000u) % 60u, n_millisec % 1000u);
}

static void put_ra_in_buff( char *buff, double ra)
{
   ra = fmod( ra, 2. * PI);
   if( ra < 0.)
      ra += PI + PI;
   show_base_60( buff, (unsigned)( 3600. * 1000. * ra * 12. / PI));
   memmove( buff, buff + 1, strlen( buff));     /* remove leading zero */
}

static void put_dec_in_buff( char *buff, const double dec)
{
   show_base_60( buff, (unsigned)( 3600. * 1000. * fabs( dec) * 180. / PI));
   *buff = (dec > 0. ? '+' : '-');
}

static double angle_between( const double *a, const double *b)
{
   const double cos_ang = dot_product( a, b) /
                              sqrt( dot_product( a, a) * dot_product( b, b));
   double rval = acose( cos_ang);

   return( rval * 180. / PI);
}

static inline bool desig_match( const tle_t *tle, const char *desig)
{
   size_t i = 0;
   bool rval = false;

   while( isdigit( desig[i]))
      i++;
   if( i == 5)
      {
      if( !desig[i])       /* desig is all digits -> it's the NORAD # */
         rval = (atoi( desig) == tle->norad_number);
      else
         {
         i = strlen( desig);
         if( i > 5 && i < 9)
            rval = !memcmp( tle->intl_desig, desig, i) && tle->intl_desig[i] <= ' ';
         }
      }
   return( rval);
}

static const char *header_text =
           "Date (UTC)  Time       R.A. (J2000)  decl   Alt   Elong   Dist(km)\n";

static int show_ephems_from( const char *path_to_tles, const ephem_t *e,
                                  const char *filename)
{
   FILE *ifile;
   char line0[100], line1[100], line2[100];
   int show_it = 1, n_lines_generated = 0;
   double jd_tle = 0., tle_range = 1e+10;

   if( verbose)
      printf( "Should examine '%s'\n", filename);
   snprintf( line0, sizeof( line0), "%s/%s", path_to_tles, filename);
   ifile = fopen( line0, "rb");
   if( !ifile)
      {
      printf( "'%s' not opened\n", line0);
      exit( 0);
      }
   *line0 = *line1 = '\0';
   while( fgets_trimmed( line2, sizeof( line2), ifile))
      {
      tle_t tle;

      if( !memcmp( line2, "# MJD ", 6))
         {
         jd_tle = atof( line2 + 6) + 2400000.5;
         tle_range = 1.;
         show_it = (jd_tle < e->jd_end && jd_tle + tle_range > e->jd_start);
         }
      else if( show_it && parse_elements( line1, line2, &tle) >= 0
                     && desig_match( &tle, e->desig))
         {
         double sat_params[N_SAT_PARAMS], jd = e->jd_start;
         size_t i, j;
         const int is_deep_type = select_ephemeris( &tle);

         if( is_deep_type)
            SDP4_init( sat_params, &tle);
         else
            SGP4_init( sat_params, &tle);
         if( verbose > 1)
            {
            printf( "Got TLEs for %f :\n", jd);
            printf( "%s\n%s\n%s\n", line0, line1, line2);
            }
         for( i = 0; i < (size_t)e->n_steps; i++, jd += e->step_size)
            if( jd >= jd_tle && jd < jd_tle + tle_range)
               {
               char buff[90], dec_buff[20], ra_buff[20];
               double pos[3], vel[3], obs_pos[3], ra, dec, dist;
               const double t_since = (jd - tle.epoch) * minutes_per_day;
               double solar_xyzr[4], topo_posn[3];

               if( !i)
                  {
                  printf( "\nEphemerides for %05d = %s%.2s-%s\n",
                              tle.norad_number,
                              (atoi( tle.intl_desig) > 57000) ? "19" : "20",
                              tle.intl_desig, tle.intl_desig + 2);
                  printf( "%s\n%s", line0, header_text);
                  }
               full_ctime( buff, jd, FULL_CTIME_YMD | FULL_CTIME_MONTHS_AS_DIGITS
                                 | FULL_CTIME_LEADING_ZEROES);
               if( is_deep_type)
                  SDP4( t_since, &tle, sat_params, pos, vel);
               else
                  SGP4( t_since, &tle, sat_params, pos, vel);
               observer_cartesian_coords( jd, e->lon, e->rho_cos_phi,
                                        e->rho_sin_phi, obs_pos);
               get_satellite_ra_dec_delta( obs_pos, pos, &ra, &dec, &dist);
               epoch_of_date_to_j2000( jd, &ra, &dec);
               put_ra_in_buff( ra_buff, ra);
               put_dec_in_buff( dec_buff, dec);
               ra_buff[10] = dec_buff[9] = '\0';
               for( j = 0; j < 3; j++)
                  topo_posn[j] = pos[j] - obs_pos[j];
               lunar_solar_position( jd, NULL, solar_xyzr);
               ecliptic_to_equatorial( solar_xyzr);
               printf( "%s  %s  %s %+05.1f %6.1f %8.0f\n", buff, ra_buff, dec_buff,
                     90. - angle_between( topo_posn, obs_pos),     /* alt */
                     angle_between( topo_posn, solar_xyzr), dist);
               n_lines_generated++;
               }
         }
      strcpy( line0, line1);
      strcpy( line1, line2);
      }
   fclose( ifile);
   return( n_lines_generated);
}

int generate_artsat_ephems( const char *path_to_tles, const ephem_t *e)
{
   FILE *ifile;
   char buff[200];
   int is_in_range = 0, id_matches = 1, ephem_lines_generated = 0;

   snprintf( buff, sizeof( buff), "%s/tle_list.txt", path_to_tles);
   ifile = fopen( buff, "rb");
   if( !ifile)
      {
      printf( "'%s' not opened\n", buff);
      exit( 0);
      }
   while( ephem_lines_generated != e->n_steps &&
                          fgets_trimmed( buff, sizeof( buff), ifile))
      {
      if( !memcmp( buff, "# Range:", 8))
         {
         buff[19] = '\0';
         if( get_time_from_string( 0., buff + 9, FULL_CTIME_YMD, NULL) < e->jd_end
                  && get_time_from_string( 0., buff + 20, FULL_CTIME_YMD, NULL) > e->jd_start)
            is_in_range = 1;
         }
      if( !memcmp( buff, "# ID:", 5))
         if( strcmp( e->desig, buff + 13) && atoi( buff + 5) != atoi( e->desig))
            id_matches = 0;
      if( !memcmp( buff, "# Include ", 10))
         {
         if( is_in_range && id_matches)
            ephem_lines_generated += show_ephems_from( path_to_tles, e, buff + 10);
         is_in_range = 0;
         id_matches = 1;
         }
      }
   fclose( ifile);
   if( ephem_lines_generated)
      printf( "%s", header_text);
   return( ephem_lines_generated);
}

static int set_location( ephem_t *e, const char *mpc_code, const char *obscode_file_name)
{
   FILE *ifile = fopen( obscode_file_name, "rb");
   char buff[200];
   int got_it = 0, planet;

   if( !ifile)
      {
      printf( "'%s' not found\n", obscode_file_name);
      exit( 0);
      }
   while( !got_it && fgets_trimmed( buff, sizeof( buff), ifile))
      if( !memcmp( mpc_code, buff, 3))
         {
         mpc_code_t c;

         planet = get_mpc_code_info( &c, buff);
         if( planet != 3)
            {
            printf( "MPC code '%s' is for planet %d\n",
                        mpc_code, planet);
            exit( 0);
            }
         e->lat = c.lat;
         e->lon = c.lon;
         e->alt = c.alt;
         e->rho_cos_phi = c.rho_cos_phi;
         e->rho_sin_phi = c.rho_sin_phi;
         got_it = 1;
         }
   fclose( ifile);
   return( got_it ? 0 : -1);
}

#ifdef ON_LINE_VERSION
   #define OBSCODES_DOT_HTML_FILENAME  "/home/projectp/public_html/cgi-bin/fo/ObsCodes.htm"
   #define PATH_TO_TLES                "/home/projectp/public_html/tles"
#else
   #define OBSCODES_DOT_HTML_FILENAME  "/home/phred/.find_orb/ObsCodes.htm"
   #define PATH_TO_TLES                "/home/phred/tles"
#endif

static const char *get_arg( const char **argv)
{
   const char *rval;

   if( argv[0][0] == '-' && argv[0][1])
      {
      if( !argv[0][2] && argv[1])
         rval = argv[1];
      else
         rval = argv[0] + 2;
      }
   else
      rval = NULL;
   if( !rval)
      {
      printf( "Can't get an argument : '%s'\n", argv[0]);
      exit( 0);
      }
   return( rval);
}

static void error_help( void)
{
   printf( "Mandatory 'sat_eph' arguments:\n"
           "   -c(MPC code) : specify location\n"
           "   -t(date/time) : starting time of ephemeris\n"
           "   -n(#) : number of ephemeris steps\n"
           "   -s(#) : ephemeris step size in days\n"
           "   -o(#) : five digit NORAD number\n"
           "Optional arguments:\n"
           "   -v(#) : level of verbosity\n");
}

int dummy_main( const int argc, const char **argv)
{
   int i;
   ephem_t e;

   if( argc < 3)
      {
      error_help( );
      return( 0);
      }
   memset( &e, 0, sizeof( ephem_t));
   for( i = 1; i < argc; i++)
      if( argv[i][0] == '-' && argv[i][1])
         {
         const char *arg = get_arg( argv + i);

         switch( argv[i][1])
            {
            case 'c':
               set_location( &e, arg, OBSCODES_DOT_HTML_FILENAME);
               break;
            case 't':
               e.jd_start = get_time_from_string( 0., arg, FULL_CTIME_YMD, NULL);
               break;
            case 'n':
               e.n_steps = atoi( arg);
               break;
            case 's':
               e.step_size = atof( arg);
               break;
            case 'o':
               /* Will handle below */
               break;
            case 'v':
               verbose = 1 + atoi( arg);
               break;
            default:
               printf( "Unrecognized option '%s'\n", argv[i]);
               error_help( );
               return( 0);
            }
         }
   e.jd_end   = e.jd_start + (double)e.n_steps * e.step_size;
   if( verbose)
      printf( "arguments parsed;  JDs %f to %f\n", e.jd_start, e.jd_end);
   for( i = 1; i < argc; i++)
      if( argv[i][0] == '-' && argv[i][1] == 'o')
         {
         e.desig = get_arg( argv + i);
         generate_artsat_ephems( PATH_TO_TLES, &e);
         }
   return( 0);
}

#ifndef ON_LINE_VERSION
int main( const int argc, const char **argv)
{
   return( dummy_main( argc, argv));
}
#else

#include <errno.h>
#ifdef __has_include
   #if __has_include(<cgi_func.h>)
       #include "cgi_func.h"
   #else
       #error   \
         'cgi_func.h' not found.  This project depends on the 'lunar'\
         library.  See www.github.com/Bill-Gray/lunar .\
         Clone that repository,  'make'  and 'make install' it.
   #endif
#else
   #include "cgi_func.h"
#endif

int main( const int unused_argc, const char **unused_argv)
{
   const char *argv[20];
   const size_t max_buff_size = 40000;       /* room for 500 obs */
   char *buff = (char *)malloc( max_buff_size);
   char field[30], time_text[80], obj_name[30];
   char num_steps[30], step_size[30], obs_code[10];
   FILE *lock_file = fopen( "lock.txt", "w");
   int cgi_status, i;
#ifndef _WIN32
   extern char **environ;

   avoid_runaway_process( 15);
#endif         /* _WIN32 */
   setbuf( lock_file, NULL);
   INTENTIONALLY_UNUSED_PARAMETER( unused_argc);
   INTENTIONALLY_UNUSED_PARAMETER( unused_argv);
   printf( "Content-type: text/html\n\n");
   printf( "<html> <body> <pre>\n");
   if( !lock_file)
      {
      printf( "<p> Server is busy.  Try again in a minute or two. </p>");
      printf( "<p> Your artsat ephemerides are very important to us! </p>");
      return( 0);
      }
   fprintf( lock_file, "We're in\n");
   *time_text = *obj_name = *num_steps = *step_size = *obs_code = '\0';
#ifndef _WIN32
   for( size_t i = 0; environ[i]; i++)
      fprintf( lock_file, "%s\n", environ[i]);
#endif
   cgi_status = initialize_cgi_reading( );
   fprintf( lock_file, "CGI status %d\n", cgi_status);
   if( cgi_status <= 0)
      {
      printf( "<p> <b> CGI data reading failed : error %d </b>", cgi_status);
      printf( "This isn't supposed to happen.</p>\n");
      return( 0);
      }
   while( !get_cgi_data( field, buff, NULL, max_buff_size))
      {
      fprintf( lock_file, "Field '%s'\n", field);
      if( !strcmp( field, "time") && strlen( buff) < sizeof( time_text))
         strcpy( time_text, buff);
      if( !strcmp( field, "obj_name") && strlen( buff) < sizeof( obj_name))
         strcpy( obj_name, buff);
      if( !strcmp( field, "num_steps") && strlen( buff) < sizeof( num_steps))
         strcpy( num_steps, buff);
      if( !strcmp( field, "step_size") && strlen( buff) < sizeof( step_size))
         {
         char *tptr = strchr( buff, 'v');

         if( tptr)
            {
            verbose = atoi( tptr + 1);
            *tptr = '\0';
            }
         strcpy( step_size, buff);
         }
      if( !strcmp( field, "obs_code") && strlen( buff) < sizeof( obs_code))
         strcpy( obs_code, buff);
      }
   fprintf( lock_file, "Fields read\n");
   argv[0] = "sat_eph";
   argv[1] = "-t";
   argv[2] = time_text;
   argv[3] = "-o";
   argv[4] = obj_name;
   argv[5] = "-n";
   argv[6] = num_steps;
   argv[7] = "-s";
   argv[8] = step_size;
   argv[9] = "-c";
   argv[10] = obs_code;
   argv[11] = NULL;
   for( i = 0; argv[i]; i++)
      fprintf( lock_file, "arg %d: '%s'\n", (int)i, argv[i]);
   dummy_main( 11, argv);
   fprintf( lock_file, "dummy_main called\n");
   free( buff);
   printf( "On-line Sat_eph compiled " __DATE__ " " __TIME__ " UTC-5h\n");
   printf( "See <a href='https://www.github.com/Bill-Gray/sat_code'>"
               "https://www.github.com/Bill-Gray/sat_code</a> for source code\n");
   printf( "</pre> </body> </html>");
   return( 0);
}
#endif
