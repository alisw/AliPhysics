/*
 * Functions related to the different time conversions in rdmc
 */

#include <math.h>
#include <time.h>
#include "rdmc.h"

#define JD2000   (2451545.0) /* Julian day for Jan 1st 2000 12:00 noon */
                /* = standard epoch to be used as reference */
#define JCENTURY  (36525.0) /* number od days per Julian century */

/************************************************************************/
/* rdmc_o_dateconv() converts a YYMMDD-date into the unix time		*/ 
/*		(sec since 1.1.70) 					*/
/************************************************************************/

time_t rdmc_o_dateconv(int date)
{
  struct tm d;

  if (date <= 0) return RDMC_NA;

  d.tm_sec = 0;               /* second (0-61, allows for leap seconds) */
  d.tm_min = 0;                                        /* minute (0-59) */
  d.tm_hour = 12;                                        /* hour (0-23) */
  d.tm_mday = date % 100;                    /* day of the month (1-31) */
  d.tm_mon = (date/100) % 100 - 1;                      /* month (0-11) */
  d.tm_year = date/10000;                           /* years since 1900 */
  d.tm_isdst = 0;        /* non-0 if daylight savings time is in effect */

  return mktime(&d);                        /* time in sec since 1.1.70 */

} /* function o_dateconv() */

/************************************************************************/
/* rdmc_o_rdateconv() converts the unix time into the YYMMDD format 	*/
/************************************************************************/
int rdmc_o_rdateconv(time_t time)
{
  struct tm *date;

  if (time <= 0) return RDMC_NA;

  date = gmtime(&time);
  return (date->tm_year * 10000 + (date->tm_mon + 1) * 100 + date->tm_mday);
} /* function o_rdateconv() */

/************************************************************************/
/* rdmc_gps_to_mjd() converts the time in GPS/UT year and day into 	*/
/*     Modified Julian Day (50084 < mjd < 51545 for 1996-99)		*/ 
/* NOTES: jdyr = (Julian date of 1200 UT 1/1/year) = GPS/UT day one     */
/************************************************************************/
int rdmc_gps_to_mjd(int gpsyear, int gpsday)

/* reference: http://www.starlink.rl.ac.uk/ function cldj.f 
 * or Zombeck p107 equation for julian date */

{
  int mjdyr, mjd;
  
  if ((gpsyear <= 0) || (gpsday <=0 ) )
    return RDMC_NA;

  mjdyr = (1461*(gpsyear+4799))/4 - (3*((gpsyear+4899)/100))/4 - 2431738.5;
  mjd = (mjdyr + gpsday) - 1;  
  return mjd;

} /* gps_to_mjd() */

/************************************************************************/
/* rdmc_mjd_to_gps() converts the time in Modified Julian Days to       */
/*              GPS/UT year and day within year                         */
/* This uses an approximation which may err high in late December.      */
/* It then checks using the exact calculation of mjdyr and gets gpsday. */
/* NOTES: jdyr = (Julian date of 1200 UT 1/1/year) = GPS/UT day one     */
/************************************************************************/
void rdmc_mjd_to_gps(int mjd, int *gpsyear, int *gpsday)
{
  int year,mjdyr;

  if ( mjd <=0  ){
    *gpsyear=RDMC_NA;
    *gpsday=RDMC_NA;
    return;
  }

  year = (4*(mjd + 678960)) / 1461; /* Approx, erring high */
  mjdyr = (1461*(year+4799))/4 - (3*((year+4899)/100))/4 - 2431738.5;
  if (mjdyr > mjd) {
    year -= 1;
    mjdyr = (1461*(year+4799))/4 - (3*((year+4899)/100))/4 - 2431738.5;
  }
  *gpsyear = year;
  *gpsday = mjd - mjdyr + 1;
  return;
} /* mjd_to_gps() */

/************************************************************************/
/* rdmc_mjd_to_tjd() converts the time in Modified Julian Days to       */
/*              Truncated Julian Days.					*/
/*		this is a trivial conversion but useful for GRB 	*/
/*		analysis with the batse catalog				*/
/************************************************************************/
void rdmc_mjd_to_tjd(int mjd, int secs, int ns, int *tjd, double *sec)
{
  *tjd=mjd-40000.0;
  *sec=(double)secs+(double)ns/1e9;
}

/************************************************************************/
/* rdmc_tjd_to_mjd() converts the time in Truncated Julian Days to      */
/*              Modified Julian Days.                                   */
/*              this is a trivial conversion but useful for GRB         */
/*              analysis with the batse catalog                         */
/************************************************************************/
void rdmc_tjd_to_mjd(int tjd, double sec, int *mjd, int *secs, int *ns)
{
  *mjd = tjd+40000.0;
  *secs = sec;
  *ns = (sec - *secs) * 1e9;
}

/************************************************************************/
/* rdmc_mjd_to_jd() converts the time in Modified Julian Days to       */
/*		Julian Days..                                  		*/
/*              this is a trivial conversion but stands here for 	*/
/*		completeness						*/
/************************************************************************/
void rdmc_mjd_to_jd(int mjd, int secs, int ns, double *jd)
{
  *jd=2400000.5+mjd+((secs+ns/1e9)/86400);
}

/************************************************************************/
/* rdmc_jd_to_mjd() converts the time in Julian Days to			*/
/*              Modified Julian Days.					*/
/*              this is a trivial conversion but stands here for        */
/*              completeness                                            */
/************************************************************************/
void rdmc_jd_to_mjd(double jd, int *mjd, int *secs, int *ns)
{
  *mjd=jd-2400000.5;
  *secs=(jd-(double)*mjd)*86400;
  *ns=((jd-(double)*mjd)*86400-(double)*secs)*1e9;
}
/************************************************************************/
/* rdmc_mjd_to_unixtime() converts the time (unix in sec since 1.1.70)  */
/*                   from mjd and secs                                  */
/************************************************************************/
time_t rdmc_mjd_to_unixtime(int mjd, int secs)
{
  time_t unix_time;

  if ((mjd <= 0) || (secs <= 0)) return RDMC_NA;  /* if no mjd->no unix time*/
  unix_time = (mjd - 40587)* 86400 + secs;     /* 40587 is the 1.1.1970 */
  return unix_time;
} /* mjd_to_unixtime() */

/************************************************************************/
/* rdmc_unixtime_to_mjd() converts UNIX time(sec since 1.1.70) into mjd */
/************************************************************************/
int rdmc_unixtime_to_mjd(time_t unix_time)
{
  if (unix_time <= 0) return RDMC_NA;
  return unix_time/86400 + 40587;
} /* unixtime_to_mjd() */


/************************************************************************/
/* rdmc_unixtime_to_mjdsecs() converts the time (unic in sec since  	*/
/*                       1.1.70) into secs after mjd                    */
/************************************************************************/
int rdmc_unixtime_to_mjdsecs(time_t unix_time)
{
  if (unix_time <= 0) return RDMC_NA;
  return unix_time % 86400;
} /* unixtime_to_mjdsecs() */


/************************************************************************/
/* rdmc_jd_to_gmst() converts mjd date and time to Greenwich mean       */
/*                       sidereal time                          	*/
/************************************************************************/
void rdmc_jd_to_gmst(double jd, double *gmst)
{
/*
 *  GMST in hours
 *  GMST from derived NOVAS package sidereal_time() function.
 *  NOVAS is online at http://aa.usno.navy.mil/aa/
 */

  double tdt, gst;
  int igst;

  tdt = (jd -  JD2000) / JCENTURY;
  gst = ((((-6.2e-6 * tdt + 0.093104) * tdt + 184.812866) * tdt
	   + 67310.54841) + 3.1644e9 * tdt);
  igst= (int)(gst/86400.0);
  *gmst = gst/3600.0 - (double)(igst*24.0);   /* gmst is smaller than 24 */
  if (*gmst < 0.0) *gmst += 24.0;     /* gmst is larger than 0 */
 
} /*  rdmc_jd_to_gmst() */



/************************************************************************/
/* rdmc_gmst_to_jd() converts Greenwich mean sidereal time  to mjd      */
/*                       sidereal time                          	*/
/************************************************************************/
void rdmc_gmst_to_jd(double gmst, int intjd, double *jd)
{
/*
 *  GMST in hours, intjd is integer part of jd
 *  GMST from derived NOVAS package sidereal_time() function.
 *  NOVAS is online at http://aa.usno.navy.mil/aa/
 *  Inverted direction according to functions GMST->UT and UT->GMST
 *  From "practical Astronomy with your Calculator", Cambridge UP 1988
 */ 


  /* Caution: Function untested in any respect */
  double tdt, gst, fracjd;

  tdt = ((double)intjd -  JD2000) / JCENTURY;
  gst = ((((-6.2e-6 * tdt + 0.093104) * tdt + 184.812866) * tdt
	   + 67310.54841) + 3.1644e9 * tdt);
  fracjd = gmst - gst;
  fracjd = fmod ((fracjd / 3600.0), 24.0); 
  if (fracjd < 0.0) fracjd += 24.0;    
  fracjd=fracjd/24.0*0.9972695663;
  *jd=fracjd+(double)intjd;

} /*  rdmc_jd_to_gmst() */

