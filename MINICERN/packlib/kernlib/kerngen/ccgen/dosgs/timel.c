/*
 * $Id$
 *
 * $Log$
 * Revision 1.1.1.1  1996/02/15 17:49:31  mclareni
 * Kernlib
 *
 */
/*>    ROUTINE TIMEL
  CERN PROGLIB# Z007    TIMEST          .VERSION KERNFOR  4.38  931108
  ORIG. 01/03/85  FCA, mod 03/11/93 GF
*/
#ifdef WIN32
#include <sys\types.h>
#else
#include <sys/types.h>
#endif
#include <time.h>
#ifndef WIN32
#include <sys/times.h>
#include <sys/param.h>
#include <sys/time.h>
#include <sys/resource.h>
#endif

#ifndef CLOCKS_PER_SEC
#define  CLOCKS_PER_SEC CLK_TCK
#endif


#ifndef RLIMIT_CPU
#define RLIMIT_CPU 0    /* For HP-UX... */
#endif
#ifndef RLIM_INFINITY
#define RLIM_INFINITY 0x7fffffff    /* For HP-UX... */
#endif

#if defined(CERNLIB_QSYSBSD)||defined(CERNLIB_QMVMI)||defined(CERNLIB_QMVAOS)
#define HZ 60.;
#endif

#ifndef HZ
#ifdef __GNUC__
#define HZ 1
#else
#define HZ 1./CLOCKS_PER_SEC
#endif
#endif

#if defined(CERNLIB_QMDOS)
#ifdef __GNUC__
#define sectim   1./timsec
#define time_t double
double timfac  = 1./CLOCKS_PER_SEC;
double timsec  = (double) 0x7fffffff/CLOCKS_PER_SEC;
#include <std.h>
#endif
#ifdef WIN32

/**************************************************************************\
*
*       For some reason 'clock()' returns elapsed time on at least the AXP
*       systems.  Substitute the routine 'ntclock()' for 'clock()'.
*
\**************************************************************************/

#define clock ntclock

#endif

#endif
#if !defined(CERNLIB_QMDOS)
struct tms tps;
#endif
static float timlim;
static time_t timstart, timlast;
static int tml_init = 1;
float deftim = 999.;

#if defined(CERNLIB_QX_SC)
#define timest timest_
#define timex  timex_
#define timed  timed_
#define timel  timel_
#endif
#if defined(CERNLIB_QXCAPT)
#define timest TIMEST
#define timex  TIMEX
#define timed  TIMED
#define timel  TIMEL
#endif

                   /*  local routine called by timst, and time_init */
static void time_st(timl)
float timl;
{
#if !defined(CERNLIB_QMDOS)
    times(&tps);
#endif
    timlim = timl;
#if !defined(CERNLIB_QMDOS)
    timstart =  tps.tms_utime+tps.tms_cutime+tps.tms_stime+tps.tms_cstime;
#endif
#if defined(CERNLIB_QMDOS)
#ifdef __GNUC__
    timstart= (long)((time(NULL)&0xfffff)*sectim) * timsec
              + (double)(clock()*timfac);
#else
    timstart= clock();
#endif

#endif
    timlast  = timstart;
    tml_init = 0;
    return;
}
                   /*  local routine to start by default  */
static void time_init()
{
#if !defined(CERNLIB_QMDOS)
        struct rlimit rlimit;
#endif
        float  maxtime;

        maxtime=deftim;

#if !defined(CERNLIB_QMDOS)
        if (getrlimit(RLIMIT_CPU, &rlimit)==0) {
                if ( rlimit.rlim_cur != RLIM_INFINITY )
                   maxtime = (float) rlimit.rlim_cur;
        }

#endif
        time_st(maxtime);
        return;
}

void timest(timl)
float *timl;
{
#if !defined(CERNLIB_QMDOS)
 struct rlimit rlimit;
#endif
 float  maxtime;

 if (tml_init != 0) {

/*  get maximum time allowed by system, and do not allow more */
    maxtime = *timl;
#if !defined(CERNLIB_QMDOS)
    if (getrlimit(RLIMIT_CPU, &rlimit)==0) {
           maxtime = (float) rlimit.rlim_cur;
           maxtime = ( maxtime > *timl ) ? *timl : maxtime;
    }
#endif
    time_st(maxtime);
 }
 return;
}
void timex(tx)
/*
C
  CERN PROGLIB# Z007    TIMEX           .VERSION KERNFOR  4.38  931108
C
*/
float *tx;
{
   time_t timnow;
   if (tml_init) {
       time_init();
       *tx = 0.;
   }
   else {
#if !defined(CERNLIB_QMDOS)
       times(&tps);
       timnow = tps.tms_utime+tps.tms_cutime+tps.tms_stime+tps.tms_cstime;
#endif
#if defined(CERNLIB_QMDOS)
#ifdef _MSDOS_
       timnow= (long)((time(NULL)&0xfffff)*sectim) * timsec
               + (double)(clock()*timfac);
#else
       timnow= clock();
#endif
#endif
       *tx = (float) (timnow - timstart) / HZ;
   }
   return;
}

void timed(td)
/*
C
  CERN PROGLIB# Z007    TIMED           .VERSION KERNFOR  4.38  931108
C
*/
float *td;
{
   time_t timnow;
   if (tml_init) {
       time_init();
       *td = timlim;
   }
   else {
#if !defined(CERNLIB_QMDOS)
       times(&tps);
       timnow = tps.tms_utime+tps.tms_cutime+tps.tms_stime+tps.tms_cstime;
#endif
#if defined(CERNLIB_QMDOS)
#ifdef _MSDOS_
       timnow= (long)((time(NULL)&0xfffff)*sectim) * timsec
               + (double)(clock()*timfac);
#else
       timnow=clock();
#endif
#endif
       *td = (float) (timnow - timlast) / HZ;
       timlast = timnow;
   }
   return;
}

void timel(tl)
/*
C
  CERN PROGLIB# Z007    TIMEL           .VERSION KERNFOR  4.38  931108
C
*/
float *tl;
{
   time_t timnow;
   if (tml_init) {
       time_init();
       *tl = timlim;
   }
   else {
#if !defined(CERNLIB_QMDOS)
       times(&tps);
       timnow = tps.tms_utime+tps.tms_cutime+tps.tms_stime+tps.tms_cstime;
#endif
#if defined(CERNLIB_QMDOS)
#ifdef _MSDOS_
       timnow= (long)((time(NULL)&0xfffff)*sectim) * timsec
               + (double)(clock()*timfac);
#else
       timnow= clock();
#endif
#endif
       *tl = timlim - (float) (timnow - timstart) / HZ;
   }
   return;
}
#ifdef __GNUC__
#undef time_t
#endif

#ifdef WIN32
#undef clock

/**************************************************************************\
*  ntclock.c -- function to return sum of user and kernel time
*
*       For some reason 'clock()' returns elapsed time on at least the AXP
*       systems.  Call native WIN32 routines to get process CPU time, and
*       return in the same (weird) units defined by 'clock()'.
*
\**************************************************************************/

#include <windows.h>
#include <stdio.h>

clock_t ntclock ()
{
    double      tTotal;
    clock_t     cTotal;
    DWORD       ret;
    FILETIME    ftKernel, ftUser, ftCreate, ftExit;

    static HANDLE hProcess = 0;

    cTotal = 0;

    if (hProcess == 0)
      hProcess = GetCurrentProcess();

    ret = GetProcessTimes (hProcess, &ftCreate, &ftExit, &ftKernel, &ftUser);
    if (ret != TRUE){
      ret = GetLastError ();
      printf ("\n* Error on GetProcessTimes in ntclock()  0x%lx", (int)ret);
      return cTotal;
    }

    /*
     * Process times are returned in a 64-bit structure, as the number of
     * 100 nanosecond ticks since 1 January 1601.  User mode and kernel mode
     * times for this process are in separate 64-bit structures.
     * To convert to floating point seconds, we will:
     *
     *          Convert sum of high 32-bit quantities to 64-bit float
     *          Multiply by 2**32
     *          Add low 32-bit tick counts
     *          Divide by 10,000,000 to convert to seconds
     *          Convert to 'clock_t' value in "clock ticks"
     */

    tTotal = ftKernel.dwHighDateTime + ftUser.dwHighDateTime;
    tTotal = tTotal * 4294967296.0;
    tTotal = tTotal + (ftKernel.dwLowDateTime + ftUser.dwLowDateTime);
    tTotal = tTotal * 1.0E-7;
    cTotal = tTotal * CLOCKS_PER_SEC;

    return cTotal;
}

#endif

/*> END <----------------------------------------------------------*/
