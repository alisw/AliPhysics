/*>    ROUTINE TIMEL
  CERN PROGLIB# Z007    TIMEST          .VERSION KERNFOR  4.38  931108
  ORIG. 01/03/85  FCA, mod 03/11/93 GF
  Version for Windows NT/Windows 95 by Valery Fine 30/05/96 (fine@vxcern.cern.ch)
*/
#include <windows.h>
#include <sys\types.h>

#ifndef gTicks
# define gTicks 1.0e-7;
#endif

#define time_t double

static float timlim;
static time_t timstart, timlast;
static HANDLE hProcess;


static int tml_init = 1;
double deftim = 999.;

#if defined(CERNLIB_QXCAPT)
#define timest type_of_call TIMEST
#define timex  type_of_call TIMEX
#define timed  type_of_call TIMED
#define timel  type_of_call TIMEL
#elif defined(CERNLIB_QX_SC)
#define timest type_of_call timest_
#define timex  type_of_call timex_
#define timed  type_of_call timed_
#define timel  type_of_call timel_
#endif


//______________________________________________________________________________
double GetRealTime(){
  union     {FILETIME ftFileTime;
             __int64  ftInt64;
            } ftRealTime; // time the process has spent in kernel mode
  SYSTEMTIME st;
  GetSystemTime(&st);
  SystemTimeToFileTime(&st,&ftRealTime.ftFileTime);
  return (double)ftRealTime.ftInt64 * gTicks;
}
 
//______________________________________________________________________________
double GetCPUTime(){
 OSVERSIONINFO OsVersionInfo;
 
//*-*         Value                      Platform
//*-*  ----------------------------------------------------
//*-*  VER_PLATFORM_WIN32s              Win32s on Windows 3.1
//*-*  VER_PLATFORM_WIN32_WINDOWS       Win32 on Windows 95
//*-*  VER_PLATFORM_WIN32_NT            Windows NT
//*-*
  OsVersionInfo.dwOSVersionInfoSize=sizeof(OSVERSIONINFO);
  GetVersionEx(&OsVersionInfo);
  if (OsVersionInfo.dwPlatformId == VER_PLATFORM_WIN32_NT) {
    DWORD       ret;
    FILETIME    ftCreate,       // when the process was created
                ftExit;         // when the process exited
 
    union     {FILETIME ftFileTime;
               __int64  ftInt64;
              } ftKernel; // time the process has spent in kernel mode
 
    union     {FILETIME ftFileTime;
               __int64  ftInt64;
              } ftUser;   // time the process has spent in user mode
 
    ret = GetProcessTimes (hProcess, &ftCreate, &ftExit,
                                     &ftKernel.ftFileTime,
                                     &ftUser.ftFileTime);
    if (ret != TRUE){
      ret = GetLastError ();
      printf("GetCPUTime", " Error on GetProcessTimes 0x%lx", (int)ret);
    }
 
    /*
     * Process times are returned in a 64-bit structure, as the number of
     * 100 nanosecond ticks since 1 January 1601.  User mode and kernel mode
     * times for this process are in separate 64-bit structures.
     * To convert to floating point seconds, we will:
     *
     *          Convert sum of high 32-bit quantities to 64-bit int
     */
 
      return (double) (ftKernel.ftInt64 + ftUser.ftInt64) * gTicks;
  }
  else
      return GetRealTime();
 
}

                   /*  local routine called by timst, and time_init */
//_______________________________________________________________
static void time_st(timl)
float timl;
{
    hProcess = GetCurrentProcess();
    timstart = GetCPUTime();
    timlast  = timstart;
    timlim   = timl;
    tml_init = 0;
    return;
}
                   /*  local routine to start by default  */
//_______________________________________________________________
static void time_init()
{
  float  maxtime;
  maxtime=deftim;
  time_st(maxtime);
  return;
}

//_______________________________________________________________
void timest(timl)
float *timl;
{
 float  maxtime;

 if (tml_init != 0) {
    maxtime = *timl;
    time_st(maxtime);
 }
 return;
}
//_______________________________________________________________
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
     timnow= GetCPUTime();
    *tx = (float) (timnow - timstart);
   }
   return;
}

//_______________________________________________________________
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
       timnow=GetCPUTime();
       *td = (float) (timnow - timlast);
       timlast = timnow;
   }
   return;
}

//_______________________________________________________________
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
       timnow= GetCPUTime();
       *tl = timlim - (float) (timnow - timstart);
   }
   return;
}
/*> END <----------------------------------------------------------*/
