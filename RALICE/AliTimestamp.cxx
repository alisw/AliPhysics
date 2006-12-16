/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// $Id$

///////////////////////////////////////////////////////////////////////////
// Class AliTimestamp
// Handling of timestamps for (astro)particle physics reserach.
//
// This class is derived from TTimeStamp and provides additional
// facilities (e.g. Julian date) which are commonly used in the
// field of (astro)particle physics.
//
// The Julian Date (JD) indicates the number of days since noon (UT) on
// 01 jan -4712 (i.e. noon 01 jan 4713 BC), being day 0 of the Julian calendar.
//
// The Modified Julian Date (MJD) indicates the number of days since midnight
// (UT) on 17-nov-1858, which corresponds to 2400000.5 days after day 0 of the
// Julian calendar.
//
// The Truncated Julian Date (TJD) corresponds to 2440000.5 days after day 0
// of the Julian calendar and consequently TJD=MJD-40000.
// This TJD date indication was used by the Vela and Batse missions in
// view of Gamma Ray Burst investigations.
//
// The Julian Epoch (JE) indicates the fractional elapsed year count since
// midnight (UT) on 01-jan at the start of the Gregorian year count.
// A year is defined to be 365.25 days, so the integer part of JE corresponds
// to the usual Gregorian year count.
// So, 01-jan-1965 00:00:00 UT corresponds to JE=1965.0
//
// Because of the fact that the Julian date indicators are all w.r.t. UT
// they provide an absolute timescale irrespective of timezone or daylight
// saving time (DST).
//
// This AliTimestamp facility allows for picosecond precision, in view
// of time of flight analyses for particle physics experiments.
// For normal date/time indication the standard nanosecond precision
// will in general be sufficient.
// Note that when the fractional JD, MJD and TJD counts are used instead
// of the integer (days,sec,ns) specification, the nanosecond precision
// may be lost due to computer accuracy w.r.t. floating point operations.
//
// The TTimeStamp EPOCH starts at 01-jan-1970 00:00:00 UTC
// which corresponds to JD=2440587.5 or the start of MJD=40587 or TJD=587.
// Using the corresponding MJD of this EPOCH allows construction of
// the yy-mm-dd hh:mm:ss:ns TTimeStamp from a given input (M/T)JD and time.
// Obviously this TTimeStamp implementation would prevent usage of values
// smaller than JD=2440587.5 or MJD=40587 or TJD=587.
// Furthermore, due to a limitation on the "seconds since the EPOCH start" count
// in TTimeStamp, the latest accessible date/time is 19-jan-2038 02:14:08 UTC.
// However, this AliTimestamp facility provides support for the full range
// of (M/T)JD values, but the setting of the corresponding TTimeStamp parameters
// is restricted to the values allowed by the TTimeStamp implementation.
// For these earlier/later (M/T)JD values, the standard TTimeStamp parameters will
// be set corresponding to the start of the TTimeStamp EPOCH.
// This implies that for these earlier/later (M/T)JD values the TTimeStamp parameters
// do not match the Julian parameters of AliTimestamp.
// As such the standard TTimeStamp parameters do not appear on the print output
// when invoking the Date() memberfunction for these earlier/later (M/T)JD values.  
//
// Examples :
// ==========
//
// Note : All TTimeStamp functionality is available as well.
//
// AliTimestamp t;
//
// t.Date();
// 
// // Retrieve Julian Date
// Int_t jd,jsec,jns;
// t.GetJD(jd,jsec,jns);
//
// // Retrieve fractional Truncated Julian Date
// Double_t tjd=t.GetTJD();
//
// // Retrieve fractional Julian Epoch
// Double_t je=t.GetJE();
//
// // Set to a specific Modified Julian Date
// Int_t mjd=50537;
// Int_t mjsec=1528;
// Int_t mjns=185643;
// t.SetMJD(mjd,mjsec,mjns);
//
// t.Date();
//
// // Time intervals for e.g. trigger or TOF analysis
// AliEvent evt;
// AliTrack* tx=evt.GetTrack(5);
// AliTimestamp* timex=tx->GetTimestamp();
// Double_t dt=evt.GetDifference(timex,"ps");
// AliTimestamp trig((AliTimestamp)evt);
// trig.Add(0,0,2,173);
// AliSignal* sx=evt.GetHit(23);
// AliTimestamp* timex=sx->GetTimestamp();
// Double_t dt=trig.GetDifference(timex,"ps");
// Int_t d,s,ns,ps;
// trig.GetDifference(timex,d,s,ns,ps);
//
// // Some practical conversion facilities
// // Note : They don't influence the actual date/time settings
// //        and as such can also be invoked as AliTimestamp::Convert(...) etc...
// Int_t y=1921;
// Int_t m=7;
// Int_t d=21;
// Int_t hh=15;
// Int_t mm=23;
// Int_t ss=47;
// Int_t ns=811743;
// Double_t jdate=t.GetJD(y,m,d,hh,mm,ss,ns);
//
// Int_t days,secs,nsecs;
// Double_t date=421.1949327;
// t.Convert(date,days,secs,nsecs);
//
// days=875;
// secs=23;
// nsecs=9118483;
// date=t.Convert(days,secs,nsecs);
//
// Double_t mjdate=40563.823744;
// Double_t epoch=t.GetJE(mjdate,"mjd");
//
//--- Author: Nick van Eijndhoven 28-jan-2005 Utrecht University.
//- Modified: NvE $Date$ Utrecht University.
///////////////////////////////////////////////////////////////////////////

#include "AliTimestamp.h"
#include "Riostream.h"

ClassImp(AliTimestamp) // Class implementation to enable ROOT I/O
 
AliTimestamp::AliTimestamp() : TTimeStamp()
{
// Default constructor
// Creation of an AliTimestamp object and initialisation of parameters.
// All attributes are initialised to the current date/time as specified
// in the docs of TTimeStamp.

 FillJulian();
 fJps=0;
}
///////////////////////////////////////////////////////////////////////////
AliTimestamp::AliTimestamp(TTimeStamp& t) : TTimeStamp(t)
{
// Creation of an AliTimestamp object and initialisation of parameters.
// All attributes are initialised to the values of the input TTimeStamp.

 FillJulian();
 fJps=0;
}
///////////////////////////////////////////////////////////////////////////
AliTimestamp::~AliTimestamp()
{
// Destructor to delete dynamically allocated memory.
}
///////////////////////////////////////////////////////////////////////////
AliTimestamp::AliTimestamp(const AliTimestamp& t) : TTimeStamp(t)
{
// Copy constructor

 fMJD=t.fMJD;
 fJsec=t.fJsec;
 fJns=t.fJns;
 fJps=t.fJps;
 fCalcs=t.fCalcs;
 fCalcns=t.fCalcns;
}
///////////////////////////////////////////////////////////////////////////
void AliTimestamp::Date(Int_t mode)
{
// Print date/time info.
//
// mode = 1 ==> Only the TTimeStamp yy-mm-dd hh:mm:ss:ns info is printed
//        2 ==> Only the Julian parameter info is printed
//        3 ==> Both the TTimeStamp and Julian parameter info is printed
//
// The default is mode=3.
//
// Note : In case the (M/T)JD falls outside the TTimeStamp range,
//        the TTimeStamp info will not be printed.

  Int_t mjd,mjsec,mjns;
  GetMJD(mjd,mjsec,mjns);
 
 if ((mode==1 || mode==3) && mjd>=40587 && (mjd<65442 || (mjd==65442 && mjsec<8047)))
 {
  cout << " " << AsString() << endl;
 }
 if (mode==2 || mode==3)
 {
  Int_t jd,jsec,jns;
  GetJD(jd,jsec,jns);
  Int_t tjd,tjsec,tjns;
  GetTJD(tjd,tjsec,tjns);
  cout << " Julian Epoch : " << setprecision(25) << GetJE() << endl;
  cout << " JD : " << jd << " sec : " << jsec << " ns : " << jns << " ps : " << fJps
       << " Fractional : " << setprecision(25) << GetJD() << endl;
  cout << " MJD : " << mjd << "  sec : " << mjsec << " ns : " << mjns << " ps : " << fJps
       << " Fractional : " << setprecision(25) << GetMJD() << endl;
  cout << " TJD : " << tjd << "  sec : " << tjsec << " ns : " << tjns << " ps : " << fJps
       << " Fractional : " << setprecision(25) << GetTJD() << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
Double_t AliTimestamp::GetJD(Int_t y,Int_t m,Int_t d,Int_t hh,Int_t mm,Int_t ss,Int_t ns) const
{
// Provide the (fractional) Julian Date (JD) corresponding to the UT date
// and time in the Gregorian calendar as specified by the input arguments.
//
// The input arguments represent the following :
// y  : year in UT (e.g. 1952, 2003 etc...)
// m  : month in UT (1=jan  2=feb etc...)
// d  : day in UT (1-31)
// hh : elapsed hours in UT (0-23) 
// mm : elapsed minutes in UT (0-59)
// ss : elapsed seconds in UT (0-59)
// ns : remaining fractional elapsed second of UT in nanosecond
//
// This algorithm is valid for all AD dates in the Gregorian calendar
// following the recipe of R.W. Sinnott Sky & Telescope 82, (aug. 1991) 183.
// See also http://scienceworld.wolfram.com/astronomy/JulianDate.html
//
// In case of invalid input, a value of -1 is returned.
//
// Note :
// ------
// This memberfunction only provides the JD corresponding to the
// UT input arguments. It does NOT set the corresponding Julian parameters
// for the current AliTimestamp instance.
// As such the TTimeStamp limitations do NOT apply to this memberfunction.
// To set the Julian parameters for the current AliTimestamp instance,
// please use the corresponding SET() memberfunctions of either AliTimestamp
// or TTimeStamp. 

 if (y<0 || m<1 || m>12 || d<1 || d>31) return -1;
 if (hh<0 || hh>23 || mm<0 || mm>59 || ss<0 || ss>59 || ns<0 || ns>1e9) return -1;

 // The UT daytime in fractional hours
 Double_t ut=double(hh)+double(mm)/60.+(double(ss)+double(ns)*1.e-9)/3600.;

 Double_t JD=0;
 
 JD=367*y-int(7*(y+int((m+9)/12))/4)
    -int(3*(int((y+(m-9)/7)/100)+1)/4)
    +int(275*m/9)+d+1721028.5+ut/24.;

 return JD;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliTimestamp::GetMJD(Int_t y,Int_t m,Int_t d,Int_t hh,Int_t mm,Int_t ss,Int_t ns) const
{
// Provide the (fractional) Modified Julian Date corresponding to the UT
// date and time in the Gregorian calendar as specified by the input arguments.
//
// The input arguments represent the following :
// y  : year in UT (e.g. 1952, 2003 etc...)
// m  : month in UT (1=jan  2=feb etc...)
// d  : day in UT (1-31)
// hh : elapsed hours in UT (0-23) 
// mm : elapsed minutes in UT (0-59)
// ss : elapsed seconds in UT (0-59)
// ns : remaining fractional elapsed second of UT in nanosecond
//
// This algorithm is valid for all AD dates in the Gregorian calendar
// following the recipe of R.W. Sinnott Sky & Telescope 82, (aug. 1991) 183.
// See also http://scienceworld.wolfram.com/astronomy/JulianDate.html
//
// In case of invalid input, a value of -1 is returned.
//
// Note :
// ------
// This memberfunction only provides the MJD corresponding to the
// UT input arguments. It does NOT set the corresponding Julian parameters
// for the current AliTimestamp instance.
// As such the TTimeStamp limitations do NOT apply to this memberfunction.
// To set the Julian parameters for the current AliTimestamp instance,
// please use the corresponding SET() memberfunctions of either AliTimestamp
// or TTimeStamp.

 Double_t JD=GetJD(y,m,d,hh,mm,ss,ns);

 if (JD<0) return JD;

 Double_t MJD=JD-2400000.5;

 return MJD; 
} 
///////////////////////////////////////////////////////////////////////////
Double_t AliTimestamp::GetTJD(Int_t y,Int_t m,Int_t d,Int_t hh,Int_t mm,Int_t ss,Int_t ns) const
{
// Provide the (fractional) Truncated Julian Date corresponding to the UT
// date and time in the Gregorian calendar as specified by the input arguments.
//
// The input arguments represent the following :
// y  : year in UT (e.g. 1952, 2003 etc...)
// m  : month in UT (1=jan  2=feb etc...)
// d  : day in UT (1-31)
// hh : elapsed hours in UT (0-23) 
// mm : elapsed minutes in UT (0-59)
// ss : elapsed seconds in UT (0-59)
// ns : remaining fractional elapsed second of UT in nanosecond
//
// This algorithm is valid for all AD dates in the Gregorian calendar
// following the recipe of R.W. Sinnott Sky & Telescope 82, (aug. 1991) 183.
// See also http://scienceworld.wolfram.com/astronomy/JulianDate.html
//
// In case of invalid input, a value of -1 is returned.
//
// Note :
// ------
// This memberfunction only provides the TJD corresponding to the
// UT input arguments. It does NOT set the corresponding Julian parameters
// for the current AliTimestamp instance.
// As such the TTimeStamp limitations do NOT apply to this memberfunction.
// To set the Julian parameters for the current AliTimestamp instance,
// please use the corresponding SET() memberfunctions of either AliTimestamp
// or TTimeStamp.

 Double_t JD=GetJD(y,m,d,hh,mm,ss,ns);

 if (JD<0) return JD;

 Double_t TJD=JD-2440000.5;

 return TJD; 
} 
///////////////////////////////////////////////////////////////////////////
Double_t AliTimestamp::GetJE(Double_t date,TString mode) const
{
// Provide the Julian Epoch (JE) corresponding to the specified date.
// The argument "mode" indicates the type of the argument "date".
//
// Available modes are :
// mode = "jd"  ==> date represents the Julian Date
//      = "mjd" ==> date represents the Modified Julian Date
//      = "tjd" ==> date represents the Truncated Julian Date
//
// The default is mode="jd".
//
// In case of invalid input, a value of -99999 is returned.
//
// Note :
// ------
// This memberfunction only provides the JE corresponding to the
// input arguments. It does NOT set the corresponding Julian parameters
// for the current AliTimestamp instance.
// As such the TTimeStamp limitations do NOT apply to this memberfunction.
// To set the Julian parameters for the current AliTimestamp instance,
// please use the corresponding SET() memberfunctions of either AliTimestamp
// or TTimeStamp.

 if ((mode != "jd") && (mode != "mjd") && (mode != "tjd")) return -99999;

 Double_t jd=date;
 if (mode=="mjd") jd=date+2400000.5;
 if (mode=="tjd") jd=date+2440000.5;

 Double_t je=2000.+(jd-2451545.)/365.25;

 return je;
}
///////////////////////////////////////////////////////////////////////////
void AliTimestamp::Convert(Double_t date,Int_t& days,Int_t& secs,Int_t& ns) const
{
// Convert date as fractional day count into integer days, secs and ns.
//
// Note : Due to computer accuracy the ns value may become inaccurate.
//
// The arguments represent the following :
// date : The input date as fractional day count
// days : Number of elapsed days
// secs : Remaining number of elapsed seconds
// ns   : Remaining fractional elapsed second in nanoseconds
//
// Note :
// ------
// This memberfunction only converts the input date into the corresponding
// integer parameters. It does NOT set the corresponding Julian parameters
// for the current AliTimestamp instance.
// As such the TTimeStamp limitations do NOT apply to this memberfunction.
// To set the Julian parameters for the current AliTimestamp instance,
// please use the corresponding SET() memberfunctions of either AliTimestamp
// or TTimeStamp.
 
 days=int(date);
 date=date-double(days);
 Int_t daysecs=24*3600;
 date=date*double(daysecs);
 secs=int(date);
 date=date-double(secs);
 ns=int(date*1.e9);
}
///////////////////////////////////////////////////////////////////////////
Double_t AliTimestamp::Convert(Int_t days,Int_t secs,Int_t ns) const
{
// Convert date in integer days, secs and ns into fractional day count. 
//
// Note : Due to computer accuracy the ns precision may be lost.
//
// The input arguments represent the following :
// days : Number of elapsed days
// secs : Remaining number of elapsed seconds
// ns   : Remaining fractional elapsed second in nanoseconds
//
// Note :
// ------
// This memberfunction only converts the input integer parameters into the
// corresponding fractional day count. It does NOT set the corresponding
// Julian parameters for the current AliTimestamp instance.
// As such the TTimeStamp limitations do NOT apply to this memberfunction.
// To set the Julian parameters for the current AliTimestamp instance,
// please use the corresponding SET() memberfunctions of either AliTimestamp
// or TTimeStamp.

 Double_t frac=double(secs)+double(ns)*1.e-9;
 Int_t daysecs=24*3600;
 frac=frac/double(daysecs);
 Double_t date=double(days)+frac;
 return date;
}
///////////////////////////////////////////////////////////////////////////
void AliTimestamp::FillJulian()
{
// Calculation and setting of the Julian date/time parameters corresponding
// to the current TTimeStamp date/time parameters.

 UInt_t y,m,d,hh,mm,ss;

 GetDate(kTRUE,0,&y,&m,&d);
 GetTime(kTRUE,0,&hh,&mm,&ss);
 Int_t ns=GetNanoSec();

 Double_t mjd=GetMJD(y,m,d,hh,mm,ss,ns);

 fMJD=int(mjd);
 fJsec=GetSec()%(24*3600); // Daytime in elapsed seconds
 fJns=ns;                  // Remaining fractional elapsed second in nanoseconds

 // Store the TTimeStamp seconds and nanoseconds values
 // for which this Julian calculation was performed.
 fCalcs=GetSec();
 fCalcns=GetNanoSec();
}
///////////////////////////////////////////////////////////////////////////
void AliTimestamp::GetMJD(Int_t& mjd,Int_t& sec,Int_t& ns)
{
// Provide the Modified Julian Date (MJD) and time corresponding to the
// currently stored AliTimestamp date/time parameters.
//
// The returned arguments represent the following :
// mjd : The modified Julian date.
// sec : The number of seconds elapsed within the MJD.
// ns  : The remaining fractional number of seconds (in ns) elapsed within the MJD.

 if (fCalcs != GetSec() || fCalcns != GetNanoSec()) FillJulian();

 mjd=fMJD;
 sec=fJsec;
 ns=fJns;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliTimestamp::GetMJD()
{
// Provide the (fractional) Modified Julian Date (MJD) corresponding to the
// currently stored AliTimestamp date/time parameters.
//
// Due to computer accuracy the ns precision may be lost.
// It is advised to use the (mjd,sec,ns) getter instead.

 Int_t mjd=0;
 Int_t sec=0;
 Int_t ns=0;
 GetMJD(mjd,sec,ns);

 Double_t date=Convert(mjd,sec,ns);

 return date;
}
///////////////////////////////////////////////////////////////////////////
void AliTimestamp::GetTJD(Int_t& tjd,Int_t& sec, Int_t& ns)
{
// Provide the Truncated Julian Date (TJD) and time corresponding to the
// currently stored AliTimestamp date/time parameters.
//
// The returned arguments represent the following :
// tjd : The modified Julian date.
// sec : The number of seconds elapsed within the MJD.
// ns  : The remaining fractional number of seconds (in ns) elapsed within the MJD.

 Int_t mjd=0;
 GetMJD(mjd,sec,ns);

 tjd=mjd-40000;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliTimestamp::GetTJD()
{
// Provide the (fractional) Truncated Julian Date (TJD) corresponding to the
// currently stored AliTimestamp date/time parameters.
//
// Due to computer accuracy the ns precision may be lost.
// It is advised to use the (mjd,sec,ns) getter instead.

 Int_t tjd=0;
 Int_t sec=0;
 Int_t ns=0;
 GetTJD(tjd,sec,ns);

 Double_t date=Convert(tjd,sec,ns);

 return date;
}
///////////////////////////////////////////////////////////////////////////
void AliTimestamp::GetJD(Int_t& jd,Int_t& sec, Int_t& ns)
{
// Provide the Julian Date (JD) and time corresponding to the currently
// stored AliTimestamp date/time parameters.
//
// The returned arguments represent the following :
// jd  : The Julian date.
// sec : The number of seconds elapsed within the JD.
// ns  : The remaining fractional number of seconds (in ns) elapsed within the JD.

 Int_t mjd=0;
 GetMJD(mjd,sec,ns);

 jd=mjd+2400000;
 sec+=12*3600;
 if (sec >= 24*3600)
 {
  sec-=24*3600;
  jd+=1;
 }
}
///////////////////////////////////////////////////////////////////////////
Double_t AliTimestamp::GetJD()
{
// Provide the (fractional) Julian Date (JD) corresponding to the currently
// stored AliTimestamp date/time parameters.
//
// Due to computer accuracy the ns precision may be lost.
// It is advised to use the (jd,sec,ns) getter instead.

 Int_t jd=0;
 Int_t sec=0;
 Int_t ns=0;
 GetJD(jd,sec,ns);

 Double_t date=Convert(jd,sec,ns);

 return date;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliTimestamp::GetJE()
{
// Provide the Julian Epoch (JE) corresponding to the currently stored
// AliTimestamp date/time parameters.

 Double_t jd=GetJD();
 Double_t je=GetJE(jd);
 return je;
}
///////////////////////////////////////////////////////////////////////////
void AliTimestamp::SetMJD(Int_t mjd,Int_t sec,Int_t ns,Int_t ps)
{
// Set the Modified Julian Date (MJD) and time and update the TTimeStamp
// parameters accordingly (if possible).
//
// Note :
// ------
// The TTimeStamp EPOCH starts at 01-jan-1970 00:00:00 UTC
// which corresponds to the start of MJD=40587.
// Using the corresponding MJD of this EPOCH allows construction of
// the yy-mm-dd hh:mm:ss:ns TTimeStamp from a given input MJD and time.
// Obviously this TTimeStamp implementation would prevent usage of MJD values
// smaller than 40587.
// Furthermore, due to a limitation on the "seconds since the EPOCH start" count
// in TTimeStamp, the latest accessible date/time is 19-jan-2038 02:14:08 UTC.
// However, this AliTimestamp facility provides support for the full range
// of (M)JD values, but the setting of the corresponding TTimeStamp parameters
// is restricted to the values allowed by the TTimeStamp implementation.
// For these earlier/later MJD values, the standard TTimeStamp parameters will
// be set corresponding to the start of the TTimeStamp EPOCH.  
// This implies that for these earlier/later MJD values the TTimeStamp parameters
// do not match the Julian parameters of AliTimestamp.  
//
// The input arguments represent the following :
// mjd : The modified Julian date.
// sec : The number of seconds elapsed within the MJD.
// ns  : The remaining fractional number of seconds (in ns) elapsed within the MJD.
// ps  : The remaining fractional number of nanoseconds (in ps) elapsed within the MJD.
//
// Note : ps=0 is the default value.

 if (sec<0 || sec>=24*3600 || ns<0 || ns>=1e9 || ps<0 || ps>=1000)
 {
  cout << " *AliTimestamp::SetMJD* Invalid input."
       << " sec : " << sec << " ns : " << ns << endl; 
  return;
 }

 fMJD=mjd;
 fJsec=sec;
 fJns=ns;
 fJps=ps;

 Int_t epoch=40587; // MJD of the start of the epoch
 Int_t limit=65442; // MJD of the latest possible TTimeStamp date/time
 
 Int_t date,time;
 if (mjd<epoch || (mjd>=limit && sec>=8047))
 {
  Set(0,kFALSE,0,kFALSE);
  date=GetDate();
  time=GetTime();
  Set(date,time,0,kTRUE,0);
 }
 else
 {
  // The elapsed time since start of EPOCH
  Int_t days=mjd-epoch;
  UInt_t secs=days*24*3600;
  secs+=sec;
  Set(secs,kFALSE,0,kFALSE);
  date=GetDate();
  time=GetTime();
  Set(date,time,ns,kTRUE,0);
 }

 // Denote that the Julian and TTimeStamp parameters are synchronised,
 // even in the case the MJD falls outside the TTimeStamp validity range.
 // The latter still allows retrieval of Julian parameters for these
 // earlier times.
 fCalcs=GetSec();
 fCalcns=GetNanoSec();
}
///////////////////////////////////////////////////////////////////////////
void AliTimestamp::SetMJD(Double_t mjd)
{
// Set the Modified Julian Date (MJD) and time and update the TTimeStamp
// parameters accordingly (if possible).
//
// Note :
// ------
// The TTimeStamp EPOCH starts at 01-jan-1970 00:00:00 UTC
// which corresponds to the start of MJD=40587.
// Using the corresponding MJD of this EPOCH allows construction of
// the yy-mm-dd hh:mm:ss:ns TTimeStamp from a given input MJD and time.
// Obviously this TTimeStamp implementation would prevent usage of MJD values
// smaller than 40587.
// Furthermore, due to a limitation on the "seconds since the EPOCH start" count
// in TTimeStamp, the latest accessible date/time is 19-jan-2038 02:14:08 UTC.
// However, this AliTimestamp facility provides support for the full range
// of (M)JD values, but the setting of the corresponding TTimeStamp parameters
// is restricted to the values allowed by the TTimeStamp implementation.
// For these earlier/later MJD values, the standard TTimeStamp parameters will
// be set corresponding to the start of the TTimeStamp EPOCH.  
// This implies that for these earlier/later MJD values the TTimeStamp parameters
// do not match the Julian parameters of AliTimestamp.  
//
// Due to computer accuracy the ns precision may be lost.
// It is advised to use the (mjd,sec,ns) setting instead.
//
// The input argument represents the following :
// mjd : The modified Julian date as fractional day count.

 Int_t days=0;
 Int_t secs=0;
 Int_t ns=0;
 Convert(mjd,days,secs,ns);
 SetMJD(days,secs,ns);
}
///////////////////////////////////////////////////////////////////////////
void AliTimestamp::SetJD(Int_t jd,Int_t sec,Int_t ns,Int_t ps)
{
// Set the Julian Date (JD) and time and update the TTimeStamp
// parameters accordingly (if possible).
//
// Note :
// ------
// The TTimeStamp EPOCH starts at 01-jan-1970 00:00:00 UTC
// which corresponds to JD=2440587.5 or the start of MJD=40587.
// Using the corresponding MJD of this EPOCH allows construction of
// the yy-mm-dd hh:mm:ss:ns TTimeStamp from a given input MJD and time.
// Obviously this TTimeStamp implementation would prevent usage of values
// smaller than JD=2440587.5.
// Furthermore, due to a limitation on the "seconds since the EPOCH start" count
// in TTimeStamp, the latest accessible date/time is 19-jan-2038 02:14:08 UTC.
// However, this AliTimestamp facility provides support for the full range
// of (M)JD values, but the setting of the corresponding TTimeStamp parameters
// is restricted to the values allowed by the TTimeStamp implementation.
// For these earlier/later JD values, the standard TTimeStamp parameters will
// be set corresponding to the start of the TTimeStamp EPOCH.  
// This implies that for these earlier/later (M)JD values the TTimeStamp parameters
// do not match the Julian parameters of AliTimestamp.  
//
// The input arguments represent the following :
// jd  : The Julian date.
// sec : The number of seconds elapsed within the JD.
// ns  : The remaining fractional number of seconds (in ns) elapsed within the JD.
// ps  : The remaining fractional number of nanoseconds (in ps) elapsed within the JD.
//
// Note : ps=0 is the default value.

 Int_t mjd=jd-2400000;
 sec-=12*3600;
 if (sec<0)
 {
  sec+=24*3600;
  mjd-=1;
 }

 SetMJD(mjd,sec,ns,ps);
}
///////////////////////////////////////////////////////////////////////////
void AliTimestamp::SetJD(Double_t jd)
{
// Set the Julian Date (JD) and time and update the TTimeStamp
// parameters accordingly (if possible).
//
// Note :
// ------
// The TTimeStamp EPOCH starts at 01-jan-1970 00:00:00 UTC
// which corresponds to JD=2440587.5 or the start of MJD=40587.
// Using the corresponding MJD of this EPOCH allows construction of
// the yy-mm-dd hh:mm:ss:ns TTimeStamp from a given input MJD and time.
// Obviously this TTimeStamp implementation would prevent usage of values
// smaller than JD=2440587.5.
// Furthermore, due to a limitation on the "seconds since the EPOCH start" count
// in TTimeStamp, the latest accessible date/time is 19-jan-2038 02:14:08 UTC.
// However, this AliTimestamp facility provides support for the full range
// of (M)JD values, but the setting of the corresponding TTimeStamp parameters
// is restricted to the values allowed by the TTimeStamp implementation.
// For these earlier/later JD values, the standard TTimeStamp parameters will
// be set corresponding to the start of the TTimeStamp EPOCH.  
// This implies that for these earlier/later (M)JD values the TTimeStamp parameters
// do not match the Julian parameters of AliTimestamp.  
//
// Due to computer accuracy the ns precision may be lost.
// It is advised to use the (jd,sec,ns) setting instead.
//
// The input argument represents the following :
// jd : The Julian date as fractional day count.

 Int_t days=0;
 Int_t secs=0;
 Int_t ns=0;
 Convert(jd,days,secs,ns);

 SetJD(days,secs,ns);
}
///////////////////////////////////////////////////////////////////////////
void AliTimestamp::SetTJD(Int_t tjd,Int_t sec,Int_t ns,Int_t ps)
{
// Set the Truncated Julian Date (TJD) and time and update the TTimeStamp
// parameters accordingly (if possible).
//
// Note :
// ------
// The TTimeStamp EPOCH starts at 01-jan-1970 00:00:00 UTC
// which corresponds to JD=2440587.5 or the start of TJD=587.
// Using the corresponding MJD of this EPOCH allows construction of
// the yy-mm-dd hh:mm:ss:ns TTimeStamp from a given input MJD and time.
// Obviously this TTimeStamp implementation would prevent usage of values
// smaller than TJD=587.
// Furthermore, due to a limitation on the "seconds since the EPOCH start" count
// in TTimeStamp, the latest accessible date/time is 19-jan-2038 02:14:08 UTC.
// However, this AliTimestamp facility provides support for the full range
// of (T)JD values, but the setting of the corresponding TTimeStamp parameters
// is restricted to the values allowed by the TTimeStamp implementation.
// For these earlier/later JD values, the standard TTimeStamp parameters will
// be set corresponding to the start of the TTimeStamp EPOCH.  
// This implies that for these earlier/later (T)JD values the TTimeStamp parameters
// do not match the Julian parameters of AliTimestamp.  
//
// The input arguments represent the following :
// tjd : The Truncated Julian date.
// sec : The number of seconds elapsed within the JD.
// ns  : The remaining fractional number of seconds (in ns) elapsed within the JD.
// ps  : The remaining fractional number of nanoseconds (in ps) elapsed within the JD.
//
// Note : ps=0 is the default value.

 Int_t mjd=tjd+40000;

 SetMJD(mjd,sec,ns,ps);
}
///////////////////////////////////////////////////////////////////////////
void AliTimestamp::SetTJD(Double_t tjd)
{
// Set the Truncated Julian Date (TJD) and time and update the TTimeStamp
// parameters accordingly (if possible).
//
// Note :
// ------
// The TTimeStamp EPOCH starts at 01-jan-1970 00:00:00 UTC
// which corresponds to JD=2440587.5 or the start of TJD=587.
// Using the corresponding MJD of this EPOCH allows construction of
// the yy-mm-dd hh:mm:ss:ns TTimeStamp from a given input MJD and time.
// Obviously this TTimeStamp implementation would prevent usage of values
// smaller than TJD=587.
// Furthermore, due to a limitation on the "seconds since the EPOCH start" count
// in TTimeStamp, the latest accessible date/time is 19-jan-2038 02:14:08 UTC.
// However, this AliTimestamp facility provides support for the full range
// of (T)JD values, but the setting of the corresponding TTimeStamp parameters
// is restricted to the values allowed by the TTimeStamp implementation.
// For these earlier/later JD values, the standard TTimeStamp parameters will
// be set corresponding to the start of the TTimeStamp EPOCH.  
// This implies that for these earlier/later (T)JD values the TTimeStamp parameters
// do not match the Julian parameters of AliTimestamp.  
//
// Due to computer accuracy the ns precision may be lost.
// It is advised to use the (jd,sec,ns) setting instead.
//
// The input argument represents the following :
// tjd : The Truncated Julian date as fractional day count.

 Int_t days=0;
 Int_t secs=0;
 Int_t ns=0;
 Convert(tjd,days,secs,ns);

 SetTJD(days,secs,ns);
}
///////////////////////////////////////////////////////////////////////////
void AliTimestamp::SetNs(Int_t ns)
{
// Set the remaining fractional number of seconds in nanosecond precision.
// Notes :
// -------
// 1) The allowed range for the argument "ns" is [0,99999999].
//    Outside that range no action is performed.
// 2) The ns fraction can also be entered directly via SetMJD() etc...
// 3) For additional accuracy see SetPs().

 if (ns>=0 && ns<=99999999) fJns=ns; 
}
///////////////////////////////////////////////////////////////////////////
Int_t AliTimestamp::GetNs() const
{
// Provide the remaining fractional number of seconds in nanosecond precision.
// This function allows trigger/timing analysis for (astro)particle physics
// experiments.
// Note : For additional accuracy see also GetPs().

 return fJns; 
}
///////////////////////////////////////////////////////////////////////////
void AliTimestamp::SetPs(Int_t ps)
{
// Set the remaining fractional number of nanoseconds in picoseconds.
// Notes :
// -------
// 1) The allowed range for the argument "ps" is [0,999].
//    Outside that range no action is performed.
// 2) The ps fraction can also be entered directly via SetMJD() etc...

 if (ps>=0 && ps<=999) fJps=ps; 
}
///////////////////////////////////////////////////////////////////////////
Int_t AliTimestamp::GetPs() const
{
// Provide remaining fractional number of nanoseconds in picoseconds.
// This function allows time of flight analysis for particle physics
// experiments.

 return fJps; 
}
///////////////////////////////////////////////////////////////////////////
void AliTimestamp::Add(Int_t d,Int_t s,Int_t ns,Int_t ps)
{
// Add (or subtract) a certain time difference to the current timestamp.
// Subtraction can be achieved by entering negative values as input arguments.
//
// The time difference is entered via the following input arguments :
//
// d  : elapsed number of days
// s  : (remaining) elapsed number of seconds
// ns : (remaining) elapsed number of nanoseconds
// ps : (remaining) elapsed number of picoseconds
//
// The specified d, s, ns and ps values will be used in an additive
// way to determine the time difference.
// So, specification of d=1, s=100, ns=0, ps=0 will result in the
// same time difference addition as d=0, s=24*3600+100, ns=0, ps=0.
// However, by making use of the latter the user should take care
// of possible integer overflow problems in the input arguments,
// which obviously will provide incorrect results. 
//
// Note : ps=0 is the default value.

 Int_t days=0;
 Int_t secs=0;
 Int_t nsec=0;
 // Use Get functions to ensure updated Julian parameters. 
 GetMJD(days,secs,nsec);
 Int_t psec=GetPs();

 psec+=ps%1000;
 nsec+=ps/1000;
 while (psec<0)
 {
  nsec-=1;
  psec+=1000;
 }
 while (psec>999)
 {
  nsec+=1;
  psec-=1000;
 }

 nsec+=ns%1000000000;
 secs+=ns/1000000000;
 while (nsec<0)
 {
  secs-=1;
  nsec+=1000000000;
 }
 while (nsec>999999999)
 {
  secs+=1;
  nsec-=1000000000;
 }

 secs+=s%(24*3600);
 days+=s/(24*3600);
 while (secs<0)
 {
  days-=1;
  secs+=24*3600;
 }
 while (secs>=24*3600)
 {
  days+=1;
  secs-=24*3600;
 }

 days+=d;

 SetMJD(days,secs,nsec,psec);
}
///////////////////////////////////////////////////////////////////////////
Int_t AliTimestamp::GetDifference(AliTimestamp* t,Int_t& d,Int_t& s,Int_t& ns,Int_t& ps)
{
// Provide the time difference w.r.t the AliTimestamp specified on the input.
// This memberfunction supports both very small (i.e. time of flight analysis
// for particle physics experiments) and very long (i.e. investigation of
// astrophysical phenomena) timescales.
//
// The time difference is returned via the following output arguments :
// d  : elapsed number of days
// s  : remaining elapsed number of seconds
// ns : remaining elapsed number of nanoseconds
// ps : remaining elapsed number of picoseconds
//
// Note :
// ------
// The calculated time difference is the absolute value of the time interval.
// This implies that the values of d, s, ns and ps are always positive or zero.
//
// The integer return argument indicates whether the AliTimestamp specified
// on the input argument occurred earlier (-1), simultaneously (0) or later (1).

 if (!t) return 0;

 // Ensure updated Julian parameters for this AliTimestamp instance 
 if (fCalcs != GetSec() || fCalcns != GetNanoSec()) FillJulian();

 // Use Get functions to ensure updated Julian parameters. 
 t->GetMJD(d,s,ns);
 ps=t->GetPs();

 d-=fMJD;
 s-=fJsec;
 ns-=fJns;
 ps-=fJps;

 if (!d && !s && !ns && !ps) return 0;

 Int_t sign=0;

 if (d>0) sign=1;
 if (d<0) sign=-1;

 if (!sign && s>0) sign=1;
 if (!sign && s<0) sign=-1;

 if (!sign && ns>0) sign=1; 
 if (!sign && ns<0) sign=-1;

 if (!sign && ps>0) sign=1; 
 if (!sign && ps<0) sign=-1;

 // In case the input stamp was earlier, take the reverse difference
 // to simplify the algebra.
 if (sign<0)
 {
  d=-d;
  s=-s;
  ns=-ns;
  ps=-ps;
 }

 // Here we always have a positive time difference
 // and can now unambiguously correct for other negative values.
 if (ps<0)
 {
  ns-=1;
  ps+=1000;
 }

 if (ns<0)
 {
  s-=1;
  ns+=1000000000;
 }

 if (s<0)
 {
  d-=1;
  s+=24*3600;
 }

 return sign;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliTimestamp::GetDifference(AliTimestamp& t,Int_t& d,Int_t& s,Int_t& ns,Int_t& ps)
{
// Provide the time difference w.r.t the AliTimestamp specified on the input.
// This memberfunction supports both very small (i.e. time of flight analysis
// for particle physics experiments) and very long (i.e. investigation of
// astrophysical phenomena) timescales.
//
// The time difference is returned via the following output arguments :
// d  : elapsed number of days
// s  : remaining elapsed number of seconds
// ns : remaining elapsed number of nanoseconds
// ps : remaining elapsed number of picoseconds
//
// Note :
// ------
// The calculated time difference is the absolute value of the time interval.
// This implies that the values of d, s, ns and ps are always positive or zero.
//
// The integer return argument indicates whether the AliTimestamp specified
// on the input argument occurred earlier (-1), simultaneously (0) or later (1).

 return GetDifference(&t,d,s,ns,ps);
}
///////////////////////////////////////////////////////////////////////////
Double_t AliTimestamp::GetDifference(AliTimestamp* t,TString u,Int_t mode)
{
// Provide the time difference w.r.t the AliTimestamp specified on the input
// argument in the units as specified by the TString argument.
// A positive return value means that the AliTimestamp specified on the input
// argument occurred later, whereas a negative return value indicates an
// earlier occurence. 
//  
// The units may be specified as :
// u = "d"  ==> Time difference returned as (fractional) day count
//     "s"  ==> Time difference returned as (fractional) second count
//     "ns" ==> Time difference returned as (fractional) nanosecond count
//     "ps" ==> Time difference returned as picosecond count
//
// It may be clear that for a time difference of several days, the picosecond
// and even the nanosecond accuracy may be lost.
// To cope with this, the "mode" argument has been introduced to allow 
// timestamp comparison on only the specified units.
//
// The following operation modes are supported :
// mode = 1 : Full time difference is returned in specified units
//        2 : Time difference is returned in specified units by
//            neglecting the elapsed time for the larger units than the
//            ones specified.
//        3 : Time difference is returned in specified units by only
//            comparing the timestamps on the level of the specified units.
//
// Example :
// ---------
// AliTimestamp t1; // Corresponding to days=3, secs=501, ns=31, ps=7 
// AliTimestamp t2; // Corresponding to days=5, secs=535, ns=12, ps=15
//
// The statement : Double_t val=t1.GetDifference(t2,....)
// would return the following values :
// val=(2*24*3600)+34-(19*1e-9)+(8*1e-12) for u="s" and mode=1
// val=34-(19*1e-9)+(8*1e-12)             for u="s" and mode=2
// val=34                                 for u="s" and mode=3
// val=-19                                for u="ns" and mode=3
//
// The default is mode=1.

 if (!t || mode<1 || mode>3) return 0;

 Double_t dt=0;

 // Ensure updated Julian parameters for this AliTimestamp instance 
 if (fCalcs != GetSec() || fCalcns != GetNanoSec()) FillJulian();

 Int_t dd=0;
 Int_t ds=0;
 Int_t dns=0;
 Int_t dps=0;

 // Use Get functions to ensure updated Julian parameters. 
 t->GetMJD(dd,ds,dns);
 dps=t->GetPs();

 dd-=fMJD;
 ds-=fJsec;
 dns-=fJns;
 dps-=fJps;

 // Time difference for the specified units only
 if (mode==3)
 {
  if (u=="d") dt=dd;
  if (u=="s") dt=ds;
  if (u=="ns") dt=dns;
  if (u=="ps") dt=dps;
  return dt;
 }

 // Suppress elapsed time for the larger units than specified
 if (mode==2)
 {
  if (u=="s") dd=0;
  if (u=="ns")
  {
   dd=0;
   ds=0;
  }
  if (u=="ps")
  {
   dd=0;
   ds=0;
   dns=0;
  }
 }

 // Compute the time difference as requested 
 if (u=="s" || u=="d")
 {
  // The time difference in (fractional) seconds
  dt=double(dd*24*3600+ds)+(double(dns)*1e-9)+(double(dps)*1e-12);
  if (u=="d") dt=dt/double(24*3600);
 }
 if (u=="ns") dt=(double(dd*24*3600+ds)*1e9)+double(dns)+(double(dps)*1e-3);
 if (u=="ps") dt=(double(dd*24*3600+ds)*1e12)+(double(dns)*1e3)+double(dps);

 return dt;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliTimestamp::GetDifference(AliTimestamp& t,TString u,Int_t mode)
{
// Provide the time difference w.r.t the AliTimestamp specified on the input
// argument in the units as specified by the TString argument.
// A positive return value means that the AliTimestamp specified on the input
// argument occurred later, whereas a negative return value indicates an
// earlier occurence. 
//  
// The units may be specified as :
// u = "d"  ==> Time difference returned as (fractional) day count
//     "s"  ==> Time difference returned as (fractional) second count
//     "ns" ==> Time difference returned as (fractional) nanosecond count
//     "ps" ==> Time difference returned as picosecond count
//
// It may be clear that for a time difference of several days, the picosecond
// and even the nanosecond accuracy may be lost.
// To cope with this, the "mode" argument has been introduced to allow 
// timestamp comparison on only the specified units.
//
// The following operation modes are supported :
// mode = 1 : Full time difference is returned in specified units
//        2 : Time difference is returned in specified units by
//            neglecting the elapsed time for the larger units than the
//            ones specified.
//        3 : Time difference is returned in specified units by only
//            comparing the timestamps on the level of the specified units.
//
// Example :
// ---------
// AliTimestamp t1; // Corresponding to days=3, secs=501, ns=31, ps=7 
// AliTimestamp t2; // Corresponding to days=5, secs=535, ns=12, ps=15
//
// The statement : Double_t val=t1.GetDifference(t2,....)
// would return the following values :
// val=(2*24*3600)+34-(19*1e-9)+(8*1e-12) for u="s" and mode=1
// val=34-(19*1e-9)+(8*1e-12)             for u="s" and mode=2
// val=34                                 for u="s" and mode=3
// val=-19                                for u="ns" and mode=3
//
// The default is mode=1.

 return GetDifference(&t,u,mode);
}
///////////////////////////////////////////////////////////////////////////
void AliTimestamp::SetUT(Int_t y,Int_t m,Int_t d,Int_t hh,Int_t mm,Int_t ss,Int_t ns,Int_t ps)
{
// Set the AliTimestamp parameters corresponding to the UT date and time
// in the Gregorian calendar as specified by the input arguments.
// This facility is exact upto picosecond precision and as such is
// for scientific observations preferable above the corresponding
// Set function(s) of TTimestamp.
// The latter has a random spread in the sub-second part, which
// might be of use in generating distinguishable timestamps while
// still keeping second precision.
//
// The input arguments represent the following :
// y  : year in UT (e.g. 1952, 2003 etc...)
// m  : month in UT (1=jan  2=feb etc...)
// d  : day in UT (1-31)
// hh : elapsed hours in UT (0-23) 
// mm : elapsed minutes in UT (0-59)
// ss : elapsed seconds in UT (0-59)
// ns : remaining fractional elapsed second of UT in nanosecond
// ps : remaining fractional elapsed nanosecond of UT in picosecond
//
// Note : ns=0 and ps=0 are the default values.
//
// This facility first determines the elapsed days, seconds etc...
// since the beginning of the specified UT year on bais of the
// input arguments. Subsequently it invokes the SetUT memberfunction
// for the elapsed timespan.
// As such this facility is valid for all AD dates in the Gregorian
// calendar with picosecond precision.

 Int_t day=GetDayOfYear(d,m,y);
 Int_t secs=hh*3600+mm*60+ss;
 SetUT(y,day-1,secs,ns,ps);
}
///////////////////////////////////////////////////////////////////////////
void AliTimestamp::SetUT(Int_t y,Int_t d,Int_t s,Int_t ns,Int_t ps)
{
// Set the AliTimestamp parameters corresponding to the specified elapsed
// timespan since the beginning of the new UT year.
// This facility is exact upto picosecond precision and as such is
// for scientific observations preferable above the corresponding
// Set function(s) of TTimestamp.
// The latter has a random spread in the sub-second part, which
// might be of use in generating distinguishable timestamps while
// still keeping second precision.
//
// The UT year and elapsed time span is entered via the following input arguments :
//
// y  : year in UT (e.g. 1952, 2003 etc...)
// d  : elapsed number of days 
// s  : (remaining) elapsed number of seconds
// ns : (remaining) elapsed number of nanoseconds
// ps : (remaining) elapsed number of picoseconds
//
// The specified d, s, ns and ps values will be used in an additive
// way to determine the elapsed timespan.
// So, specification of d=1, s=100, ns=0, ps=0 will result in the
// same elapsed time span as d=0, s=24*3600+100, ns=0, ps=0.
// However, by making use of the latter the user should take care
// of possible integer overflow problems in the input arguments,
// which obviously will provide incorrect results. 
//
// Note : ns=0 and ps=0 are the default values.
//
// This facility first sets the (M)JD corresponding to the start (01-jan 00:00:00)
// of the specified UT year following the recipe of R.W. Sinnott
// Sky & Telescope 82, (aug. 1991) 183.
// Subsequently the day and (sub)second parts are added to the AliTimestamp.
// As such this facility is valid for all AD dates in the Gregorian calendar.

 Double_t jd=GetJD(y,1,1,0,0,0,0);
 SetJD(jd);

 Int_t mjd,sec,nsec;
 GetMJD(mjd,sec,nsec);
 SetMJD(mjd,0,0,0);
 Add(d,s,ns,ps);
}
///////////////////////////////////////////////////////////////////////////
