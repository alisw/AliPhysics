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
// This AliTimestamp facility allows for nanosecond precision.
// However, when the fractional JD, MJD and TJD counts are used instead
// of the integer (days,sec,ns) specification, the nanosecond precision
// may be lost due to computer accuracy w.r.t. floating point operations.
//
// The TTimeStamp EPOCH starts at 01-jan-1970 00:00:00 UTC
// which corresponds to JD=2440587.5 or the start of MJD=40587 or TJD=587.
// Using the corresponding MJD of this EPOCH allows construction of
// the yy-mm-dd hh:mm:ss:ns TTimeStamp from a given input (M/T)JD and time.
// Obviously this TTimeStamp implementation would prevent usage of values
// smaller than JD=2440587.5 or MJD=40587 or TJD=587.
// However, this AliTimestamp facility provides support for the full range
// of (M/T)JD values, but the setting of the corresponding TTimeStamp parameters
// is restricted to the values allowed by the TTimeStamp implementation.
// For these earlier (M/T)JD values, the standard TTimeStamp parameters will
// be set corresponding to the start of the TTimeStamp EPOCH.
// This implies that for these earlier (M/T)JD values the TTimeStamp parameters
// do not match the Julian parameters of AliTimestamp.
// As such the standard TTimeStamp parameters do not appear on the print output
// when invoking the Date() memberfunction for these earlier (M/T)JD values.  
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
}
///////////////////////////////////////////////////////////////////////////
AliTimestamp::AliTimestamp(TTimeStamp& t) : TTimeStamp(t)
{
// Creation of an AliTimestamp object and initialisation of parameters.
// All attributes are initialised to the values of the input TTimeStamp.

 FillJulian();
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

 if ((mode==1 || mode==3) && mjd>=40587) cout << " " << AsString() << endl;
 if (mode==2 || mode==3)
 {
  Int_t jd,jsec,jns;
  GetJD(jd,jsec,jns);
  Int_t tjd,tjsec,tjns;
  GetTJD(tjd,tjsec,tjns);
  cout << " Julian Epoch : " << setprecision(25) << GetJE() << endl;
  cout << " JD : " << jd << " sec : " << jsec << " ns : " << jns
       << " Fractional : " << setprecision(25) << GetJD() << endl;
  cout << " MJD : " << mjd << "  sec : " << mjsec << " ns : " << mjns
       << " Fractional : " << setprecision(25) << GetMJD() << endl;
  cout << " TJD : " << tjd << "  sec : " << tjsec << " ns : " << tjns
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
void AliTimestamp::GetMJD(Int_t& mjd,Int_t& sec, Int_t& ns)
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
void AliTimestamp::SetMJD(Int_t mjd,Int_t sec,Int_t ns)
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
// However, this AliTimestamp facility provides support for the full range
// of (M)JD values, but the setting of the corresponding TTimeStamp parameters
// is restricted to the values allowed by the TTimeStamp implementation.
// For these earlier MJD values, the standard TTimeStamp parameters will
// be set corresponding to the start of the TTimeStamp EPOCH.  
// This implies that for these earlier MJD values the TTimeStamp parameters
// do not match the Julian parameters of AliTimestamp.  
//
// The input arguments represent the following :
// mjd : The modified Julian date.
// sec : The number of seconds elapsed within the MJD.
// ns  : The remaining fractional number of seconds (in ns) elapsed within the MJD.

 if (sec<0 || ns<0)
 {
  cout << " *AliTimestamp::SetMJD* Invalid input."
       << " sec : " << sec << " ns : " << ns << endl; 
  return;
 }

 fMJD=mjd;
 fJsec=sec;
 fJns=ns;

 Int_t epoch=40587;
 
 if (mjd<epoch)
 {
  Set(0,kFALSE,0,kFALSE);
 }
 else
 {
  // The elapsed time since start of EPOCH
  Int_t days=mjd-epoch;
  UInt_t secs=days*24*3600;
  secs+=sec;
  Set(secs,kFALSE,0,kFALSE);
  Int_t date=GetDate();
  Int_t time=GetTime();
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
// However, this AliTimestamp facility provides support for the full range
// of (M)JD values, but the setting of the corresponding TTimeStamp parameters
// is restricted to the values allowed by the TTimeStamp implementation.
// For these earlier MJD values, the standard TTimeStamp parameters will
// be set corresponding to the start of the TTimeStamp EPOCH.  
// This implies that for these earlier MJD values the TTimeStamp parameters
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
void AliTimestamp::SetJD(Int_t jd,Int_t sec,Int_t ns)
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
// However, this AliTimestamp facility provides support for the full range
// of (M)JD values, but the setting of the corresponding TTimeStamp parameters
// is restricted to the values allowed by the TTimeStamp implementation.
// For these earlier JD values, the standard TTimeStamp parameters will
// be set corresponding to the start of the TTimeStamp EPOCH.  
// This implies that for these earlier (M)JD values the TTimeStamp parameters
// do not match the Julian parameters of AliTimestamp.  
//
// The input arguments represent the following :
// jd  : The Julian date.
// sec : The number of seconds elapsed within the JD.
// ns  : The remaining fractional number of seconds (in ns) elapsed within the JD.

 Int_t mjd=jd-2400000;
 sec-=12*3600;
 if (sec<0)
 {
  sec+=24*3600;
  mjd-=1;
 }

 SetMJD(mjd,sec,ns);
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
// However, this AliTimestamp facility provides support for the full range
// of (M)JD values, but the setting of the corresponding TTimeStamp parameters
// is restricted to the values allowed by the TTimeStamp implementation.
// For these earlier JD values, the standard TTimeStamp parameters will
// be set corresponding to the start of the TTimeStamp EPOCH.  
// This implies that for these earlier (M)JD values the TTimeStamp parameters
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
void AliTimestamp::SetTJD(Int_t tjd,Int_t sec,Int_t ns)
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
// However, this AliTimestamp facility provides support for the full range
// of (T)JD values, but the setting of the corresponding TTimeStamp parameters
// is restricted to the values allowed by the TTimeStamp implementation.
// For these earlier JD values, the standard TTimeStamp parameters will
// be set corresponding to the start of the TTimeStamp EPOCH.  
// This implies that for these earlier (T)JD values the TTimeStamp parameters
// do not match the Julian parameters of AliTimestamp.  
//
// The input arguments represent the following :
// tjd : The Truncated Julian date.
// sec : The number of seconds elapsed within the JD.
// ns  : The remaining fractional number of seconds (in ns) elapsed within the JD.

 Int_t mjd=tjd+40000;

 SetMJD(mjd,sec,ns);
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
// However, this AliTimestamp facility provides support for the full range
// of (T)JD values, but the setting of the corresponding TTimeStamp parameters
// is restricted to the values allowed by the TTimeStamp implementation.
// For these earlier JD values, the standard TTimeStamp parameters will
// be set corresponding to the start of the TTimeStamp EPOCH.  
// This implies that for these earlier (T)JD values the TTimeStamp parameters
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
