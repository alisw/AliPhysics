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
// Class AliAstrolab
// Virtual lab to correlate measurements with astrophysical phenomena.
//
// This class is derived from TTask, but the only reason for this
// is to enable this class to serve as a base class for other TTask
// derived classes (e.g. AliEventSelector) without the need for
// multiple virtual inheritance.
// So, this AliAstrolab class itself does not provide any TTask
// related functionality.
//
// The lab can be given a terrestrial location via the usual longitude
// and latitude specifications.
// Since this class is derived from AliTimestamp, a lab can also be
// given a specific timestamp. Together with the terrestrial location
// this provides access to local (sidereal) times etc...
// In addition to the usual astronomical reference frames, a local
// lab reference frame can also be specified. Together with the lab's
// timestamp this uniquely defines all the coordinate transformations
// between the various reference frames.
// These lab characteristics allow a space-time correlation of lab
// observations with external (astrophysical) phenomena.
//
// Observations are entered as generic signals containing a position,
// reference frame specification and a timestamp.
// These observations can then be analysed in various reference frames
// via the available GET functions.
//
// Various external (astrophysical) phenomena may be entered as
// so-called reference signals.
// This class provides facilities (e.g. MatchRefSignal) to check
// correlations of the stored measurement with these reference signals.
// 
// Coding example :
// ----------------
// gSystem->Load("ralice");
//
// AliAstrolab lab("IceCube","The South Pole Neutrino Observatory");
// lab.SetLabPosition(0,-90,"deg"); // South Pole
//
// lab.SetLocalFrame(90,180,90,270,0,0); // Local frame has X-axis to the North
//
// lab.Data(1,"dms"); // Print laboratory parameters
//
// // Enter some observed event to be investigated
// AliTimestamp ts;
// ts.SetUT(1989,7,30,8,14,23,738504,0);
// Float_t vec[3]={1,23.8,118.65};
// Ali3Vector r;
// r.SetVector(vec,"sph","deg");
// lab.SetSignal(&r,"loc","M",&ts,0,"Event10372");
//
// // Enter some reference signals
// Float_t alpha=194818.0;
// Float_t delta=84400.;
// lab.SetSignal(alpha,delta,"B",1950,"M",-1,"Altair");
// alpha=124900.0;
// delta=272400.;
// lab.SetSignal(alpha,delta,"B",1950,"M",-1,"NGP");
// alpha=64508.917;
// delta=-164258.02;
// lab.SetSignal(alpha,delta,"J",2000,"M",-1,"Sirius");
// alpha=23149.08;
// delta=891550.8;
// lab.SetSignal(alpha,delta,"J",2000,"M",-1,"Polaris");
// alpha=43600.;
// delta=163100.;
// lab.SetSignal(alpha,delta,"J",2000,"M",-1,"Aldebaran");
// Float_t l=327.531;
// Float_t b=-35.8903;
// Float_t pos[3]={1,90.-b,l};
// r.SetVector(pos,"sph","deg");
// lab.SetUT(1989,7,30,8,14,16,0,0);
// lab.SetSignal(&r,"gal","M",0,-1,"GRB890730");
//
// // List all stored objects
// lab.ListSignals("equ","M",5);
//
// // Determine minimal space and time differences with reference signals
// Double_t da,dt;
// Int_t ia,it;
// da=lab.GetDifference(0,"deg",dt,"s",1,&ia,&it);
// cout << " Minimal differences damin (deg) : " << da << " dtmin (s) : " << dt
//      << " damin-index : " << ia << " dtmin-index : " << it << endl;
// cout << " damin for "; lab->PrintSignal("equ","T",&ts,5,ia); cout << endl;
// cout << " dtmin for "; lab->PrintSignal("equ","T",&ts,5,it); cout << endl;
//
// // Search for space and time match with the reference signals
// da=5;
// dt=10;
// TArrayI* arr=lab.MatchRefSignal(da,"deg",dt,"s");
// Int_t index=0;
// if (arr)
// {
//  for (Int_t i=0; i<arr->GetSize(); i++)
//  {
//   index=arr->At(i);
//   cout << " Match found for index : " << index << endl;
//   cout << " Corresponding ref. object "; lab->PrintSignal("equ","T",&ts,5,index); cout << endl;
//  }
// }
//
//
//--- Author: Nick van Eijndhoven 15-mar-2007 Utrecht University
//- Modified: NvE $Date$ Utrecht University
///////////////////////////////////////////////////////////////////////////

#include "AliAstrolab.h"
#include "Riostream.h"
 
ClassImp(AliAstrolab) // Class implementation to enable ROOT I/O
 
AliAstrolab::AliAstrolab(const char* name,const char* title) : TTask(name,title),AliTimestamp()
{
// Default constructor

 fToffset=0;
 fXsig=0;
 fRefs=0;
 fBias=0;
 fGal=0;
 fIndices=0;
}
///////////////////////////////////////////////////////////////////////////
AliAstrolab::~AliAstrolab()
{
// Destructor to delete all allocated memory.

 if (fXsig)
 {
  delete fXsig;
  fXsig=0;
 }
 if (fRefs)
 {
  delete fRefs;
  fRefs=0;
 }
 if (fIndices)
 {
  delete fIndices;
  fIndices=0;
 }
}
///////////////////////////////////////////////////////////////////////////
AliAstrolab::AliAstrolab(const AliAstrolab& t) : TTask(t),AliTimestamp(t)
{
// Copy constructor

 fToffset=t.fToffset;
 fLabPos=t.fLabPos;
 fXsig=0;
 if (t.fXsig) fXsig=new AliSignal(*(t.fXsig));
 fRefs=0;
 if (t.fRefs)
 {
  Int_t size=t.fRefs->GetSize();
  fRefs=new TObjArray(size);
  for (Int_t i=0; i<size; i++)
  {
   AliSignal* sx=(AliSignal*)t.fRefs->At(i);
   if (sx) fRefs->AddAt(sx->Clone(),i);
  }
 }
 fBias=0;
 fGal=0;
 fIndices=0;
}
///////////////////////////////////////////////////////////////////////////
void AliAstrolab::Data(Int_t mode,TString u)
{
// Provide lab information.
//
// "mode" indicates the mode of the timestamp info (see AliTimestamp::Date).
//
// The string argument "u" allows to choose between different angular units
// in case e.g. a spherical frame is selected.
// u = "rad" : angles provided in radians
//     "deg" : angles provided in degrees
//     "dms" : angles provided in ddd:mm:ss.sss
//     "hms" : angles provided in hh:mm:ss.sss
//
// The defaults are mode=1 and u="deg".
 
 const char* name=GetName();
 const char* title=GetTitle();
 cout << " *" << ClassName() << "::Data*";
 if (strlen(name))  cout << " Name : " << GetName();
 if (strlen(title)) cout << " Title : " << GetTitle();
 cout << endl;

 Double_t l,b;
 GetLabPosition(l,b,"deg");
 cout << " Lab position longitude : "; PrintAngle(l,"deg",u,2);
 cout << " latitude : "; PrintAngle(b,"deg",u,2);
 cout << endl;
 cout << " Lab time offset w.r.t. UT : "; PrintTime(fToffset,12); cout << endl;

 // UT and Local time info
 Date(mode,fToffset);
} 
///////////////////////////////////////////////////////////////////////////
void AliAstrolab::SetLabPosition(Ali3Vector& p)
{
// Set the lab position in the terrestrial coordinates.
// The right handed reference frame is defined such that the North Pole
// corresponds to a polar angle theta=0 and the Greenwich meridian corresponds
// to an azimuth angle phi=0, with phi increasing eastwards.

 fLabPos.SetPosition(p);

 // Determine local time offset in fractional hours w.r.t. UT.
 Double_t vec[3];
 p.GetVector(vec,"sph","deg");
 Double_t l=vec[2];
 fToffset=l/15.;
}
///////////////////////////////////////////////////////////////////////////
void AliAstrolab::SetLabPosition(Double_t l,Double_t b,TString u)
{
// Set the lab position in the terrestrial longitude (l) and latitude (b).
// Positions north of the equator have b>0, whereas b<0 indicates
// locations south of the equator.
// Positions east of Greenwich have l>0, whereas l<0 indicates
// locations west of Greenwich.
//
// The string argument "u" allows to choose between different angular units
// u = "rad" : angles provided in radians
//     "deg" : angles provided in degrees
//     "dms" : angles provided in dddmmss.sss
//     "hms" : angles provided in hhmmss.sss
//
// The default is u="deg".

 Double_t r=1,theta=0,phi=0;

 l=ConvertAngle(l,u,"deg");
 b=ConvertAngle(b,u,"deg");

 Double_t offset=90.;

 theta=offset-b;
 phi=l;

 Double_t p[3]={r,theta,phi};
 fLabPos.SetPosition(p,"sph","deg");

 // Local time offset in fractional hours w.r.t. UT.
 fToffset=l/15.;
}
///////////////////////////////////////////////////////////////////////////
AliPosition AliAstrolab::GetLabPosition() const
{
// Provide the lab position in the terrestrial coordinates.
// The right handed reference frame is defined such that the North Pole
// corresponds to a polar angle theta=0 and the Greenwich meridian corresponds
// to an azimuth angle phi=0, with phi increasing eastwards.

 return fLabPos;
}
///////////////////////////////////////////////////////////////////////////
void AliAstrolab::GetLabPosition(Double_t& l,Double_t& b,TString u) const
{
// Provide the lab position in the terrestrial longitude (l) and latitude (b).
// Positions north of the equator have b>0, whereas b<0 indicates
// locations south of the equator.
// Positions east of Greenwich have l>0, whereas l<0 indicates
// locations west of Greenwich.
//
// The string argument "u" allows to choose between different angular units
// u = "rad" : angles provided in radians
//     "deg" : angles provided in degrees
//
// The default is u="deg".

 Double_t pi=acos(-1.);

 Double_t offset=90.;
 if (u=="rad") offset=pi/2.;

 Double_t p[3];
 fLabPos.GetPosition(p,"sph",u);
 b=offset-p[1];
 l=p[2];
}
///////////////////////////////////////////////////////////////////////////
Double_t AliAstrolab::GetLT()
{
// Provide the Lab's local time in fractional hours.
// A mean solar day lasts 24h (i.e. 86400s).
//
// In case a hh:mm:ss format is needed, please use the Convert() facility. 
 
 Double_t h=GetLT(fToffset);
 return h;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliAstrolab::GetLMST()
{
// Provide the Lab's Local Mean Sidereal Time (LMST) in fractional hours.
// A sidereal day corresponds to 23h 56m 04.09s (i.e. 86164.09s) mean solar time.
// The definition of GMST is such that a sidereal clock corresponds with
// 24 sidereal hours per revolution of the Earth.
// As such, local time offsets w.r.t. UT and GMST can be treated similarly. 
//
// In case a hh:mm:ss format is needed, please use the Convert() facility. 

 Double_t h=GetLMST(fToffset);
 return h;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliAstrolab::GetLAST()
{
// Provide the Lab's Local Apparent Sidereal Time (LAST) in fractional hours.
// A sidereal day corresponds to 23h 56m 04.09s (i.e. 86164.09s) mean solar time.
// The definition of GMST and GAST is such that a sidereal clock corresponds with
// 24 sidereal hours per revolution of the Earth.
// As such, local time offsets w.r.t. UT, GMST and GAST can be treated similarly. 
//
// In case a hh:mm:ss format is needed, please use the Convert() facility. 

 Double_t h=GetLAST(fToffset);
 return h;
}
///////////////////////////////////////////////////////////////////////////
void AliAstrolab::PrintAngle(Double_t a,TString in,TString out,Int_t ndig) const
{
// Printing of angles in various formats.
//
// The input argument "a" denotes the angle to be printed. 
// The string arguments "in" and "out" specify the angular I/O formats.
//
// in = "rad" : input angle provided in radians
//      "deg" : input angle provided in degrees
//      "dms" : input angle provided in dddmmss.sss
//      "hms" : input angle provided in hhmmss.sss
//
// out = "rad" : output angle provided in radians
//       "deg" : output angle provided in degrees
//       "dms" : output angle provided in dddmmss.sss
//       "hms" : output angle provided in hhmmss.sss
//
// The argument "ndig" specifies the number of digits for the fractional
// part (e.g. ndig=6 for "dms" corresponds to micro-arcsecond precision).
// No rounding will be performed, so an arcsecond count of 3.473 with ndig=1
// will appear as 03.4 on the output.
// Due to computer accuracy, precision on the pico-arcsecond level may get lost.
//
// The default is ndig=1.
//
// Note : The angle info is printed without additional spaces or "endline".
//        This allows the print to be included in various composite output formats.

 Double_t b=ConvertAngle(a,in,out);

 if (out=="deg" || out=="rad")
 {
  cout.setf(ios::fixed,ios::floatfield);
  cout << setprecision(ndig) << b << " " << out.Data();
  cout.unsetf(ios::fixed);
  return; 
 }

 Double_t epsilon=1.e-12; // Accuracy in (arc)seconds
 Int_t word=0,ddd=0,hh=0,mm=0,ss=0;
 Double_t s;
 ULong64_t sfrac=0;

 if (out=="dms")
 {
  word=Int_t(b);
  word=abs(word);
  ddd=word/10000;
  word=word%10000;
  mm=word/100;
  ss=word%100;
  s=fabs(b)-Double_t(ddd*10000+mm*100+ss);
  if (s>(1.-epsilon))
  {
   s=0.;
   ss++;
  }
  while (ss>=60)
  {
   ss-=60;
   mm++;
  }
  while (mm>=60)
  {
   mm-=60;
   ddd++;
  }
  while (ddd>=360)
  {
   ddd-=360;
  }
  s*=pow(10.,ndig);
  sfrac=ULong64_t(s);
  if (b<0) cout << "-";
  cout << ddd << "d " << mm << "' " << ss << "."
       << setfill('0') << setw(ndig) << sfrac << "\"";
  return;
 }

 if (out=="hms")
 {
  word=Int_t(b);
  word=abs(word);
  hh=word/10000;
  word=word%10000;
  mm=word/100;
  ss=word%100;
  s=fabs(b)-Double_t(hh*10000+mm*100+ss);
  if (s>(1.-epsilon))
  {
   s=0.;
   ss++;
  }
  while (ss>=60)
  {
   ss-=60;
   mm++;
  }
  while (mm>=60)
  {
   mm-=60;
   hh++;
  }
  while (hh>=24)
  {
   hh-=24;
  }
  s*=pow(10.,ndig);
  sfrac=ULong64_t(s);
  if (b<0) cout << "-";
  cout << hh << "h " << mm << "m " << ss << "."
       << setfill('0') << setw(ndig) << sfrac << "s";
  return;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliAstrolab::SetSignal(Ali3Vector* r,TString frame,TString mode,AliTimestamp* ts,Int_t jref,TString name)
{
// Store a signal as specified by the position r and the timestamp ts.
// The position is stored in International Celestial Reference System (ICRS) coordinates.
// The ICRS is a fixed, time independent frame and as such provides a unique reference
// frame without the need of specifying any epoch etc...
// The ICRS coordinate definitions match within 20 mas with the mean ones of the J2000.0
// equatorial system. Nevertheless, to obtain the highest accuracy, the slight
// coordinate correction between J2000 and ICRS is performed here via the
// so-called frame bias matrix.
// For further details see the U.S. Naval Observatory (USNO) circular 179 (2005),
// which is available on http://aa.usno,navy.mil/publications/docs/Circular_179.pdf.
//
// The input parameter "frame" allows the user to specify the frame to which
// the components of r refer. Available options are :
//
//  frame = "equ" ==> Equatorial coordinates with right ascension (a) and declination (d),
//                    where the "sph" components of r correspond to theta=(pi/2)-d and phi=a.
//          "gal" ==> Galactic coordinates with longitude (l) and lattitude (b).
//                    where the "sph" components of r correspond to theta=(pi/2)-b and phi=l.
//          "ecl" ==> Ecliptic coordinates with longitude (l) and lattitude (b),
//                    where the "sph" components of r correspond to theta=(pi/2)-b and phi=l.
//          "hor" ==> Horizontal coordinates at the AliAstrolab location, where the "sph"
//                    components of r correspond to theta=zenith angle and phi=pi-azimuth.
//          "icr" ==> ICRS coordinates with longitude (l) and lattitude (b),
//                    where the "sph" components of r correspond to theta=(pi/2)-b and phi=l.
//          "loc" ==> Local coordinates at the AliAstrolab location, where the "sph"
//                    components of r correspond to the usual theta and phi angles.
//
// In case the coordinates are the equatorial right ascension and declination (a,d),
// they can represent so-called "mean" and "true" values.
// The distinction between these two representations is the following :
//
// mean values : (a,d) are only corrected for precession and not for nutation
// true values : (a,d) are corrected for both precession and nutation
//
// The input parameter "mode" allows the user to specifiy either "mean" or "true"
// values for the input in case of equatorial (a,d) coordinates.
//
// mode = "M" --> Input coordinates are the mean values 
//        "T" --> Input coordinates are the true values 
//
// The input parameter "jref" allows the user to store so-called "reference" signals.
// These reference signals may be used to check space-time event coincidences with the
// stored measurement (e.g. coincidence of the measurement with transient phenomena).
//
// jref = 0 --> Storage of the measurement
//        j --> Storage of a reference signal at the j-th position (j=1 is first)
//      < 0 --> Add a reference signal at the next available position
//
// Via the input argument "name" the user can give the stored signal also a name.
//
// The default values are jref=0 and name="".
//
// Note : In case ts=0 the current timestamp of the lab will be taken.

 if (!r)
 {
  if (!jref && fXsig)
  {
   delete fXsig;
   fXsig=0;
  }
  return;
 }

 if (frame!="equ" && frame!="gal" && frame!="ecl" && frame!="hor" && frame!="icr" && frame!="loc")
 {
  if (!jref && fXsig)
  {
   delete fXsig;
   fXsig=0;
  }
  return;
 }

 if (frame=="equ" && mode!="M" && mode!="m" && mode!="T" && mode!="t")
 {
  if (!jref && fXsig)
  {
   delete fXsig;
   fXsig=0;
  }
  return;
 }

 if (!ts) ts=(AliTimestamp*)this;

 Double_t vec[3]={1,0,0};
 vec[1]=r->GetX(2,"sph","rad");
 vec[2]=r->GetX(3,"sph","rad");
 Ali3Vector q;
 q.SetVector(vec,"sph","rad"); 

 AliSignal* sxref=0;

 if (!jref) // Storage of the measurement
 {
  if (!fXsig)
  {
   fXsig=new AliSignal();
  }
  else
  {
   fXsig->Reset(1);
  }
  if (name != "") fXsig->SetName(name);
  fXsig->SetTitle("Event in ICRS coordinates");
  fXsig->SetTimestamp(*ts);
 }
 else // Storage of a reference signal
 {
  if (!fRefs) 
  {
   fRefs=new TObjArray();
   fRefs->SetOwner();
  }
  // Expand array size if needed
  if (jref>0 && jref>=fRefs->GetSize()) fRefs->Expand(jref+1);
  sxref=new AliSignal();
  if (name != "") sxref->SetName(name);
  sxref->SetTitle("Reference event in ICRS coordinates");
  sxref->SetTimestamp(*ts);
 }

 if (frame=="loc")
 {
  // Convert to horizontal coordinates
  q=q.GetUnprimed(&fL);

  // Store the signal
  SetSignal(&q,"hor",mode,ts,jref);
  return;
 }

 if (frame=="equ")
 {
  // Convert to "mean" values at specified epoch
  if (mode=="T" || mode=="t")
  {
   SetNmatrix(ts);
   q=q.GetUnprimed(&fN);
  }

  // Convert to "mean" values at J2000
  SetPmatrix(ts);
  q=q.GetUnprimed(&fP);

  // Convert to ICRS values
  if (!fBias) SetBmatrix(); 
  q=q.GetUnprimed(&fB);
 }

 if (frame=="gal")
 {
  // Convert to J2000 equatorial mean coordinates
  if (fGal != 2) SetGmatrix("J");
  q=q.GetUnprimed(&fG);

  // Convert to ICRS values
  if (!fBias) SetBmatrix(); 
  q=q.GetUnprimed(&fB);
 }

 if (frame=="ecl")
 {
  // Convert to mean equatorial values at specified epoch
  SetEmatrix(ts);
  q=q.GetUnprimed(&fE);

  // Convert to "mean" values at J2000
  SetPmatrix(ts);
  q=q.GetUnprimed(&fP);

  // Convert to ICRS values
  if (!fBias) SetBmatrix(); 
  q=q.GetUnprimed(&fB);
 }

 if (frame=="hor")
 {
  // Convert to "true" equatorial values at the specified timestamp
  SetHmatrix(ts);
  q=q.GetUnprimed(&fH);

  // Convert to "mean" values at specified timestamp
  SetNmatrix(ts);
  q=q.GetUnprimed(&fN);

  // Convert to "mean" values at J2000
  SetPmatrix(ts);
  q=q.GetUnprimed(&fP);

  // Convert to ICRS values
  if (!fBias) SetBmatrix(); 
  q=q.GetUnprimed(&fB);
 }

 // Store the signal in ICRS coordinates
 if (!jref) // Storage of a regular signal
 {
  fXsig->SetPosition(q);
 }
 else // Storage of a reference signal
 {
  sxref->SetPosition(q);
  if (jref<0)
  {
   fRefs->Add(sxref);
  }
  else
  {
   fRefs->AddAt(sxref,jref-1);
  }
 }
}
///////////////////////////////////////////////////////////////////////////
void AliAstrolab::SetSignal(Double_t a,Double_t d,TString s,Double_t e,TString mode,Int_t jref,TString name)
{
// Store a signal with right ascension (a) and declination (d) given for epoch e.
// The position is stored in International Celestial Reference System (ICRS) coordinates.
// The ICRS is a fixed, time independent frame and as such provides a unique reference
// frame without the need of specifying any epoch etc...
// The ICRS coordinate definitions match within 20 mas the mean ones of the J2000.0
// equatorial system. Nevertheless, to obtain the highest accuracy, the slight
// coordinate correction between J2000 and ICRS is performed here via the
// so-called frame bias matrix.
// For further details see the U.S. Naval Observatory (USNO) circular 179 (2005),
// which is available on http://aa.usno,navy.mil/publications/docs/Circular_179.pdf.
//
// The coordinates (a,d) can represent so-called "mean" and "true" values.
// The distinction between these two representations is the following :
//
// mean values : (a,d) are only corrected for precession and not for nutation
// true values : (a,d) are corrected for both precession and nutation
//
// The input parameter "mode" allows the user to specifiy either "mean" or "true"
// values for the input (a,d) coordinates.
//
// a    : Right ascension in hhmmss.sss
// d    : Declination in dddmmss.sss
// s    = "B" --> Besselian reference epoch.
//        "J" --> Julian reference epoch.
// e    : Reference epoch for the input coordinates (e.g. 1900, 1950, 2000,...)
// mode = "M" --> Input coordinates are the mean values 
//        "T" --> Input coordinates are the true values 
//
// The input parameter "jref" allows the user to store so-called "reference" signals.
// These reference signals may be used to check space-time event coincidences with the
// stored measurement (e.g. coincidence of the measurement with transient phenomena).
//
// jref = 0 --> Storage of the measurement
//        j --> Storage of a reference signal at the j-th position (j=1 is first)
//      < 0 --> Add a reference signal at the next available position
//
// Via the input argument "name" the user can give the stored signal also a name.
//
// The default values are jref=0 and name="".

 if (s!="B" && s!="b" && s!="J" && s!="j")
 {
  if (!jref && fXsig)
  {
   delete fXsig;
   fXsig=0;
  }
  return;
 }

 if (mode!="M" && mode!="m" && mode!="T" && mode!="t")
 {
  if (!jref && fXsig)
  {
   delete fXsig;
   fXsig=0;
  }
  return;
 }

 // Convert coordinates to fractional degrees.
 a=ConvertAngle(a,"hms","deg");
 d=ConvertAngle(d,"dms","deg");


 AliTimestamp tx;
 tx.SetEpoch(e,s);

 Ali3Vector r;
 Double_t vec[3]={1.,90.-d,a};
 r.SetVector(vec,"sph","deg");

 SetSignal(&r,"equ",mode,&tx,jref,name);
}
///////////////////////////////////////////////////////////////////////////
void AliAstrolab::SetSignal(Double_t a,Double_t d,TString mode,AliTimestamp* ts,Int_t jref,TString name)
{
// Store a signal with right ascension (a) and declination (d) given for timestamp ts.
// The position is stored in International Celestial Reference System (ICRS) coordinates.
// The ICRS is a fixed, time independent frame and as such provides a unique reference
// frame without the need of specifying any epoch etc...
// The ICRS coordinate definitions match within 20 mas the mean ones of the J2000.0
// equatorial system. Nevertheless, to obtain the highest accuracy, the slight
// coordinate correction between J2000 and ICRS is performed here via the
// so-called frame bias matrix.
// For further details see the U.S. Naval Observatory (USNO) circular 179 (2005),
// which is available on http://aa.usno,navy.mil/publications/docs/Circular_179.pdf.
//
// The coordinates (a,d) can represent so-called "mean" and "true" values.
// The distinction between these two representations is the following :
//
// mean values : (a,d) are only corrected for precession and not for nutation
// true values : (a,d) are corrected for both precession and nutation
//
// The input parameter "mode" allows the user to specifiy either "mean" or "true"
// values for the input (a,d) coordinates.
//
// a    : Right ascension in hhmmss.sss
// d    : Declination in dddmmss.sss
// mode = "M" --> Input coordinates are the mean values 
//        "T" --> Input coordinates are the true values 
//
// The input parameter "jref" allows the user to store so-called "reference" signals.
// These reference signals may be used to check space-time event coincidences with the
// stored measurement (e.g. coincidence of the measurement with transient phenomena).
//
// jref = 0 --> Storage of the measurement
//        j --> Storage of a reference signal at the j-th position (j=1 is first)
//      < 0 --> Add a reference signal at the next available position
//
// Via the input argument "name" the user can give the stored signal also a name.
//
// The default values are jref=0 and name="".
//
// Note : In case ts=0 the current timestamp of the lab will be taken.

 // Convert coordinates to fractional degrees.
 a=ConvertAngle(a,"hms","deg");
 d=ConvertAngle(d,"dms","deg");

 Ali3Vector r;
 Double_t vec[3]={1.,90.-d,a};
 r.SetVector(vec,"sph","deg");

 SetSignal(&r,"equ",mode,ts,jref,name);
}
///////////////////////////////////////////////////////////////////////////
AliSignal* AliAstrolab::GetSignal(Ali3Vector& r,TString frame,TString mode,AliTimestamp* ts,Int_t jref)
{
// Provide the user specified coordinates of a signal at the specific timestamp ts.
// The coordinates are returned via the vector argument "r".
// In addition also a pointer to the stored signal object is provided.
// In case no stored signal was available or one of the input arguments was
// invalid, the returned pointer will be 0.
//
// The input parameter "frame" allows the user to specify the frame to which
// the components of r refer. Available options are :
//
//  frame = "equ" ==> Equatorial coordinates with right ascension (a) and declination (d),
//                    where the "sph" components of r correspond to theta=(pi/2)-d and phi=a.
//          "gal" ==> Galactic coordinates with longitude (l) and lattitude (b).
//                    where the "sph" components of r correspond to theta=(pi/2)-b and phi=l.
//          "ecl" ==> Ecliptic coordinates with longitude (l) and lattitude (b),
//                    where the "sph" components of r correspond to theta=(pi/2)-b and phi=l.
//          "hor" ==> Horizontal coordinates at the AliAstrolab location, where the "sph"
//                    components of r correspond to theta=zenith angle and phi=pi-azimuth.
//          "icr" ==> ICRS coordinates with longitude (l) and lattitude (b),
//                    where the "sph" components of r correspond to theta=(pi/2)-b and phi=l.
//          "loc" ==> Local coordinates at the AliAstrolab location, where the "sph"
//                    components of r correspond to the usual theta and phi angles.
//
// In case the coordinates are the equatorial right ascension and declination (a,d),
// they can represent so-called "mean" and "true" values.
// The distinction between these two representations is the following :
//
// mean values : (a,d) are only corrected for precession and not for nutation
// true values : (a,d) are corrected for both precession and nutation
//
// The input parameter "mode" allows the user to specifiy either "mean" or "true"
// values for the input in case of equatorial (a,d) coordinates.
//
// mode = "M" --> Input coordinates are the mean values 
//        "T" --> Input coordinates are the true values 
//
// The input parameter "jref" allows the user to access so-called "reference" signals.
// These reference signals may be used to check space-time event coincidences with the
// stored measurement (e.g. coincidence of the measurement with transient phenomena).
//
// jref = 0 --> Access to the measurement
//        j --> Access to the reference signal at the j-th position (j=1 is first)
//
// Default value is jref=0.
//
// Note : In case ts=0 the current timestamp of the lab will be taken.

 r.SetZero();

 if (frame!="equ" && frame!="gal" && frame!="ecl" && frame!="hor" && frame!="icr" && frame!="loc") return 0;

 if (frame=="equ" && mode!="M" && mode!="m" && mode!="T" && mode!="t") return 0;

 AliSignal* sx=GetSignal(jref);

 if (!sx) return 0;

 if (!ts) ts=(AliTimestamp*)this;

 Double_t vec[3];
 sx->GetPosition(vec,"sph","rad");
 Ali3Vector q;
 q.SetVector(vec,"sph","rad");

 if (frame=="icr")
 {
  r.Load(q);
  return sx;
 }

 // Convert from ICRS to equatorial J2000 coordinates
 if (!fBias) SetBmatrix();
 q=q.GetPrimed(&fB);

 if (frame=="equ")
 {
  // Precess to specified timestamp
  AliTimestamp ts1;
  ts1.SetEpoch(2000,"J");
  Precess(q,&ts1,ts);

  // Nutation correction if requested
  if (mode=="T" || mode=="t") Nutate(q,ts);
 }

 if (frame=="gal")
 {
  // Convert from equatorial J2000 to galactic
  if (fGal != 2) SetGmatrix("J");
  q=q.GetPrimed(&fG);
 }

 if (frame=="ecl")
 {
  // Precess to specified timestamp
  AliTimestamp ts1;
  ts1.SetEpoch(2000,"J");
  Precess(q,&ts1,ts);

  // Convert from equatorial to ecliptic coordinates
  SetEmatrix(ts);
  q=q.GetPrimed(&fE);
 }

 if (frame=="hor")
 {
  // Precess to specified timestamp
  AliTimestamp ts1;
  ts1.SetEpoch(2000,"J");
  Precess(q,&ts1,ts);

  // Nutation correction
  Nutate(q,ts);

  // Convert from equatorial to horizontal coordinates
  SetHmatrix(ts);
  q=q.GetPrimed(&fH);
 }

 if (frame=="loc")
 {
  // Get the signal in horizontal coordinates
  GetSignal(q,"hor",mode,ts);

  // Convert from horizontal local-frame coordinates
  q=q.GetPrimed(&fL);
 }

 r.Load(q);
 return sx;
}
///////////////////////////////////////////////////////////////////////////
AliSignal* AliAstrolab::GetSignal(Ali3Vector& r,TString frame,TString mode,AliTimestamp* ts,TString name)
{
// Provide the user specified coordinates of the signal with the specified
// name at the specific timestamp ts.
// The coordinates are returned via the vector argument "r".
// In addition also a pointer to the stored signal object is provided.
// In case no such stored signal was available or one of the input arguments was
// invalid, the returned pointer will be 0.
//
// The input parameter "frame" allows the user to specify the frame to which
// the components of r refer. Available options are :
//
//  frame = "equ" ==> Equatorial coordinates with right ascension (a) and declination (d),
//                    where the "sph" components of r correspond to theta=(pi/2)-d and phi=a.
//          "gal" ==> Galactic coordinates with longitude (l) and lattitude (b).
//                    where the "sph" components of r correspond to theta=(pi/2)-b and phi=l.
//          "ecl" ==> Ecliptic coordinates with longitude (l) and lattitude (b),
//                    where the "sph" components of r correspond to theta=(pi/2)-b and phi=l.
//          "hor" ==> Horizontal coordinates at the AliAstrolab location, where the "sph"
//                    components of r correspond to theta=zenith angle and phi=pi-azimuth.
//          "icr" ==> ICRS coordinates with longitude (l) and lattitude (b),
//                    where the "sph" components of r correspond to theta=(pi/2)-b and phi=l.
//          "loc" ==> Local coordinates at the AliAstrolab location, where the "sph"
//                    components of r correspond to the usual theta and phi angles.
//
// In case the coordinates are the equatorial right ascension and declination (a,d),
// they can represent so-called "mean" and "true" values.
// The distinction between these two representations is the following :
//
// mean values : (a,d) are only corrected for precession and not for nutation
// true values : (a,d) are corrected for both precession and nutation
//
// The input parameter "mode" allows the user to specifiy either "mean" or "true"
// values for the input in case of equatorial (a,d) coordinates.
//
// mode = "M" --> Input coordinates are the mean values 
//        "T" --> Input coordinates are the true values 
//
// Note : In case ts=0 the current timestamp of the lab will be taken.

 AliSignal* sx=0;
 Int_t j=GetSignalIndex(name);
 if (j>=0) sx=GetSignal(r,frame,mode,ts,j);
 return sx;
}
///////////////////////////////////////////////////////////////////////////
AliSignal* AliAstrolab::GetSignal(Double_t& a,Double_t& d,TString mode,AliTimestamp* ts,Int_t jref)
{
// Provide precession (and nutation) corrected right ascension (a) and
// declination (d) of the stored signal object at the specified timestamp.
// In addition also a pointer to the stored signal object is provided.
// In case no stored signal was available or one of the input arguments was
// invalid, the returned pointer will be 0.
//
// The coordinates (a,d) can represent so-called "mean" and "true" values.
// The distinction between these two representations is the following :
//
// mean values : (a,d) are only corrected for precession and not for nutation
// true values : (a,d) are corrected for both precession and nutation
//
// The input parameter "mode" allows the user to select either
// "mean" or "true" values for (a,d).
//
// The correction methods used are the new IAU 2000 ones as described in the 
// U.S. Naval Observatory (USNO) circular 179 (2005), which is available on
// http://aa.usno,navy.mil/publications/docs/Circular_179.pdf.
//
// a    : Right ascension in hhmmss.sss
// d    : Declination in dddmmss.sss
// mode = "M" --> Output coordinates are the mean values 
//        "T" --> Output coordinates are the true values 
// ts   : Timestamp at which the corrected coordinate values are requested.
//
// The input parameter "jref" allows the user to access so-called "reference" signals.
// These reference signals may be used to check space-time event coincidences with the
// stored regular signal (e.g. coincidence of the stored signal with transient phenomena).
//
// jref = 0 --> Access to the measurement
//        j --> Access to the reference signal at the j-th position (j=1 is first)
//
// Default value is jref=0.
//
// Note : In case ts=0 the current timestamp of the lab will be taken.

 a=0;
 d=0;

 Ali3Vector r;
 AliSignal* sx=GetSignal(r,"equ",mode,ts,jref);

 if (!sx) return 0;

 // Retrieve the requested (a,d) values
 Double_t vec[3];
 r.GetVector(vec,"sph","deg");
 d=90.-vec[1];
 a=vec[2];

 while (a<-360.)
 {
  a+=360.;
 }
 while (a>360.)
 {
  a-=360.;
 }
 while (d<-90.)
 {
  d+=90.;
 }
 while (d>90.)
 {
  d-=90.;
 }

 // Convert coordinates to appropriate format
 a=ConvertAngle(a,"deg","hms");
 d=ConvertAngle(d,"deg","dms");

 return sx;
}
///////////////////////////////////////////////////////////////////////////
AliSignal* AliAstrolab::GetSignal(Double_t& a,Double_t& d,TString mode,AliTimestamp* ts,TString name)
{
// Provide precession (and nutation) corrected right ascension (a) and
// declination (d) of the stored signal object with the specified name
// at the specific timestamp ts.
// In addition also a pointer to the stored signal object is provided.
// In case no stored signal was available or one of the input arguments was
// invalid, the returned pointer will be 0.
//
// The coordinates (a,d) can represent so-called "mean" and "true" values.
// The distinction between these two representations is the following :
//
// mean values : (a,d) are only corrected for precession and not for nutation
// true values : (a,d) are corrected for both precession and nutation
//
// The input parameter "mode" allows the user to select either
// "mean" or "true" values for (a,d).
//
// The correction methods used are the new IAU 2000 ones as described in the 
// U.S. Naval Observatory (USNO) circular 179 (2005), which is available on
// http://aa.usno,navy.mil/publications/docs/Circular_179.pdf.
//
// a    : Right ascension in hhmmss.sss
// d    : Declination in dddmmss.sss
// mode = "M" --> Output coordinates are the mean values 
//        "T" --> Output coordinates are the true values 
// ts   : Timestamp at which the corrected coordinate values are requested.
// name : Name of the requested signal object
//
// Note : In case ts=0 the current timestamp of the lab will be taken.

 AliSignal* sx=0;
 Int_t j=GetSignalIndex(name);
 if (j>=0) sx=GetSignal(a,d,mode,ts,j);
 return sx;
}
///////////////////////////////////////////////////////////////////////////
AliSignal* AliAstrolab::GetSignal(Double_t& a,Double_t& d,TString s,Double_t e,TString mode,Int_t jref)
{
// Provide precession (and nutation) corrected right ascension (a) and
// declination (d) of the stored signal object at the specified epoch e.
// In addition also a pointer to the stored signal object is provided.
// In case no stored signal was available or one of the input arguments was
// invalid, the returned pointer will be 0.
//
// The coordinates (a,d) can represent so-called "mean" and "true" values.
// The distinction between these two representations is the following :
//
// mean values : (a,d) are only corrected for precession and not for nutation
// true values : (a,d) are corrected for both precession and nutation
//
// The input parameter "mode" allows the user to specifiy either "mean" or "true"
// values for the input (a,d) coordinates.
//
// a    : Right ascension in hhmmss.sss
// d    : Declination in dddmmss.sss
// s    = "B" --> Besselian reference epoch.
//        "J" --> Julian reference epoch.
// e    : Reference epoch for the input coordinates (e.g. 1900, 1950, 2000,...)
// mode = "M" --> Input coordinates are the mean values 
//        "T" --> Input coordinates are the true values 
//
// The input parameter "jref" allows the user to access so-called "reference" signals.
// These reference signals may be used to check space-time event coincidences with the
// stored measurement (e.g. coincidence of the measurement with transient phenomena).
//
// jref = 0 --> Access to the measurement
//        j --> Access to the reference signal at the j-th position (j=1 is first)
//
// Default value is jref=0.

 a=0;
 d=0;

 if (s!="B" && s!="b" && s!="J" && s!="j") return 0;

 if (mode!="M" && mode!="m" && mode!="T" && mode!="t") return 0;

 // Convert coordinates to fractional degrees.
 a=ConvertAngle(a,"hms","deg");
 d=ConvertAngle(d,"dms","deg");


 AliTimestamp tx;
 tx.SetEpoch(e,s);

 AliSignal* sx=GetSignal(a,d,mode,&tx,jref);
 return sx;
}
///////////////////////////////////////////////////////////////////////////
AliSignal* AliAstrolab::GetSignal(Double_t& a,Double_t& d,TString s,Double_t e,TString mode,TString name)
{
// Provide precession (and nutation) corrected right ascension (a) and
// declination (d) of the stored signal object with the specified name
// at the specific epoch e.
// In addition also a pointer to the stored signal object is provided.
// In case no stored signal was available or one of the input arguments was
// invalid, the returned pointer will be 0.
//
// The coordinates (a,d) can represent so-called "mean" and "true" values.
// The distinction between these two representations is the following :
//
// mean values : (a,d) are only corrected for precession and not for nutation
// true values : (a,d) are corrected for both precession and nutation
//
// The input parameter "mode" allows the user to specifiy either "mean" or "true"
// values for the input (a,d) coordinates.
//
// a    : Right ascension in hhmmss.sss
// d    : Declination in dddmmss.sss
// s    = "B" --> Besselian reference epoch.
//        "J" --> Julian reference epoch.
// e    : Reference epoch for the input coordinates (e.g. 1900, 1950, 2000,...)
// mode = "M" --> Input coordinates are the mean values 
//        "T" --> Input coordinates are the true values 
// name : Name of the requested signal object

 AliSignal* sx=0;
 Int_t j=GetSignalIndex(name);
 if (j>=0) sx=GetSignal(a,d,s,e,mode,j);
 return sx;
}
///////////////////////////////////////////////////////////////////////////
AliSignal* AliAstrolab::GetSignal(Int_t jref)
{
// Provide the pointer to a stored signal object.
//
// The input parameter "jref" allows the user to access so-called "reference" signals.
// These reference signals may be used to check space-time event coincidences with the
// stored measurement (e.g. coincidence of the measurement with transient phenomena).
//
// jref = 0 --> Access to the measurement
//        j --> Access to the reference signal at the j-th position (j=1 is first)
//
// Default value is jref=0.

 AliSignal* sx=0;
 if (!jref)
 {
  sx=fXsig;
 }
 else
 {
  if (jref>0 && jref<fRefs->GetSize()) sx=(AliSignal*)fRefs->At(jref-1);
 }
 return sx;
}
///////////////////////////////////////////////////////////////////////////
AliSignal* AliAstrolab::GetSignal(TString name)
{
// Provide the pointer to the stored signal object with the specified name.

 AliSignal* sx=0;
 Int_t j=GetSignalIndex(name);
 if (j>=0) sx=GetSignal(j);
 return sx;
}
///////////////////////////////////////////////////////////////////////////
void AliAstrolab::RemoveRefSignal(Int_t j,Int_t compress)
{
// Remove the reference signal which was stored at the j-th position (j=1 is first).
// Note : j=0 means that all stored reference signals will be removed.
//        j<0 allows array compression (see below) without removing any signals. 
//
// The "compress" parameter allows compression of the ref. signal storage array.
//
// compress = 1 --> Array will be compressed
//            0 --> Array will not be compressed
//
// Note : Compression of the storage array means that the indices of the
//        reference signals in the storage array will change.

 if (!fRefs) return;

 // Clearing of the complete storage
 if (!j)
 {
  delete fRefs;
  fRefs=0;
  return;
 }

 // Removing a specific reference signal
 if (j>0 && j<fRefs->GetSize())
 {
  TObject* obj=fRefs->RemoveAt(j-1);
  if (obj) delete obj;
 }

 // Compression of the storage array
 if (compress) fRefs->Compress();
}
///////////////////////////////////////////////////////////////////////////
void AliAstrolab::RemoveRefSignal(TString name,Int_t compress)
{
// Remove the reference signal with the specified name.
//
// The "compress" parameter allows compression of the ref. signal storage array.
//
// compress = 1 --> Array will be compressed
//            0 --> Array will not be compressed
//
// Note : Compression of the storage array means that the indices of the
//        reference signals in the storage array will change.

 Int_t j=GetSignalIndex(name);
 if (j>0) RemoveRefSignal(j,compress);
}
///////////////////////////////////////////////////////////////////////////
Int_t AliAstrolab::GetSignalIndex(TString name)
{
// Provide storage index of the signal with the specified name.
// In case the name matches with the stored measurement,
// the value 0 is returned.
// In case no signal with the specified name was found, the value -1 is returned.

 Int_t index=-1;
 
 if (fXsig)
 {
  if (name==fXsig->GetName()) return 0;
 }

 if (!fRefs) return -1;

 for (Int_t i=0; i<fRefs->GetSize(); i++)
 {
  AliSignal* sx=(AliSignal*)fRefs->At(i);
  if (!sx) continue;

  if (name==sx->GetName())
  {
   index=i+1;
   break;
  }
 }

 return index;
}
///////////////////////////////////////////////////////////////////////////
void AliAstrolab::PrintSignal(TString frame,TString mode,AliTimestamp* ts,Int_t ndig,Int_t jref)
{
// Print data of a stored signal in user specified coordinates at the specific timestamp ts.
// In case no stored signal was available or one of the input arguments was
// invalid, no printout is produced.
//
// The argument "ndig" specifies the number of digits for the fractional
// part (e.g. ndig=6 for "dms" corresponds to micro-arcsecond precision).
// No rounding will be performed, so an arcsecond count of 3.473 with ndig=1
// will appear as 03.4 on the output.
// Due to computer accuracy, precision on the pico-arcsecond level may get lost.
//
// Note : The angle info is printed without additional spaces or "endline".
//        This allows the print to be included in various composite output formats.
//
// The input parameter "frame" allows the user to specify the frame to which
// the coordinates refer. Available options are :
//
//  frame = "equ" ==> Equatorial coordinates with right ascension (a) and declination (d).
//
//          "gal" ==> Galactic coordinates with longitude (l) and lattitude (b).
//
//          "ecl" ==> Ecliptic coordinates with longitude (l) and lattitude (b).
//
//          "hor" ==> Horizontal azimuth and altitude coordinates at the AliAstrolab location.
//
//          "icr" ==> ICRS coordinates with longitude (l) and lattitude (b).
//
//          "loc" ==> Local spherical angles theta and phi at the AliAstrolab location.
//
// In case the coordinates are the equatorial right ascension and declination (a,d),
// they can represent so-called "mean" and "true" values.
// The distinction between these two representations is the following :
//
// mean values : (a,d) are only corrected for precession and not for nutation
// true values : (a,d) are corrected for both precession and nutation
//
// The input parameter "mode" allows the user to specifiy either "mean" or "true"
// values for the input in case of equatorial (a,d) coordinates.
//
// mode = "M" --> Input coordinates are the mean values 
//        "T" --> Input coordinates are the true values 
//
// The input parameter "mode" also determines which type of time and
// local hour angle will appear in the printout.
//
// mode = "M" --> Mean Sidereal Time (MST) and Local Mean Hour Angle (LMHA)
//        "T" --> Apparent Sidereal Time (AST) and Local Apparent Hour Angle (LAHA)
//
// The input parameter "jref" allows printing of a so-called "reference" signal.
// These reference signals may serve to check space-time event coincidences with the
// stored measurement (e.g. coincidence of the measurement with transient phenomena).
//
// jref = 0 --> Printing of the measurement
//        j --> Printing of the j-th reference signal
//
// Default value is jref=0.
//
// Note : In case ts=0 the current timestamp of the lab will be taken.

 Ali3Vector r;
 AliSignal* sx=GetSignal(r,frame,mode,ts,jref);

 if (!sx) return;

 // Local Hour Angle of the signal
 Double_t lha=GetHourAngle("A",ts,jref);
 TString slha="LAHA";
 if (mode=="M" || mode=="m")
 {
  lha=GetHourAngle("M",ts,jref);
  slha="LMHA";
 }

 TString name=sx->GetName();
 if (name != "") cout << name.Data() << " ";

 if (frame=="equ")
 {
  Double_t a,d;
  d=90.-r.GetX(2,"sph","deg");
  a=r.GetX(3,"sph","rad");
  cout << "Equatorial (" << mode.Data() <<") a : "; PrintAngle(a,"rad","hms",ndig);
  cout << " d : "; PrintAngle(d,"deg","dms",ndig);
  cout << " " << slha.Data() << " : "; PrintAngle(lha,"deg","hms",ndig);
  return;
 }

 if (frame=="gal")
 {
  Double_t l,b;
  b=90.-r.GetX(2,"sph","deg");
  l=r.GetX(3,"sph","deg");
  cout << "Galactic l : "; PrintAngle(l,"deg","deg",ndig);
  cout << " b : "; PrintAngle(b,"deg","deg",ndig); 
  cout << " " << slha.Data() << " : "; PrintAngle(lha,"deg","hms",ndig);
  return;
 }

 if (frame=="icr")
 {
  Double_t a,d;
  d=90.-r.GetX(2,"sph","deg");
  a=r.GetX(3,"sph","rad");
  cout << "ICRS l : "; PrintAngle(a,"rad","hms",ndig);
  cout << " b : "; PrintAngle(d,"deg","dms",ndig);
  cout << " " << slha.Data() << " : "; PrintAngle(lha,"deg","hms",ndig);
  return;
 }

 if (frame=="ecl")
 {
  Double_t a,d;
  d=90.-r.GetX(2,"sph","deg");
  a=r.GetX(3,"sph","deg");
  cout << "Ecliptic l : "; PrintAngle(a,"deg","deg",ndig);
  cout << " b : "; PrintAngle(d,"deg","deg",ndig);
  cout << " " << slha.Data() << " : "; PrintAngle(lha,"deg","hms",ndig);
  return;
 }

 if (frame=="hor")
 {
  Double_t alt=90.-r.GetX(2,"sph","deg");
  Double_t azi=180.-r.GetX(3,"sph","deg");
  while (azi>360)
  {
   azi-=360.;
  }
  while (azi<0)
  {
   azi+=360.;
  }
  cout << "Horizontal azi : "; PrintAngle(azi,"deg","deg",ndig);
  cout << " alt : "; PrintAngle(alt,"deg","deg",ndig);
  cout << " " << slha.Data() << " : "; PrintAngle(lha,"deg","hms",ndig);
  return;
 }

 if (frame=="loc")
 {
  Double_t theta=r.GetX(2,"sph","deg");
  Double_t phi=r.GetX(3,"sph","deg");
  cout << "Local-frame phi : "; PrintAngle(phi,"deg","deg",ndig);
  cout << " theta : "; PrintAngle(theta,"deg","deg",ndig);
  cout << " " << slha.Data() << " : "; PrintAngle(lha,"deg","hms",ndig);
  return;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliAstrolab::PrintSignal(TString frame,TString mode,AliTimestamp* ts,Int_t ndig,TString name)
{
// Print data of the stored signal with the specified name in user specified coordinates
// at the specific timestamp ts.
// In case such stored signal was available or one of the input arguments was
// invalid, no printout is produced.
//
// The argument "ndig" specifies the number of digits for the fractional
// part (e.g. ndig=6 for "dms" corresponds to micro-arcsecond precision).
// No rounding will be performed, so an arcsecond count of 3.473 with ndig=1
// will appear as 03.4 on the output.
// Due to computer accuracy, precision on the pico-arcsecond level may get lost.
//
// Note : The angle info is printed without additional spaces or "endline".
//        This allows the print to be included in various composite output formats.
//
// The input parameter "frame" allows the user to specify the frame to which
// the coordinates refer. Available options are :
//
//  frame = "equ" ==> Equatorial coordinates with right ascension (a) and declination (d).
//
//          "gal" ==> Galactic coordinates with longitude (l) and lattitude (b).
//
//          "ecl" ==> Ecliptic coordinates with longitude (l) and lattitude (b).
//
//          "hor" ==> Horizontal azimuth and altitude coordinates at the AliAstrolab location.
//
//          "icr" ==> ICRS coordinates with longitude (l) and lattitude (b).
//
//          "loc" ==> Local spherical angles theta and phi at the AliAstrolab location.
//
// In case the coordinates are the equatorial right ascension and declination (a,d),
// they can represent so-called "mean" and "true" values.
// The distinction between these two representations is the following :
//
// mean values : (a,d) are only corrected for precession and not for nutation
// true values : (a,d) are corrected for both precession and nutation
//
// The input parameter "mode" allows the user to specifiy either "mean" or "true"
// values for the input in case of equatorial (a,d) coordinates.
//
// mode = "M" --> Input coordinates are the mean values 
//        "T" --> Input coordinates are the true values 
//
// The input parameter "mode" also determines which type of time and
// local hour angle will appear in the printout.
//
// mode = "M" --> Mean Sidereal Time (MST) and Local Mean Hour Angle (LMHA)
//        "T" --> Apparent Sidereal Time (AST) and Local Apparent Hour Angle (LAHA)
//
// Note : In case ts=0 the current timestamp of the lab will be taken.

 Int_t j=GetSignalIndex(name);
 if (j>=0) PrintSignal(frame,mode,ts,ndig,j);
}
///////////////////////////////////////////////////////////////////////////
void AliAstrolab::ListSignals(TString frame,TString mode,Int_t ndig)
{
// List all stored signals in user specified coordinates at the timestamp
// of the actual recording of the stored measurement under investigation.
// In case no (timestamp of the) actual measurement is available,
// the current timestamp of the lab will be taken.
// In case no stored signal is available or one of the input arguments is
// invalid, no printout is produced.
//
// The argument "ndig" specifies the number of digits for the fractional
// part (e.g. ndig=6 for "dms" corresponds to micro-arcsecond precision).
// No rounding will be performed, so an arcsecond count of 3.473 with ndig=1
// will appear as 03.4 on the output.
// Due to computer accuracy, precision on the pico-arcsecond level may get lost.
//
// The default value is ndig=1.
//
// Note : The angle info is printed without additional spaces or "endline".
//        This allows the print to be included in various composite output formats.
//
// The input parameter "frame" allows the user to specify the frame to which
// the coordinates refer. Available options are :
//
//  frame = "equ" ==> Equatorial coordinates with right ascension (a) and declination (d).
//
//          "gal" ==> Galactic coordinates with longitude (l) and lattitude (b).
//
//          "ecl" ==> Ecliptic coordinates with longitude (l) and lattitude (b).
//
//          "hor" ==> Horizontal azimuth and altitude coordinates at the AliAstrolab location.
//
//          "icr" ==> ICRS coordinates with longitude (l) and lattitude (b).
//
//          "loc" ==> Local spherical angles theta and phi at the AliAstrolab location.
//
// In case the coordinates are the equatorial right ascension and declination (a,d),
// they can represent so-called "mean" and "true" values.
// The distinction between these two representations is the following :
//
// mean values : (a,d) are only corrected for precession and not for nutation
// true values : (a,d) are corrected for both precession and nutation
//
// The input parameter "mode" allows the user to specifiy either "mean" or "true"
// values for the input in case of equatorial (a,d) coordinates.
//
// mode = "M" --> Input coordinates are the mean values 
//        "T" --> Input coordinates are the true values 
//
// The input parameter "mode" also determines which type of time and
// local hour angle will appear in the listing.
//
// mode = "M" --> Mean Sidereal Time (MST) and Local Mean Hour Angle (LMHA)
//        "T" --> Apparent Sidereal Time (AST) and Local Apparent Hour Angle (LAHA)

 Int_t iprint=0;

 AliTimestamp* tx=0;

 Int_t dform=1;
 if (mode=="T" || mode=="t") dform=-1;

 if (fXsig)
 {
  tx=fXsig->GetTimestamp();
  if (!tx) tx=(AliTimestamp*)this;
  cout << " *AliAstrolab::ListSignals* List of all stored signals." << endl;
  cout << " === The measurement under investigation ===" << endl;
  cout << " Timestamp of the actual observation" << endl;
  tx->Date(dform,fToffset);
  cout << " Location of the actual observation" << endl;
  cout << " "; PrintSignal(frame,mode,tx,ndig); cout << endl;
  iprint=1;
 }

 if (!fRefs) return;

 for (Int_t i=1; i<=fRefs->GetSize(); i++)
 {
  AliSignal* sx=GetSignal(i);
  if (!sx) continue;

  if (!iprint)
  {
   cout << " *AliAstrolab::ListRefSignals* List of all stored signals." << endl;
   tx=(AliTimestamp*)this;
   cout << " Current timestamp of the laboratory" << endl;
   tx->Date(dform,fToffset);
   iprint=1;
  }
  if (iprint==1)
  {
   cout << " === All stored reference signals according to the above timestamp ===" << endl;
   iprint=2;
  }
  cout << " Index : " << i << " "; PrintSignal(frame,mode,tx,ndig,i); cout << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliAstrolab::Precess(Ali3Vector& r,AliTimestamp* ts1,AliTimestamp* ts2)
{
// Correct mean right ascension and declination given for Julian date "jd"
// for the earth's precession corresponding to the specified timestamp.
// The results are the so-called "mean" (i.e. precession corrected) values,
// corresponding to the specified timestamp.
// The method used is the new IAU 2000 one as described in the 
// U.S. Naval Observatory (USNO) circular 179 (2005), which is available on
// http://aa.usno,navy.mil/publications/docs/Circular_179.pdf.
// Since the standard reference epoch is J2000, this implies that all
// input (a,d) coordinates will be first internally converted to the
// corresponding J2000 values before the precession correction w.r.t. the
// specified lab timestamp will be applied.
//
// r  : Input vector containing the right ascension and declination information
//      in the form of standard Ali3Vector coordinates.
//      In spherical coordinates the phi corresponds to the right ascension,
//      whereas the declination corresponds to (pi/2)-theta.
// jd : Julian date corresponding to the input coordinate values.
// ts : Timestamp corresponding to the requested corrected coordinate values.
//
// Note : In case ts=0 the current timestamp of the lab will be taken.

 // Convert back to J2000 values
 Ali3Vector r0;
 SetPmatrix(ts1);
 r0=r.GetUnprimed(&fP);

 // Precess to the specified timestamp
 if (!ts2) ts2=(AliTimestamp*)this;
 SetPmatrix(ts2);
 r=r0.GetPrimed(&fP);
}
///////////////////////////////////////////////////////////////////////////
void AliAstrolab::Nutate(Ali3Vector& r,AliTimestamp* ts)
{
// Correct mean right ascension and declination for the earth's nutation
// corresponding to the specified timestamp.
// The results are the so-called "true" (i.e. nutation corrected) values,
// corresponding to the specified timestamp.
// The method used is the new IAU 2000 one as described in the 
// U.S. Naval Observatory (USNO) circular 179 (2005), which is available on
// http://aa.usno,navy.mil/publications/docs/Circular_179.pdf.
//
// r  : Input vector containing the right ascension and declination information
//      in the form of standard Ali3Vector coordinates.
//      In spherical coordinates the phi corresponds to the right ascension,
//      whereas the declination corresponds to (pi/2)-theta.
// ts : Timestamp for which the corrected coordinate values are requested.
//
// Note : In case ts=0 the current timestamp of the lab will be taken.

 // Nutation correction for the specified timestamp
 if (!ts) ts=(AliTimestamp*)this;
 SetNmatrix(ts);
 r=r.GetPrimed(&fN);
}
///////////////////////////////////////////////////////////////////////////
void AliAstrolab::SetBmatrix()
{
// Set the frame bias matrix elements.
// The formulas and numerical constants used are the ones from the 
// U.S. Naval Observatory (USNO) circular 179 (2005), which is available on
// http://aa.usno,navy.mil/publications/docs/Circular_179.pdf.

 Double_t pi=acos(-1.);
 
 // Parameters in mas
 Double_t a=-14.6;
 Double_t x=-16.6170;
 Double_t e=-6.8192;

 // Convert to radians
 a*=pi/(180.*3600.*1000.);
 x*=pi/(180.*3600.*1000.);
 e*=pi/(180.*3600.*1000.);

 Double_t mat[9];
 mat[0]=1.-0.5*(a*a+x*x);
 mat[1]=a;
 mat[2]=-x;
 mat[3]=-a-e*x;
 mat[4]=1.-0.5*(a*a+e*e);
 mat[5]=-e;
 mat[6]=x-e*a;
 mat[7]=e+x*a;
 mat[8]=1.-0.5*(e*e+x*x);

 fB.SetMatrix(mat);
 fBias=1;
}
///////////////////////////////////////////////////////////////////////////
void AliAstrolab::SetPmatrix(AliTimestamp* ts)
{
// Set precession matrix elements for Julian date jd w.r.t. J2000.
// The formulas and numerical constants used are the ones from the 
// U.S. Naval Observatory (USNO) circular 179 (2005), which is available on
// http://aa.usno,navy.mil/publications/docs/Circular_179.pdf.
// All numerical constants refer to the standard reference epoch J2000.

 Double_t mat[9]={0,0,0,0,0,0,0,0,0};
 if (!ts)
 {
  fP.SetMatrix(mat);
  return;
 }

 Double_t pi=acos(-1.);

 Double_t t=(ts->GetJD()-2451545.0)/36525.; // Julian centuries since J2000.0

 // Parameters for the precession matrix in arcseconds
 Double_t eps0=84381.406; // Mean ecliptic obliquity at J2000.0
 Double_t psi=5038.481507*t-1.0790069*pow(t,2)-0.00114045*pow(t,3)+0.000132851*pow(t,4)
                -0.0000000951*pow(t,4);
 Double_t om=eps0-0.025754*t+0.0512623*pow(t,2)-0.00772503*pow(t,3)-0.000000467*pow(t,4)
                 +0.0000003337*pow(t,5);
 Double_t chi=10.556403*t-2.3814292*pow(t,2)-0.00121197*pow(t,3)+0.000170663*pow(t,4)
              -0.0000000560*pow(t,5);

 // Convert to radians
 eps0*=pi/(180.*3600.);
 psi*=pi/(180.*3600.);
 om*=pi/(180.*3600.);
 chi*=pi/(180.*3600.);

 Double_t s1=sin(eps0);
 Double_t s2=sin(-psi);
 Double_t s3=sin(-om);
 Double_t s4=sin(chi);
 Double_t c1=cos(eps0);
 Double_t c2=cos(-psi);
 Double_t c3=cos(-om);
 Double_t c4=cos(chi);

 mat[0]=c4*c2-s2*s4*c3;
 mat[1]=c4*s2*c1+s4*c3*c2*c1-s1*s4*s3;
 mat[2]=c4*s2*s1+s4*c3*c2*s1+c1*s4*s3;
 mat[3]=-s4*c2-s2*c4*c3;
 mat[4]=-s4*s2*c1+c4*c3*c2*c1-s1*c4*s3;
 mat[5]=-s4*s2*s1+c4*c3*c2*s1+c1*c4*s3;
 mat[6]=s2*s3;
 mat[7]=-s3*c2*c1-s1*c3;
 mat[8]=-s3*c2*s1+c3*c1;

 fP.SetMatrix(mat);
}
///////////////////////////////////////////////////////////////////////////
void AliAstrolab::SetNmatrix(AliTimestamp* ts)
{
// Set nutation matrix elements for the specified Julian date jd.
// The formulas and numerical constants used are the ones from the 
// U.S. Naval Observatory (USNO) circular 179 (2005), which is available on
// http://aa.usno,navy.mil/publications/docs/Circular_179.pdf.

 Double_t mat[9]={0,0,0,0,0,0,0,0,0};
 if (!ts)
 {
  fN.SetMatrix(mat);
  return;
 }

 Double_t pi=acos(-1.);
 
 Double_t dpsi,deps,eps;
 ts->Almanac(&dpsi,&deps,&eps);

 // Convert to radians
 dpsi*=pi/(180.*3600.);
 deps*=pi/(180.*3600.);
 eps*=pi/(180.*3600.);

 Double_t s1=sin(eps);
 Double_t s2=sin(-dpsi);
 Double_t s3=sin(-(eps+deps));
 Double_t c1=cos(eps);
 Double_t c2=cos(-dpsi);
 Double_t c3=cos(-(eps+deps));

 mat[0]=c2;
 mat[1]=s2*c1;
 mat[2]=s2*s1;
 mat[3]=-s2*c3;
 mat[4]=c3*c2*c1-s1*s3;
 mat[5]=c3*c2*s1+c1*s3;
 mat[6]=s2*s3;
 mat[7]=-s3*c2*c1-s1*c3;
 mat[8]=-s3*c2*s1+c3*c1;

 fN.SetMatrix(mat);
}
///////////////////////////////////////////////////////////////////////////
void AliAstrolab::SetGmatrix(TString mode)
{
// Set the mean equatorial to galactic coordinate conversion matrix.
// The B1950 parameters were taken from section 22.3 of the book
// "An Introduction to Modern Astrophysics" by Carrol and Ostlie (1996).
// The J2000 parameters are obtained by precession of the B1950 values.
//
// Via the input argument "mode" the required epoch can be selected
// mode = "B" ==> B1950
//        "J" ==> J2000

 Ali3Vector x; // The Galactic x-axis in the equatorial frame
 Ali3Vector y; // The Galactic y-axis in the equatorial frame
 Ali3Vector z; // The Galactic z-axis in the equatorial frame

 Double_t a,d;
 Double_t vec[3]={1,0,0};

 fGal=1; // Set flag to indicate B1950 matrix values

 // B1950 equatorial coordinates of the North Galactic Pole (NGP)
 a=124900.;
 d=272400.;
 a=ConvertAngle(a,"hms","deg");
 d=ConvertAngle(d,"dms","deg");
 vec[1]=90.-d;
 vec[2]=a;
 z.SetVector(vec,"sph","deg");

 // B1950 equatorial coordinates of the Galactic l=b=0 point
 a=174224.;
 d=-285500.;
 a=ConvertAngle(a,"hms","deg");
 d=ConvertAngle(d,"dms","deg");
 vec[1]=90.-d;
 vec[2]=a;
 x.SetVector(vec,"sph","deg");

 // Precess to the corresponding J2000 values if requested
 if (mode=="J")
 {
  fGal=2; // Set flag to indicate J2000 matrix values
  AliTimestamp t1;
  t1.SetEpoch(1950,"B");
  AliTimestamp t2;
  t2.SetEpoch(2000,"J");
  Precess(z,&t1,&t2);
  Precess(x,&t1,&t2);
 }

 // The Galactic y-axis is determined for the right handed frame
 y=z.Cross(x);

 fG.SetAngles(x.GetX(2,"sph","deg"),x.GetX(3,"sph","deg"),
              y.GetX(2,"sph","deg"),y.GetX(3,"sph","deg"),
              z.GetX(2,"sph","deg"),z.GetX(3,"sph","deg"));
}
///////////////////////////////////////////////////////////////////////////
void AliAstrolab::SetEmatrix(AliTimestamp* ts)
{
// Set the mean equatorial to ecliptic coordinate conversion matrix
// for the specified timestamp.
// A nice sketch and explanation of the two frames can be found
// in chapter 3 of the book "Astronomy Methods" by Hale Bradt (2004).

 Double_t dpsi,deps,eps;
 ts->Almanac(&dpsi,&deps,&eps);

 // Convert to degrees
 eps/=3600.;

 // Positions of the ecliptic axes w.r.t. the equatorial ones
 // at the moment of the specified timestamp 
 Double_t theta1=90; // Ecliptic x-axis
 Double_t phi1=0;
 Double_t theta2=90.-eps; //Ecliptic y-axis
 Double_t phi2=90;
 Double_t theta3=eps; // Ecliptic z-axis
 Double_t phi3=270;

 fE.SetAngles(theta1,phi1,theta2,phi2,theta3,phi3);
}
///////////////////////////////////////////////////////////////////////////
void AliAstrolab::SetHmatrix(AliTimestamp* ts)
{
// Set the mean equatorial to horizontal coordinate conversion matrix
// for the specified timestamp.
// A nice sketch and explanation of the two frames can be found
// in chapter 3 of the book "Astronomy Methods" by Hale Bradt (2004).
//
// Note : In order to simplify the calculations, we use here a
//        right-handed horizontal frame.

 Ali3Vector x; // The (South pointing) horizontal x-axis in the equatorial frame
 Ali3Vector y; // The (East pointing) horizontal y-axis in the equatorial frame
 Ali3Vector z; // The (Zenith pointing) horizontal z-axis in the equatorial frame

 Double_t l,b;
 GetLabPosition(l,b,"deg");

 Double_t a;
 Double_t vec[3]={1,0,0};

 // Equatorial coordinates of the horizontal z-axis
 // at the moment of the specified timestamp 
 a=ts->GetLAST(fToffset);
 a*=15.; // Convert fractional hours to degrees 
 vec[1]=90.-b;
 vec[2]=a;
 z.SetVector(vec,"sph","deg");

 // Equatorial coordinates of the horizontal x-axis
 // at the moment of the specified timestamp 
 vec[1]=180.-b;
 vec[2]=a;
 x.SetVector(vec,"sph","deg");

 // The horizontal y-axis is determined for the right handed frame
 y=z.Cross(x);

 fH.SetAngles(x.GetX(2,"sph","deg"),x.GetX(3,"sph","deg"),
              y.GetX(2,"sph","deg"),y.GetX(3,"sph","deg"),
              z.GetX(2,"sph","deg"),z.GetX(3,"sph","deg"));
}
///////////////////////////////////////////////////////////////////////////
void AliAstrolab::SetLocalFrame(Double_t t1,Double_t p1,Double_t t2,Double_t p2,Double_t t3,Double_t p3)
{
// Specification of the orientations of the local-frame axes.
// The input arguments represent the angles (in degrees) of the local-frame axes
// w.r.t. a so called Master Reference Frame (MRF), with the same convention
// as the input arguments of TRrotMatix::SetAngles.
//
// The right handed Master Reference Frame is defined as follows :
//  Z-axis : Points to the local Zenith
//  X-axis : Makes an angle of 90 degrees with the Z-axis and points South
//  Y-axis : Makes an angle of 90 degrees with the Z-axis and points East
//
// Once the user has specified the local reference frame, any observed event
// can be related to astronomical space-time locations via the SetSignal
// and GetSignal memberfunctions.

 // Set the matrix for the conversion of our reference frame coordinates
 // into the local-frame ones.

 fL.SetAngles(t1,p1,t2,p2,t3,p3);
}
///////////////////////////////////////////////////////////////////////////
Double_t AliAstrolab::ConvertAngle(Double_t a,TString in,TString out) const
{
// Conversion of various angular formats.
//
// The input argument "a" denotes the angle to be converted. 
// The string arguments "in" and "out" specify the angular I/O formats.
//
// in = "rad" : input angle provided in radians
//      "deg" : input angle provided in degrees
//      "dms" : input angle provided in dddmmss.sss
//      "hms" : input angle provided in hhmmss.sss
//
// out = "rad" : output angle provided in radians
//       "deg" : output angle provided in degrees
//       "dms" : output angle provided in dddmmss.sss
//       "hms" : output angle provided in hhmmss.sss
 
 if (in==out) return a;

 // Convert input to its absolute value in (fractional) degrees. 
 Double_t pi=acos(-1.);
 Double_t epsilon=1.e-12; // Accuracy in (arc)seconds
 Int_t word=0,ddd=0,hh=0,mm=0,ss=0;
 Double_t s=0;

 Double_t b=fabs(a);

 if (in=="rad") b*=180./pi;

 if (in=="dms")
 {
  word=Int_t(b);
  ddd=word/10000;
  word=word%10000;
  mm=word/100;
  ss=word%100;
  s=b-Double_t(ddd*10000+mm*100+ss);
  b=Double_t(ddd)+Double_t(mm)/60.+(Double_t(ss)+s)/3600.;
 }

 if (in=="hms")
 {
  word=Int_t(b);
  hh=word/10000;
  word=word%10000;
  mm=word/100;
  ss=word%100;
  s=b-Double_t(hh*10000+mm*100+ss);
  b=15.*(Double_t(hh)+Double_t(mm)/60.+(Double_t(ss)+s)/3600.);
 }

 while (b>360)
 {
  b-=360.;
 }

 if (out=="rad") b*=pi/180.;

 if (out=="dms")
 {
  ddd=Int_t(b);
  b=b-Double_t(ddd);
  b*=60.;
  mm=Int_t(b);
  b=b-Double_t(mm);
  b*=60.;
  ss=Int_t(b);
  s=b-Double_t(ss);
  if (s>(1.-epsilon))
  {
   s=0.;
   ss++;
  }
  while (ss>=60)
  {
   ss-=60;
   mm++;
  }
  while (mm>=60)
  {
   mm-=60;
   ddd++;
  }
  while (ddd>=360)
  {
   ddd-=360;
  }
  b=Double_t(10000*ddd+100*mm+ss)+s;
 }

 if (out=="hms")
 {
  b/=15.;
  hh=Int_t(b);
  b=b-Double_t(hh);
  b*=60.;
  mm=Int_t(b);
  b=b-Double_t(mm);
  b*=60.;
  ss=Int_t(b);
  s=b-Double_t(ss);
  if (s>(1.-epsilon))
  {
   s=0.;
   ss++;
  }
  while (ss>=60)
  {
   ss-=60;
   mm++;
  }
  while (mm>=60)
  {
   mm-=60;
   hh++;
  }
  while (hh>=24)
  {
   hh-=24;
  }
  b=Double_t(10000*hh+100*mm+ss)+s;
 }

 if (a<0) b=-b;

 return b;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliAstrolab::GetHourAngle(TString mode,AliTimestamp* ts,Int_t jref)
{
// Provide the Local Hour Angle (in fractional degrees) of a stored signal
// object at the specified timestamp.
//
// The input parameter "mode" allows the user to select either the
// "mean" or "apparent" value for the returned Hour Angle.
//
// mode = "M" --> Output is the Mean Hour Angle
//        "A" --> Output is the Apparent Hour Angle
// ts   : Timestamp at which the hour angle is requested.
//
// The input parameter "jref" allows the user to specify so-called "reference" signals.
// These reference signals may be used to check space-time event coincidences with the
// stored measurement (e.g. coincidence of the measurement with transient phenomena).
//
// jref = 0 --> Use the stored measurement
//        j --> Use the reference signal at the j-th position (j=1 is first)
//
// Default value is jref=0.
//
// Note : In case ts=0 the current timestamp of the lab will be taken.

 if (!ts) ts=(AliTimestamp*)this;

 // Get corrected right ascension and declination for the specified timestamp.
 Double_t a,d;
 if (mode=="M" || mode=="m") GetSignal(a,d,"M",ts,jref);
 if (mode=="A" || mode=="a") GetSignal(a,d,"T",ts,jref);

 // Convert coordinates to fractional degrees.
 a=ConvertAngle(a,"hms","deg");
 d=ConvertAngle(d,"dms","deg");

 a/=15.; // Convert a to fractional hours
 Double_t ha=0;
 if (mode=="M" || mode=="m") ha=ts->GetLMST(fToffset)-a;
 if (mode=="A" || mode=="a") ha=ts->GetLAST(fToffset)-a;
 ha*=15.; // Convert to (fractional) degrees

 return ha;
}
///////////////////////////////////////////////////////////////////////////
void AliAstrolab::SetLT(Int_t y,Int_t m,Int_t d,Int_t hh,Int_t mm,Int_t ss,Int_t ns,Int_t ps)
{
// Set the AliTimestamp parameters corresponding to the LT date and time
// in the Gregorian calendar as specified by the input arguments.
//
// The input arguments represent the following :
// y  : year in LT (e.g. 1952, 2003 etc...)
// m  : month in LT (1=jan  2=feb etc...)
// d  : day in LT (1-31)
// hh : elapsed hours in LT (0-23) 
// mm : elapsed minutes in LT (0-59)
// ss : elapsed seconds in LT (0-59)
// ns : remaining fractional elapsed second of LT in nanosecond
// ps : remaining fractional elapsed nanosecond of LT in picosecond
//
// Note : ns=0 and ps=0 are the default values.
//

 SetLT(fToffset,y,m,d,hh,mm,ss,ns,ps);
}
///////////////////////////////////////////////////////////////////////////
void AliAstrolab::SetLT(Int_t y,Int_t d,Int_t s,Int_t ns,Int_t ps)
{
// Set the AliTimestamp parameters corresponding to the specified elapsed
// timespan since the beginning of the new LT year.
//
// The LT year and elapsed time span is entered via the following input arguments :
//
// y  : year in LT (e.g. 1952, 2003 etc...)
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

 SetLT(fToffset,y,d,s,ns,ps);
}
///////////////////////////////////////////////////////////////////////////
Double_t AliAstrolab::GetDifference(Int_t j,TString au,Double_t& dt,TString tu,Int_t mode,Int_t* ia,Int_t* it)
{
// Provide space and time difference between the j-th reference signal
// (j=1 indicates first) and the stored measurement.
// 
// The return value of this memberfunction provides the positional angular
// difference, whereas the output argument "dt" provides the time difference.
//
// The units of the angular difference can be specified via the the "au"
// input argument, where
//
// au = "rad" --> Angular difference in (fractional) radians
//      "deg" --> Angular difference in (fractional) degrees
//
// The units of the time difference can be specified via the "tu" and "mode"
// input arguments. For details please refer to AliTimestamp::GetDifference().
// Also here mode=1 is the default value.
//
// For the time difference the reference signal is used as the standard.
// This means that in case of a positive time difference, the stored
// measurement occurred later than the reference signal.
//
// In case j=0, the stored measurement will be compared with each
// reference signal and the returned angular and time differences are
// the minimal differences which were encountered.
// In this case the user may obtain the indices of the two stored reference signals
// which had the minimal angular and minimal time difference via the output
// arguments "ia" and "it" as follows :
//
// ia = Index of the stored reference signal with minimial angular difference
// it = Index of the stored reference signal with minimial time difference
//
// In case these indices are the same, there obviously was 1 single reference signal
// which showed both the minimal angular and time difference.
//
// The default values are mode=1, ia=0 and it=0;

 Double_t dang=999;
 dt=1.e20;
 if (ia) *ia=0;
 if (it) *it=0;

 if (!fRefs) return dang;

 Ali3Vector rx; // Position of the measurement
 Ali3Vector r0; // Position of the reference signal

 AliSignal* sx=GetSignal(rx,"icr","M",0);
 if (!sx) return dang;

 AliTimestamp* tx=sx->GetTimestamp();
 if (!tx) return dang;

 // Space and time difference w.r.t. a specific reference signal
 if (j>0)
 {
  AliSignal* s0=GetSignal(r0,"icr","M",0,j);
  if (!s0) return dang;
  AliTimestamp* t0=s0->GetTimestamp();

  if (!t0) return dang;

  dang=r0.GetOpeningAngle(rx,au);
  dt=t0->GetDifference(tx,tu,mode);
  return dang;
 }

 // Minimal space and time difference encountered over all reference signals
 Double_t dangmin=dang;
 Double_t dtmin=dt;
 for (Int_t i=1; i<=fRefs->GetSize(); i++)
 {
  AliSignal* s0=GetSignal(r0,"icr","M",0,i);
  if (!s0) continue;
  AliTimestamp* t0=s0->GetTimestamp();
  if (!t0) continue;
  dang=r0.GetOpeningAngle(rx,au);
  dt=t0->GetDifference(tx,tu,mode);
  if (fabs(dang)<dangmin)
  {
   dangmin=fabs(dang);
   *ia=i;
  }
  if (fabs(dt)<dtmin)
  {
   dtmin=fabs(dt);
   *it=i;
  }
 }

 dt=dtmin;
 return dangmin;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliAstrolab::GetDifference(TString name,TString au,Double_t& dt,TString tu,Int_t mode)
{
// Provide space and time difference between the reference signal
// with the specified name and the stored measurement.
// 
// The return value of this memberfunction provides the positional angular
// difference, whereas the output argument "dt" provides the time difference.
//
// The units of the angular difference can be specified via the the "au"
// input argument, where
//
// au = "rad" --> Angular difference in (fractional) radians
//      "deg" --> Angular difference in (fractional) degrees
//
// The units of the time difference can be specified via the "tu" and "mode"
// input arguments. For details please refer to AliTimestamp::GetDifference().
// Also here mode=1 is the default value.
//
// For the time difference the reference signal is used as the standard.
// This means that in case of a positive time difference, the stored
// measurement occurred later than the reference signal.

 Double_t dang=999;
 dt=1.e20;

 Int_t j=GetSignalIndex(name);
 if (j>0) dang=GetDifference(j,au,dt,tu,mode);
 return dang;
}
///////////////////////////////////////////////////////////////////////////
TArrayI* AliAstrolab::MatchRefSignal(Double_t da,TString au,Double_t dt,TString tu,Int_t mode)
{
// Provide the storage indices of the reference signals which match in space
// and time with the stored measurement.
// The indices are returned via a pointer to a TArrayI object.
// In case no matches were found, the null pointer is returned.
// A reference signal is regarded as matching with the stored measurement
// if the positional angular difference doesn't exceed "da" and the absolute
// value of the time difference doesn't exceed "dt".
//
// The units of the angular difference "da" can be specified via the the "au"
// input argument, where
//
// au = "rad" --> Angular difference in (fractional) radians
//      "deg" --> Angular difference in (fractional) degrees
//
// The units of the time difference "dt" can be specified via the "tu" and "mode"
// input arguments. For details please refer to AliTimestamp::GetDifference().
// Also here mode=1 is the default value.

 if (!fXsig || !fRefs) return 0;

 Int_t nrefs=fRefs->GetEntries();

 if (fIndices) delete fIndices;
 fIndices=new TArrayI(nrefs);

 Double_t dang,dtime;
 Int_t jfill=0; 
 for (Int_t i=1; i<=fRefs->GetSize(); i++)
 {
  dang=GetDifference(i,au,dtime,tu,mode);
  if (fabs(dang)<=da && fabs(dtime)<=dt)
  {
   fIndices->AddAt(i,jfill);
   jfill++;
  }
 }

 fIndices->Set(jfill);
 return fIndices;
}
///////////////////////////////////////////////////////////////////////////
TObject* AliAstrolab::Clone(const char* name) const
{
// Make a deep copy of the current object and provide the pointer to the copy.
// This memberfunction enables automatic creation of new objects of the
// correct type depending on the object type, a feature which may be very useful
// for containers when adding objects in case the container owns the objects.

 AliAstrolab* lab=new AliAstrolab(*this);
 if (name)
 {
  if (strlen(name)) lab->SetName(name);
 }
 return lab;
}
///////////////////////////////////////////////////////////////////////////
