//-----------------------------------------------------------------------
// File and Version Information:
//
// Copyright Information: See EvtGen/COPYRIGHT
//
//
// Description:
//   DFN model:
//      F(k+) = N (1-x)^a exp((1+a)x) ,x=k+/(mB-mb) 
//      the fermi motion distribution according to
//      hep-ph/9905351 v2
//   BLNP model:
//      F(wtilde,Lambda,b) = pow(_b,_b)/(tgamma(_b)*_Lambda)*pow(wtilde/_Lambda,_b-1)*
//                           exp(-_b*wtilde/Lambda);
//      the leading order shape function (exp) (hep-ph/0504071)
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Sven Menke (DFN model)
//      Alexei Volk (BLNP model)
//-----------------------------------------------------------------------
//-----------------------
// This Class's Header --
//-----------------------
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenModels/EvtPFermi.hh"
#include "EvtGenBase/EvtReport.hh"
//---------------
// C Headers --
//---------------
#include <math.h>
#include <stdlib.h>

//----------------
// Constructors --
//----------------

//for DFN model
EvtPFermi::EvtPFermi(const double &a, const double &mB, const double &mb)
{
  _a = a;
  _mb = mb;
  _mB = mB;
}

// for BLNP modell
EvtPFermi::EvtPFermi(const double &Lambda, const double &b)
{
  _Lambda = Lambda;
  _b = b;
}


//--------------
// Destructor --
//--------------

EvtPFermi::~EvtPFermi( )
{
}

//-----------
// Methods --
//-----------

double EvtPFermi::getFPFermi(const double &kplus)
{
  double FKplus;
  double x = kplus/(_mB-_mb);

  if ( x      >= 1)   return 0;
  if ( kplus <= -_mb) return 0; 

  FKplus = pow(1-x,_a)*exp((1+_a)*x);

  return FKplus;
}

// get value for the leading order exponential SF 
double EvtPFermi::getSFBLNP(const double &what)
{
  double SF;
  double massB = 5.2792; 
  

  if ( what      > massB )   return 0;
  if ( what < 0 ) return 0; 

#if defined(__SUNPRO_CC)
  report(Severity::Error,"EvtGen") << "The tgamma function is not available on this platform\n";
  report(Severity::Error,"EvtGen") <<"Presumably, you are getting the wrong answer, so I abort..";
  ::abort();
#else
  SF = pow(_b,_b)/(tgamma(_b)*_Lambda)*pow(what/_Lambda,_b-1)*exp(-_b*what/_Lambda); 
#endif
  
  return SF;
}

