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
//      F(what,Lambda,b) = pow(_b,_b)/(tgamma(_b)*_Lambda)*pow(what/_Lambda,_b-1)*
//                           exp(-_b*what/Lambda);
//      the leading order shape function (exp) (hep-ph/0504071)
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Sven Menke (DFN model)
//      Alexei Volk (BLNP model)
//-----------------------------------------------------------------------

#ifndef EVTPFERMI_HH
#define EVTPFERMI_HH

class EvtPFermi {

public:
  
  // Constructors

  EvtPFermi(const double &a, const double &mB, const double &mb);
  EvtPFermi(const double &Lambda, const double &b);
  
  // Destructor

  virtual ~EvtPFermi( );

  // Operators

  // Selectors 

  // Modifiers

  // Methods

  double getFPFermi(const double &kplus);
  double getSFBLNP(const double &what);
  
protected:
  
  // Helper functions

private:

  // Friends
  
  // Data members

  double _a;
  double _mb;
  double _mB;
  double _Lambda;
  double _b;
};


#endif // EVTPFERMI_HH


