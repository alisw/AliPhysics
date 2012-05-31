/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 * Author:  Denis Dujmic, SLAC
 *
 * Copyright (C) 2005 SLAC
 *******************************************************************************/

// Flatte propagator: S.M.Flatte, Phys. Lett. B63, 224 (1976)

#ifndef EVT_PROP_FLATTE_HH
#define EVT_PROP_FLATTE_HH

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtPropagator.hh"


class EvtPropFlatte :  public EvtPropagator {
public:

  EvtPropFlatte(double m0, 
		double g0, double m0a, double m0b, 
		double g1, double m1a, double m1b); 
  EvtPropFlatte(const EvtPropFlatte& other); 
  ~EvtPropFlatte(); 

  EvtAmplitude<EvtPoint1D>* clone() const;

protected:

  EvtComplex amplitude(const EvtPoint1D& x) const;

  double _m0a;
  double _m0b;

  double _g1;
  double _m1a;
  double _m1b;

};

#endif

