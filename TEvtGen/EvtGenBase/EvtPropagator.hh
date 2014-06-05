/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtPropagator.hh,v 1.2 2009-03-16 16:42:03 robbep Exp $
 *  Author: Alexei Dvoretskii, dvoretsk@slac.stanford.edu, 2001-2002
 *
 * Copyright (C) 2002 Caltech
 *******************************************************************************/

// Defines propagator as a function of mass and width

#ifndef EVT_PROPAGATOR_HH
#define EVT_PROPAGATOR_HH

#include <assert.h>
#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtAmplitude.hh"
#include "EvtGenBase/EvtPoint1D.hh"

class EvtPropagator : public EvtAmplitude<EvtPoint1D> {
public:
  
  EvtPropagator(double m0, double g0)
    : _m0(m0), _g0(g0)
  {
    assert(m0 > 0);
    assert(g0 >= 0);
  }
  EvtPropagator(const EvtPropagator& other)
    : EvtAmplitude<EvtPoint1D>(other), _m0(other._m0), _g0(other._g0)
  {}
  virtual ~EvtPropagator()
  {}
  
  // Accessors 

  inline double m0() const  { return _m0; }
  inline double g0() const  { return _g0; }
  
  // Modifiers (can be useful e.g. for fitting!)

  inline void set_m0(double m0)  { assert(m0>0);  _m0 = m0; }
  inline void set_g0(double g0)  { assert(g0>=0); _g0 = g0; }

protected:

  double _m0;
  double _g0;

}; 

#endif

