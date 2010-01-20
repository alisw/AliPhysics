/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtFlatAmp.hh,v 1.4 2009/02/18 03:31:38 ryd Exp $
 *  Author: Alexei Dvoretskii, dvoretsk@slac.stanford.edu, 2001-2002
 *
 * Copyright (C) 2002 Caltech
 *******************************************************************************/

// Flat amplitude

#ifndef EVT_FLAT_AMP_HH
#define EVT_FLAT_AMP_HH

#include "EvtGenBase/EvtAmplitude.hh"

template <class T>
class EvtFlatAmp : public EvtAmplitude<T> {  
public:

  EvtFlatAmp() {}
  EvtFlatAmp(const EvtFlatAmp<T>& other) : EvtAmplitude<T>(other) {}
  virtual ~EvtFlatAmp() {}

  virtual EvtAmplitude<T>* clone() const { return new EvtFlatAmp<T>(*this); }
  virtual EvtComplex amplitude(const T& ) const { return EvtComplex(1.,0.); }
}; 

#endif
