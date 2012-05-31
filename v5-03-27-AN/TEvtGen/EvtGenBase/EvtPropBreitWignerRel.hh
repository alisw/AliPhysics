/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtPropBreitWignerRel.hh,v 1.3 2003/06/20 17:20:11 dvoretsk Exp $
 *  Author: Alexei Dvoretskii, dvoretsk@slac.stanford.edu, 2001-2002
 *
 * Copyright (C) 2002 Caltech
 *******************************************************************************/

// Relativistic Breit-Wigner Propagator

#ifndef EVT_PROP_BREIT_WIGNER_REL_HH
#define EVT_PROP_BREIT_WIGNER_REL_HH

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtPropagator.hh"

class EvtPropBreitWignerRel : public EvtPropagator {
public:

  EvtPropBreitWignerRel(double m0, double g0); 
  EvtPropBreitWignerRel(const EvtPropBreitWignerRel& other); 
  ~EvtPropBreitWignerRel(); 

  EvtAmplitude<EvtPoint1D>* clone() const;

protected:

  EvtComplex amplitude(const EvtPoint1D& x) const;
};

#endif

