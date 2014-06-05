/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtAmpAmpPdf.hh,v 1.2 2009-03-16 16:43:40 robbep Exp $
 *  Author: Alexei Dvoretskii, dvoretsk@slac.stanford.edu, 2001-2002
 *
 * Copyright (C) 2002 Caltech
 *******************************************************************************/

#ifndef EVT_AMP_AMP_PDF_HH
#define EVT_AMP_AMP_PDF_HH

// From the product A1A2* four PDF terms can be constructed, by taking the positive
// and the negative parts or the real and imaginary part of the product.

#include <assert.h>
#include "EvtGenBase/EvtMacros.hh"
#include "EvtGenBase/EvtAmplitude.hh"
#include "EvtGenBase/EvtPdf.hh"

enum {POSRE=0,NEGRE,POSIM,NEGIM};

template <class T>
class EvtAmpAmpPdf : public EvtPdf<T> {  
public:

  EvtAmpAmpPdf() {}
  EvtAmpAmpPdf(int type, const EvtAmplitude<T>& amp1, const EvtAmplitude<T>& amp2) 
    : EvtPdf<T>(), _type(type), _amp1(amp1.clone()), _amp2(amp2.clone()) 
  {}
  EvtAmpAmpPdf(const EvtAmpAmpPdf<T>& other) 
    : EvtPdf<T>(other), _type(other._type), COPY_PTR(_amp1), COPY_PTR(_amp2) 
  {}
  virtual ~EvtAmpAmpPdf() 
  { 
    delete _amp1; 
    delete _amp2; 
  }
    
  virtual EvtAmpAmpPdf<T>* clone() const { return new EvtAmpAmpPdf(*this); }

  virtual double pdf(const T& p) const
  {
    EvtComplex amp1 = _amp1->evaluate(p);
    EvtComplex amp2 = _amp2->evaluate(p);
    EvtComplex pr = amp1 * conj(amp2);

    if(_type == POSRE) return real(pr) > 0 ? real(pr) : 0.;
    if(_type == NEGRE) return real(pr) < 0 ? -real(pr) : 0.;
    if(_type == POSIM) return imag(pr) > 0 ? imag(pr) : 0.;
    if(_type == NEGIM) return imag(pr) < 0 ? -imag(pr) : 0.;
    
    assert(0);
  }
  
private:
  
  int _type;
  EvtAmplitude<T>* _amp1;
  EvtAmplitude<T>* _amp2;
};

#endif

