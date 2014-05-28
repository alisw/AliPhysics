/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtAmplitudeSum.hh,v 1.2 2009-03-16 16:42:46 robbep Exp $
 *  Author: Alexei Dvoretskii, dvoretsk@slac.stanford.edu, 2001-2002
 *
 * Copyright (C) 2002 Caltech
 *******************************************************************************/

#ifndef EVT_AMPLITUDE_SUM_HH
#define EVT_AMPLITUDE_SUM_HH

#include <stdio.h>
#include <assert.h>
#include <vector>
#include "EvtGenBase/EvtAmplitude.hh"

template <class T>
class EvtAmplitudeSum : public EvtAmplitude<T> {
  
public:
  
  EvtAmplitudeSum() {}
  EvtAmplitudeSum(const EvtAmplitudeSum<T>& other)
    : EvtAmplitude<T>(other)
  {
    int i;
    for(i=0;i<other.nTerms();i++) {
      
      EvtComplex c = other.c(i);      
      _c.push_back(c);
      EvtAmplitude<T>* amp = other.getTerm(i);
      assert(amp);
      EvtAmplitude<T>* amp1 = amp->clone();
      assert(amp1);
      _term.push_back(amp1);
    }
  }
  
  virtual ~EvtAmplitudeSum()
  {
    for(size_t i=0;i<_term.size();i++) {
      
      delete _term[i];
    }
  }  
  
  virtual EvtAmplitude<T>* clone() const
  {
    return new EvtAmplitudeSum<T>(*this);
  }
  
  
  void addTerm(EvtComplex c,const EvtAmplitude<T>& amp)
  {
    _c.push_back(c);
    _term.push_back(amp.clone());
  }
  
  void addOwnedTerm(EvtComplex c, EvtAmplitude<T>* amp)
  {
    assert(amp);
    _c.push_back(c);
    _term.push_back(amp);
  }
  
  int nTerms() const { return _term.size(); }  // number of terms

  void print() const {

    int N = nTerms();
    printf("Amplitude has %d terms\n",N);
    int i;
    for(i=0;i<N;i++) {
      printf("c%d = (%f,%f)\n",i,real(_c[i]),imag(_c[i]));
      assert(_term[i]);
    }
  }


  inline EvtComplex   c(int i) const { return _c[i]; }
  inline EvtAmplitude<T>*  getTerm(int i) const { return _term[i]; }

protected:

  virtual EvtComplex amplitude(const T& p) const
  {
    if(_term.size() == 0) 
      printf("Warning: amplitude sum has zero terms\n");
    
    EvtComplex value = 0.;

    for(size_t i=0;i<_term.size();i++) {    
      value+=_c[i]*_term[i]->evaluate(p);    
    }    
    return value;
  }

private:
  
  std::vector<EvtComplex> _c;              // coefficients
  std::vector<EvtAmplitude<T>*> _term;  // pointers to amplitudes
}; 


#endif


