/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtPdfSum.hh,v 1.3 2009-03-16 16:40:16 robbep Exp $
 *  Author: Alexei Dvoretskii, dvoretsk@slac.stanford.edu, 2001-2002
 *
 * Copyright (C) 2002 Caltech
 *******************************************************************************/

// Sum of PDF functions. 

#ifndef EVT_PDF_SUM_HH
#define EVT_PDF_SUM_HH

#include <stdio.h>
#include <vector>
using std::vector;
#include "EvtGenBase/EvtPdf.hh"

template <class T>
class EvtPdfSum : public EvtPdf<T> {
public:

  EvtPdfSum() {}
  EvtPdfSum(const EvtPdfSum<T>& other);
  virtual ~EvtPdfSum();
  virtual EvtPdf<T>* clone() const { return new EvtPdfSum(*this); }


  // Manipulate terms and coefficients
  
  void addTerm(double c,const EvtPdf<T>& pdf)
  { assert(c >= 0.); _c.push_back(c); _term.push_back(pdf.clone()); }

  void addOwnedTerm(double c, EvtPdf<T>* pdf)
  { _c.push_back(c); _term.push_back(pdf); }
  
  size_t nTerms() const { return _term.size(); }  // number of terms
  
  inline double   c(int i) const { return _c[i]; }
  inline EvtPdf<T>* getPdf(int i) const { return _term[i]; }


  // Integrals

  virtual EvtValError compute_integral() const;
  virtual EvtValError compute_integral(int N) const;
  virtual T randomPoint();
  
protected:
  
  virtual double pdf(const T& p) const;
  
  vector<double> _c;                     // coefficients
  vector<EvtPdf<T>*> _term;       // pointers to pdfs
}; 


template <class T>
EvtPdfSum<T>::EvtPdfSum(const EvtPdfSum<T>& other)
  : EvtPdf<T>(other)
{
  for(size_t i = 0; i < other.nTerms(); i++) {
    _c.push_back(other._c[i]);
    _term.push_back(other._term[i]->clone());
  }
}

template <class T>
EvtPdfSum<T>::~EvtPdfSum()
{
  for(size_t i = 0; i < _c.size(); i++) {
    delete _term[i];
  }
}


template <class T>
double EvtPdfSum<T>::pdf(const T& p) const
{
  double ret = 0.;
  for(size_t i=0; i < _c.size(); i++) {
    ret += _c[i] * _term[i]->evaluate(p);
  }
  return ret;
}

/*
 * Compute the sum integral by summing all term integrals.
 */

template <class T>
EvtValError EvtPdfSum<T>::compute_integral() const 
{
  EvtValError itg(0.0,0.0);
  for(size_t i=0;i<nTerms();i++) {
    itg += _c[i]*_term[i]->getItg();
  }
  return itg;
}

template <class T>
EvtValError EvtPdfSum<T>::compute_integral(int N) const
{
  EvtValError itg(0.0,0.0);
  for(size_t i=0;i<nTerms();i++) itg += _c[i]*_term[i]->getItg(N);
  return itg;
}


/*
 * Sample points randomly according to the sum of PDFs. First throw a random number uniformly
 * between zero and the value of the sum integral. Using this random number select one
 * of the PDFs. The generate a random point according to that PDF.
 */

template <class T>
T EvtPdfSum<T>::randomPoint()
{
  if(!this->_itg.valueKnown()) this->_itg = compute_integral();      
  
  double max = this->_itg.value();
  double rnd = EvtRandom::Flat(0,max);
  
  double sum = 0.;
  size_t i;
  for(i = 0; i < nTerms(); i++) {
    double itg = _term[i]->getItg().value();
    sum += _c[i] * itg;
    if(sum > rnd) break;
  }
  
  return _term[i]->randomPoint(); 
}

#endif


