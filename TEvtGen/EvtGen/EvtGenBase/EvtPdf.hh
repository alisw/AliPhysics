/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtPdf.hh,v 1.2 2009-03-16 16:40:15 robbep Exp $
 *  Author: Alexei Dvoretskii, dvoretsk@slac.stanford.edu, 2001-2002
 *
 * Copyright (C) 2002 Caltech
 *******************************************************************************/

/*
 *  All classes are templated on the point type T
 *
 * EvtPdf:
 *
 * Probability density function defined on an interval of phase-space.
 * Integral over the interval can be calculated by Monte Carlo integration.
 * Some (but not all) PDFs are analytic in the sense that they can be integrated
 * by numeric quadrature and distributions can be generated according to them.
 *
 * EvtPdfGen:
 * 
 * Generator adaptor. Can be used to generate random points
 * distributed according to the PDF for analytic PDFs.
 *
 * EvtPdfPred:
 *
 * Predicate adaptor for PDFs. Can be used for generating random points distributed
 * according to the PDF for any PDF using rejection method. (See "Numerical Recipes").
 *
 * EvtPdfUnary:
 *
 * Adapter for generic algorithms. Evaluates the PDF and returns the value
 *
 * EvtPdfDiv:
 *
 * PDF obtained by division of one PDF by another. Because the two PDFs are 
 * arbitrary this PDF is not analytic. When importance sampling is used the 
 * original PDF is divided by the analytic comparison function. EvtPdfDiv is 
 * used to represent the modified PDF.
 */

#ifndef EVT_PDF_HH
#define EVT_PDF_HH

#include <assert.h>
#include <stdio.h>
#include "EvtGenBase/EvtValError.hh"
#include "EvtGenBase/EvtPredGen.hh"
#include "EvtGenBase/EvtStreamInputIterator.hh"
#include "EvtGenBase/EvtPdfMax.hh"
#include "EvtGenBase/EvtMacros.hh"
#include "EvtGenBase/EvtRandom.hh"

template <class T> class EvtPdfPred;
template <class T> class EvtPdfGen;

template <class T> class EvtPdf {
public:
  
  EvtPdf() {}
  EvtPdf(const EvtPdf& other) : _itg(other._itg) {}
  virtual ~EvtPdf() {}
  virtual EvtPdf<T>* clone() const = 0;
  
  double evaluate(const T& p) const { 
    if(p.isValid()) return pdf(p); 
    else return 0.;
  }

  // Find PDF maximum. Points are sampled according to pc

  EvtPdfMax<T> findMax(const EvtPdf<T>& pc, int N);

  // Find generation efficiency.

  EvtValError findGenEff(const EvtPdf<T>& pc, int N, int nFindMax);
  
  // Analytic integration. Calls cascade down until an overridden
  // method is called.

  void setItg(EvtValError itg) {_itg = itg; }
  
  EvtValError getItg() const {
    if(!_itg.valueKnown()) _itg = compute_integral();
    return _itg;
  }
  EvtValError getItg(int N) const {
    if(!_itg.valueKnown()) _itg = compute_integral(N);
    return _itg;
  }
  
  virtual EvtValError compute_integral() const
    //make sun happy - return something
  { printf("Analytic integration of PDF is not defined\n"); assert(0); return compute_integral();}
  virtual EvtValError compute_integral(int) const { return compute_integral(); }

  //  Monte Carlo integration.

  EvtValError compute_mc_integral(const EvtPdf<T>& pc, int N);
  
  // Generation. Create predicate accept-reject generators.
  // nMax iterations will be used to find the maximum of the accept-reject predicate
  
  EvtPredGen<EvtPdfGen<T>,EvtPdfPred<T> >  accRejGen(const EvtPdf<T>& pc, int nMax, double factor = 1.);
  
  virtual T randomPoint();

protected:

  virtual double pdf(const T&) const = 0;
  mutable EvtValError _itg;
};


template <class T> class EvtPdfGen {
public:
  typedef T result_type;

  EvtPdfGen() : _pdf(0) {}
  EvtPdfGen(const EvtPdfGen<T>& other) :
    _pdf(other._pdf ? other._pdf->clone() : 0)
  {}
  EvtPdfGen(const EvtPdf<T>& pdf) : 
    _pdf(pdf.clone())
  {}
  ~EvtPdfGen() { delete _pdf;}
  
  result_type operator()() {return _pdf->randomPoint();}
  
private:
  
  EvtPdf<T>* _pdf;
};


template <class T> class EvtPdfPred {
public:
  typedef T    argument_type;
  typedef bool result_type;
  
  EvtPdfPred() {}
  EvtPdfPred(const EvtPdf<T>& thePdf) : itsPdf(thePdf.clone()) {}
  EvtPdfPred(const EvtPdfPred& other) : COPY_PTR(itsPdf), COPY_MEM(itsPdfMax) {}
  ~EvtPdfPred() { delete itsPdf; }
  
  result_type operator()(argument_type p)
  {
    assert(itsPdf);
    assert(itsPdfMax.valueKnown());
    
    double random = EvtRandom::Flat(0.,itsPdfMax.value());
    return (random <= itsPdf->evaluate(p));     
  }
  
  EvtPdfMax<T> getMax() const { return itsPdfMax; }  
  void setMax(const EvtPdfMax<T>& max) { itsPdfMax = max; }
  template <class InputIterator> void compute_max(InputIterator it, InputIterator end,
						  double factor = 1.)
  {
    T p = *it++;
    itsPdfMax = EvtPdfMax<T>(p,itsPdf->evaluate(p)*factor);
    
    while(!(it == end)) {      
      T p = *it++;
      double val = itsPdf->evaluate(p)*factor;
      if(val > itsPdfMax.value()) itsPdfMax = EvtPdfMax<T>(p,val);
    }
  }
  
private:
  
  EvtPdf<T>*   itsPdf;
  EvtPdfMax<T> itsPdfMax;
};


template <class T> class EvtPdfUnary {
public:
  typedef double result_type;
  typedef T      argument_type;

  EvtPdfUnary() {}
  EvtPdfUnary(const EvtPdf<T>& thePdf) : itsPdf(thePdf.clone()) {}
  EvtPdfUnary(const EvtPdfUnary& other) : COPY_PTR(itsPdf) {}
  ~EvtPdfUnary() { delete itsPdf; }

  result_type operator()(argument_type p)
  {
    assert(itsPdf);
    double ret = itsPdf->evaluate(p);
    return ret;    
  }
  
private:
  
  EvtPdf<T>* itsPdf;
};


template <class T> class EvtPdfDiv : public EvtPdf<T> {
public:

  EvtPdfDiv() : itsNum(0), itsDen(0) {}
  EvtPdfDiv(const EvtPdf<T>& theNum, const EvtPdf<T>& theDen)
    : EvtPdf<T>(), itsNum(theNum.clone()), itsDen(theDen.clone())
  {}
  EvtPdfDiv(const EvtPdfDiv<T>& other)
    : EvtPdf<T>(other), COPY_PTR(itsNum), COPY_PTR(itsDen)
  {} 
  virtual ~EvtPdfDiv() { delete itsNum; delete itsDen; }
  virtual EvtPdf<T>* clone() const
  { return new EvtPdfDiv(*this); }
  
  virtual double pdf(const T& p) const
  {
    double num = itsNum->evaluate(p);
    double den = itsDen->evaluate(p);
    assert(den != 0);
    return num/den;
  }
  
private:
  
  EvtPdf<T>* itsNum; // numerator
  EvtPdf<T>* itsDen; // denominator
};  


template <class T>
EvtPdfMax<T> EvtPdf<T>::findMax(const EvtPdf<T>& pc, int N)
{
  EvtPdfPred<T> pred(*this);
  EvtPdfGen<T> gen(pc);
  pred.compute_max(iter(gen,N),iter(gen));
  EvtPdfMax<T> p = pred.getMax();
  return p;
}


template <class T>
EvtValError EvtPdf<T>::findGenEff(const EvtPdf<T>& pc, int N, int nFindMax)
{
  assert(N > 0 || nFindMax > 0);
  EvtPredGen<EvtPdfGen<T>,EvtPdfPred<T> > gen = accRejGen(pc,nFindMax);
  int i;
  for(i=0;i<N;i++) gen();
  double eff = double(gen.getPassed())/double(gen.getTried());
  double err = sqrt(double(gen.getPassed()))/double(gen.getTried());
  return EvtValError(eff,err);
}

template <class T>
EvtValError EvtPdf<T>::compute_mc_integral(const EvtPdf<T>& pc, int N)
{
  assert(N > 0);

  EvtValError otherItg = pc.getItg();
  EvtPdfDiv<T> pdfdiv(*this,pc);
  EvtPdfUnary<T> unary(pdfdiv);  
  
  EvtPdfGen<T> gen(pc);    
  EvtStreamInputIterator<T> begin = iter(gen,N);
  EvtStreamInputIterator<T> end;

  double sum = 0.;
  double sum2 = 0.;
  while(!(begin == end)) {
    
    double value = pdfdiv.evaluate(*begin++);
    sum += value;
    sum2 += value*value;
  }
  
  EvtValError x;
  if(N > 0) {
    double av = sum/((double) N);
    if(N > 1) {
      double dev2 = (sum2 - av*av*N)/((double) (N - 1));
      // Due to numerical precision dev2 may sometimes be negative
      if(dev2 < 0.) dev2 = 0.;
      double error = sqrt(dev2/((double) N));
      x = EvtValError(av,error);
    }
    else x = EvtValError(av);
  }
  _itg = x * pc.getItg();
  return _itg;
}

template <class T>
T EvtPdf<T>::randomPoint()
{
  printf("Function defined for analytic PDFs only\n");
  assert(0);
  T temp;
  return temp;
}

template <class T>
EvtPredGen<EvtPdfGen<T>,EvtPdfPred<T> > 
EvtPdf<T>::accRejGen(const EvtPdf<T>& pc, int nMax, double factor)
{
  EvtPdfGen<T> gen(pc);
  EvtPdfDiv<T> pdfdiv(*this,pc);
  EvtPdfPred<T> pred(pdfdiv);
  pred.compute_max(iter(gen,nMax),iter(gen),factor);
  return EvtPredGen<EvtPdfGen<T>,EvtPdfPred<T> >(gen,pred);
}

#endif





