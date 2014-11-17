/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtDalitzResPdf.hh,v 1.2 2009-03-16 16:42:46 robbep Exp $
 *  Author: Alexei Dvoretskii, dvoretsk@slac.stanford.edu, 2001-2002
 *
 * Copyright (C) 2002 Caltech
 *******************************************************************************/

/*
 * Pole compensating function for terms that exibit a resonant structure
 * in one dimension only. 
 * 
 * f =  1    g*m0
 *     --  ------------------
 *     pi  (q-q0)^2 + g^2m0^2
 *
 * m is the mass of the resonance, g is its width. The approximation works well for a narrow
 * resonance. It is also readily integrable over the Dalitz plot coordinate to produce
 * 
 * Int = 1/pi atan((q-q0)/(g*m0))
 */

#ifndef EVT_DALITZ_RES_PDF_HH
#define EVT_DALITZ_RES_PDF_HH

#include "EvtGenBase/EvtPdf.hh"
#include "EvtGenBase/EvtDalitzPoint.hh"
#include "EvtGenBase/EvtCyclic3.hh"

class EvtDalitzResPdf : public EvtPdf<EvtDalitzPoint> {

public:

  EvtDalitzResPdf(const EvtDalitzPlot& dp,double m0, double g0, EvtCyclic3::Pair pairRes);
  EvtDalitzResPdf(const EvtDalitzResPdf& other);
  virtual ~EvtDalitzResPdf();
  
  
  EvtPdf<EvtDalitzPoint>* clone() const { return new EvtDalitzResPdf(*this); }
  
  virtual EvtValError compute_integral(int N) const;
  virtual EvtDalitzPoint randomPoint();
  double pdfMaxValue() const;

protected:

  virtual double pdf(const EvtDalitzPoint&) const; 

private:

  EvtDalitzPlot _dp;
  double _m0;                  // mass
  double _g0;                  // width
  EvtCyclic3::Pair _pair;      // resonant pair
};

#endif

