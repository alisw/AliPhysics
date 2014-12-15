/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtIntegPdf1D.hh,v 1.2 2009-03-16 16:42:03 robbep Exp $
 *  Author: Alexei Dvoretskii, dvoretsk@slac.stanford.edu, 2001-2002
 *
 * Copyright (C) 2002 Caltech
 *******************************************************************************/

// Analytically integrable one dimensional PDF.

#ifndef EVT_INTEG_PDF_1D_HH
#define EVT_INTEG_PDF_1D_HH

#include "EvtGenBase/EvtPdf.hh"
#include "EvtGenBase/EvtPoint1D.hh"

class EvtIntegPdf1D : public EvtPdf<EvtPoint1D> {

public:
  
  EvtIntegPdf1D(double min, double max);
  EvtIntegPdf1D(const EvtIntegPdf1D&);
  virtual ~EvtIntegPdf1D();

  // Pdf integral function and its inverse to be defined in subclasses  

  virtual double pdfIntegral(double x) const = 0;
  virtual double pdfIntegralInverse(double x) const = 0;

  virtual EvtValError compute_integral() const;
  virtual EvtPoint1D randomPoint();

protected: 
  
  double _min;
  double _max;
};


#endif


