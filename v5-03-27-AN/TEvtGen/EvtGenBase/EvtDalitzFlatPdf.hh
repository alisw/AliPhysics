/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtDalitzFlatPdf.hh,v 1.1 2003/06/20 17:20:03 dvoretsk Exp $
 *  Author: Alexei Dvoretskii, dvoretsk@slac.stanford.edu, 2001-2002
 *
 * Copyright (C) 2002 Caltech
 *******************************************************************************/

/*
 * Uniform PDF defined on a Dalitz plot.
 */

#ifndef EVT_DALITZ_FLAT_PDF_HH
#define EVT_DALITZ_FLAT_PDF_HH

#include <assert.h>
#include "EvtGenBase/EvtPdf.hh"
#include "EvtGenBase/EvtDalitzPlot.hh"
#include "EvtGenBase/EvtDalitzPoint.hh"

class EvtDalitzFlatPdf : public EvtPdf<EvtDalitzPoint> {
public:
  
  EvtDalitzFlatPdf(const EvtDalitzPlot& dp);
  EvtDalitzFlatPdf(const EvtDalitzFlatPdf& other);
  virtual ~EvtDalitzFlatPdf();
  virtual EvtPdf<EvtDalitzPoint>* clone() const;
  
  virtual EvtValError compute_integral(int N) const;
  virtual EvtDalitzPoint randomPoint();
  
protected:

  virtual double pdf(const EvtDalitzPoint&) const;

  EvtDalitzPlot _dp;
};

#endif
