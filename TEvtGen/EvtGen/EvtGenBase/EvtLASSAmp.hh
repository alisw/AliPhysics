/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *  Authors: Denis Dujmic, ddujmic@slac.stanford.edu
 *           Andrew Wagner, apwagner@slac.stanford.edu
 *
 * Copyright (C) 2005 SLAC
 *******************************************************************************/

#ifndef EVT_LASS_AMP_HH
#define EVT_LASS_AMP_HH

#include <string>

#include "EvtGenBase/EvtCyclic3.hh"
#include "EvtGenBase/EvtDalitzPoint.hh"
#include "EvtGenBase/EvtDalitzPlot.hh"
#include "EvtGenBase/EvtAmplitude.hh"

class EvtComplex;

class EvtLASSAmp : public EvtAmplitude<EvtDalitzPoint> {
  
public:

  EvtLASSAmp( EvtDalitzPlot *dp, 
	      EvtCyclic3::Pair pair,
	      double m0, double g0,
	      double a, double r, double cutoff, std::string subtype = "LASS" );

  EvtLASSAmp(const EvtLASSAmp& other);

  ~EvtLASSAmp();

  virtual EvtComplex amplitude(const EvtDalitzPoint& p) const;

  virtual EvtAmplitude<EvtDalitzPoint>* clone() const
  { return new EvtLASSAmp(*this); }

private:

  EvtDalitzPlot *_dalitzSpace;

  EvtCyclic3::Pair      _pair;   

  double _m0;
  double _g0;
  double _q0;
  double _r;
  double _a;
  double _cutoff;
  std::string _subtype;

};

#endif








