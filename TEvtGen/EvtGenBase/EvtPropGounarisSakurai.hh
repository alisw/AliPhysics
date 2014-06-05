/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtPropGounarisSakurai.hh,v 1.1 2009-03-16 16:50:50 robbep Exp $
 *  Author: Matt Graham 
 *  modified from EvtPropBreitWignerRel...this should be used for rho's
 *******************************************************************************/

// Relativistic Breit-Wigner Propagator

#ifndef EVT_PROP_GOUNARIS_SAKURAI_HH
#define EVT_PROP_GOUNARIS_SAKURAI_HH

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtPropagator.hh"
#include "EvtGenBase/EvtCyclic3.hh"
#include "EvtGenBase/EvtDalitzPoint.hh"
#include "EvtGenBase/EvtDalitzPlot.hh"

class EvtPropGounarisSakurai : public EvtPropagator {
public:

  EvtPropGounarisSakurai(EvtDalitzPlot *dp, 
                       EvtCyclic3::Pair pair, double m0, double g0); 
  EvtPropGounarisSakurai(const EvtPropGounarisSakurai& other); 
  ~EvtPropGounarisSakurai(); 

  EvtAmplitude<EvtPoint1D>* clone() const;

protected:

  EvtComplex amplitude(const EvtPoint1D& x) const;

private:
  EvtDalitzPlot *_dalitzSpace;

  EvtCyclic3::Pair      _pair;   
  double _gbase;
  double _m1;
  double _m2;
  double _dfun;
  double   dFun      ( double s  ) const;
  double   dh_dsFun  ( double s  ) const;
  double   hFun      ( double s  ) const;
  double   kFun      ( double s  ) const;
  double   fsFun     ( double s  ) const;

};

#endif

