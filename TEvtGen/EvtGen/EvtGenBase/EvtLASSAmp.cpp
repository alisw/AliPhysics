/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory 
 * Package: EvtGenBase 
 *  Author: Denis Dujmic, ddujmic@slac.stanford.edu               
 *  
 * Copyright (C) 2005 SLAC   
 *******************************************************************************/

#include <assert.h>
#include <math.h>
#include <iostream>
#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtLASSAmp.hh"
#include "EvtGenBase/EvtDalitzCoord.hh"
#include "EvtGenBase/EvtdFunction.hh"
#include "EvtGenBase/EvtCyclic3.hh"
using std::endl;
using EvtCyclic3::Index;
using EvtCyclic3::Pair;


EvtLASSAmp::EvtLASSAmp(EvtDalitzPlot *dp, 
		       EvtCyclic3::Pair pair,
		       double m0, double g0,
		       double a, double r, double cutoff, std::string subtype) :
  EvtAmplitude<EvtDalitzPoint>(),
  _pair(pair),
  _m0(m0),
  _g0(g0),
  _r(r),
  _a(a),
  _cutoff(cutoff),
  _subtype(subtype)
{
  _dalitzSpace = dp;
  double ma = dp->m( first(pair) );
  double mb = dp->m( second(pair) );
  double E0a  = 0.5 * (_m0*_m0 + ma*ma - mb*mb) / _m0;
  _q0 = E0a*E0a - ma*ma;
  assert(_q0>0);
  _q0 = sqrt(_q0);
}



EvtLASSAmp::EvtLASSAmp(const EvtLASSAmp& other) :
  EvtAmplitude<EvtDalitzPoint>(other),
  _pair(other._pair),
  _m0(other._m0),
  _g0(other._g0),
  _q0(other._q0),
  _r(other._r),
  _a(other._a),
  _cutoff(other._cutoff),
  _subtype(other._subtype)

{  
  _dalitzSpace = other._dalitzSpace;
}


EvtLASSAmp::~EvtLASSAmp() {}



EvtComplex 
EvtLASSAmp::amplitude(const EvtDalitzPoint &dalitzPoint) const {
  
  /*

    Parameterization of Kpi S-wave using LASS scattering data.
    - Nucl.Phys.B296, 493 (1988)
    - W.Dunwoodie,http://www.slac.stanford.edu/~wmd/kpi_swave/kpi_swave_fit.note

            m                                     m0^2*Gamma0/q0
    ----------------- + exp(2*i*delta) * --------------------------------
    q*cot(delta)-i*q                     m0^2-m^2 - i*m0*Gamma0*q/m*m0/q0


    where q = momentum of K or pi in Kpi system
          
          q*cot(delta) = 1/ a   + 1/2 * [ r * q**2 ]

	  a = scattering length

	  r = effective range

  */

  double s = dalitzPoint.q(_pair);
  double m = sqrt(s);
  double q = dalitzPoint.p(first(_pair), _pair);

  // elastic scattering 
  double qcotd = 1./_a + 0.5*_r*q*q;
  EvtComplex lass_elastic = m<_cutoff ? m / ( qcotd - EvtComplex(0,q) ) : 0;

  // relative phase
  double cosd=1;
  double sind=0;
  if (q>0) {
    cosd = qcotd*qcotd/(q*q);
    cosd = sqrt( cosd/(1+cosd) );
    sind = sqrt( 1-cosd*cosd );
  }
  EvtComplex lass_phase( cosd, sind);
  lass_phase *= lass_phase;

  // K*(1430)
  double gamma = _g0 * q/m * _m0/_q0;
  EvtComplex lass_Kstar = (_m0*_m0)*(_g0/_q0)/(_m0*_m0-m*m-EvtComplex(0.,_m0*gamma));   

  EvtComplex theAmplitude(0.0, 0.0);

  if (_subtype == "LASS_ELASTIC") {

    theAmplitude = lass_elastic;

  } else if (_subtype == "LASS_RESONANT") {

    theAmplitude = lass_phase * lass_Kstar;

  } else {

    theAmplitude = lass_phase * lass_Kstar + lass_elastic;

  }

  return theAmplitude;
 
}
