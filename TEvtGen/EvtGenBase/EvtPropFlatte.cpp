#include "EvtGenBase/EvtPatches.hh"
/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 * Author : D. Dujmic, J. Thompson
 *
 * Copyright (C) 2005 SLAC
 *******************************************************************************/

#include <math.h>
#include "EvtGenBase/EvtPropFlatte.hh"

#include <iostream>
using std::cout;
using std::endl;

EvtPropFlatte::EvtPropFlatte(double m0, 
			     double g0, double m0a, double m0b, 
			     double g1, double m1a, double m1b) :
  EvtPropagator( m0, g0),
  _m0a(m0a),
  _m0b(m0b),
  _g1 (g1),
  _m1a(m1a),
  _m1b(m1b)
{}


EvtPropFlatte::EvtPropFlatte(const EvtPropFlatte& other) : 
  EvtPropagator(other),
  _m0a (other._m0a),
  _m0b (other._m0b),
  _g1  (other._g1),
  _m1a (other._m1a),
  _m1b (other._m1b)
{}


EvtPropFlatte::~EvtPropFlatte() 
{}
  

EvtAmplitude<EvtPoint1D>* EvtPropFlatte::clone() const
{ 
  return new EvtPropFlatte(*this); 
}



EvtComplex EvtPropFlatte::amplitude(const EvtPoint1D& x) const
{

  /*

  Use BES parameterization:

                        1.
      -----------------------------------------
       m0^2 - m^2 - i*m0*( g1*rho1 + g2*rho2 )

       
  Resonance mass: m0
  Channel1: m0a, m0b, g0
  Channel2: m1a, m1b, g1

  where breakup momenta q's are:

          E0a = (m^2 + m0a^2 - m0b^2) / 2m
          q0  = sqrt( E0a^2 - m0a^2 )

          E1a = (m^2 + m1a^2 - m1b^2) / 2m
          q1  = sqrt( E1a^2 - m1a^2 )


  */



  double s = x.value()*x.value();
  double m = x.value();
  
  double E0a  = 0.5 * (s + _m0a*_m0a - _m0b*_m0b) / m;
  double qSq0 = E0a*E0a - _m0a*_m0a;

  double E1a  = 0.5 * (s + _m1a*_m1a - _m1b*_m1b) / m;
  double qSq1 = E1a*E1a - _m1a*_m1a;

  EvtComplex gamma0 = qSq0 >= 0 ?  EvtComplex(  _g0 * sqrt(qSq0), 0)  : EvtComplex( 0, _g0 * sqrt(-qSq0) );
  EvtComplex gamma1 = qSq1 >= 0 ?  EvtComplex(  _g1 * sqrt(qSq1), 0)  : EvtComplex( 0, _g1 * sqrt(-qSq1) );

  EvtComplex gamma = gamma0 + gamma1;
  
  EvtComplex a = 1.0/( _m0*_m0 - s - EvtComplex(0.0,2*_m0/m)*gamma  );

  return a;
}

