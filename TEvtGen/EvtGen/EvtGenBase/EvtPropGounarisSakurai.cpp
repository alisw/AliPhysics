#include "EvtGenBase/EvtPatches.hh"
/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtPropGounarisSakurai.cpp,v 1.1 2009-03-16 16:47:03 robbep Exp $
 *  Author: Matt Graham 
 *  modified from EvtPropBreitWignerRel...this should be used for rho's
 *******************************************************************************/

#include <math.h>
#include "EvtGenBase/EvtPropGounarisSakurai.hh"


EvtPropGounarisSakurai::EvtPropGounarisSakurai(EvtDalitzPlot *dp, 
              EvtCyclic3::Pair pair, double m0, double g0) 
  : EvtPropagator(m0,g0),_pair(pair),_gbase(g0)
{
  _dalitzSpace = dp;
  _m1= dp->m(EvtCyclic3::first(_pair));
  _m2= dp->m(EvtCyclic3::second(_pair));
}


EvtPropGounarisSakurai::EvtPropGounarisSakurai(const EvtPropGounarisSakurai& other) 
  : EvtPropagator(other), _pair(other._pair), _gbase(other._gbase), 
    _m1(other._m1), _m2(other._m2)
{
_dalitzSpace = other._dalitzSpace;
}


EvtPropGounarisSakurai::~EvtPropGounarisSakurai() 
{}
  

EvtAmplitude<EvtPoint1D>* EvtPropGounarisSakurai::clone() const
{ 
  return new EvtPropGounarisSakurai(*this); 
}


EvtComplex EvtPropGounarisSakurai::amplitude(const EvtPoint1D& x) const
{
  double m = x.value();
  double s = m*m;
  double m2=_m0*_m0;
  double _width=_gbase;
  double _mass=_m0;

  double A = ( 1 + dFun( m2 )*_width/_mass );
  double B = s - m2 - fsFun( s );
  //  double C = sqrt(s)*_g0;//wrong!
  double C = sqrt(m2)*_g0;//correct!
  double D = B*B + C*C;

  EvtComplex rpt( A*B/D, - A*C/D );
  return rpt;


}

//  adapted from RhoPiTools
double EvtPropGounarisSakurai::fsFun( double s ) const
{
  double m2=_m0*_m0;

  EvtTwoBodyKine vd(_m1,_m2,sqrt(s));
  EvtTwoBodyKine vR(_m1,_m2,_m0);
  double k_s   = vd.p();
  double k_Am2 = vR.p();
  //  
  double f     = _gbase * m2 / pow( k_Am2, 3 )
                   * ( 
                       pow( k_s, 2 ) * (hFun( s ) - hFun( m2 ))
                       + (m2 - s) * pow( k_Am2, 2 ) * dh_dsFun( m2 )
                     );

  return f;
}


double EvtPropGounarisSakurai::hFun( double s ) const
{
  double sm    = _m1 + _m2;
  double SQRTs = sqrt(s);
  EvtTwoBodyKine vd(_m1,_m2,sqrt(s));
  double k_s   = vd.p();

  return 2/EvtConst::pi * (k_s/SQRTs) * log( (SQRTs + 2*k_s)/(sm) );
}

double EvtPropGounarisSakurai::dh_dsFun( double s ) const
{  
  EvtTwoBodyKine vd(_m1,_m2,sqrt(s));
  double k_s   = vd.p();
 
  return hFun(s) * ( 1/(8*pow( k_s, 2)) - 1/(2*s) ) + 1/(2*EvtConst::pi*s);
}

double EvtPropGounarisSakurai::dFun( double s ) const
{
  double sm   = _m1 + _m2;
  double sm24 = sm*sm/4;
  double m    = sqrt(s);
  EvtTwoBodyKine vd(_m1,_m2,sqrt(s));
  double k_m2   = vd.p();
  double _pi=EvtConst::pi;

  return 3.0/_pi * sm24/pow( k_m2, 2 ) * log( (m + 2*k_m2)/sm ) 
         + m/(2*_pi*k_m2) - sm24*m/(_pi * pow( k_m2, 3 ));
}

