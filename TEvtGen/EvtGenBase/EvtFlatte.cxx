//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 1998      Caltech, UCSB
//
// Module: EvtFlatte.cc
//
// Description: resonance-defining class 
//
// Modification history:
//
//    ponyisi  18 Feb 2008  created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <math.h>
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtKine.hh"
#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtFlatte.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtConst.hh"

EvtFlatte::~EvtFlatte(){}

//operator

EvtFlatte& EvtFlatte::operator = ( const EvtFlatte  &n)
{
  if ( &n == this ) return *this;
  _p4_p = n._p4_p;
  _p4_d1 = n._p4_d1;
  _p4_d2 = n._p4_d2;
  _ampl = n._ampl;
  _theta = n._theta;
  _mass = n._mass;
  _params = n._params;
  //  _m1a = n._m1a;
  //  _m1b = n._m1b;
  //  _g1 = n._g1;
  //  _m2a = n._m2a;
  //  _m2b = n._m2b;
  //  _g2 = n._g2;
  return  *this;
}

//constructor
 
EvtFlatte::EvtFlatte(const EvtVector4R& p4_p, const EvtVector4R& p4_d1,
                     const  EvtVector4R& p4_d2, double ampl,
                     double theta, double mass, 
                     vector<EvtFlatteParam>& params
		     //                   double m1a, double m1b, double g1,
		     //                   double m2a, double m2b, double g2
                     ): 
  _p4_p(p4_p),_p4_d1(p4_d1), _p4_d2(p4_d2), _ampl(ampl), _theta(theta), 
  _mass(mass),
  _params(params)
  //  _m1a(m1a), _m1b(m1b), _g1(g1),
  //  _m2a(m2a), _m2b(m2b), _g2(g2)
{}

//amplitude function

EvtComplex EvtFlatte::resAmpl() {
 
  double pi180inv = 1.0/EvtConst::radToDegrees;

  //   EvtComplex ampl(cos(_theta*pi180inv), sin(_theta*pi180inv));
  //   ampl *= _ampl;

  // SCALARS ONLY
  double mR = (_p4_d1+_p4_d2).mass();

  EvtComplex w;

  for (vector<EvtFlatteParam>::const_iterator param = _params.begin();
       param != _params.end();
       ++param) {

    double m1 = (*param).m1(); double m2 = (*param).m2(); 
    double g = (*param).g();
    w += (g*g
	  *sqrtCplx((1-((m1-m2)*(m1-m2))/(mR*mR))*
		    (1-((m1+m2)*(m1+m2))/(mR*mR))));
    //     cout << m1 << " " << mR << " " << w << endl;
  }

  EvtComplex denom = _mass*_mass - mR*mR - EvtComplex(0,1)*w;
  EvtComplex ampl = _ampl*EvtComplex(cos(_theta*pi180inv), sin(_theta*pi180inv))/denom;
  //  cout << abs(1/denom) << endl;
  return ampl;
}


