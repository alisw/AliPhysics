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
// Module: EvtResonance.cc
//
// Description: resonance-defining class 
//
// Modification history:
//
//    NK        September 4, 1997      Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <math.h>
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtKine.hh"
#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtResonance.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtConst.hh"
using std::endl;

EvtResonance::~EvtResonance(){}


EvtResonance& EvtResonance::operator = ( const EvtResonance  &n)
{
  if ( &n == this ) return *this;
  _p4_p = n._p4_p;
  _p4_d1 = n._p4_d1;
  _p4_d2 = n._p4_d2;
  _ampl = n._ampl;
  _theta = n._theta;
  _gamma = n._gamma;
  _spin = n._spin;
  _bwm = n._bwm;
   return  *this;
}

 
EvtResonance::EvtResonance(const EvtVector4R& p4_p, const EvtVector4R& p4_d1,
			   const  EvtVector4R& p4_d2, double ampl, 
			   double theta, double gamma, double bwm, int spin): 
  _p4_p(p4_p),_p4_d1(p4_d1), _p4_d2(p4_d2),_ampl(ampl), _theta(theta), 
  _gamma(gamma), _bwm(bwm), _spin(spin) {}

EvtComplex EvtResonance::resAmpl() {
 
  double pi180inv = 1.0/EvtConst::radToDegrees;

  EvtComplex ampl;
  //EvtVector4R  _p4_d3 = _p4_p-_p4_d1-_p4_d2;

  //get cos of the angle between the daughters from their 4-momenta
  //and the 4-momentum of the parent

  //in general, EvtDecayAngle(parent, part1+part2, part1) gives the angle
  //the missing particle (not listed in the arguments) makes
  //with part2 in the rest frame of both
  //listed particles (12)
 
  //angle 3 makes with 2 in rest frame of 12 (CS3)  
  double cos_phi_0 = EvtDecayAngle(_p4_p, _p4_d1+_p4_d2, _p4_d1);
  //angle 3 makes with 1 in 12 is, of course, -cos_phi_0

  switch (_spin) {

  case 0 : 
    ampl=(_ampl*EvtComplex(cos(_theta*pi180inv),sin(_theta*pi180inv))*
	  sqrt(_gamma/EvtConst::twoPi)*
	  (1.0/((_p4_d1+_p4_d2).mass()-_bwm-EvtComplex(0.0,0.5*_gamma)))); 
    break;

  case 1 : 
    ampl=(_ampl*EvtComplex(cos(_theta*pi180inv),sin(_theta*pi180inv))*
	  sqrt(_gamma/EvtConst::twoPi)*
	  (cos_phi_0/((_p4_d1+_p4_d2).mass()-_bwm-EvtComplex(0.0,0.5*_gamma))));
    break;

  case 2: 
    ampl=(_ampl*EvtComplex(cos(_theta*pi180inv),sin(_theta*pi180inv))*
	  sqrt(_gamma/EvtConst::twoPi)*
	  ((1.5*cos_phi_0*cos_phi_0-0.5)/((_p4_d1+_p4_d2).mass()-_bwm-EvtComplex(0.0, 0.5*_gamma))));
    break;
             
  case 3:  
    ampl=(_ampl*EvtComplex(cos(_theta*pi180inv),sin(_theta*pi180inv))*
	  sqrt(_gamma/EvtConst::twoPi)*
	  ((2.5*cos_phi_0*cos_phi_0*cos_phi_0-1.5*cos_phi_0)/((_p4_d1+_p4_d2).mass()-_bwm-EvtComplex(0.0, 0.5*_gamma))));
    break;

  default:
    report(Severity::Debug,"EvtGen") << "EvtGen: wrong spin in EvtResonance" << endl;
    ampl = EvtComplex(0.0);
    break;         

  }

  return ampl;
}

EvtComplex EvtResonance::relBrWig(int i) {
    
//this function returns relativistic Breit-Wigner amplitude
//for a given resonance (for P-wave decays of scalars only at the moment!)

    EvtComplex BW;
    EvtVector4R  _p4_d3 = _p4_p-_p4_d1-_p4_d2;
    EvtVector4R _p4_12 = _p4_d1 + _p4_d2;

    double msq13 = (_p4_d1 + _p4_d3).mass2();
    double msq23 = (_p4_d2 + _p4_d3).mass2();
    double msqParent = _p4_p.mass2();
    double msq1 = _p4_d1.mass2();
    double msq2 = _p4_d2.mass2();
    double msq3 = _p4_d3.mass2();  

    double M;

    double p2 = sqrt((_p4_12.mass2() - (_p4_d1.mass() + _p4_d2.mass())*(_p4_d1.mass() + _p4_d2.mass()))*(_p4_12.mass2() - (_p4_d1.mass() - _p4_d2.mass())*(_p4_d1.mass() - _p4_d2.mass())))/(2.0*_p4_12.mass());
    
    double p2R = sqrt((_bwm*_bwm - (_p4_d1.mass() + _p4_d2.mass())*(_p4_d1.mass() + _p4_d2.mass()))*(_bwm*_bwm - (_p4_d1.mass() - _p4_d2.mass())*(_p4_d1.mass() - _p4_d2.mass())))/(2.0*_bwm);

    double gam, R;

    if (i == 1) {

	R = 2.0/(0.197);

    }
    else R = 5.0/(0.197);
  
    gam = _gamma*(_bwm/_p4_12.mass())*(p2/p2R)*(p2/p2R)*(p2/p2R)*((1 + R*R*p2R*p2R)/(1 + R*R*p2*p2));
    M = (msq13 - msq23 - (msqParent - msq3)*(msq1 - msq2)/(_bwm*_bwm))*sqrt((1 + R*R*p2R*p2R)/(1 + R*R*p2*p2)); 
    
    BW = sqrt(_gamma)*M/((_bwm*_bwm - _p4_12.mass2()) - EvtComplex(0.0,1.0)*gam*_bwm);
    
    return BW;

}


