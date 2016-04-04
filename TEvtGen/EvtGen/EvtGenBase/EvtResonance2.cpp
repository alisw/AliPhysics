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
// Module: EvtResonance2.cc
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
#include "EvtGenBase/EvtResonance2.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtConst.hh"

EvtResonance2::~EvtResonance2(){}

EvtResonance2& EvtResonance2::operator = ( const EvtResonance2  &n)
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
  _invmass_angdenom = n._invmass_angdenom;
   return  *this;
}

 
EvtResonance2::EvtResonance2(const EvtVector4R& p4_p, const EvtVector4R& p4_d1,
			     const  EvtVector4R& p4_d2, double ampl, 
			     double theta, double gamma, double bwm, int spin,
			     bool invmass_angdenom): 
  _p4_p(p4_p),_p4_d1(p4_d1), _p4_d2(p4_d2),_ampl(ampl), _theta(theta), 
  _gamma(gamma), _bwm(bwm), _spin(spin), _invmass_angdenom(invmass_angdenom) {}


EvtComplex EvtResonance2::resAmpl() {
 
  double pi180inv = 1.0/EvtConst::radToDegrees;

  EvtComplex ampl;
  EvtVector4R  p4_d3 = _p4_p-_p4_d1-_p4_d2;

  //get cos of the angle between the daughters from their 4-momenta
  //and the 4-momentum of the parent

  //in general, EvtDecayAngle(parent, part1+part2, part1) gives the angle
  //the missing particle (not listed in the arguments) makes
  //with part2 in the rest frame of both
  //listed particles (12)
 
  //angle 3 makes with 2 in rest frame of 12 (CS3)  
  //double cos_phi_0 = EvtDecayAngle(_p4_p, _p4_d1+_p4_d2, _p4_d1);
  //angle 3 makes with 1 in 12 is, of course, -cos_phi_0

  //first compute several quantities...follow CLEO preprint 00-23

  double mAB=(_p4_d1+_p4_d2).mass();
  double mBC=(_p4_d2+p4_d3).mass();
  double mAC=(_p4_d1+p4_d3).mass();
  double mA=_p4_d1.mass(); 
  double mB=_p4_d2.mass(); 
  double mD=_p4_p.mass();
  double mC=p4_d3.mass();
  
  double mR=_bwm;
  double gammaR=_gamma;
  double mdenom = _invmass_angdenom ? mAB : mR;
  double pAB=sqrt( (((mAB*mAB-mA*mA-mB*mB)*(mAB*mAB-mA*mA-mB*mB)/4.0) -
		    mA*mA*mB*mB)/(mAB*mAB));
  double pR=sqrt( (((mR*mR-mA*mA-mB*mB)*(mR*mR-mA*mA-mB*mB)/4.0) -
		   mA*mA*mB*mB)/(mR*mR));

  double pD= (((mD*mD-mR*mR-mC*mC)*(mD*mD-mR*mR-mC*mC)/4.0) -
		   mR*mR*mC*mC)/(mD*mD);
  if ( pD>0 ) { pD=sqrt(pD); } else {pD=0;}
  double pDAB=sqrt( (((mD*mD-mAB*mAB-mC*mC)*(mD*mD-mAB*mAB-mC*mC)/4.0) -
		   mAB*mAB*mC*mC)/(mD*mD));



  double fR=1;
  double fD=1;
  int power=0;
  switch (_spin) {
  case 0:
    fR=1.0;
    fD=1.0;
    power=1;
    break;
  case 1:
    fR=sqrt(1.0+1.5*1.5*pR*pR)/sqrt(1.0+1.5*1.5*pAB*pAB);
    fD=sqrt(1.0+5.0*5.0*pD*pD)/sqrt(1.0+5.0*5.0*pDAB*pDAB);
    power=3;
    break;
  case 2:
    fR = sqrt( (9+3*pow((1.5*pR),2)+pow((1.5*pR),4))/(9+3*pow((1.5*pAB),2)+pow((1.5*pAB),4)) );
    fD = sqrt( (9+3*pow((5.0*pD),2)+pow((5.0*pD),4))/(9+3*pow((5.0*pDAB),2)+pow((5.0*pDAB),4)) );
    power=5;
    break;
  default:
    report(Severity::Info,"EvtGen") << "Incorrect spin in EvtResonance22.cc\n";
  }
  
  double gammaAB= gammaR*pow(pAB/pR,power)*(mR/mAB)*fR*fR;
  switch (_spin) {
  case 0:
    ampl=_ampl*EvtComplex(cos(_theta*pi180inv),sin(_theta*pi180inv))*
          fR*fD/(mR*mR-mAB*mAB-EvtComplex(0.0,mR*gammaAB));
    break;
  case 1:
    ampl=_ampl*EvtComplex(cos(_theta*pi180inv),sin(_theta*pi180inv))*
      (fR*fD*(mAC*mAC-mBC*mBC+((mD*mD-mC*mC)*(mB*mB-mA*mA)/(mdenom*mdenom)))/
       (mR*mR-mAB*mAB-EvtComplex(0.0,mR*gammaAB)));
    break;
  case 2:
    ampl=_ampl*EvtComplex(cos(_theta*pi180inv),sin(_theta*pi180inv))*
      fR*fD/(mR*mR-mAB*mAB-EvtComplex(0.0,mR*gammaAB))*
      (pow((mBC*mBC-mAC*mAC+(mD*mD-mC*mC)*(mA*mA-mB*mB)/(mdenom*mdenom)),2)-
       (1.0/3.0)*(mAB*mAB-2*mD*mD-2*mC*mC+pow((mD*mD- mC*mC)/mdenom, 2))*
       (mAB*mAB-2*mA*mA-2*mB*mB+pow((mA*mA-mB*mB)/mdenom,2))); 
  break;

  default:
    report(Severity::Info,"EvtGen") << "Incorrect spin in EvtResonance22.cc\n";
  }

  return ampl;
}



