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
// Module: EvtGen/EvtYmSToYnSpipiCLEO.hh
//
// Description: This model is based on matrix element method used by
//              CLEO in Phys.Rev.D76:072001,2007 to model the dipion mass
//              and helicity angle distribution in the decays Y(mS) -> pi pi Y(nS),
//              where m,n are integers and m>n and m<4.
//              This model has two parameters, Re(B/A) and Im(B/A), which
//              are coefficients of the matrix element's terms determined by
//              the CLEO fits.
//
// Example:
//
// Decay  Upsilon(3S)
//  1.0000    Upsilon      pi+  pi-     YMSTOYNSPIPICLEO -2.523 1.189;
// Enddecay
// Decay  Upsilon(3S)
//  1.0000    Upsilon(2S)  pi+  pi-     YMSTOYNSPIPICLEO -0.395 0.001;
// Enddecay
// Decay  Upsilon(2S)
//  1.0000    Upsilon      pi+  pi-     YMSTOYNSPIPICLEO -0.753 0.000;
// Enddecay
//
//   --> the order of parameters is: Re(B/A) Im(B/A)
//
// Modification history:
//
//    SEKULA  Jan. 28, 2008         Module created
//
//------------------------------------------------------------------------


#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenModels/EvtYmSToYnSpipiCLEO.hh"
#include <string>
using std::endl;

EvtYmSToYnSpipiCLEO::~EvtYmSToYnSpipiCLEO() {}

std::string EvtYmSToYnSpipiCLEO::getName(){

  return "YMSTOYNSPIPICLEO";     

}


EvtDecayBase* EvtYmSToYnSpipiCLEO::clone(){

  return new EvtYmSToYnSpipiCLEO;

}

void EvtYmSToYnSpipiCLEO::init(){

  static EvtId PIP=EvtPDL::getId("pi+");
  static EvtId PIM=EvtPDL::getId("pi-");
  static EvtId PI0=EvtPDL::getId("pi0");

  // check that there are 2 arguments
  checkNArg(2);
  checkNDaug(3);

  checkSpinParent(EvtSpinType::VECTOR);
  checkSpinDaughter(0,EvtSpinType::VECTOR);



  if ((!(getDaug(1)==PIP&&getDaug(2)==PIM))&&
      (!(getDaug(1)==PI0&&getDaug(2)==PI0))) {
    report(Severity::Error,"EvtGen") << "EvtYmSToYnSpipiCLEO generator expected "
                           << " pi+ and pi- (or pi0 and pi0) "
			   << "as 2nd and 3rd daughter. "<<endl;
    report(Severity::Error,"EvtGen") << "Will terminate execution!"<<endl;
    ::abort();
  }

}

void EvtYmSToYnSpipiCLEO::initProbMax() {
  setProbMax(2.0);
}      

void EvtYmSToYnSpipiCLEO::decay( EvtParticle *p){


  // We want to simulate the following process:
  //
  // Y(mS) -> Y(nS) X, X -> pi+ pi- (pi0 pi0)
  //
  // The CLEO analysis assumed such an intermediate process
  // were occurring, and wrote down the matrix element
  // and its components according to this assumption. 
  //
  //


  double ReB_over_A = getArg(0);
  double ImB_over_A = getArg(1);

  p->makeDaughters(getNDaug(),getDaugs());
  EvtParticle *v,*s1,*s2;
  v=p->getDaug(0);
  s1=p->getDaug(1);
  s2=p->getDaug(2);

  double m_pi = s1->getP4().mass();
  double M_mS = p->getP4().mass();
  double M_nS = v->getP4().mass();

// //   report(Severity::Info,"EvtYmSToYnSpipiCLEO")  << "M_nS = " << v->getP4().mass() << endl;

  EvtVector4R P_nS;
  EvtVector4R P_pi1;
  EvtVector4R P_pi2;

  // Perform a simple accept/reject until we get a configuration of the
  // dipion system that passes
  bool acceptX = false;

  while( false == acceptX ) 
    {

      // Begin by generating a random X mass between the kinematic
      // boundaries, 2*m_pi and M(mS) - M(nS)

      double mX = EvtRandom::Flat(2.0 * m_pi, M_mS-M_nS);

      //   report(Severity::Info,"EvtYmSToYnSpipiCLEO")  << "m_X = " << mX << endl;

      // Now create a two-body decay from the Y(mS) in its rest frame
      // of Y(mS) -> Y(nS) + X

      double masses[2];
      masses[0] = M_nS;
      masses[1] = mX;

      EvtVector4R p4[2];

      EvtGenKine::PhaseSpace( 2, masses, p4, M_mS );

      P_nS = p4[0];
      EvtVector4R P_X  = p4[1];

      // Now create the four-vectors for the two pions in the X
      // rest frame, X -> pi pi
  
      masses[0] = s1->mass();
      masses[1] = s2->mass();

      EvtGenKine::PhaseSpace( 2, masses, p4, P_X.mass() );

      // compute cos(theta), the polar helicity angle between a pi+ and
      // the direction opposite the Y(mS) in the X rest frame. If the pions are pi0s, then
      // choose the one where cos(theta) = [0:1].

      EvtVector4R P_YmS_X = boostTo(p->getP4(), P_X);
      double costheta = - p4[0].dot(P_YmS_X)/(p4[0].d3mag()*P_YmS_X.d3mag());
      if (EvtPDL::name(s1->getId()) == "pi0") {
	if (costheta < 0) {
	  costheta = - p4[1].dot(P_YmS_X)/(p4[1].d3mag()*P_YmS_X.d3mag());
	}
      }
      if (EvtPDL::name(s1->getId()) == "pi-") {
	costheta = - p4[1].dot(P_YmS_X)/(p4[1].d3mag()*P_YmS_X.d3mag());
      }
  
      // //   report(Severity::Info,"EvtYmSToYnSpipiCLEO")  << "cos(theta) = " << costheta << endl;
    
    

      // Now boost the pion four vectors into the Y(mS) rest frame
      P_pi1 = boostTo(p4[0],P_YmS_X);
      P_pi2 = boostTo(p4[1],P_YmS_X);

      // Use a simple accept-reject to test this dipion system

      // Now compute the components of the matrix-element squared
      //
      // M(x,y)^2 = Q(x,y)^2 + |B/A|^2 * E1E2(x,y)^2 + 2*Re(B/A)*Q(x,y)*E1E2(x,y)
      //
      // x=m_pipi^2 and y = cos(theta), and where 
      //
      //   Q(x,y) = (x^2 + 2*m_pi^2)
      //  
      //   E1E2(x,y) = (1/4) * ( (E1 + E2)^2 - (E2 - E1)^2_max * cos(theta)^2 )
      //
      // and E1 + E2 = M_mS - M_nS and (E2 - E1)_max is the maximal difference
      // in the energy of the two pions allowed for a given mX value.
      //

      double Q    = (mX*mX - 2.0 * m_pi * m_pi);

      double deltaEmax = 
	- 2.0 * 
	sqrt( P_nS.get(0)*P_nS.get(0) - M_nS*M_nS ) *
	sqrt( 0.25 - pow(m_pi/mX,2.0));

      double sumE = (M_mS*M_mS - M_nS*M_nS + mX*mX)/(2.0 * M_mS);

      double E1E2 = 0.25 * ( pow(sumE, 2.0) - pow( deltaEmax * costheta, 2.0) );

      double M2 = Q*Q + (pow(ReB_over_A,2.0) + pow(ImB_over_A,2.0)) * E1E2*E1E2 + 2.0 * ReB_over_A * Q * E1E2;

      // phase space factor
      //
      // this is given as d(PS) = C * p(*)_X * p(X)_{pi+} * d(cosTheta) * d(m_X)
      // 
      // where C is a normalization constant, p(*)_X is the X momentum magnitude in the
      // Y(mS) rest frame, and p(X)_{pi+} is the pi+/pi0 momentum in the X rest frame
      //
  
      double dPS = 
	sqrt( (M_mS*M_mS - pow(M_nS + mX,2.0)) * (M_mS*M_mS - pow(M_nS - mX,2.0)) ) * // p(*)_X
	sqrt(mX*mX - 4*m_pi*m_pi); // p(X)_{pi}

      // the double-differential decay rate dG/(dcostheta dmX)
      double dG = M2 * dPS;

      // Throw a uniform random number from 0 --> probMax and do accept/reject on this
      
      double rnd = EvtRandom::Flat(0.0,getProbMax(0.0));

      if (rnd < dG)
	acceptX = true;

    }


  // initialize the daughters
  v->init(  getDaugs()[0], P_nS);
  s1->init( getDaugs()[1], P_pi1);
  s2->init( getDaugs()[2], P_pi2); 

//   report(Severity::Info,"EvtYmSToYnSpipiCLEO")  << "M_nS = " << v->getP4().mass() << endl;
//   report(Severity::Info,"EvtYmSToYnSpipiCLEO")  << "m_pi = " << s1->getP4().mass() << endl;
//   report(Severity::Info,"EvtYmSToYnSpipiCLEO")  << "m_pi = " << s2->getP4().mass() << endl;
//   report(Severity::Info,"EvtYmSToYnSpipiCLEO")  << "M2 = "   << M2 << endl;
  
  // Pass the polarization of the parent Upsilon
  EvtVector4C ep0,ep1,ep2;  
  
  ep0=p->eps(0);
  ep1=p->eps(1);
  ep2=p->eps(2);


  vertex(0,0,(ep0*v->epsParent(0).conj()));
  vertex(0,1,(ep0*v->epsParent(1).conj()));
  vertex(0,2,(ep0*v->epsParent(2).conj()));
  
  vertex(1,0,(ep1*v->epsParent(0).conj()));
  vertex(1,1,(ep1*v->epsParent(1).conj()));
  vertex(1,2,(ep1*v->epsParent(2).conj()));
  
  vertex(2,0,(ep2*v->epsParent(0).conj()));
  vertex(2,1,(ep2*v->epsParent(1).conj()));
  vertex(2,2,(ep2*v->epsParent(2).conj()));


  return ;

}




