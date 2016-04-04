//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information:
//      Copyright (C) 1998      Caltech, UCSB
//
// Module: EvtBaryonVminusAFF.cc
//
// Description: Routine to implement semileptonic form factors
//              according to the model BaryonVminusA
//
// Modification history:
//
//    R.J. Tesarek     May 28, 2004     Module created
//    Karen Gibson     1/20/2006        Module updated for 1/2+->1/2+,
//                                      1/2+->1/2-, 1/2+->3/2- Lambda decays
//
//--------------------------------------------------------------------------
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtBaryonPCRFF.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtIdSet.hh"
#include "EvtGenBase/EvtConst.hh"
#include <string>
#include <math.h>
#include <stdlib.h>
using std::endl;

void EvtBaryonPCRFF::getdiracff(EvtId parent, EvtId daught,
				double q2, double /* mass */ , 
				double *f1, double *f2, double *f3, 
				double *g1, double *g2, double *g3 ) {

  // Baryons (partial list 5/28/04)
  static EvtId LAMCP=EvtPDL::getId("Lambda_c+");
  static EvtId LAMCM=EvtPDL::getId("anti-Lambda_c-");
  static EvtId LAMC1P=EvtPDL::getId("Lambda_c(2593)+");
  static EvtId LAMC1M=EvtPDL::getId("anti-Lambda_c(2593)-");
  static EvtId LAMB=EvtPDL::getId("Lambda_b0");
  static EvtId LAMBB=EvtPDL::getId("anti-Lambda_b0");

  double F1, F2, F3, G1, G2, G3;

  if( parent==LAMB || parent==LAMBB ) {
    // Implement constituent quark model form factors predicted 
    // by M. Pervin, W. Roberst, and S. Capstick, Phys. Rev. C72, 035201 (2005)

    if( daught==LAMCP|| daught==LAMCM ) {
      
      //  Parameters needed in the calculation;
      double mQ = 5.28; 
      double mq = 1.89;
      double md = 0.40;
      double MLamB = EvtPDL::getMass(parent);
      double MLamC = EvtPDL::getMass(daught);
      
      double aL  = 0.59;
      double aLp = 0.55;
      
      double aL2  = aL*aL;
      double aLp2 = aLp*aLp;
      double aLLp2 = 0.5*(aL2+aLp2);
      
      // relativistic correction factor
      double k2 = 1.0;
      double rho2 = 3.*md*md/(2.*k2*aLLp2);  
      
      // w = scalar product of the 4 velocities of the Lb and Lc.
      double w = 0.5*(MLamB*MLamB + MLamC*MLamC - q2)/MLamB/MLamC;

      double I = pow(aL*aLp/aLLp2, 1.5)*exp(-rho2*(w*w-1.));

      // Calculate the form factors
      F1 = I*( 1.0 + (md/aLLp2)*( (aLp2/mq)+(aL2/mQ) ) );
      F2 = -I*( (md/mq)*(aLp2/aLLp2) - aL2*aLp2/(4.*aLLp2*mq*mQ) );
      F3 = -I*md*aL2/(mQ*aLLp2);

      G1 = I*( 1.0 - (aL2*aLp2)/(12.*aLLp2*mq*mQ) );
      G2 = -I*( md*aLp2/(mq*aLLp2) 
		+ (aL2*aLp2)/(12.*aLLp2*mq*mQ)*(1.+12.*md*md/aLLp2) );
      G3 = I*( md*aL2/(mQ*aLLp2) + md*md*aL2*aLp2/(mq*mQ*aLLp2*aLLp2) );
       
      // Set form factors to be passed to the amplitude calc.
      *f1 = F1;
      *f2 = F2;
      *f3 = F3;
      *g1 = G1;
      *g2 = G2;
      *g3 = G3;

    }

    else if( daught==LAMC1P || daught==LAMC1M ) {

      double mQ = 5.28; 
      double mq = 1.89;
      double md = 0.40;
      double MLamB = EvtPDL::getMass(parent);
      double MLamC = EvtPDL::getMass(daught);
      
      double aL  = 0.59;
      double aLp = 0.47;
      
      double aL2  = aL*aL;
      double aLp2 = aLp*aLp;
      double aLLp2 = 0.5*(aL2+aLp2);
      
      // relativistic correction factor
      double k2 = 1.0;
      double rho2 = 3.*md*md/(2.*k2*aLLp2);  

      // w = scalar product of the 4 velocities of the Lb and Lc.
      double w = 0.5*(MLamB*MLamB + MLamC*MLamC - q2)/MLamB/MLamC;
      
      double I = pow(aL*aLp/aLLp2, 2.5)*exp(-rho2*(w*w-1.));

      // Calculate the form factors
      F1 = I*aL/6.0*( 3.0/mq - 1.0/mQ );
      F2 = -I*( 2.0*md/aL - aL/(2.0*mq) + 2.*md*md*aL/(mQ*aLLp2) 
		- (md*aL/(6.*mq*mQ*aLLp2))*( 3.*aL2 - 2.*aLp2));
      F3 = I*2.*md*md*aL/(mQ*aLLp2);

      G1 = I*( 2.0*md/aL - aL/(6.*mQ) 
	       + (md*aL/(6.*mq*mQ*aLLp2))*( 3.*aL2 - 2.*aLp2));
      G2 = I*( -2.*md/aL + aL/(2.*mq) + aL/(3.*mQ) );
      G3 = I*aL/(3.*mQ)*( 1.0 - (md/(2.*mq*aLLp2))*( 3.*aL2 - 2.*aLp2));
       
      // Set form factors to be passed to the amplitude calc.
      *f1 = F1;
      *f2 = F2;
      *f3 = F3;
      *g1 = G1;
      *g2 = G2;
      *g3 = G3;
    }
  }

  else {
    *f1 = 1.0;
    *f2 = 1.0;
    *f3 = 0.0;
    *g1 = 1.0;
    *g2 = 1.0;
    *g3 = 0.0;
  }

  return ;
}


void EvtBaryonPCRFF::getraritaff( EvtId parent, EvtId daught,
				  double q2, double /* mass */, 
				  double *f1, double *f2, double *f3, double *f4, 
				  double *g1, double *g2, double *g3, double *g4 ) {

  // Baryons (partial list 5/28/04)
  static EvtId LAMB=EvtPDL::getId("Lambda_b0");
  static EvtId LAMBB=EvtPDL::getId("anti-Lambda_b0");
  static EvtId LAMC2P=EvtPDL::getId("Lambda_c(2625)+");
  static EvtId LAMC2M=EvtPDL::getId("anti-Lambda_c(2625)-");

  double F1, F2, F3, F4, G1, G2, G3, G4;

  if( parent==LAMB || parent==LAMBB ) {
    // Implement constituent quark model form factors predicted 
    // by M. Pervin, W. Roberst, and S. Capstick, Phys. Rev. C72, 035201 (2005)

    if( daught==LAMC2P|| daught==LAMC2M ) {
      
      double mQ = 5.28; 
      double mq = 1.89;
      double md = 0.40;
      double MLamB = EvtPDL::getMass(parent);
      double MLamC = EvtPDL::getMass(daught);
      
      double aL  = 0.59;
      double aLp = 0.47;
      
      double aL2  = aL*aL;
      double aLp2 = aLp*aLp;
      double aLLp2 = 0.5*(aL2+aLp2);
      
      // relativistic correction factor
      double k2 = 1.0;
      double rho2 = 3.*md*md/(2.*k2*aLLp2);  

      // w = scalar product of the 4 velocities of the Lb and Lc.
      double w = 0.5*(MLamB*MLamB + MLamC*MLamC - q2)/MLamB/MLamC;
      
      double I = -(1./sqrt(3.))*pow(aL*aLp/aLLp2, 2.5)*exp(-rho2*(w*w-1.));
    
      // Calculate the form factors
      F1 = I*3.0*md/aL*( 1.0 + (md/aLLp2)*( (aLp2/mq)+(aL2/mQ) ) );
      F2 = -I*( (3.*md*md/mq)*(aLp2/(aLLp2*aL2)) - 5.*aL*aLp2*md/(4.*aLLp2*mq*mQ) );
      F3 = -I*( 3.*md*md*aL/(mQ*aLLp2) + aL/(2.*mQ) );
      F4 = I*aL/mQ;

      G1 = I*( 3.0*md/aL - (aL/(2.*mQ))*(1. + 3.*md*aLp2/(2.*aLLp2*mq) ) );
      G2 = -I*( (3.*md*md/mq)*(aLp2/(aLLp2*aL)) + aL*aLp2*md/(4.*aLLp2*aLLp2*mq*mQ)*(aLLp2+12.*md*md) );
      G3 = I*aL/(mQ*aLLp2)*( aLLp2/2. + 3.*md*md + aLp2*md/(mq*aLLp2)*(aLLp2+6.*md*md) );
      G4 = -I*( aL/mQ + md/(mq*mQ)*aLp2*aL/aLLp2 );
       
      // Set form factors to be passed to the amplitude calc.
      *f1 = F1;
      *f2 = F2;
      *f3 = F3;
      *f4 = F4;
      *g1 = G1;
      *g2 = G2;
      *g3 = G3;
      *g4 = G4;
    }
  }
  
  else {
    *f1 = 1.0;
    *f2 = 1.0;
    *f3 = 0.0;
    *f4 = 0.0;
    *g1 = 1.0;
    *g2 = 1.0;
    *g3 = 0.0;
    *g4 = 0.0;

  }

  return ;


}

void EvtBaryonPCRFF::getscalarff(EvtId, EvtId, double, double, double*, double*) {

  report(Severity::Error,"EvtGen") << "Not implemented :getscalarff in EvtBaryonPCRFF.\n";  
  ::abort();

}

void EvtBaryonPCRFF::getvectorff(EvtId, EvtId, double, double, double*, double*,
				 double*, double*) {

  report(Severity::Error,"EvtGen") << "Not implemented :getvectorff in EvtBaryonPCRFF.\n";  
  ::abort();

}

void EvtBaryonPCRFF::gettensorff(EvtId, EvtId, double, double, double*, double*,
				 double*, double*) {

  report(Severity::Error,"EvtGen") << "Not implemented :gettensorff in EvtBaryonPCRFF.\n";  
  ::abort();

}

void EvtBaryonPCRFF::getbaryonff(EvtId, EvtId, double, double, double*, 
				 double*, double*, double*){
  
  report(Severity::Error,"EvtGen") << "Not implemented :getbaryonff in EvtBaryonPCRFF.\n";  
  ::abort();

}
