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
// Module: EvtBCVFF.cc
//
// Description: form factors for B->Vlnu 
//
// Modification history:
//
//    AVL Jul 6, Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtBCVFF.hh"

#include "EvtGenBase/EvtId.hh"
#include <string>
#include <math.h>
#include "EvtGenBase/EvtPDL.hh"
#include <stdlib.h>
#include <iostream>

using namespace std;

EvtBCVFF::EvtBCVFF(int idV, int fit) {

  idVector = idV;
  whichfit = fit;
  //cout<<"==== EvtBCVFF:: idVector="<<idVector<<" whichfit="<<whichfit<<endl;
  return;
}

void EvtBCVFF::getvectorff(EvtId,EvtId,
			   double t, double, double *a1f,
			   double *a2f, double *vf, double *a0f ){
  double q2=t;

  if(whichfit == 0) {
    *vf = 0;
    *a0f = 0;
    *a1f = 1;
    *a2f=0;
    
    return;
  };

  if( idVector == EvtPDL::getId("J/psi").getId() ) { // Bc -> J/psi
      if(whichfit == 1) { // SR form factor set from [Kiselev, hep-ph/0211021]
	double Mbc = 6.277, Mpsi=3.0967; // Experimental values
	double Mpole2 = 4.5*4.5, den = 1./(1.-q2/Mpole2);
	double FV = 0.11*den,
	FAp = -0.074*den,
	FA0 = 5.9*den,
	FAm = 0.12*den;
	*vf = (Mbc + Mpsi)*FV;
	*a2f = -(Mbc+Mpsi)*FAp;
	*a1f = FA0/(Mbc+Mpsi);
	*a0f = (q2*FAm + (Mbc+Mpsi)*(*a1f)-(Mbc-Mpsi)*(*a2f))/(2*Mpsi);    
	return;
      }
      else if(whichfit == 2) {  // form factor set from  [Ebert, hep-ph/0306306] 
	*vf = (0.49077824756158533 - 0.0012925655191347828*q2)/(1 - 0.06292520325875656*q2);
	*a0f = (0.4160345034630221 - 0.0024720095310225023*q2)/(1 - 0.061603451915567785*q2);
	*a1f = (0.4970212860605933 - 0.0067519730024654745*q2)/(1 - 0.050487026667172176*q2);
	*a2f = (0.7315284919705497 + 0.0014263826220727142*q2 -  0.0006946090066269195*q2*q2)/(1 - 0.04885587273651653*q2);
	return;
      };
  }
  else if(idVector == EvtPDL::getId("psi(2S)").getId()) { // Bc -> psi((2S)
      if(whichfit == 1) {
	////cout<<"BC2:: psi2S, Kiselev, q2="<<q2<<endl;
	double Mbc = 6.277, Mpsi=3.0967, Mpsi2S = 3.686, kappa = Mpsi/Mpsi2S; // Experimental values
	double Mpole2 = 4.5*4.5, den=1./(1.-q2/Mpole2);
	double FV = 0.11*den*kappa/3.1,
	FAp = -0.074*den*kappa/4.9,
	FA0 = 5.9*den*kappa/3.5,
	FAm = 0.12*den*kappa/2.3;
	*vf = (Mbc + Mpsi2S)*FV;
	*a2f = -(Mbc+Mpsi2S)*FAp;
	*a1f = FA0/(Mbc+Mpsi2S);
	*a0f = (q2*FAm + (Mbc+Mpsi2S)*(*a1f)-(Mbc-Mpsi2S)*(*a2f))/(2*Mpsi2S);  
	return;
      }
      else if(whichfit == 2) {
	////cout<<"BC2:: psi2S, Ebert, q2="<<q2<<endl;
	*vf  =  (0.24177223968739653 - 0.053589051007278135*q2)/(1 - 0.0977848994260899*q2);
	*a0f = (0.23996026570086615 - 0.03530198514007337*q2)/(1 - 0.09371162519983989*q2);
	*a1f = (0.17418379258849329 - 0.004129699022085851*q2*q2)/(1 + 0.06607665248402918*q2);
	*a2f = (0.1352376939112041 - 0.040361722565209444*q2 + 0.003343515369431853*q2*q2)/(1 - 0.1463698128333418*q2);
	return;
      };
  }
  else {
  report(Severity::Error,"EvtGen") << "Not implemented :getbaryonff in EvtBCVFF.\n";  
    ::abort();
  };
}


void EvtBCVFF::getscalarff(EvtId, EvtId, double, double, double*, 
			       double*){
  
  report(Severity::Error,"EvtGen") << "Not implemented :getbaryonff in EvtBCVFF.\n";  
  ::abort();

}



void EvtBCVFF::gettensorff(EvtId, EvtId, double, double, double*, 
			       double*, double*, double*){
  
  report(Severity::Error,"EvtGen") << "Not implemented :getbaryonff in EvtBCVFF.\n";  
  ::abort();

}



void EvtBCVFF::getbaryonff(EvtId, EvtId, double, double, double*, 
			       double*, double*, double*){
  
  report(Severity::Error,"EvtGen") << "Not implemented :getbaryonff in EvtBCVFF.\n";  
  ::abort();

}

void EvtBCVFF::getdiracff(EvtId, EvtId, double, double, double*, double*,
			       double*, double*, double*, double*) {
  
  report(Severity::Error,"EvtGen") << "Not implemented :getdiracff in EvtBCVFF.\n";
  ::abort();

}

void EvtBCVFF::getraritaff(EvtId, EvtId, double, double, double*, double*, 
				double*, double*, double*, double*, double*, double*) {
  
  report(Severity::Error,"EvtGen") << "Not implemented :getraritaff in EvtBCVFF.\n";
  ::abort();

}

