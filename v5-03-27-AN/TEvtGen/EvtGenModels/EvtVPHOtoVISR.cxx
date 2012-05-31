//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2004      Cornell
//
// Module: EvtVPHOtoVISR.cc
//
// Description: Routine to decay vpho -> vector ISR photon
//
// Modification history:
//
//    Ryd       March 20, 2004       Module created
//
//------------------------------------------------------------------------
// 
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenModels/EvtVPHOtoVISR.hh"
#include <string>

EvtVPHOtoVISR::~EvtVPHOtoVISR() {}

std::string EvtVPHOtoVISR::getName(){

  return "VPHOTOVISR"; 
    
}


EvtDecayBase* EvtVPHOtoVISR::clone(){

  return new EvtVPHOtoVISR;

}

void EvtVPHOtoVISR::init(){

  // check that there are 0 or 2 arguments
  checkNArg(0,2);

  // check that there are 2 daughters
  checkNDaug(2);

  // check the parent and daughter spins
  checkSpinParent(EvtSpinType::VECTOR);
  checkSpinDaughter(0,EvtSpinType::VECTOR);
  checkSpinDaughter(1,EvtSpinType::PHOTON);
}

void EvtVPHOtoVISR::initProbMax() {

  //setProbMax(100000.0);

}      

void EvtVPHOtoVISR::decay( EvtParticle *p){

  //take photon along z-axis, either forward or backward.
  //Implement this as generating the photon momentum along 
  //the z-axis uniformly 

  double w=p->mass();
  double s=w*w;

  double L=2.0*log(w/0.000511);
  double alpha=1/137.0;
  double beta=(L-1)*2.0*alpha/EvtConst::pi;

  //This uses the fact that there is a daughter of the 
  //psi(3770)
  assert(p->getDaug(0)->getDaug(0)!=0);
  double md=EvtPDL::getMeanMass(p->getDaug(0)->getDaug(0)->getId());

  static double mD0=EvtPDL::getMeanMass(EvtPDL::getId("D0"));
  static double mDp=EvtPDL::getMeanMass(EvtPDL::getId("D+"));

  double pgmax=(s-4.0*md*md)/(2.0*w);

  assert(pgmax>0.0);

  double pgz=0.99*pgmax*exp(log(EvtRandom::Flat(1.0))/beta);

  if (EvtRandom::Flat(1.0)<0.5) pgz=-pgz;

  double k=fabs(pgz);

  EvtVector4R p4g(k,0.0,0.0,pgz);

  EvtVector4R p4res=p->getP4Restframe()-p4g;

  double mres=p4res.mass();

  double ed=mres/2.0;

  assert(ed>md);

  double pd=sqrt(ed*ed-md*md);


  //std::cout << "k, mres, w, md, ed, pd:"<<k<<" "<<mres<<" "<<w<<" "<<md<<" "<<ed<<" "<<pd<<std::endl;

  p->getDaug(0)->init(getDaug(0),p4res);
  p->getDaug(1)->init(getDaug(1),p4g);


  double sigma=beta*pow(2/w,beta)*(1+alpha*(1.5*L-2.0+EvtConst::pi*EvtConst::pi/3.0)/EvtConst::pi);

  double m=EvtPDL::getMeanMass(p->getDaug(0)->getId());
  double Gamma=EvtPDL::getWidth(p->getDaug(0)->getId());

  //mres is the energy of the psi(3770)

  double p0=0.0;
  if (ed>mD0) p0=sqrt(ed*ed-mD0*mD0);
  double pp=0.0;
  if (ed>mDp) pp=sqrt(ed*ed-mDp*mDp);

  double p0norm=sqrt(0.25*m*m-mD0*mD0);
  double ppnorm=sqrt(0.25*m*m-mDp*mDp);

  double r0=12.7;
  double rp=12.7;

  if (getNArg()==2){
    r0=getArg(0);
    rp=getArg(1);
  }

  double GammaTot=Gamma*(pp*pp*pp/(1+pp*pp*rp*rp)+p0*p0*p0/(1+p0*p0*r0*r0))/
    (ppnorm*ppnorm*ppnorm/(1+ppnorm*ppnorm*rp*rp)+
     p0norm*p0norm*p0norm/(1+p0norm*p0norm*r0*r0));
  

  sigma*=pd*pd*pd/((mres-m)*(mres-m)+0.25*GammaTot*GammaTot);

  assert(sigma>0.0);

  static double sigmax=sigma;

  if (sigma>sigmax){
    sigmax=sigma;
  }



  static int count=0;

  count++;

  //if (count%10000==0){
  //  std::cout << "sigma :"<<sigma<<std::endl;
  //  std::cout << "sigmax:"<<sigmax<<std::endl;
  //}

  double norm=sqrt(sigma);

  vertex(0,0,0,norm*p->eps(0)*p->epsParent(0).conj());
  vertex(1,0,0,norm*p->eps(1)*p->epsParent(0).conj());
  vertex(2,0,0,norm*p->eps(2)*p->epsParent(0).conj());

  vertex(0,1,0,norm*p->eps(0)*p->epsParent(1).conj());
  vertex(1,1,0,norm*p->eps(1)*p->epsParent(1).conj());
  vertex(2,1,0,norm*p->eps(2)*p->epsParent(1).conj());

  vertex(0,2,0,norm*p->eps(0)*p->epsParent(2).conj());
  vertex(1,2,0,norm*p->eps(1)*p->epsParent(2).conj());
  vertex(2,2,0,norm*p->eps(2)*p->epsParent(2).conj());

  vertex(0,0,1,norm*p->eps(0)*p->epsParent(0).conj());
  vertex(1,0,1,norm*p->eps(1)*p->epsParent(0).conj());
  vertex(2,0,1,norm*p->eps(2)*p->epsParent(0).conj());

  vertex(0,1,1,norm*p->eps(0)*p->epsParent(1).conj());
  vertex(1,1,1,norm*p->eps(1)*p->epsParent(1).conj());
  vertex(2,1,1,norm*p->eps(2)*p->epsParent(1).conj());

  vertex(0,2,1,norm*p->eps(0)*p->epsParent(2).conj());
  vertex(1,2,1,norm*p->eps(1)*p->epsParent(2).conj());
  vertex(2,2,1,norm*p->eps(2)*p->epsParent(2).conj());

  return;
}

