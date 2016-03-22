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
// Description: Routine to decay vpho -> (DDx) + ISR photon from 3.9 to 4.3 GeV, using CLEO-c data (Brian Lang)
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
#include "EvtGenModels/EvtVPHOtoVISRHi.hh"
#include <string>

using std::endl;

EvtVPHOtoVISRHi::~EvtVPHOtoVISRHi() {}

std::string EvtVPHOtoVISRHi::getName(){

  return "VPHOTOVISRHI"; 
    
}


EvtDecayBase* EvtVPHOtoVISRHi::clone(){

  return new EvtVPHOtoVISRHi;

}

void EvtVPHOtoVISRHi::init(){

  // check that there are 0 or 1 arguments
  checkNArg(0,1);

  // check that there are 2 daughters
  checkNDaug(2);

  // check the parent and daughter spins
  checkSpinParent(EvtSpinType::VECTOR);
  checkSpinDaughter(0,EvtSpinType::VECTOR);
  checkSpinDaughter(1,EvtSpinType::PHOTON);
}

void EvtVPHOtoVISRHi::initProbMax() {

   setProbMax(20.0);

}      

void EvtVPHOtoVISRHi::decay( EvtParticle *p){
  //take photon along z-axis, either forward or backward.
  //Implement this as generating the photon momentum along 
  //the z-axis uniformly 
   double power=1;
   if (getNArg()==1) power=getArg(0);
   // define particle names
  static EvtId D0=EvtPDL::getId("D0");
  static EvtId D0B=EvtPDL::getId("anti-D0");
  static EvtId DP=EvtPDL::getId("D+");
  static EvtId DM=EvtPDL::getId("D-");
  static EvtId DSM=EvtPDL::getId("D_s-");
  static EvtId DSP=EvtPDL::getId("D_s+");
  static EvtId DSMS=EvtPDL::getId("D_s*-");
  static EvtId DSPS=EvtPDL::getId("D_s*+");
  static EvtId D0S=EvtPDL::getId("D*0");
  static EvtId D0BS=EvtPDL::getId("anti-D*0");
  static EvtId DPS=EvtPDL::getId("D*+");
  static EvtId DMS=EvtPDL::getId("D*-");
  // setup some parameters
  double w=p->mass();
  double s=w*w;
  double L=2.0*log(w/0.000511);
  double alpha=1/137.0;
  double beta=(L-1)*2.0*alpha/EvtConst::pi;
  // make sure only 2 or 3 body are present
  assert (p->getDaug(0)->getNDaug() == 2 || p->getDaug(0)->getNDaug() == 3);

  // determine minimum rest mass of parent
  double md1 = EvtPDL::getMeanMass(p->getDaug(0)->getDaug(0)->getId());
  double md2 = EvtPDL::getMeanMass(p->getDaug(0)->getDaug(1)->getId());
  double minResMass = md1+md2;
  if (p->getDaug(0)->getNDaug() == 3) {
     double md3 = EvtPDL::getMeanMass(p->getDaug(0)->getDaug(2)->getId());
     minResMass = minResMass + md3;
  }
  
  // calculate the maximum energy of the ISR photon
  double pgmax=(s-minResMass*minResMass)/(2.0*w);
  double pgz=0.99*pgmax*exp(log(EvtRandom::Flat(1.0))/(beta*power));
  if (EvtRandom::Flat(1.0)<0.5) pgz=-pgz;
  
  double k=fabs(pgz);
  // print of ISR energy 
  // std::cout << "Energy ISR :"<< k <<std::endl;
  EvtVector4R p4g(k,0.0,0.0,pgz);

  EvtVector4R p4res=p->getP4Restframe()-p4g;

  double mres=p4res.mass();

  // set masses
  p->getDaug(0)->init(getDaug(0),p4res);
  p->getDaug(1)->init(getDaug(1),p4g);

  
  // determine XS - langbw
  // very crude way of determining XS just a simple straight line Approx.
  // this was determined by eye.
  // lots of cout statements to make plots to check that things are working as expected
  double sigma=9.0;
  if (mres<=3.9) sigma = 0.00001;

  bool sigmacomputed(false);

  // DETERMINE XS FOR D*D*
  if (p->getDaug(0)->getNDaug() == 2 
      &&((p->getDaug(0)->getDaug(0)->getId()==D0S 
          && p->getDaug(0)->getDaug(1)->getId()==D0BS)
         ||(p->getDaug(0)->getDaug(0)->getId()==DPS 
            && p->getDaug(0)->getDaug(1)->getId()==DMS))){
     if(mres>4.18) {
        sigma*=5./9.*(1.-1.*sqrt((4.18-mres)*(4.18-mres))/(4.3-4.18));
     }  
     else if(mres>4.07 && mres<=4.18) {
     sigma*=5./9.;
     }  
     else if (mres<=4.07&&mres>4.03)
     {
        sigma*=(5./9. - 1.5/9.*sqrt((4.07-mres)*(4.07-mres))/(4.07-4.03)); 
     }
     else if (mres<=4.03&& mres>=4.013)
     {
        sigma*=(3.5/9. - 3.5/9.*sqrt((4.03-mres)*(4.03-mres))/(4.03-4.013)); 
     }
     else{     
        sigma=0.00001; 
     }
     sigmacomputed = true;
//     std::cout << "DSDSXS "<<sigma<< " " <<  mres<<std::endl;
  }
  
  // DETERMINE XS FOR D*D
  if(p->getDaug(0)->getNDaug() == 2 && ((p->getDaug(0)->getDaug(0)->getId()==D0S 
                                         && p->getDaug(0)->getDaug(1)->getId()==D0B)
                                        ||(p->getDaug(0)->getDaug(0)->getId()==DPS 
                                           && p->getDaug(0)->getDaug(1)->getId()==DM) 
                                        ||(p->getDaug(0)->getDaug(0)->getId()==D0BS
                                           && p->getDaug(0)->getDaug(1)->getId()==D0)
                                        ||(p->getDaug(0)->getDaug(0)->getId()==DMS 
                                           && p->getDaug(0)->getDaug(1)->getId()==DP)) )
  {
     if(mres>=4.2){
        sigma*=1.5/9.;
     }
     else if( mres>4.06 && mres<4.2){
        sigma*=((1.5/9.+2.5/9.*sqrt((4.2-mres)*(4.2-mres))/(4.2-4.06)));
     }  
     else if(mres>=4.015 && mres<4.06){
        sigma*=((4./9.+3./9.*sqrt((4.06-mres)*(4.06-mres))/(4.06-4.015)));
     }  
     else if (mres<4.015 && mres>=3.9){
        sigma*=((7./9.-7/9.*sqrt((4.015-mres)*(4.015-mres))/(4.015-3.9))); 
     } 
     else { 
        sigma = 0.00001; 
     }
     sigmacomputed = true;
//     std::cout << "DSDXS "<<sigma<< " " <<  mres<<std::endl;
  }
     
  // DETERMINE XS FOR Ds*Ds*
  if (((p->getDaug(0)->getDaug(0)->getId()==DSPS && p->getDaug(0)->getDaug(1)->getId()==DSMS)))
  {
     if(mres>(2.112+2.112)){
      sigma=0.4; 
     }
     else  {
//      sigma=0.4; 
//      sigma = 0 surely below Ds*Ds* threshold? - ponyisi
        sigma=0.00001;
     }
     sigmacomputed = true;
//     std::cout << "DsSDsSXS "<<sigma<< " " <<  mres<<std::endl;
  }

  // DETERMINE XS FOR Ds*Ds
  if (p->getDaug(0)->getNDaug() == 2 && ((p->getDaug(0)->getDaug(0)->getId()==DSPS 
                                          && p->getDaug(0)->getDaug(1)->getId()==DSM)
                                         || (p->getDaug(0)->getDaug(0)->getId()==DSMS
                                             && p->getDaug(0)->getDaug(1)->getId()==DSP)))
  {
     if(mres>4.26){
        sigma=0.05; 
     } 
     else if (mres>4.18 && mres<=4.26){
        sigma*=1./9.*(0.05+0.95*sqrt((4.26-mres)*(4.26-mres))/(4.26-4.18));
     } 
     else if (mres>4.16 && mres<=4.18){
        sigma*=1/9.; 
     } 
     else if (mres<=4.16 && mres>4.08){
        sigma*=1/9.*(1-sqrt((4.16-mres)*(4.16-mres))/(4.16-4.08)); 
     }
     else if (mres<=(4.08)){
        sigma=0.00001; 
     }
     sigmacomputed = true;
//     std::cout << "DsSDsXS "<<sigma<< " " <<  mres<<std::endl;
  }

  // DETERMINE XS FOR DD
  if (p->getDaug(0)->getNDaug() == 2 && ((p->getDaug(0)->getDaug(0)->getId()==D0 
                                          && p->getDaug(0)->getDaug(1)->getId()==D0B)
                                         ||(p->getDaug(0)->getDaug(0)->getId()==DP 
                                            && p->getDaug(0)->getDaug(1)->getId()==DM))){ 
     sigma*=0.4/9.;  
     sigmacomputed = true;
//     std::cout << "DDXS "<<sigma<< " " <<  mres<<std::endl;
  } 
  
  // DETERMINE XS FOR DsDs
  if (p->getDaug(0)->getNDaug() == 2 && ((p->getDaug(0)->getDaug(0)->getId()==DSP && p->getDaug(0)->getDaug(1)->getId()==DSM))){
     sigma*=0.2/9.;
     sigmacomputed = true;
//     std::cout << "DsDsXS "<<sigma<< " " <<  mres<<std::endl;
  } 

  // DETERMINE XS FOR MULTIBODY
  if (p->getDaug(0)->getNDaug() == 3){
     if(mres>4.03){
        sigma*=0.5/9.;
     }
     else {
        sigma=0.00001; 
     }
     sigmacomputed = true;
//     std::cout << "DSDpiXS "<<sigma<< " " <<  mres<<std::endl;
  }

  if (! sigmacomputed) {
    report(Severity::Error,"EvtGen") << "VPHOTOVISRHI: This model requires daughters to be listed in a particular order." << endl
                           << "The following are acceptable:" << endl
                           << "D0 anti-D0" << endl
                           << "D+ D-" << endl
                           << "D*0 anti-D0" << endl
                           << "anti-D*0 D0" << endl
                           << "D*+ D-" << endl
                           << "D*- D+" << endl
                           << "D*0 anti-D*0" << endl
                           << "D*+ D*-" << endl
                           << "D_s+ D_s-" << endl
                           << "D_s*+ D_s-" << endl
                           << "D_s*- D_s+" << endl
                           << "D_s*+ D_s*-" << endl
                           << "(D* D pi can be in any order)" << endl
                           << "Aborting..." << endl;
    assert(0);
  }

  if(sigma<0) sigma = 0.0;

//   static double sigmax=sigma;
//   if (sigma>sigmax){
//      sigmax=sigma;
//   }
  
  static int count=0;
  
  count++;
  
//   if (count%10000==0){
//      std::cout << "sigma :"<<sigma<<std::endl;
//      std::cout << "sigmax:"<<sigmax<<std::endl;
//   }
  
  double norm=sqrt(sigma);
  
//  EvtParticle* d=p->getDaug(0);
  
  
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
