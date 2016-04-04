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
// Module: EvtGoityRoberts.cc
//
// Description: Routine to decay vector-> scalar scalar
//
// Modification history:
//
//    RYD     November 24, 1996        Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtGoityRoberts.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"
#include <string>
#include "EvtGenBase/EvtVector4C.hh"

EvtGoityRoberts::~EvtGoityRoberts() {}

std::string EvtGoityRoberts::getName(){

  return "GOITY_ROBERTS";     

}


EvtDecayBase* EvtGoityRoberts::clone(){

  return new EvtGoityRoberts;

}

void EvtGoityRoberts::init(){

  // check that there are 0 arguments
  checkNArg(0);
  checkNDaug(4);

  checkSpinParent(EvtSpinType::SCALAR);
  checkSpinDaughter(1,EvtSpinType::SCALAR);
  checkSpinDaughter(2,EvtSpinType::DIRAC);
  checkSpinDaughter(3,EvtSpinType::NEUTRINO);

}


void EvtGoityRoberts::initProbMax() {

   setProbMax( 3000.0);
}      

void EvtGoityRoberts::decay( EvtParticle *p){

  //added by Lange Jan4,2000
  static EvtId DST0=EvtPDL::getId("D*0");
  static EvtId DSTB=EvtPDL::getId("anti-D*0");
  static EvtId DSTP=EvtPDL::getId("D*+");
  static EvtId DSTM=EvtPDL::getId("D*-");
  static EvtId D0=EvtPDL::getId("D0");
  static EvtId D0B=EvtPDL::getId("anti-D0");
  static EvtId DP=EvtPDL::getId("D+");
  static EvtId DM=EvtPDL::getId("D-");



  EvtId meson=getDaug(0);

  if (meson==DST0||meson==DSTP||meson==DSTM||meson==DSTB) {
    DecayBDstarpilnuGR(p,getDaug(0),getDaug(2),getDaug(3));
  }
  else{
    if (meson==D0||meson==DP||meson==DM||meson==D0B) {
      DecayBDpilnuGR(p,getDaug(0),getDaug(2),getDaug(3));
    }
    else{
      report(Severity::Error,"EvtGen") << "Wrong daugther in EvtGoityRoberts!\n";
    }
  }
  return ;
}

void EvtGoityRoberts::DecayBDstarpilnuGR(EvtParticle *pb,EvtId ndstar,
					 EvtId nlep, EvtId /*nnu*/)
{

  pb->initializePhaseSpace(getNDaug(),getDaugs());

  //added by Lange Jan4,2000
  static EvtId EM=EvtPDL::getId("e-");
  static EvtId EP=EvtPDL::getId("e+");
  static EvtId MUM=EvtPDL::getId("mu-");
  static EvtId MUP=EvtPDL::getId("mu+");

  EvtParticle *dstar, *pion, *lepton, *neutrino;
  
  // pb->makeDaughters(getNDaug(),getDaugs());
  dstar=pb->getDaug(0);
  pion=pb->getDaug(1);
  lepton=pb->getDaug(2);
  neutrino=pb->getDaug(3);

  EvtVector4C l1, l2, et0, et1, et2;
  
  EvtVector4R v,vp,p4_pi;
  double w;
  
  v.set(1.0,0.0,0.0,0.0);       //4-velocity of B meson
  vp=(1.0/dstar->getP4().mass())*dstar->getP4();  //4-velocity of D*
  p4_pi=pion->getP4();

  w=v*vp;                       //four velocity transfere.

  EvtTensor4C omega;

  double mb=EvtPDL::getMeanMass(pb->getId());     //B mass
  double md=EvtPDL::getMeanMass(ndstar);   //D* mass

  EvtComplex dmb(0.0460,-0.5*0.00001);   // B*-B mass splitting ?
  EvtComplex dmd(0.1421,-0.5*0.00006);
                             // The last two sets of numbers should
                             // be correctly calculated from the
                             // dstar and pion charges.
  double g = 0.5;         // EvtAmplitude proportional to these coupling constants
  double alpha3 =  0.690; // See table I in G&R's paper
  double alpha1 = -1.430;
  double alpha2 = -0.140;
  double f0 = 0.093;      // The pion decay constants set to 93 MeV

  EvtComplex dmt3(0.563,-0.5*0.191); // Mass splitting = dmt - iGamma/2  
  EvtComplex dmt1(0.392,-0.5*1.040);
  EvtComplex dmt2(0.709,-0.5*0.405);
                   
  double betas=0.285;      // magic number for meson wave function ground state
  double betap=0.280;      // magic number for meson wave function state "1"
  double betad=0.260;      // magic number for meson wave function state "2"
  double betasp=betas*betas+betap*betap;
  double betasd=betas*betas+betad*betad;

  double lambdabar=0.750;  //M(0-,1-) - mQ From Goity&Roberts's code

// Isgur&Wise fct
  double xi = exp(lambdabar*lambdabar*(1.0-w*w)/(4*betas*betas));
  double xi1= -1.0*sqrt(2.0/3.0)*(
       lambdabar*lambdabar*(w*w-1.0)/(4*betas*betas))*
       exp(lambdabar*lambdabar*(1.0-w*w)/(4*betas*betas));
  double rho1= sqrt(1.0/2.0)*(lambdabar/betas)*
               pow((2*betas*betap/(betasp)),2.5)*
               exp(lambdabar*lambdabar*(1.0-w*w)/(2*betasp));
  double rho2= sqrt(1.0/8.0)*(lambdabar*lambdabar/(betas*betas))*
               pow((2*betas*betad/(betasd)),3.5)*
               exp(lambdabar*lambdabar*(1.0-w*w)/(2*betasd));

  //report(Severity::Info,"EvtGen") <<"rho's:"<<rho1<<rho2<<endl;

  EvtComplex h1nr,h2nr,h3nr,f1nr,f2nr;
  EvtComplex f3nr,f4nr,f5nr,f6nr,knr,g1nr,g2nr,g3nr,g4nr,g5nr;
  EvtComplex h1r,h2r,h3r,f1r,f2r,f3r,f4r,f5r,f6r,kr,g1r,g2r,g3r,g4r,g5r;
  EvtComplex h1,h2,h3,f1,f2,f3,f4,f5,f6,k,g1,g2,g3,g4,g5;

  // Non-resonance part
  h1nr = -g*xi*(p4_pi*v)/(f0*mb*md*(EvtComplex(p4_pi*v,0.0)+dmb));
  h2nr = -g*xi/(f0*mb*(EvtComplex(p4_pi*v,0.0)+dmb));
  h3nr = -(g*xi/(f0*md))*(1.0/(EvtComplex(p4_pi*v,0.0)+dmb)
                         -EvtComplex((1.0+w)/(p4_pi*vp),0.0));

  f1nr = -(g*xi/(2*f0*mb))*(1.0/(EvtComplex(p4_pi*v,0.0)+dmb) -
         1.0/(EvtComplex(p4_pi*vp,0.0)+dmd));
  f2nr = f1nr*mb/md;
  f3nr = EvtComplex(0.0);
  f4nr = EvtComplex(0.0);
  f5nr = (g*xi/(2*f0*mb*md))*(EvtComplex(1.0,0.0)
                 +(p4_pi*v)/(EvtComplex(p4_pi*v,0.0)+dmb));
  f6nr = (g*xi/(2*f0*mb))*(1.0/(EvtComplex(p4_pi*v,0.0)+dmb)
                 -EvtComplex(1.0/(p4_pi*vp),0.0));

  knr = (g*xi/(2*f0))*((p4_pi*(vp-w*v))/(EvtComplex(p4_pi*v,0.0)+dmb) +
                 EvtComplex((p4_pi*(v-w*vp))/(p4_pi*vp),0.0));
  
  g1nr = EvtComplex(0.0);
  g2nr = EvtComplex(0.0);
  g3nr = EvtComplex(0.0);
  g4nr = (g*xi)/(f0*md*EvtComplex(p4_pi*vp));
  g5nr = EvtComplex(0.0);

  // Resonance part (D** removed by hand - alainb)
  h1r = -alpha1*rho1*(p4_pi*v)/(f0*mb*md*(EvtComplex(p4_pi*v,0.0)+dmt1)) +
        alpha2*rho2*(p4_pi*(v+2.0*w*v-vp))
        /(3*f0*mb*md*(EvtComplex(p4_pi*v,0.0)+dmt2)) -
        alpha3*xi1*(p4_pi*v)/(f0*mb*md*EvtComplex(p4_pi*v,0.0)+dmt3); 
  h2r = -alpha2*(1+w)*rho2/(3*f0*mb*(EvtComplex(p4_pi*v,0.0)+dmt2)) -
        alpha3*xi1/(f0*mb*(EvtComplex(p4_pi*v,0.0)+dmt3));
  h3r = alpha2*rho2*(1+w)/(3*f0*md*(EvtComplex(p4_pi*v,0.0)+dmt2)) -
        alpha3*xi1/(f0*md*(EvtComplex(p4_pi*v,0.0)+dmt3));

  f1r = -alpha2*rho2*(w-1.0)/(6*f0*mb*(EvtComplex(p4_pi*v,0.0)+dmt2)) -
        alpha3*xi1/(2*f0*mb*(EvtComplex(p4_pi*v,0.0)+dmt3));
  f2r = f1r*mb/md;
  f3r = EvtComplex(0.0);
  f4r = EvtComplex(0.0);
  f5r = alpha1*rho1*(p4_pi*v)/(2*f0*mb*md*(EvtComplex(p4_pi*v,0.0)+dmt1)) +
        alpha2*rho2*(p4_pi*(vp-v/3.0-2.0/3.0*w*v))/
        (2*f0*mb*md*(EvtComplex(p4_pi*v,0.0)+dmt2)) +
        alpha3*xi1*(p4_pi*v)/(2*f0*mb*md*(EvtComplex(p4_pi*v,0.0)+dmt3));
  f6r = alpha2*rho2*(w-1.0)/(6*f0*mb*(EvtComplex(p4_pi*v,0.0)+dmt2)) +
        alpha3*xi1/(2*f0*mb*(EvtComplex(p4_pi*v,0.0)+dmt3));

  kr = -alpha1*rho1*(w-1.0)*(p4_pi*v)/(2*f0*(EvtComplex(p4_pi*v,0.0)+dmt1)) -
       alpha2*rho2*(w-1.0)*(p4_pi*(vp-w*v))
       /(3*f0*(EvtComplex(p4_pi*v,0.0)+dmt2)) +
       alpha3*xi1*(p4_pi*(vp-w*v))/(2*f0*(EvtComplex(p4_pi*v,0.0)+dmt3));
  
  g1r = EvtComplex(0.0);
  g2r = EvtComplex(0.0);
  g3r = -g2r;
  g4r = 2.0*alpha2*rho2/(3*f0*md*(EvtComplex(p4_pi*v,0.0)+dmt2));
  g5r = EvtComplex(0.0);

  //Sum
  h1 = h1nr + h1r;
  h2 = h2nr + h2r;
  h3 = h3nr + h3r;

  f1 = f1nr + f1r;
  f2 = f2nr + f2r;
  f3 = f3nr + f3r;
  f4 = f4nr + f4r;
  f5 = f5nr + f5r;
  f6 = f6nr + f6r;

  k = knr+kr;
  
  g1 = g1nr + g1r;
  g2 = g2nr + g2r;
  g3 = g3nr + g3r;
  g4 = g4nr + g4r;
  g5 = g5nr + g5r;

  EvtTensor4C g_metric;
  g_metric.setdiag(1.0,-1.0,-1.0,-1.0);

  if (nlep==EM||nlep==MUM){ 
    omega=EvtComplex(0.0,0.5)*dual(h1*mb*md*EvtGenFunctions::directProd(v,vp)+
                             h2*mb*EvtGenFunctions::directProd(v,p4_pi)+
                             h3*md*EvtGenFunctions::directProd(vp,p4_pi))+
        f1*mb*EvtGenFunctions::directProd(v,p4_pi)+f2*md*EvtGenFunctions::directProd(vp,p4_pi)+
                       f3*EvtGenFunctions::directProd(p4_pi,p4_pi)+f4*mb*mb*EvtGenFunctions::directProd(v,v)+
        f5*mb*md*EvtGenFunctions::directProd(vp,v)+f6*mb*EvtGenFunctions::directProd(p4_pi,v)+k*g_metric+
        EvtComplex(0.0,0.5)*EvtGenFunctions::directProd(dual(EvtGenFunctions::directProd(vp,p4_pi)).cont2(v),
                              (g1*p4_pi+g2*mb*v))+
        EvtComplex(0.0,0.5)*EvtGenFunctions::directProd((g3*mb*v+g4*md*vp+g5*p4_pi),
                             dual(EvtGenFunctions::directProd(vp,p4_pi)).cont2(v));

   l1=EvtLeptonVACurrent(lepton->spParent(0),neutrino->spParentNeutrino());
   l2=EvtLeptonVACurrent(lepton->spParent(1),neutrino->spParentNeutrino());
  }
  else{
    if (nlep==EP||nlep==MUP){ 
      omega=EvtComplex(0.0,-0.5)*dual(h1*mb*md*EvtGenFunctions::directProd(v,vp)+
                             h2*mb*EvtGenFunctions::directProd(v,p4_pi)+
                                      h3*md*EvtGenFunctions::directProd(vp,p4_pi))+
        f1*mb*EvtGenFunctions::directProd(v,p4_pi)+f2*md*EvtGenFunctions::directProd(vp,p4_pi)+
                       f3*EvtGenFunctions::directProd(p4_pi,p4_pi)+f4*mb*mb*EvtGenFunctions::directProd(v,v)+
        f5*mb*md*EvtGenFunctions::directProd(vp,v)+f6*mb*EvtGenFunctions::directProd(p4_pi,v)+k*g_metric+
        EvtComplex(0.0,-0.5)*EvtGenFunctions::directProd(dual(EvtGenFunctions::directProd(vp,p4_pi)).cont2(v),
                              (g1*p4_pi+g2*mb*v))+
        EvtComplex(0.0,-0.5)*EvtGenFunctions::directProd((g3*mb*v+g4*md*vp+g5*p4_pi),
                             dual(EvtGenFunctions::directProd(vp,p4_pi)).cont2(v));

   l1=EvtLeptonVACurrent(neutrino->spParentNeutrino(),lepton->spParent(0));
   l2=EvtLeptonVACurrent(neutrino->spParentNeutrino(),lepton->spParent(1));
    }
    else{
   report(Severity::Debug,"EvtGen") << "42387dfs878w wrong lepton number\n";
    }
 }

  et0=omega.cont2( dstar->epsParent(0).conj() );
  et1=omega.cont2( dstar->epsParent(1).conj() );
  et2=omega.cont2( dstar->epsParent(2).conj() );

  vertex(0,0,l1.cont(et0));
  vertex(0,1,l2.cont(et0));

  vertex(1,0,l1.cont(et1));
  vertex(1,1,l2.cont(et1));

  vertex(2,0,l1.cont(et2));
  vertex(2,1,l2.cont(et2));

  return;

}

void EvtGoityRoberts::DecayBDpilnuGR(EvtParticle *pb,EvtId nd,
				     EvtId nlep, EvtId /*nnu*/)

{
  //added by Lange Jan4,2000
  static EvtId EM=EvtPDL::getId("e-");
  static EvtId EP=EvtPDL::getId("e+");
  static EvtId MUM=EvtPDL::getId("mu-");
  static EvtId MUP=EvtPDL::getId("mu+");

  EvtParticle *d, *pion, *lepton, *neutrino;

  pb->initializePhaseSpace(getNDaug(),getDaugs());
  d=pb->getDaug(0);
  pion=pb->getDaug(1);
  lepton=pb->getDaug(2);
  neutrino=pb->getDaug(3);

  EvtVector4C l1, l2, et0, et1, et2;
 
  EvtVector4R v,vp,p4_pi;
  double w;
  
  v.set(1.0,0.0,0.0,0.0);       //4-velocity of B meson
  vp=(1.0/d->getP4().mass())*d->getP4();  //4-velocity of D
  p4_pi=pion->getP4();                  //4-momentum of pion
  w=v*vp;                       //four velocity transfer.
  
  double mb=EvtPDL::getMeanMass(pb->getId());     //B mass
  double md=EvtPDL::getMeanMass(nd);   //D* mass
  EvtComplex dmb(0.0460,-0.5*0.00001);   //B mass splitting ?
                      //The last two numbers should be
                      //correctly calculated from the
                      //dstar and pion particle number.

  double g = 0.5;         // Amplitude proportional to these coupling constants
  double alpha3 =  0.690; // See table I in G&R's paper
  double alpha1 = -1.430;
  double alpha2 = -0.140;
  double f0=0.093;        // The pion decay constant set to 93 MeV

  EvtComplex dmt3(0.563,-0.5*0.191); // Mass splitting = dmt - iGamma/2  
  EvtComplex dmt1(0.392,-0.5*1.040);
  EvtComplex dmt2(0.709,-0.5*0.405);
                   
  double betas=0.285;      // magic number for meson wave function ground state
  double betap=0.280;      // magic number for meson wave function state "1"
  double betad=0.260;      // magic number for meson wave function state "2"
  double betasp=betas*betas+betap*betap;
  double betasd=betas*betas+betad*betad;

  double lambdabar=0.750;  //M(0-,1-) - mQ From Goity&Roberts's code

  // Isgur&Wise fct
  double xi = exp(lambdabar*lambdabar*(1.0-w*w)/(4*betas*betas));
  double xi1= -1.0*sqrt(2.0/3.0)*(lambdabar*lambdabar*(w*w-1.0)/(4*betas*betas))*
              exp(lambdabar*lambdabar*(1.0-w*w)/(4*betas*betas));
  double rho1= sqrt(1.0/2.0)*(lambdabar/betas)*
               pow((2*betas*betap/(betasp)),2.5)*
               exp(lambdabar*lambdabar*(1.0-w*w)/(2*betasp));
  double rho2= sqrt(1.0/8.0)*(lambdabar*lambdabar/(betas*betas))*
               pow((2*betas*betad/(betasd)),3.5)*
               exp(lambdabar*lambdabar*(1.0-w*w)/(2*betasd));

  EvtComplex h,a1,a2,a3;
  EvtComplex hnr,a1nr,a2nr,a3nr;
  EvtComplex hr,a1r,a2r,a3r;

// Non-resonance part (D* and D** removed by hand - alainb)
  hnr = g*xi*(1.0/(EvtComplex(p4_pi*v,0.0)+dmb))/(2*f0*mb*md);
  a1nr= -1.0*g*xi*(1+w)*(1.0/(EvtComplex(p4_pi*v,0.0)+dmb))/(2*f0);
  a2nr= g*xi*((p4_pi*(v+vp))/(EvtComplex(p4_pi*v,0.0)+dmb))/(2*f0*mb);
  a3nr=EvtComplex(0.0,0.0);

// Resonance part (D** remove by hand - alainb)
  hr = alpha2*rho2*(w-1)*(1.0/(EvtComplex(p4_pi*v,0.0)+dmt2))/(6*f0*mb*md) +
       alpha3*xi1*(1.0/(EvtComplex(p4_pi*v,0.0)+dmt3))/(2*f0*mb*md);
  a1r= -1.0*alpha2*rho2*(w*w-1)*(1.0/(EvtComplex(p4_pi*v,0.0)+dmt2))/(6*f0) -
       alpha3*xi1*(1+w)*(1.0/(EvtComplex(p4_pi*v,0.0)+dmt3))/(2*f0);
  a2r= alpha1*rho1*((p4_pi*v)/(EvtComplex(p4_pi*v,0.0)+dmt1))/(2*f0*mb) +
       alpha2*rho2*(0.5*p4_pi*(w*vp-v)+p4_pi*(vp-w*v))/
                  (3*f0*mb*(EvtComplex(p4_pi*v,0.0)+dmt2)) +
       alpha3*xi1*((p4_pi*(v+vp))/(EvtComplex(p4_pi*v,0.0)+dmt3))/(2*f0*mb);
  a3r= -1.0*alpha1*rho1*((p4_pi*v)/(EvtComplex(p4_pi*v,0.0)+dmt1))/(2*f0*md) -
       alpha2*rho2*((p4_pi*(vp-w*v))/(EvtComplex(p4_pi*v,0.0)+dmt2))/(2*f0*md);

// Sum
  h=hnr+hr;
  a1=a1nr+a1r;
  a2=a2nr+a2r;
  a3=a3nr+a3r;

  EvtVector4C omega;

  if ( nlep==EM|| nlep==MUM ) {
    omega=EvtComplex(0.0,-1.0)*h*mb*md*dual(EvtGenFunctions::directProd(vp,p4_pi)).cont2(v)+
                 a1*p4_pi+a2*mb*v+a3*md*vp;
    l1=EvtLeptonVACurrent(
             lepton->spParent(0),neutrino->spParentNeutrino());
    l2=EvtLeptonVACurrent(
             lepton->spParent(1),neutrino->spParentNeutrino());
  }
  else{
    if ( nlep==EP|| nlep==MUP ) {
     omega=EvtComplex(0.0,1.0)*h*mb*md*dual(EvtGenFunctions::directProd(vp,p4_pi)).cont2(v)+
                 a1*p4_pi+a2*mb*v+a3*md*vp;
     l1=EvtLeptonVACurrent(
              neutrino->spParentNeutrino(),lepton->spParent(0));
     l2=EvtLeptonVACurrent(
              neutrino->spParentNeutrino(),lepton->spParent(1));
    }
    else{
     report(Severity::Error,"EvtGen") << "42387dfs878w wrong lepton number\n";
    }
  }

  vertex(0,l1*omega);
  vertex(1,l2*omega);

return;

}

