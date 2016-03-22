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
// Module: EvtVectorIsr.cc
//
// Description: 
//   This is a special decay model to generate e+e- -> phi gamma + soft gammas
//   using soft collinear ISR calculation from AfkQed
//   This is implemented as a decay of the VPHO.
//
// Modification history:
//
//    Joe Izen        Oct, 2005             Soft Colinear Photons (secondary ISR) ported from AfkQed
//    Joe Izen        Dec  16, 2002         Fix cos_theta distribution - prevents boom at cos_theta=+/-1 
//    RYD/Adriano     June 16, 1998         Module created
//
//------------------------------------------------------------------------
//
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>

#include <math.h>
#include <iostream>
#include <iomanip>
#include <sstream>


#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPhotonParticle.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtAbsLineShape.hh"
#include "EvtGenModels/EvtVectorIsr.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtAbsLineShape.hh"
#include <string>
#include "EvtGenBase/EvtVector4C.hh"

EvtVectorIsr::~EvtVectorIsr() {}

std::string EvtVectorIsr::getName(){

  return "VECTORISR";     
}

EvtDecayBase* EvtVectorIsr::clone(){

  return new EvtVectorIsr;
}

void EvtVectorIsr::init(){

  // check that there are 2 arguments
  
  checkNDaug(2);
  
  checkSpinParent(EvtSpinType::VECTOR);
  checkSpinDaughter(0,EvtSpinType::VECTOR);
  checkSpinDaughter(1,EvtSpinType::PHOTON);

  int narg = getNArg();
  if ( narg > 4 ) checkNArg(4);

  csfrmn=1.;
  csbkmn=1.;
  fmax=1.2;
  firstorder=false;

  if ( narg > 0 ) csfrmn=getArg(0);
  if ( narg > 1 ) csbkmn=getArg(1);
  if ( narg > 2 ) fmax=getArg(2);
  if ( narg > 3 ) firstorder=true;
}


void EvtVectorIsr::initProbMax(){

  noProbMax();
}

void EvtVectorIsr::decay( EvtParticle *p ){

  //the elctron mass
  double electMass=EvtPDL::getMeanMass(EvtPDL::getId("e-"));

  static EvtId gammaId=EvtPDL::getId("gamma");

  EvtParticle *phi;
  EvtParticle *gamma;

  //4-mom of the two colinear photons to the decay of the vphoton
  EvtVector4R p4softg1(0.,0.,0.,0.);
  EvtVector4R p4softg2(0.,0.,0.,0.);


  //get pointers to the daughters set
  //get masses/initial phase space - will overwrite the
  //p4s below to get the kinematic distributions correct
  p->initializePhaseSpace(getNDaug(),getDaugs());
  phi=p->getDaug(0);
  gamma=p->getDaug(1);

  //Generate soft colinear photons and the electron and positron energies after emission.
  //based on method of AfkQed and notes of Vladimir Druzhinin.
  //
  //function ckhrad(eb,q2m,r1,r2,e01,e02,f_col)
  //eb:      energy of incoming electrons in CM frame
  //q2m:     minimum invariant mass of the virtual photon after soft colinear photon emission
  //returned arguments
  //e01,e02: energies of e+ and e- after soft colinear photon emission
  //fcol:    weighting factor for Born cross section for use in an accept/reject test.


  double wcm=p->mass();
  double eb=0.5*wcm;

  //TO guarantee the collinear photons are softer than the ISR photon, require q2m > m*wcm
  double q2m=phi->mass()*wcm;
  double f_col(0.);
  double e01(0.);
  double e02(0.);
  double ebeam=eb;
  double wcm_new = wcm;
  double s_new = wcm*wcm;

  double fran = 1.;
  double f = 0;
  int m = 0;
  double largest_f=0;//only used when determining max weight for this vector particle mass
    
  if (!firstorder){
    while (fran > f){
      m++;    
    
      int n=0;
      while (f_col == 0.){
	n++;
	ckhrad(eb,q2m,e01,e02,f_col);
	if (n > 10000){
	  report(Severity::Info,"EvtGen") << "EvtVectorIsr is having problems. Called ckhrad 10000 times.\n";
	  assert(0);
	}
      }
    
      //Effective beam energy after soft photon emission (neglecting electron mass)
      ebeam = sqrt(e01*e02);
      wcm_new = 2*ebeam;
      s_new = wcm_new*wcm_new;
    
      //The Vector mass should never be greater than wcm_new
      if (phi->mass() > wcm_new){
	report(Severity::Info,"EvtGen") << "EvtVectorIsr finds Vector mass="<<phi->mass()<<" > Weff=" << wcm_new<<".  Should not happen\n";
	assert(0);
      }
 
      //Determine Born cross section @ wcm_new for e+e- -> gamma V.  We aren't interested in the absolute normalization
      //Just the functional dependence. Assuming a narrow resonance when determining cs_Born
      double cs_Born = 1.;
      if (EvtPDL::getMaxRange(phi->getId()) > 0.) {
	double x0 = 1 - EvtPDL::getMeanMass(phi->getId())*EvtPDL::getMeanMass(phi->getId())/s_new;
      
	//L = log(s/(electMass*electMass)  
	double L = 2.*log(wcm_new/electMass);
      
	// W(x0) is actually 2*alpha/pi times the following
	double W = (L-1.)*(1. - x0 +0.5*x0*x0);
      
	//Born cross section is actually 12*pi*pi*Gammaee/EvtPDL::getMeanMass(phi->getId()) times the following
	//(we'd need the full W(x0) as well)
	cs_Born = W/s_new;
      }
    
      f = cs_Born*f_col;

      //if fmax was set properly, f should NEVER be larger than fmax
      if (f > fmax && fmax > 0.){
	  report(Severity::Info,"EvtGen") << "EvtVectorIsr finds a problem with fmax, the maximum weight setting\n"
	     << "fmax is the third decay argument in the .dec file. VectorIsr attempts to set it reasonably if it wasn't provided\n"
	     << "To determine a more appropriate value, build GeneratorQAApp, and set the third argument for this decay <0.\n"
	     << "If you haven't been providing the first 2 arguments, set them to be 1. 1.). The program will report\n"
	     << "the largest weight it finds.  You should set fmax to be slightly larger.\n"
	     << "Alternatively try the following values for various vector particles: "
	     << "phi->1.15   J/psi-psi(4415)->0.105\n"
	     << "The current value of f and fmax for " << EvtPDL::name(phi->getId()) << " are " << f << "  " << fmax << "\n"
	     << "Will now assert\n";
	assert(0);
      }
 

      if (fmax > 0.) {
	fran = fmax*EvtRandom::Flat(0.0,1.0);
      }
    
      else {
	//determine max weight for this vector particle mass
	if (f>largest_f) {
	  largest_f = f;
	  report(Severity::Info,"EvtGen")  << m << " " <<  EvtPDL::name(phi->getId()) << " "
	       << "vector_mass " 
	       << " " << EvtPDL::getMeanMass(phi->getId()) << "  fmax should be at least " << largest_f 
	       << ".        f_col cs_B = " << f_col << " " << cs_Born 
	       << std::endl;
	}
	if (m%10000 == 0) {  
	  report(Severity::Info,"EvtGen") << m << " " <<  EvtPDL::name(phi->getId()) << " "
	       << "vector_mass " 
	       << " " << EvtPDL::getMeanMass(phi->getId()) << "  fmax should be at least " << largest_f 
	       << ".        f_col cs_B = " << f_col << " " << cs_Born 
	       << std::endl;
	}
      
	f_col = 0.;
	f = 0.;
	//determine max weight for this vector particle mass
      }
    
      if (m > 100000){
      
	if (fmax > 0.) report(Severity::Info,"EvtGen") << "EvtVectorIsr is having problems. Check the fmax value - the 3rd argument in the .dec file\n"
					     << "Recommended values for various vector particles: "
					     << "phi->1.15   J/psi-psi(4415)->0.105   "
					     << "Upsilon(1S,2S,3S)->0.14\n";
	assert(0);
      }
    }//while (fran > f)
  
  }//if (firstorder)
  
  //Compute parameters for boost to/from the system after colinear radiation

  double bet_l;
  double gam_l;
  double betgam_l;
  
  double csfrmn_new;
  double csbkmn_new;

  if (firstorder){
    bet_l = 0.;
    gam_l = 1.;
    betgam_l = 0.;
    csfrmn_new = csfrmn;
    csbkmn_new = csbkmn;
  } else {  
    double xx       = e02/e01;
    double sq_xx    = sqrt(xx);
    bet_l    = (1.-xx)/(1.+xx);
    gam_l    = (1.+xx)/(2.*sq_xx);
    betgam_l = (1.-xx)/(2.*sq_xx);
  
    //Boost photon cos_theta limits in lab to limits in the system after colinear rad
    csfrmn_new=(csfrmn - bet_l)/(1. - bet_l*csfrmn);
    csbkmn_new=(csbkmn - bet_l)/(1. - bet_l*csbkmn);
  }
 
//    //generate kinematics according to Bonneau-Martin article
//    //Nucl. Phys. B27 (1971) 381-397

  // For backward compatibility with .dec files before SP5, the backward cos limit for
  //the ISR photon is actually given as *minus* the actual limit. Sorry, this wouldn't be
  //my choice.  -Joe

   //gamma momentum in the vpho restframe *after* soft colinear radiation
  double pg = (s_new - phi->mass()*phi->mass())/(2.*wcm_new);


  //calculate the beta of incoming electrons after  colinear rad in the frame where e= and e- have equal momentum
  double beta=electMass/ebeam; //electMass/Ebeam = 1/gamma
  beta=sqrt(1. - beta*beta);   //sqrt (1 - (1/gamma)**2)

  double ymax=log((1.+beta*csfrmn_new)/(1.-beta*csfrmn_new));
  double ymin=log((1.-beta*csbkmn_new)/(1.+beta*csbkmn_new));

  // photon theta distributed as  2*beta/(1-beta**2*cos(theta)**2)
  double y=(ymax-ymin)*EvtRandom::Flat(0.0,1.0) + ymin;
  double cs=exp(y);
  cs=(cs - 1.)/(cs + 1.)/beta;
  double sn=sqrt(1-cs*cs);

  double fi=EvtRandom::Flat(EvtConst::twoPi);

  //four-vector for the phi
  double phi_p0 = sqrt(phi->mass()*phi->mass()+pg*pg);
  double phi_p3 = -pg*cs;


  //boost back to frame before colinear radiation.
  EvtVector4R p4phi(gam_l*phi_p0 + betgam_l*phi_p3,
		    -pg*sn*cos(fi),
		    -pg*sn*sin(fi),
		    betgam_l*phi_p0 + gam_l*phi_p3);

  double isr_p0 = pg;
  double isr_p3 = -phi_p3;
  EvtVector4R p4gamma(gam_l*isr_p0 + betgam_l*isr_p3,
		      -p4phi.get(1),
		      -p4phi.get(2),
		      betgam_l*isr_p0 + gam_l*isr_p3);

  
  //four-vectors of the collinear photons
  if (!firstorder) {
    p4softg1.set(0, eb-e02);    p4softg1.set(3, e02-eb);
    p4softg2.set(0, eb-e01);    p4softg2.set(3, eb-e01);
  }
  
  //save momenta for particles
  phi->init( getDaug(0),p4phi);
  gamma->init( getDaug(1),p4gamma);


  //add the two colinear photons as vphoton daughters
  EvtPhotonParticle *softg1=new EvtPhotonParticle;;
  EvtPhotonParticle *softg2=new EvtPhotonParticle;;
  softg1->init(gammaId,p4softg1);
  softg2->init(gammaId,p4softg2);
  softg1->addDaug(p);
  softg2->addDaug(p);

  //try setting the spin density matrix of the phi
  //get polarization vector for phi in its parents restframe.
  EvtVector4C phi0=phi->epsParent(0);
  EvtVector4C phi1=phi->epsParent(1);
  EvtVector4C phi2=phi->epsParent(2);

  //get polarization vector for a photon in its parents restframe.
  EvtVector4C gamma0=gamma->epsParentPhoton(0);
  EvtVector4C gamma1=gamma->epsParentPhoton(1);

  EvtComplex r1p=phi0*gamma0;
  EvtComplex r2p=phi1*gamma0;
  EvtComplex r3p=phi2*gamma0;


  EvtComplex r1m=phi0*gamma1;
  EvtComplex r2m=phi1*gamma1;
  EvtComplex r3m=phi2*gamma1;

  EvtComplex rho33=r3p*conj(r3p)+r3m*conj(r3m);
  EvtComplex rho22=r2p*conj(r2p)+r2m*conj(r2m);
  EvtComplex rho11=r1p*conj(r1p)+r1m*conj(r1m);

  EvtComplex rho13=r3p*conj(r1p)+r3m*conj(r1m);
  EvtComplex rho12=r2p*conj(r1p)+r2m*conj(r1m);
  EvtComplex rho23=r3p*conj(r2p)+r3m*conj(r2m);

  EvtComplex rho31=conj(rho13);
  EvtComplex rho32=conj(rho23);
  EvtComplex rho21=conj(rho12);


  EvtSpinDensity rho;
  rho.setDim(3);

  rho.set(0,0,rho11);
  rho.set(0,1,rho12);
  rho.set(0,2,rho13);
  rho.set(1,0,rho21);
  rho.set(1,1,rho22);
  rho.set(1,2,rho23);
  rho.set(2,0,rho31);
  rho.set(2,1,rho32);
  rho.set(2,2,rho33);

  setDaughterSpinDensity(0);
  phi->setSpinDensityForward(rho);

  return ;
}

double EvtVectorIsr::ckhrad1(double xx, double a, double b){
  //port of AfkQed/ckhrad.F function ckhrad1
  double yy = xx*xx; 
  double zz = 1. - 2*xx + yy; 
  return  0.5* (1. + yy + zz/(a-1.) + 0.25*b*( -0.5*(1. + 3*yy)*log(xx)) - zz  );
}

void EvtVectorIsr::ckhrad(const double& e_beam,const double& q2_min,double& e01,double& e02,double& f){
  //port of AfkQed/ckhrad.F subroutine ckhrad
  const double adp   = 1. / 137.0359895 / EvtConst::pi;
  const double pi2   = EvtConst::pi*EvtConst::pi;
  //  const double dme   = 0.00051099906;
  const double dme   = EvtPDL::getMeanMass(EvtPDL::getId("e-"));

  double r1=EvtRandom::Flat();//Generates Flat from 0 - 1
  double r2=EvtRandom::Flat();

  double sss    = 4.*e_beam*e_beam;
  double biglog = log(sss/(dme*dme));
  double beta   = 2.*adp*(biglog - 1.);
  double betae_lab = beta;
  double p3     = adp*(pi2/3. - 0.5);
  double p12    = adp*adp * (11./8. - 2.*pi2/3.);
  double coefener =  1. + 0.75*betae_lab + p3;
  double coef1 = coefener + 0.125*pi2*beta*beta;
  double coef2 = p12* biglog*biglog;
  double facts  = coef1 + coef2; 
  
  double y1_min = 0;
  double e1min  = 0.25 * q2_min/e_beam; 
  double y1_max = pow( 1. - e1min/e_beam, 0.5*beta );
  double y1     = y1_min +r1 *(y1_max - y1_min);
  e01           = e_beam *(1. - pow(y1, 2./beta) );
  
  double y2_min = 0.;
  double e2min  = 0.25 * q2_min/e01; 
  double y2_max = pow( 1. - e2min/e_beam, 0.5*beta);
  double y2     = y2_min +r2 *(y2_max - y2_min);
  e02           = e_beam *(1. - pow(y2, 2./beta) );
  

  double xx1 = e01/e_beam;
  double xx2 = e02/e_beam;

  f = y1_max * y2_max * ckhrad1(xx1,biglog,betae_lab) * ckhrad1(xx2,biglog,betae_lab) * facts;

  return;
 }

