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
// Module: EvtSemiLeptonicAmp.cc
//
// Description: Base class for semileptonic decays
//
// Modification history:
//
//    DJL       April 17,1998       Module created
//
//------------------------------------------------------------------------
//
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtSemiLeptonicScalarAmp.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtAmp.hh"
#include "EvtGenBase/EvtScalarParticle.hh"
#include "EvtGenBase/EvtVectorParticle.hh"
#include "EvtGenBase/EvtTensorParticle.hh"

double EvtSemiLeptonicAmp::CalcMaxProb( EvtId parent, EvtId meson, 
					EvtId lepton, EvtId nudaug,
                     EvtSemiLeptonicFF *FormFactors ) {

  //This routine takes the arguements parent, meson, and lepton
  //number, and a form factor model, and returns a maximum
  //probability for this semileptonic form factor model.  A
  //brute force method is used.  The 2D cos theta lepton and
  //q2 phase space is probed.

  //Start by declaring a particle at rest.

  //It only makes sense to have a scalar parent.  For now. 
  //This should be generalized later.

  EvtScalarParticle *scalar_part;
  EvtParticle *root_part;

  scalar_part=new EvtScalarParticle;

  //cludge to avoid generating random numbers!
  scalar_part->noLifeTime();

  EvtVector4R p_init;
  
  p_init.set(EvtPDL::getMass(parent),0.0,0.0,0.0);
  scalar_part->init(parent,p_init);
  root_part=(EvtParticle *)scalar_part;
  root_part->setDiagonalSpinDensity();      

  EvtParticle *daughter, *lep, *trino;
  
  EvtAmp amp;

  EvtId listdaug[3];
  listdaug[0] = meson;
  listdaug[1] = lepton;
  listdaug[2] = nudaug;

  amp.init(parent,3,listdaug);

  root_part->makeDaughters(3,listdaug);
  daughter=root_part->getDaug(0);
  lep=root_part->getDaug(1);
  trino=root_part->getDaug(2);

  //cludge to avoid generating random numbers!
  daughter->noLifeTime();
  lep->noLifeTime();
  trino->noLifeTime();


  //Initial particle is unpolarized, well it is a scalar so it is 
  //trivial
  EvtSpinDensity rho;
  rho.setDiag(root_part->getSpinStates());
  
  double mass[3];
  
  double m = root_part->mass();
  
  EvtVector4R p4meson, p4lepton, p4nu, p4w;
  double q2max;

  double q2, elepton, plepton;
  int i,j;
  double erho,prho,costl;

  double maxfoundprob = 0.0;
  double prob = -10.0;
  int massiter;

  for (massiter=0;massiter<3;massiter++){

    mass[0] = EvtPDL::getMeanMass(meson);
    mass[1] = EvtPDL::getMeanMass(lepton);
    mass[2] = EvtPDL::getMeanMass(nudaug);
    if ( massiter==1 ) {
      mass[0] = EvtPDL::getMinMass(meson);
    }
    if ( massiter==2 ) {
      mass[0] = EvtPDL::getMaxMass(meson);
      if ( (mass[0]+mass[1]+mass[2])>m) mass[0]=m-mass[1]-mass[2]-0.00001; 
    }

    q2max = (m-mass[0])*(m-mass[0]);
    
    //loop over q2

    for (i=0;i<25;i++) {
      q2 = ((i+0.5)*q2max)/25.0;
      
      erho = ( m*m + mass[0]*mass[0] - q2 )/(2.0*m);
      
      prho = sqrt(erho*erho-mass[0]*mass[0]);
      
      p4meson.set(erho,0.0,0.0,-1.0*prho);
      p4w.set(m-erho,0.0,0.0,prho);
      
      //This is in the W rest frame
      elepton = (q2+mass[1]*mass[1])/(2.0*sqrt(q2));
      plepton = sqrt(elepton*elepton-mass[1]*mass[1]);
      
      double probctl[3];

      for (j=0;j<3;j++) {
	
        costl = 0.99*(j - 1.0);
	
	//These are in the W rest frame. Need to boost out into
	//the B frame.
        p4lepton.set(elepton,0.0,
		  plepton*sqrt(1.0-costl*costl),plepton*costl);
        p4nu.set(plepton,0.0,
		 -1.0*plepton*sqrt(1.0-costl*costl),-1.0*plepton*costl);

	EvtVector4R boost((m-erho),0.0,0.0,1.0*prho);
        p4lepton=boostTo(p4lepton,boost);
        p4nu=boostTo(p4nu,boost);

	//Now initialize the daughters...

        daughter->init(meson,p4meson);
        lep->init(lepton,p4lepton);
        trino->init(nudaug,p4nu);

        CalcAmp(root_part,amp,FormFactors);

	//Now find the probability at this q2 and cos theta lepton point
        //and compare to maxfoundprob.

	//Do a little magic to get the probability!!
	prob = rho.normalizedProb(amp.getSpinDensity());

	probctl[j]=prob;
      }

      //probclt contains prob at ctl=-1,0,1.
      //prob=a+b*ctl+c*ctl^2

      double a=probctl[1];
      double b=0.5*(probctl[2]-probctl[0]);
      double c=0.5*(probctl[2]+probctl[0])-probctl[1];

      prob=probctl[0];
      if (probctl[1]>prob) prob=probctl[1];
      if (probctl[2]>prob) prob=probctl[2];

      if (fabs(c)>1e-20){
	double ctlx=-0.5*b/c;
	if (fabs(ctlx)<1.0){
	  double probtmp=a+b*ctlx+c*ctlx*ctlx;
	  if (probtmp>prob) prob=probtmp;
	} 

      }

      if ( prob > maxfoundprob ) {
	maxfoundprob = prob; 
      }

    }
    if ( EvtPDL::getWidth(meson) <= 0.0 ) {
      //if the particle is narrow dont bother with changing the mass.
      massiter = 4;
    }

  }
  root_part->deleteTree();  

  maxfoundprob *=1.1;
  return maxfoundprob;
  
}

