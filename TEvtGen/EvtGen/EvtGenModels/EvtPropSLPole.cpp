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
// Module: EvtPropSLPole.cc
//
// Description: Routine to implement semileptonic decays according
//              to light cone sum rules
//
// Modification history:
//
//    DJL       April 23, 1998       Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtPropSLPole.hh"
#include "EvtGenModels/EvtSLPoleFF.hh"
#include "EvtGenBase/EvtSemiLeptonicScalarAmp.hh"
#include "EvtGenBase/EvtSemiLeptonicVectorAmp.hh"
#include "EvtGenBase/EvtSemiLeptonicTensorAmp.hh"
#include "EvtGenBase/EvtIntervalFlatPdf.hh"
#include "EvtGenBase/EvtScalarParticle.hh"
#include "EvtGenBase/EvtVectorParticle.hh"
#include "EvtGenBase/EvtTensorParticle.hh"
#include "EvtGenBase/EvtTwoBodyVertex.hh"
#include "EvtGenBase/EvtPropBreitWignerRel.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtAmpPdf.hh"
#include "EvtGenBase/EvtMassAmp.hh"
#include "EvtGenBase/EvtSpinType.hh"
#include "EvtGenBase/EvtDecayTable.hh"
#include <string>

EvtPropSLPole::~EvtPropSLPole() {}

std::string EvtPropSLPole::getName(){

  return "PROPSLPOLE";     

}


EvtDecayBase* EvtPropSLPole::clone(){

  return new EvtPropSLPole;

}

void EvtPropSLPole::decay( EvtParticle *p ){

  if(! _isProbMaxSet){

     EvtId parnum,mesnum,lnum,nunum;

     parnum = getParentId();
     mesnum = getDaug(0);
     lnum = getDaug(1);
     nunum = getDaug(2);

     double mymaxprob = calcMaxProb(parnum,mesnum,
                           lnum,nunum,SLPoleffmodel);

     setProbMax(mymaxprob);

     _isProbMaxSet = true;

  }

  double minKstMass = EvtPDL::getMinMass(p->getDaug(0)->getId());
  double maxKstMass = EvtPDL::getMaxMass(p->getDaug(0)->getId());

  EvtIntervalFlatPdf flat(minKstMass, maxKstMass);
  EvtPdfGen<EvtPoint1D> gen(flat);
  EvtPoint1D point = gen(); 
 
  double massKst = point.value();

  p->getDaug(0)->setMass(massKst);
  p->initializePhaseSpace(getNDaug(),getDaugs());

//  EvtVector4R p4meson = p->getDaug(0)->getP4();

  calcamp->CalcAmp(p,_amp2,SLPoleffmodel); 

  EvtParticle *mesonPart = p->getDaug(0);
  
  double meson_BWAmp = calBreitWigner(mesonPart, point);  

  int list[2];
  list[0]=0; list[1]=0;
  _amp2.vertex(0,0,_amp2.getAmp(list)*meson_BWAmp);
  list[0]=0; list[1]=1;
  _amp2.vertex(0,1,_amp2.getAmp(list)*meson_BWAmp);

  list[0]=1; list[1]=0;
  _amp2.vertex(1,0,_amp2.getAmp(list)*meson_BWAmp);
  list[0]=1; list[1]=1;
  _amp2.vertex(1,1,_amp2.getAmp(list)*meson_BWAmp);

  list[0]=2; list[1]=0;
  _amp2.vertex(2,0,_amp2.getAmp(list)*meson_BWAmp);
  list[0]=2; list[1]=1;
  _amp2.vertex(2,1,_amp2.getAmp(list)*meson_BWAmp);
     
  
  return;

}

void EvtPropSLPole::initProbMax(){

  _isProbMaxSet = false;

  return;

}


void EvtPropSLPole::init(){
  
  checkNDaug(3);

  //We expect the parent to be a scalar 
  //and the daughters to be X lepton neutrino

  checkSpinParent(EvtSpinType::SCALAR);
  checkSpinDaughter(1,EvtSpinType::DIRAC);
  checkSpinDaughter(2,EvtSpinType::NEUTRINO);

  EvtSpinType::spintype mesontype=EvtPDL::getSpinType(getDaug(0));

  SLPoleffmodel = new EvtSLPoleFF(getNArg(),getArgs());
  
  if ( mesontype==EvtSpinType::SCALAR ) { 
    calcamp = new EvtSemiLeptonicScalarAmp; 
  }
  if ( mesontype==EvtSpinType::VECTOR ) { 
    calcamp = new EvtSemiLeptonicVectorAmp; 
  }
  if ( mesontype==EvtSpinType::TENSOR ) { 
    calcamp = new EvtSemiLeptonicTensorAmp; 
  }

}


double EvtPropSLPole::calBreitWignerBasic(double maxMass){

  if ( _width< 0.0001) return 1.0;
  //its not flat - but generated according to a BW

  double mMin=_massMin;
  double mMax=_massMax;
  if ( maxMass>-0.5 && maxMass< mMax) mMax=maxMass;

  double massGood = EvtRandom::Flat(mMin, mMax);

  double ampVal = sqrt(1.0/(pow(massGood-_mass, 2.0) + pow(_width, 2.0)/4.0));

  return ampVal;

}


double EvtPropSLPole::calBreitWigner(EvtParticle *pmeson, EvtPoint1D point){

  EvtId mesnum = pmeson->getId(); 
  double _mass = EvtPDL::getMeanMass(mesnum);
  double _width = EvtPDL::getWidth(mesnum);
  double _maxRange = EvtPDL::getMaxRange(mesnum);
  EvtSpinType::spintype mesontype=EvtPDL::getSpinType(mesnum);
  _includeDecayFact=true;
  _includeBirthFact=true;
  _spin = mesontype;
  _blatt = 3.0;

  double maxdelta = 15.0*_width;

  if ( _maxRange > 0.00001 ) {
    _massMax=_mass+maxdelta;
    _massMin=_mass-_maxRange;
  }
  else{
    _massMax=_mass+maxdelta;
    _massMin=_mass-15.0*_width;
  }

  _massMax=_mass+maxdelta;
  if ( _massMin< 0. ) _massMin=0.;


  EvtParticle* par=pmeson->getParent();
  double maxMass=-1.;
  if ( par != 0 ) {
    if ( par->hasValidP4() ) maxMass=par->mass();
    for ( size_t i=0;i<par->getNDaug();i++) {
      EvtParticle *tDaug=par->getDaug(i);
      if ( pmeson != tDaug )
        maxMass-=EvtPDL::getMinMass(tDaug->getId());
    }
  }

  EvtId *dauId=0;
  double *dauMasses=0;
  size_t nDaug = pmeson->getNDaug();
  if ( nDaug > 0) {
     dauId=new EvtId[nDaug];
     dauMasses=new double[nDaug];
     for (size_t j=0;j<nDaug;j++) {
       dauId[j]=pmeson->getDaug(j)->getId();
       dauMasses[j]=pmeson->getDaug(j)->mass();
     }
   }
   EvtId *parId=0;
   EvtId *othDaugId=0;
   EvtParticle *tempPar=pmeson->getParent();
   if (tempPar) {
     parId=new EvtId(tempPar->getId());
     if ( tempPar->getNDaug()==2 ) {
       if ( tempPar->getDaug(0) == pmeson ) othDaugId=new EvtId(tempPar->getDaug(1)->getId());
       else othDaugId=new EvtId(tempPar->getDaug(0)->getId());
     }
   }

  if ( nDaug!=2) return calBreitWignerBasic(maxMass);

  if ( _width< 0.00001) return 1.0;

  //first figure out L - take the lowest allowed.

  EvtSpinType::spintype spinD1=EvtPDL::getSpinType(dauId[0]);
  EvtSpinType::spintype spinD2=EvtPDL::getSpinType(dauId[1]);

  int t1=EvtSpinType::getSpin2(spinD1);
  int t2=EvtSpinType::getSpin2(spinD2);
  int t3=EvtSpinType::getSpin2(_spin);

  int Lmin=-10;

  // allow for special cases.
  if (Lmin<-1 ) {

    //There are some things I don't know how to deal with
    if ( t3>4) return calBreitWignerBasic(maxMass);
    if ( t1>4) return calBreitWignerBasic(maxMass);
    if ( t2>4) return calBreitWignerBasic(maxMass);

    //figure the min and max allowwed "spins" for the daughters state
    Lmin=std::max(t3-t2-t1,std::max(t2-t3-t1,t1-t3-t2));
    if (Lmin<0) Lmin=0;
    assert(Lmin==0||Lmin==2||Lmin==4);
  }

  //double massD1=EvtPDL::getMeanMass(dauId[0]);
  //double massD2=EvtPDL::getMeanMass(dauId[1]);
  double massD1=dauMasses[0];
  double massD2=dauMasses[1];

  // I'm not sure how to define the vertex factor here - so retreat to nonRel code.
  if ( (massD1+massD2)> _mass ) return  calBreitWignerBasic(maxMass);

  //parent vertex factor not yet implemented
  double massOthD=-10.;
  double massParent=-10.;
  int birthl=-10;
  if ( othDaugId) {
    EvtSpinType::spintype spinOth=EvtPDL::getSpinType(*othDaugId);
    EvtSpinType::spintype spinPar=EvtPDL::getSpinType(*parId);

    int tt1=EvtSpinType::getSpin2(spinOth);
    int tt2=EvtSpinType::getSpin2(spinPar);
    int tt3=EvtSpinType::getSpin2(_spin);

    //figure the min and max allowwed "spins" for the daughters state
    if ( (tt1<=4) && ( tt2<=4) ) {
      birthl=std::max(tt3-tt2-tt1,std::max(tt2-tt3-tt1,tt1-tt3-tt2));
      if (birthl<0) birthl=0;

      massOthD=EvtPDL::getMeanMass(*othDaugId);
      massParent=EvtPDL::getMeanMass(*parId);

    }

  }
  double massM=_massMax;
  if ( (maxMass > -0.5) && (maxMass < massM) ) massM=maxMass;

  //special case... if the parent mass is _fixed_ we can do a little better
  //and only for a two body decay as that seems to be where we have problems

  // Define relativistic propagator amplitude

  EvtTwoBodyVertex vd(massD1,massD2,_mass,Lmin/2);
  vd.set_f(_blatt);
  EvtPropBreitWignerRel bw(_mass,_width);
  EvtMassAmp amp(bw,vd);
//  if ( _fixMassForMax) amp.fixUpMassForMax();
//  else std::cout << "problem problem\n";
  if ( _includeDecayFact) {
    amp.addDeathFact();
    amp.addDeathFactFF();
  }
  if ( massParent>-1.) {
    if ( _includeBirthFact ) {

      EvtTwoBodyVertex vb(_mass,massOthD,massParent,birthl/2);
      amp.setBirthVtx(vb);
      amp.addBirthFact();
      amp.addBirthFactFF();
    }
  }

  EvtAmpPdf<EvtPoint1D> pdf(amp);

  double ampVal = sqrt(pdf.evaluate(point));

  if ( parId) delete parId;
  if ( othDaugId) delete othDaugId;
  if ( dauId) delete [] dauId;
  if ( dauMasses) delete [] dauMasses;

  return ampVal;
 
}


double EvtPropSLPole::calcMaxProb( EvtId parent, EvtId meson,
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
//  root_part->set_type(EvtSpinType::SCALAR);
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

  EvtDecayBase *decayer;
  decayer = EvtDecayTable::getInstance()->getDecayFunc(daughter);
  if ( decayer ) {
    daughter->makeDaughters(decayer->nRealDaughters(),decayer->getDaugs());
    for(int ii=0; ii<decayer->nRealDaughters(); ii++){
      daughter->getDaug(ii)->setMass(EvtPDL::getMeanMass(daughter->getDaug(ii)->getId()));
    }
  }  

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

        calcamp->CalcAmp(root_part,amp,FormFactors);

        EvtPoint1D *point = new EvtPoint1D(mass[0]);

        double meson_BWAmp = calBreitWigner(daughter, *point);

        int list[2];
        list[0]=0; list[1]=0;
        amp.vertex(0,0,amp.getAmp(list)*meson_BWAmp);
        list[0]=0; list[1]=1;
        amp.vertex(0,1,amp.getAmp(list)*meson_BWAmp);

        list[0]=1; list[1]=0;
        amp.vertex(1,0,amp.getAmp(list)*meson_BWAmp);
        list[0]=1; list[1]=1;
        amp.vertex(1,1,amp.getAmp(list)*meson_BWAmp);

        list[0]=2; list[1]=0;
        amp.vertex(2,0,amp.getAmp(list)*meson_BWAmp);
        list[0]=2; list[1]=1;
        amp.vertex(2,1,amp.getAmp(list)*meson_BWAmp);

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

      //report(Severity::Debug,"EvtGen") << "prob,probctl:"<<prob<<" "
      //                            << probctl[0]<<" "
      //                            << probctl[1]<<" "
      //                            << probctl[2]<<endl;

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

