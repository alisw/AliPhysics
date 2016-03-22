//--------------------------------------------------------------------------
//
//
// Copyright Information: See EvtGen/COPYRIGHT
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Module: EvtBtoXsEtap.cc
//
// Description: Routine to perform two-body non-resonant B->Xs,gluon decays.
// It generates an X_s mass spectrum based on a parameterisation of the
// b->s,gluon  spectrum of Atwood-Soni. The resultant X_s particles may
// be decayed by JETSET.
//
// Modification history:
//
//    Adlene Hicheur       January 10, 2001       Module created
//------------------------------------------------------------------------
//
#include "EvtGenBase/EvtPatches.hh"

#include <stdlib.h>
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtBtoXsEtap.hh"
#include <string>
#include "EvtGenBase/EvtConst.hh"
using std::endl;

EvtBtoXsEtap::~EvtBtoXsEtap() {}

std::string EvtBtoXsEtap::getName(){

  return "BTOXSETAP";     

}

EvtDecayBase* EvtBtoXsEtap::clone(){

  return new EvtBtoXsEtap;

}


void EvtBtoXsEtap::init(){

  // check that there are no arguments

  checkNArg(0);
}

void EvtBtoXsEtap::initProbMax(){

  noProbMax();

}

void EvtBtoXsEtap::decay( EvtParticle *p ){

  // useless
  //  if ( p->getNDaug() != 0 ) {
  //  //Will end up here because maxrate multiplies by 1.2
  //  report(Severity::Debug,"EvtGen") << "In EvtBtoXsEtap: X_s daughters should not be here!"<<endl;
  //  return;
  //}

  double m_b;
  int i;
  p->makeDaughters(getNDaug(),getDaugs());
  EvtParticle *pdaug[MAX_DAUG];

  for(i=0;i<getNDaug();i++){
     pdaug[i]=p->getDaug(i);   
  }

  static EvtVector4R p4[MAX_DAUG];
  static double mass[MAX_DAUG];

  m_b = p->mass();

  // Prepare for phase space routine.

  mass[1] = EvtPDL::getMass(getDaug(1));

  double xbox, ybox, min, max,hichfit;
  min=0.493;
  max=4.3;
  const double TwoPi = EvtConst::twoPi;
  int Xscode = EvtPDL::getStdHep(getDaug(0));

  // A five parameters fit, the shape is taken from Atwood & Soni

  //  double par[18];
  double par[6];
  if ((Xscode == 30343) || (Xscode == -30343) || 
      (Xscode == 30353) || (Xscode == -30353)) { // Xsu or Xsd
    min=0.6373; //  Just above K pi threshold for Xsd/u
    //min=0.6333; //  K pi threshold for neutral Xsd
    //    par[0]=-2057.2380371094;
    par[0]=2.36816;
    //    par[1]=2502.2556152344;
    par[1]=0.62325725;
    //    par[2]=1151.5632324219;
    par[2]=2.2;
    //    par[3]=0.82431584596634;
    par[3]=-0.2109375;
    //    par[4]=-4110.5234375000;
    par[4]=2.7;
    //    par[5]=8445.6757812500;
    par[5]=0.54;
    //    par[6]=-3034.1894531250;
    //    par[7]=1.1557708978653;
    //    par[8]=1765.9311523438;
    //    par[9]=1.3730158805847;
    //    par[10]=0.51371538639069;
    //    par[11]=2.0056934356689;
    //    par[12]=37144.097656250;
    //    par[13]=-50296.781250000;
    //    par[14]=27319.095703125;
    //    par[15]=-7408.0678710938;
    //    par[16]=1000.8093261719;
    //    par[17]=-53.834449768066;
  } else {
    report(Severity::Debug,"EvtGen") << "In EvtBtoXsEtap: Particle with id " << Xscode << " is not a Xsd/u particle"<<endl;
    return;
  }

  double boxheight=par[5];
  double boxwidth=max-min;

  mass[0]=0.0;
  while ((mass[0] > max) || (mass[0] < min)){
    xbox = EvtRandom::Flat(boxwidth)+min;
    ybox=EvtRandom::Flat(boxheight);
    if (xbox<par[2]) {

      hichfit=(1/sqrt(TwoPi*par[1]))*exp(-0.5*pow((xbox-par[0])/par[1],2));      
      //      alifit=par[0]+par[1]*xbox+par[2]*pow(xbox,2);
      //    } else if (xbox<par[7]) {
      //      alifit=par[4]+par[5]*xbox+par[6]*pow(xbox,2);
      //    } else if (xbox<par[11]) {
      //      alifit=par[8]*exp(-0.5*pow((xbox-par[9])/par[10],2));
    } else {
      hichfit=par[3]*pow((xbox-par[4]),2)+par[5];
      //      alifit=par[12]+par[13]*xbox+par[14]*pow(xbox,2)+par[15]*pow(xbox,3)+par[16]*pow(xbox,4)+par[17]*pow(xbox,5);
    }
    if (ybox>hichfit) {
      mass[0]=0.0;
    } else {
      mass[0]=xbox;
    }
  }

  // debug stuff:  report(Severity::Info,"EvtGen") << "Xscode " << Xscode << " daughter 1 mass " << mass[0] << " daughter 2 mass " << mass[1] << endl;

  EvtGenKine::PhaseSpace( getNDaug(), mass, p4, m_b );

  for(i=0;i<getNDaug();i++){
     pdaug[i]->init( getDaugs()[i], p4[i] );
  }

  return ;
}


