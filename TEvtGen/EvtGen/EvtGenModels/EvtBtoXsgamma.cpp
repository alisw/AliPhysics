//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Module: EvtBtoXsgamma.cc
//
// Description: Routine to perform two-body non-resonant B->Xs,gamma decays.
//              Set the first input parameter to 1 to use the Ali-Greub model,
//              or 2 to use the Kagan-Neubert model.
//
// Modification history:
//
//    Mark Ian Williams       July 20, 2000       Module created
//    Mark Ian Williams       July 21, 2000       Module works
//    Mark Ian Williams       July 25, 2000       Works for all Xs modes
//    Mark Ian Williams       Aug  09, 2000       New values for mass minima
//    Mark Ian Williams       Sept 06, 2000       14 parameter M_Xs function
//    Mark Ian Williams       Sept 07, 2000       18 parameter M_Xs function
//    Mark Ian Williams       Sept 07, 2000       Tidied up the code
//    Mark Ian Williams       Sept 10, 2000       Updated parameters
//    Mark Ian Williams       Sept 11, 2000       Finalised code
//    Jane Tinslay            March 21, 2001      Re-worked so that you can choose
//                                                between the Ali-Greub and Kagan-Neubert
//                                                Modules.                          
//------------------------------------------------------------------------
//

#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtBtoXsgamma.hh"
#include <string>
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenModels/EvtBtoXsgammaAliGreub.hh"
#include "EvtGenModels/EvtBtoXsgammaKagan.hh"
#include "EvtGenModels/EvtBtoXsgammaFixedMass.hh"
#include "EvtGenModels/EvtBtoXsgammaFlatEnergy.hh"
using std::endl;

EvtBtoXsgamma::~EvtBtoXsgamma() {

  delete _model; _model=0;

}

std::string EvtBtoXsgamma::getName(){

  return "BTOXSGAMMA";     

}

EvtDecayBase* EvtBtoXsgamma::clone(){

  return new EvtBtoXsgamma;

}

void EvtBtoXsgamma::init(){
  //Arguments:
  // 0: Ali-Greub model = 1, Kagan model = 2
 //No more arguments for Ali-Greub model
  // 1:
  // 2:
  // 3:

  // check that at least one b->sg model has been selected
  if (getNArg() == 0) {
    
    report(Severity::Error,"EvtGen") << "EvtBtoXsgamma generator expected "
                           << " at least 1 argument but found: "<<getNArg()<<endl;
    report(Severity::Error,"EvtGen") << "Will terminate execution!"<<endl;
    ::abort();
  }    
}

void EvtBtoXsgamma::initProbMax(){

  noProbMax();

}

void EvtBtoXsgamma::decay( EvtParticle *p ){

  //initialize here. -- its too damn slow otherwise.

  if ( _model == 0 ) {
    
    if (getArg(0) == 1) _model = new EvtBtoXsgammaAliGreub();
    else if (getArg(0) == 2) _model = new EvtBtoXsgammaKagan();
    else if (getArg(0) == 3) _model = new EvtBtoXsgammaFixedMass();
    else if (getArg(0) == 4) _model = new EvtBtoXsgammaFlatEnergy();
    else{
      report(Severity::Error,"EvtGen") << "No valid EvtBtoXsgamma generator model selected "
			     << "Set arg(0) to 1 for Ali-Greub model or 2 for "
			     <<" Kagan model or 3 for a fixed mass"<<endl;
      report(Severity::Error,"EvtGen") << "Will terminate execution!"<<endl;
      ::abort();
      
    }
    _model->init(getNArg(),getArgs());
  }


  //  if ( p->getNDaug() != 0 ) {
    //Will end up here because maxrate multiplies by 1.2
  //  report(Severity::Debug,"EvtGen") << "In EvtBtoXsgamma: X_s daughters should not be here!"<<endl;
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

  mass[1] = EvtPDL::getMass(getDaug(1));
 
  int Xscode = EvtPDL::getStdHep(getDaug(0));
   
  mass[0] = _model->GetMass(Xscode);

  EvtGenKine::PhaseSpace( getNDaug(), mass, p4, m_b );

  for(i=0;i<getNDaug();i++){
     pdaug[i]->init( getDaugs()[i], p4[i] );
  }

}

