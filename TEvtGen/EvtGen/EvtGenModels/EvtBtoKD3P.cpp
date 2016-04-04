//-----------------------------------------------------------------------
// File and Version Information: 
//      $Id: EvtBtoKD3P.cpp,v 1.2 2009-04-02 15:22:28 robbep Exp $
// 
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations. If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information:
//      Copyright (C) 2003, Colorado State University
//
// Module creator:
//      Abi soffer, CSU, 2003
//-----------------------------------------------------------------------
#include "EvtGenBase/EvtPatches.hh"

// Decay model that does the decay B+->D0K, D0->3 psudoscalars

#include <assert.h>

#include "EvtGenModels/EvtBtoKD3P.hh"
#include "EvtGenBase/EvtDecayTable.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenModels/EvtPto3P.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtDalitzPoint.hh"
#include "EvtGenBase/EvtCyclic3.hh"
using std::endl;

//------------------------------------------------------------------
EvtBtoKD3P::EvtBtoKD3P() :
  _model1(0),
  _model2(0),
  _decayedOnce(false)
{
}

//------------------------------------------------------------------
EvtBtoKD3P::EvtBtoKD3P(const EvtBtoKD3P & other) : 
    EvtDecayAmp( other ){
}

//------------------------------------------------------------------
EvtBtoKD3P::~EvtBtoKD3P(){
}

//------------------------------------------------------------------
EvtDecayBase * EvtBtoKD3P::clone(){ 
  return new EvtBtoKD3P(); 
} 

//------------------------------------------------------------------
std::string EvtBtoKD3P::getName(){
  return "BTOKD3P";     
}

//------------------------------------------------------------------
void EvtBtoKD3P::init(){
  checkNArg(2); // r, phase
  checkNDaug(3); // K, D0(allowed), D0(suppressed). 
                 // The last two daughters are really one particle

  // check that the mother and all daughters are scalars:
  checkSpinParent  (  EvtSpinType::SCALAR);
  checkSpinDaughter(0,EvtSpinType::SCALAR);
  checkSpinDaughter(1,EvtSpinType::SCALAR);
  checkSpinDaughter(2,EvtSpinType::SCALAR);

  // Check that the B dtr types are K D D:

  // get the parameters:
  _r = getArg(0);
  double phase = getArg(1);
  _exp = EvtComplex(cos(phase), sin(phase));
}

//------------------------------------------------------------------
void EvtBtoKD3P::initProbMax(){
  setProbMax(1); // this is later changed in decay()
}

//------------------------------------------------------------------
void EvtBtoKD3P::decay(EvtParticle *p){
  // tell the subclass that we decay the daughter:
  _daugsDecayedByParentModel = true;

  // the K is the 1st daughter of the B EvtParticle.
  // The decay mode of the allowed D (the one produced in b->c decay) is 2nd
  // The decay mode of the suppressed D (the one produced in b->u decay) is 3rd
  const int KIND = 0;
  const int D1IND = 1;
  const int D2IND = 2;

  // generate kinematics of daughters (K and D): 
  EvtId tempDaug[2] = {getDaug(KIND), getDaug(D1IND)};
  p->initializePhaseSpace(2, tempDaug);  

  // Get the D daughter particle and the decay models of the allowed
  // and suppressed D modes:
  EvtParticle * theD = p->getDaug(D1IND); 
  EvtPto3P * model1 = (EvtPto3P*)(EvtDecayTable::getInstance()->getDecayFunc(theD));

  // for the suppressed mode, re-initialize theD as the suppressed D alias:
  theD->init(getDaug(D2IND), theD->getP4());
  EvtPto3P * model2 = (EvtPto3P*)(EvtDecayTable::getInstance()->getDecayFunc(theD));

  // on the first call: 
  if (false == _decayedOnce) {
    _decayedOnce = true;

    // store the D decay model pointers:
    _model1 = model1;
    _model2 = model2;

    // check the decay models of the first 2 daughters and that they
    // have the same final states:
    std::string name1=model1->getName();
    std::string name2=model2->getName();

    if (name1 != "PTO3P") {
      report(Severity::Error,"EvtGen") 
	<< "D daughters of EvtBtoKD3P decay must decay via the \"PTO3P\" model"
	<< endl
	<< "    but found to decay via " << name1.c_str() 
	<< " or " << name2.c_str() 
	<< ". Will terminate execution!" << endl;
      assert(0);
    }

    EvtId * daugs1 = model1->getDaugs();
    EvtId * daugs2 = model2->getDaugs();
    
    bool idMatch = true;
    int d;
    for (d = 0; d < 2; ++d) {
      if (daugs1[d] != daugs2[d]) {
	idMatch = false;
      }
    }
    if (false == idMatch) {
      report(Severity::Error,"EvtGen") 
	<< "D daughters of EvtBtoKD3P decay must decay to the same final state"
	<< endl
	<< "   particles in the same order (not CP-conjugate order)," << endl
	<< "   but they were found to decay to" << endl;
      for (d = 0; d < model1->getNDaug(); ++d) {
	report(Severity::Error,"") << "   " << EvtPDL::name(daugs1[d]).c_str() << " ";
      }
      report(Severity::Error,"") << endl;
      for (d = 0; d < model1->getNDaug(); ++d) {
	report(Severity::Error,"") << "   "  << EvtPDL::name(daugs2[d]).c_str() << " ";
      }
      report(Severity::Error,"") << endl << ". Will terminate execution!" << endl;
      assert(0);
    }      

    // estimate the probmax. Need to know the probmax's of the 2
    // models for this:
    setProbMax(model1->getProbMax(0) 
	       + _r * _r * model2->getProbMax(0)
	       + 2 * _r * sqrt(model1->getProbMax(0) * model2->getProbMax(0)));

  } // end of things to do on the first call
  
  // make sure the models haven't changed since the first call:
  if (_model1 != model1 || _model2 != model2) {
    report(Severity::Error,"EvtGen") 
      << "D daughters of EvtBtoKD3P decay should have only 1 decay modes, "
      << endl
      << "    but a new decay mode was found after the first call" << endl
      << "    Will terminate execution!" << endl;
    assert(0);
  }

  // get the cover function for each of the models and add them up.
  // They are summed with coefficients 1 because we are willing to
  // take a small inefficiency (~50%) in order to ensure that the
  // cover function is large enough without getting into complications
  // associated with the smallness of _r:
  EvtPdfSum<EvtDalitzPoint> * pc1 = model1->getPC();
  EvtPdfSum<EvtDalitzPoint> * pc2 = model2->getPC();
  EvtPdfSum<EvtDalitzPoint> pc;
  pc.addTerm(1.0, *pc1);
  pc.addTerm(1.0, *pc2);

  // from this combined cover function, generate the Dalitz point:
  EvtDalitzPoint x = pc.randomPoint();

  // get the aptitude for each of the models on this point and add them up:
  EvtComplex amp1 = model1->amplNonCP(x);
  EvtComplex amp2 = model2->amplNonCP(x);
  EvtComplex amp = amp1 + amp2 * _r * _exp;

  // get the value of the cover function for this point and set the
  // relative amplitude for this decay:

  double comp = sqrt(pc.evaluate (x));
  vertex (amp/comp);
  
  // Make the daughters of theD:
  bool massTreeOK = theD->generateMassTree();
  if (massTreeOK == false) {return;}

  // Now generate the p4's of the daughters of theD:
  std::vector<EvtVector4R> v = model2->initDaughters(x);  
  
  if(v.size() != theD->getNDaug()) {    
    report(Severity::Error,"EvtGen") 
      << "Number of daughters " << theD->getNDaug() 
      << " != " << "Momentum vector size " << v.size() 
      << endl
      << "     Terminating execution." << endl;
    assert(0);
  }
  
  // Apply the new p4's to the daughters:
  for(unsigned int i=0; i<theD->getNDaug(); ++i){
    theD->getDaug(i)->init(model2->getDaugs()[i], v[i]);
  }    
}

