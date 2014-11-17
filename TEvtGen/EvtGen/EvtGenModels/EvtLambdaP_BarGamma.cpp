//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2003      Caltech
//
// Module: EvtGen/EvtRadiativeBaryonicPenguins.hh
//
// Description:Implementation of the decay B- -> lambda p_bar gamma according to
// Cheng, Yang; hep-ph/0201015
//
// Modification history:
//
//    JFS     December 16th, 2003         Module created
//
//------------------------------------------------------------------------
#include "EvtGenBase/EvtPatches.hh"

#include "EvtGenModels/EvtLambdaP_BarGamma.hh"
#include "EvtGenBase/EvtGammaMatrix.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtSpinType.hh"
#include "EvtGenBase/EvtDiracParticle.hh"
#include "EvtGenBase/EvtPhotonParticle.hh"
#include <stdlib.h>
using std::cout;
using std::endl;

EvtLambdaP_BarGamma::EvtLambdaP_BarGamma() :
  _mLambdab   ( 5.624),            // Lambda_b mass
  _mLambda0    ( 1.115684),         // Lambda0 mass
  _c7Eff       ( -0.31),            // Wilson coefficient                                      
  _mb          (  4.4),             // running b mass                                          
  _mV          (  5.42),            // pole mass vector current                                
  _mA          (  5.86),            // pole mass axial current                                 
  _GF          (  1.166E-5),        // Fermi constant                                          
  _gLambdab   (  16),               // coupling constant Lambda_b -> B- p
  _e0          (  1),               // electromagnetic coupling (+1)                           
  _g1          (  0.64),            // heavy-light form factors at q_mSqare                    
  _g2          ( -0.10),     
  _f1          (  0.64),
  _f2          ( -0.31),
  _VtbVtsStar ( 0.038)            // |V_tb V_ts^*|
{
}



std::string EvtLambdaP_BarGamma::getName(){
    return "B_TO_LAMBDA_PBAR_GAMMA";
}

EvtDecayBase* EvtLambdaP_BarGamma::clone(){
    return new EvtLambdaP_BarGamma;
}

void EvtLambdaP_BarGamma::init() {
    // no arguments, daughter lambda p_bar gamma
    checkNArg(0);
    checkNDaug(3);
    
    checkSpinParent(EvtSpinType::SCALAR);
    checkSpinDaughter(0, EvtSpinType::DIRAC);
    checkSpinDaughter(1, EvtSpinType::DIRAC);
    checkSpinDaughter(2, EvtSpinType::PHOTON);    
}


// initialize phasespace and calculate the amplitude
void EvtLambdaP_BarGamma::decay(EvtParticle* p) {
    EvtComplex I(0, 1);
    
    p->initializePhaseSpace(getNDaug(), getDaugs());
    
    EvtDiracParticle* theLambda = static_cast<EvtDiracParticle*>(p->getDaug(0));
    EvtVector4R lambdaMomentum = theLambda->getP4Lab();
    
    EvtDiracParticle* theAntiP = static_cast<EvtDiracParticle*>(p->getDaug(1));
    
    EvtPhotonParticle* thePhoton = static_cast<EvtPhotonParticle*>(p->getDaug(2));
    EvtVector4R photonMomentum = thePhoton->getP4Lab();     // get momentum in the same frame
        
    // loop over all possible spin states
    for (int i=0; i<2; ++i) {
      EvtDiracSpinor lambdaPol = theLambda->spParent(i);
      for (int j=0; j<2; ++j)  {
	EvtDiracSpinor antiP_Pol = theAntiP->spParent(j);
	for (int k=0; k<2; ++k) {
	  EvtVector4C photonPol = thePhoton->epsParentPhoton(k); // one of two possible polarization states
	  EvtGammaMatrix photonGamma; // sigma[mu][nu] * epsilon[mu] * k[nu] (watch lower indices)
	  for (int mu=0; mu<4; ++mu)
	    for (int nu=0; nu<4; ++nu)
	      photonGamma += EvtGammaMatrix::sigmaLower(mu, nu) * photonPol.get(mu) * photonMomentum.get(nu);
	  
	  EvtComplex amp = 
	    -I*_gLambdab * lambdaPol.adjoint() * 
      ((constA()*EvtGammaMatrix::id() + constB()*EvtGammaMatrix::g5())
       * photonGamma * (EvtGenFunctions::slash(lambdaMomentum) + 
                        EvtGenFunctions::slash(photonMomentum) + 
                        _mLambdab*EvtGammaMatrix::id())
       / ((lambdaMomentum + photonMomentum)*(lambdaMomentum + photonMomentum) - _mLambdab*_mLambdab)
       * EvtGammaMatrix::g5() * antiP_Pol);
	  // use of parentheses so I do not have to define EvtDiracSpinor*EvtGammaMatrix, which shouldn't be defined to prevent errors in indexing

	  vertex(i, j, k, amp);
	}
      }
    }
}

void EvtLambdaP_BarGamma::initProbMax()
{
    // setProbMax(1);
    setProbMax(9.0000E-13); // found by trial and error
}

// form factors at 0
double EvtLambdaP_BarGamma::f0(double fqm, int n) const {
    return fqm * pow(1 - pow(_mLambdab - _mLambda0, 2) / (_mV * _mV), n);
}

double EvtLambdaP_BarGamma::g0(double gqm, int n) const {
    return gqm * pow(1 - pow(_mLambdab - _mLambda0, 2) / (_mA * _mA), n);
}

double EvtLambdaP_BarGamma::constA() const {
    return _GF/sqrt(2.) * _e0 / (8 * EvtConst::pi*EvtConst::pi) * 2 * _c7Eff * _mb * _VtbVtsStar
        * (f0(_f1) - f0(_f2));
}

double EvtLambdaP_BarGamma::constB() const {
    return _GF/sqrt(2.) * _e0 / (8 * EvtConst::pi*EvtConst::pi) * 2 * _c7Eff * _mb * _VtbVtsStar
        * (g0(_g1) - (_mLambdab - _mLambda0) / (_mLambdab + _mLambda0) * g0(_g2));
}
