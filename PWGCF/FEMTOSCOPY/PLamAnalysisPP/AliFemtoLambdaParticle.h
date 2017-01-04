#ifndef ALIFEMTOLambdaParticle_H
#define ALIFEMTOLambdaParticle_H
//
//Class AliFemtoLambdaParticle, AliFemtoLambdaEvent, AliFemtoLambdaEventCollection
//
//AliFemtoLambdaParticle, AliFemtoLambdaEvent, AliFemtoLambdaEventCollection
//


#include <iostream>
#include <string>
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TBits.h"
#include "TObject.h"
#include "TVector2.h"
#include "AliESDtrack.h"
#include "TVector3.h"

using namespace std;

class AliFemtoLambdaParticle // Reconstructed Lambdas parameters needed for correlations
{
 public:
  
  AliFemtoLambdaParticle();
  virtual ~AliFemtoLambdaParticle();
  //AliFemtoLambdaParticle(const AliFemtoLambdaParticle &obj);
  AliFemtoLambdaParticle &operator=(const AliFemtoLambdaParticle &obj);
  
  TVector3 fMomentum;  //v0 momentum
  TVector3 fMomentumMC;  //v0 Monte Carlo momentum
  TVector3 fMomentumMCMother;  //v0 Monte Carlo momentum of mother particle
  int fPDGCode;
  int fPDGCodeMother;
  double fPt;           //v0 transverse momentum
  double fMass;         //v0 reconstructed mass
  short fDaughterID1;   //Daughter (proton) AODtrack ID
  short fDaughterID2;   //Daughter (pion) AODtrack ID
  bool fV0tag;
  double fPointing;
  bool fReal;

  //stuff related to daughter tracks:
  TVector3 fMomentumPosDaughter;//momentum of positive daughter particle
  TVector3 fMomentumNegDaughter;//momentum of negative daughter particle
  double fPhiPosdaughter;
  double fEtaPosdaughter;
  double fPhiNegdaughter;
  double fEtaNegdaughter;
  double fPosDaughPosTPC[9][3];
  double fNegDaughPosTPC[9][3];
  
#ifdef __ROOT__
  ClassDef(AliFemtoLambdaParticle, 1);
#endif
};
#endif
