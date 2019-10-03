#ifndef ALIFEMTOProtonParticle_H
#define ALIFEMTOProtonParticle_H
//
//Class AliFemtoLambdaParticle, AliFemtoLambdaEvent, AliFemtoLambdaEventCollection
//
//AliFemtoLambdaParticle, AliFemtoLambdaEvent, AliFemtoLambdaEventCollection
//authors: 
//        Dhevan Gangadharan (dhevan.raja.gangadharan@cern.ch)
//        Matthew Steinpreis (matthew.steinpreis@cern.ch)
//


#include <iostream>
#include <string>
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TBits.h"
#include "TObject.h"
#include "AliESDtrack.h"
#include "TVector3.h"


using namespace std;

//class TVector3;

class AliFemtoProtonParticle  // Proton parameters are stored in this class
{
 public:
  
  AliFemtoProtonParticle();
  virtual ~AliFemtoProtonParticle();
//  AliFemtoProtonParticle(const AliFemtoProtonParticle &obj);
  AliFemtoProtonParticle &operator=(const AliFemtoProtonParticle &obj);
 
  TVector3 fMomentum;  //proton momentum
  TVector3 fMomentumMC;  //proton monte carlo momentum
  TVector3 fMomentumMCMother;  //momentum of mother particle of proton
  TVector3 fMomentumMCMotherParton;
  TVector3 fPositionTPC[9];

  int fPDGCode;
  int fPDGCodeMother; // PDG code of the weakly decaying mother resonance in the case of feed-down
  int fPDGCodePartonMother; //PDG code of the parton that created the proton
  int fPartonMotherLabel; //Label of the particle of the corresponding mother parton defined a line above
  float fPt;           //proton transverse momentum
  short fID;   //Daughter (pion) AODtrack ID
  float fPhi;
  float fPhistar[9];
  float fEta;
  //double fPrimPosTPC[9][3];
  Bool_t fReal;
  Bool_t fProtonTag;

  ClassDef(AliFemtoProtonParticle, 3)
};


#endif
