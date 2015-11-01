////////////////////////////////////////////////////////////////////////////////
//
//  This class provides storage for event and track information which 
//  are used for same-event as well as mixed-event analyses in AliAnalysisTaskhCascadeFemto 
//
//  author: maria.nicassio@cern.ch (derived and adapted from D. Gangadharan PWGCF/FEMTOSCOPY/Chaoticity/AliChaoticityEventCollection
//                                  and J. Salzwedel PWGCF/FEMTOSCOPY/V0LamAnalysis/AliAnalysisV0LamEventCollection) 
//
////////////////////////////////////////////////////////////////////////////////



#include "AliAnalysishCascadeEventCollection.h"


AliAnalysishCascadeEventCollection::~AliAnalysishCascadeEventCollection(){

  for(int i = 0; i < fifo; i++){
     if((fEvt + i)->fReconstructedXi != NULL){
       delete [] (fEvt + i)->fReconstructedXi;
     }
     if((fEvt + i)->fReconstructedProton != NULL){
       delete [] (fEvt + i)->fReconstructedProton;
     }
   }
   delete [] fEvt;
 }
 //_____________________________________________________________________________


AliAnalysishCascadeEventCollection::AliAnalysishCascadeEventCollection() : fEvt(0x0), fifo(0) {


}

//_____________________________________________________________________________

AliAnalysishCascadeEventCollection::AliAnalysishCascadeEventCollection(short eventBuffSize, int maxXiMult, int maxProtonMult) : fEvt(0x0), fifo(0) {

  SetBuffSize(eventBuffSize);

  fEvt = new AliAnalysishCascadeEvent[fifo];  //allocate pointer array of AliAnalysishCascadeEvents

  for(int ii = 0; ii < fifo; ii++){ //Initialize particle table pointers to NULL

    (fEvt + ii)->fReconstructedXi = NULL;

    (fEvt + ii)->fNumberCandidateXi = 0;

    (fEvt + ii)->fReconstructedXi = new AliReconstructedXi[maxXiMult];

    for(int j=0; j < 3; j++){

      (fEvt + ii)->fPrimaryVertex[j] = 0.;

    }

    (fEvt + ii)->fReconstructedProton = NULL;

    (fEvt + ii)->fNumberCandidateProton = 0;

    (fEvt + ii)->fReconstructedProton = new AliReconstructedProton[maxProtonMult];


  }

}


//_____________________________________________________________________________
void AliAnalysishCascadeEventCollection::FifoShift() { //Shift elements in FIFO by one and clear last element in FIFO 

  for(unsigned short i=fifo-1 ; i > 0; i--) {

    for(int j=0; j<(fEvt + i-1)->fNumberCandidateXi; j++){

      (fEvt + i)->fReconstructedXi[j] = (fEvt + i-1)->fReconstructedXi[j];

    }

    (fEvt + i)->fNumberCandidateXi = (fEvt + i-1)->fNumberCandidateXi;

    for(int j=0; j<(fEvt + i-1)->fNumberCandidateProton; j++){

      (fEvt + i)->fReconstructedProton[j] = (fEvt + i-1)->fReconstructedProton[j];

    }

    (fEvt + i)->fNumberCandidateProton = (fEvt + i-1)->fNumberCandidateProton;

    for(int j=0; j<3; j++){

      (fEvt + i)->fPrimaryVertex[j] = (fEvt + i-1)->fPrimaryVertex[j];

    }

  }

   (fEvt)->fNumberCandidateXi=0;
   (fEvt)->fNumberCandidateProton=0;

   for(int j=0; j<3; j++) {

     (fEvt)->fPrimaryVertex[j] = 0.;

   }

}



//_____________________________________________________________________________

AliReconstructedXi::AliReconstructedXi() :

   xiPt(0.),
   xiEta(0.),
   xiPhi(0.),
   xiRap(0.),
   xiMass(0.),
   isXiMass(kFALSE),
   isXi(kFALSE),
   isaXi(kFALSE),
   indexB(0),
   indexP(0),
   indexN(0),
   doSkipOver(kFALSE),
   doPickOne(kFALSE),


   mcXiOriginType(),

   daughterBacEtaS(0),
   daughterPosEtaS(0),
   daughterNegEtaS(0),
   daugherBacPhiS(0),
   daugherPosPhiS(0),
   daugherNegPhiS(0)


 {

  std::fill(xiMomentum,xiMomentum+3,0.);
  std::fill(xiMomentumTruth,xiMomentumTruth+3,0.);
  std::fill(daughterPosShiftedGlobalPosition,daughterPosShiftedGlobalPosition+3,0.);
  std::fill(daughterNegShiftedGlobalPosition,daughterNegShiftedGlobalPosition+3,0.);
  std::fill(daughterBacShiftedGlobalPosition,daughterBacShiftedGlobalPosition+3,0.);


 }

//_____________________________________________________________________________

AliReconstructedXi::~AliReconstructedXi()

 {

 }
//_____________________________________________________________________________

AliReconstructedProton::AliReconstructedProton() :

   pPt(0.),
   pEta(0.),
   pPhi(0.),
   pRap(0.),
   pCharge(0.),
   nSigmaProtonTPC(0),
   nSigmaProtonTOF(0),
   isTOFmismatch(0),
   isP(0),
   isaP(0),
   index(0),
   mcProtonOriginType(kUnassigned),
   doSkipOver(kFALSE),


   pEtaS(0.),
   pPhiS(0.)


 {
  std::fill(pMomentum,pMomentum+3,0.);
  std::fill(pMomentumTruth,pMomentumTruth+3,0.);
  std::fill(pShiftedGlobalPosition,pShiftedGlobalPosition+3,0.);
  std::fill(iptoPV,iptoPV+2,0.);
 }
//_____________________________________________________________________________

AliReconstructedProton::~AliReconstructedProton()

 {

 }

