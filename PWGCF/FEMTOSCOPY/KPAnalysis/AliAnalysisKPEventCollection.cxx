////////////////////////////////////////////////////////////////////////////////
//
//  This class provides storage for event and track information which 
//  are used for same-event as well as mixed-event analyses in AliAnalysisTaskKPFemto 
//  Author: ramona.lea@cern.ch 
//  was: maria.nicassio@cern.ch (derived and adapted from D. Gangadharan PWGCF/FEMTOSCOPY/Chaoticity/AliChaoticityEventCollection
//                                  and J. Salzwedel PWGCF/FEMTOSCOPY/V0LamAnalysis/AliAnalysisV0LamEventCollection) 
//
////////////////////////////////////////////////////////////////////////////////



#include "AliAnalysisKPEventCollection.h"

//_____________________________________________________________________________
// Default constructor 
AliReconstructedFirst::AliReconstructedFirst() :
  fPt(0),
  fEta(0),
  fTheta(0),
  fPhi(0),
  fRap(0),
  fCharge(0),
  fDCAxy(0),
  fDCAz(0),
  isTOFmismatch(kFALSE),
  isMCptc(kFALSE),
  fMCcode(0),
  fPDGcode(0),
  fMCmumIdx(0),
  fMCmumPDG(0),
  fMCgrandmumIdx(0),
  fMCgrandmumPDG(0),
  index(0),
  mcFirstOriginType(kUnassigned),
  doSkipOver(kFALSE),
  fEtaS(0),
  fPhiS(0),
  isP(0),
  isaP(0)
{
  // std::fill(fMomentum,fMomentum+3,0.);
  // std::fill(fMomentumTruth,fMomentumTruth+3,0.);
  // std::fill(fShiftedGlobalPosition,fShiftedGlobalPosition+3,0.);
  // std::fill(iptoPV,iptoPV+2,0.);
  // std::fill(nSigmaFirstTPC,nSigmaFirstTPC+5,0.);
  // std::fill(nSigmaFirstTOF,nSigmaFirstTOF+5,0.);

  // default constructor constructor
}
// //_____________________________________________________________________________
// AliReconstructedFirst::AliReconstructedFirst(const AliReconstructedFirst &obj) :
//   fPt(0),
//   fEta(0),
//   fTheta(0),
//   fPhi(0),
//   fRap(0),
//   fCharge(0),
//   fDCAxy(0),
//   fDCAz(0),
//   isTOFmismatch(kFALSE),
//   isMCptc(kFALSE),
//   fMCcode(0),
//   fPDGcode(0),
//   fMCmumIdx(0),
//   fMCmumPDG(0),
//   fMCgrandmumIdx(0),
//   fMCgrandmumPDG(0),
//   index(0),
//   mcFirstOriginType(kUnassigned),
//   doSkipOver(kFALSE),
//   fEtaS(0),
//   fPhiS(0),
//   isP(0),
//   isaP(0)
// {
//   // copy constructor
// }
// //_____________________________________________________________________________
// AliReconstructedFirst &AliReconstructedFirst::operator=(const AliReconstructedFirst &obj)
// {
//   //Assignment operator
//   if(this == &obj) return *this;
//   fPt            = obj.fPt;
//   fEta           = obj.fEta;
//   fTheta         = obj.fTheta;
//   fPhi           = obj.fPhi;
//   fRap           = obj.fRap;
//   fCharge        = obj.fCharge;
//   fDCAxy         = obj.fDCAxy;
//   fDCAz          = obj.fDCAz;
//   isTOFmismatch  = obj.isTOFmismatch;
//   isMCptc        = obj.isMCptc;
//   fMCcode        = obj.fMCcode;
//   fPDGcode       = obj.fPDGcode;
//   fMCmumIdx      = obj.fMCmumIdx;
//   fMCmumPDG      = obj.fMCmumPDG;
//   fMCgrandmumIdx = obj.fMCgrandmumIdx;
//   fMCgrandmumPDG = obj.fMCgrandmumPDG;
//   index          = obj.index;
//   mcFirstOriginType = obj.mcFirstOriginType;
//   doSkipOver     = obj.doSkipOver;
//   fEtaS          = obj.fEtaS;
//   fPhiS          = obj.fPhiS;
//   isP            = obj.isP;
//   isaP           = obj.isaP;

//   return (*this);

// }

//_____________________________________________________________________________

AliReconstructedFirst::~AliReconstructedFirst()

{

}

//_____________________________________________________________________________

AliReconstructedSecond::AliReconstructedSecond() :
  sPt(0),
  sEta(0),
  sTheta(0),
  sPhi(0),
  sRap(0),
  sCharge(0),
  sDCAxy(0),
  sDCAz(0),
  isTOFmismatch(kFALSE),
  isMCptc(kFALSE),
  sMCcode(0),
  sPDGcode(0),
  sMCmumIdx(0),
  sMCmumPDG(0),
  sMCgrandmumIdx(0),
  sMCgrandmumPDG(0),
  index(0),
  mcSecondOriginType(kUnassigned),
  doSkipOver(kFALSE),
  sEtaS(0),
  sPhiS(0),
  isP(0),
  isaP(0)
{
  // std::fill(sMomentum,sMomentum+3,0.);
  // std::fill(sMomentumTruth,sMomentumTruth+3,0.);
  // std::fill(sShiftedGlobalPosition,sShiftedGlobalPosition+3,0.);
  // std::fill(iptoPV,iptoPV+2,0.);
  // std::fill(nSigmaSecondTPC,nSigmaSecondTPC+5,0.);
  // std::fill(nSigmaSecondTOF,nSigmaSecondTOF+5,0.);
  
  // Default constructor

}
// //_____________________________________________________________________________
// AliReconstructedSecond::AliReconstructedSecond(const AliReconstructedSecond &obj) :
//   sPt(0),
//   sEta(0),
//   sTheta(0),
//   sPhi(0),
//   sRap(0),
//   sCharge(0),
//   sDCAxy(0),
//   sDCAz(0),
//   isTOFmismatch(kFALSE),
//   isMCptc(kFALSE),
//   sMCcode(0),
//   sPDGcode(0),
//   sMCmumIdx(0),
//   sMCmumPDG(0),
//   sMCgrandmumIdx(0),
//   sMCgrandmumPDG(0),
//   index(0),
//   mcSecondOriginType(kUnassigned),
//   doSkipOver(kFALSE),
//   sEtaS(0),
//   sPhiS(0),
//   isP(0),
//   isaP(0)
// {
//   // copy constructor
// }
// //_____________________________________________________________________________
// AliReconstructedSecond &AliReconstructedSecond::operator=(const AliReconstructedSecond &obj)
// {
//   //Assignment operator
//   if(this == &obj) return *this;
//   sPt            = obj.sPt;
//   sEta           = obj.sEta;
//   sTheta         = obj.sTheta;
//   sPhi           = obj.sPhi;
//   sRap           = obj.sRap;
//   sCharge        = obj.sCharge;
//   sDCAxy         = obj.sDCAxy;
//   sDCAz          = obj.sDCAz;
//   isTOFmismatch  = obj.isTOFmismatch;
//   isMCptc        = obj.isMCptc;
//   sMCcode        = obj.sMCcode;
//   sPDGcode       = obj.sPDGcode;
//   sMCmumIdx      = obj.sMCmumIdx;
//   sMCmumPDG      = obj.sMCmumPDG;
//   sMCgrandmumIdx = obj.sMCgrandmumIdx;
//   sMCgrandmumPDG = obj.sMCgrandmumPDG;
//   index          = obj.index;
//   mcSecondOriginType = obj.mcSecondOriginType;
//   doSkipOver     = obj.doSkipOver;
//   sEtaS          = obj.sEtaS;
//   sPhiS          = obj.sPhiS;
//   isP            = obj.isP;
//   isaP           = obj.isaP;

//   return (*this);

// }

//_____________________________________________________________________________

AliReconstructedSecond::~AliReconstructedSecond()

 {

 }

//_____________________________________________________________________________

AliAnalysisKPEvent::AliAnalysisKPEvent():
  fNumberCandidateFirst(0),
  fNumberCandidateSecond(0),
  fReconstructedFirst(0x0),
  fReconstructedSecond(0x0)
{
  //Default constructor
}
//_____________________________________________________________________________
AliAnalysisKPEvent::~AliAnalysisKPEvent()
{
  //Destructor
  
  if(fReconstructedFirst){
    delete fReconstructedFirst;
    fReconstructedFirst= NULL;
  }
  if(fReconstructedSecond){
    delete fReconstructedSecond;
    fReconstructedSecond= NULL;
  }
  
}
//_____________________________________________________________________________
AliAnalysisKPEventCollection::AliAnalysisKPEventCollection() : 
  fEvt(0x0), 
  fifo(0) 
{
  
}

//______________________________________________________________________________

AliAnalysisKPEventCollection::~AliAnalysisKPEventCollection(){

  for(Int_t i = 0; i < fifo; i++){
     if((fEvt + i)->fReconstructedFirst != NULL){
       delete [] (fEvt + i)->fReconstructedFirst;
     }
     if((fEvt + i)->fReconstructedSecond != NULL){
       delete [] (fEvt + i)->fReconstructedSecond;
     }
   }
   delete [] fEvt;fEvt = NULL;
 }


//_____________________________________________________________________________

AliAnalysisKPEventCollection::AliAnalysisKPEventCollection(Short_t eventBuffSize, Int_t maxFirstMult, Int_t maxSecondMult) : fEvt(0x0), fifo(0) {

  SetBuffSize(eventBuffSize);

  fEvt = new AliAnalysisKPEvent[fifo];  //allocate pointer array of AliAnalysisKPEvents

  for(Int_t ii = 0; ii < fifo; ii++){ //Initialize particle table pointers to NULL

    (fEvt + ii)->fReconstructedFirst = NULL;

    (fEvt + ii)->fNumberCandidateFirst = 0;

    (fEvt + ii)->fReconstructedFirst = new AliReconstructedFirst[maxFirstMult];


    (fEvt + ii)->fReconstructedSecond = NULL;

    (fEvt + ii)->fNumberCandidateSecond = 0;

    (fEvt + ii)->fReconstructedSecond = new AliReconstructedSecond[maxSecondMult];


  }

}


//_____________________________________________________________________________
void AliAnalysisKPEventCollection::FifoShift() { //Shift elements in FIFO by one and clear last element in FIFO 

  for(unsigned short i=fifo-1 ; i > 0; i--) {

    for(Int_t j=0; j<(fEvt + i-1)->fNumberCandidateFirst; j++){

      (fEvt + i)->fReconstructedFirst[j] = (fEvt + i-1)->fReconstructedFirst[j];

    }

    (fEvt + i)->fNumberCandidateFirst = (fEvt + i-1)->fNumberCandidateFirst;

    for(Int_t j=0; j<(fEvt + i-1)->fNumberCandidateSecond; j++){

      (fEvt + i)->fReconstructedSecond[j] = (fEvt + i-1)->fReconstructedSecond[j];

    }

    (fEvt + i)->fNumberCandidateSecond = (fEvt + i-1)->fNumberCandidateSecond;

    for(Int_t j=0; j<3; j++){

      (fEvt + i)->fPrimaryVertex[j] = (fEvt + i-1)->fPrimaryVertex[j];

    }

  }

   (fEvt)->fNumberCandidateFirst=0;
   (fEvt)->fNumberCandidateSecond=0;

   for(Int_t j=0; j<3; j++) {

     (fEvt)->fPrimaryVertex[j] = 0.;

   }

}




