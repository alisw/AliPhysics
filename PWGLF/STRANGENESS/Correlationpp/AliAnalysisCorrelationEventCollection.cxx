////////////////////////////////////////////////////////////////////////////////
//
//  This class provides storage for event and track information which 
//  are used for same-event as well as mixed-event analyses in AliAnalysisTaskCorrelationhhK0s
//  Author: chiara.de.martin@cern.ch
//  was: ramona.lea@cern.ch 
//  was: maria.nicassio@cern.ch (derived and adapted from D. Gangadharan PWGCF/FEMTOSCOPY/Chaoticity/AliChaoticityEventCollection
//                                  and J. Salzwedel PWGCF/FEMTOSCOPY/V0LamAnalysis/AliAnalysisV0LamEventCollection) 
//
////////////////////////////////////////////////////////////////////////////////

#include "AliAnalysisCorrelationEventCollection.h"

//_____________________________________________________________________________
// Default constructor 
AliReconstructedFirstC::AliReconstructedFirstC() :
  fPt(0),
  fEta(0),
  fTheta(0),
  fPhi(0),
  fRap(0),
  fCharge(0),
  fDCAxy(0),
  fDCAz(0),
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
  fMultiplicity(0),
  fZvertex(0)
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

AliReconstructedFirstC::~AliReconstructedFirstC()

{

}

//_____________________________________________________________________________

AliReconstructedSecondC::AliReconstructedSecondC() :
  sPt(0),
  sEta(0),
  sTheta(0),
  sPhi(0),
  sRap(0),
  sCharge(0),
  sDCAxy(0),
  sDCAz(0),
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
  sDcaPosV0(0),
  sDcaNegV0(0),
  sPtArmV0(0),
  sAlphaV0(0),
  sInvMassK0s(0),
  sInvMassLambda(0),
  sInvMassAntiLambda(0),
  sCosPointingAngle(0),
  sDcaV0ToPV(0),
  sMultiplicity(0),
  sZvertex(0),
  sctau(0), 
  sLabelMotherPos(0),
  sLabelPos(0),
  sLabelNeg(0),

//variables used for cascades as associated particles
  cLabelMotherBach(0),
  cisPrimCasc(0),
  cInvMassLambda(0),
  cInvMassXi(0),
  cInvMassOmega(0),
  cCosPointingAngleXi(0),
  cCosPointingAngleV0ToXi(0),
  cDCAXiDaughters(0),
  cRapCasc(0),
  cPt(0),
  cctau(0),
  cEta(0),
  cTheta(0),
  cPhi(0),
  cCharge(0),
  cAssocOrNot(0)

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

AliReconstructedSecondC::~AliReconstructedSecondC()

 {

 }

//_____________________________________________________________________________

/*AliReconstructedCascC::AliReconstructedCascC() :
  
  cLabelMotherBach(0),
  cisPrimXi(0),
  cInvMassLambda(0),
  cInvMassXi(0),
  cInvMassOmega(0),
  cCosPointingAngleXi(0),
  cCosPointingAngleV0ToX(0),
  cDCAXiDaughters(0),
  cRapCasc(0),
  cPt(0),
  cctau(0),
  cEta(0),
  cTheta(0),
  cPhi(0),
  cCharge(0),
{
  // std::fill(sMomentum,sMomentum+3,0.);
  // std::fill(sMomentumTruth,sMomentumTruth+3,0.);
  // std::fill(sShiftedGlobalPosition,sShiftedGlobalPosition+3,0.);
  // std::fill(iptoPV,iptoPV+2,0.);
  // std::fill(nSigmaSecondTPC,nSigmaSecondTPC+5,0.);
  // std::fill(nSigmaSecondTOF,nSigmaSecondTOF+5,0.);
  
  // Default constructor

}

//_____________________________________________________________________________

AliReconstructedCascC::~AliReconstructedCascC()

 {

 }

*/
//_____________________________________________________________________________

AliAnalysisCorrelationEvent::AliAnalysisCorrelationEvent():
  fNumberCandidateFirst(0),
  fNumberCandidateSecond(0),
  fReconstructedFirst(0x0),
  fReconstructedSecond(0x0)
{
  //Default constructor
}
//_____________________________________________________________________________
AliAnalysisCorrelationEvent::~AliAnalysisCorrelationEvent()
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
  /*
  if(fReconstructedCasc){
    delete fReconstructedCasc;
    fReconstructedCasc= NULL;
  }
  */
  
}
//_____________________________________________________________________________
AliAnalysisCorrelationEventCollection::AliAnalysisCorrelationEventCollection() : 
  fEvt(0x0), 
  fifo(0) 
{
  
}

//______________________________________________________________________________

AliAnalysisCorrelationEventCollection::~AliAnalysisCorrelationEventCollection(){

  for(Int_t i = 0; i < fifo; i++){
     if((fEvt + i)->fReconstructedFirst != NULL){
       delete [] (fEvt + i)->fReconstructedFirst;
     }
     if((fEvt + i)->fReconstructedSecond != NULL){
       delete [] (fEvt + i)->fReconstructedSecond;
     }
     /*
     if((fEvt + i)->fReconstructedCasc != NULL){
       delete [] (fEvt + i)->fReconstructedCasc;
     }
     */
   }
   delete [] fEvt;fEvt = NULL;
 }


//_____________________________________________________________________________

AliAnalysisCorrelationEventCollection::AliAnalysisCorrelationEventCollection(Short_t eventBuffSize, Int_t maxFirstMult, Int_t maxSecondMult) : fEvt(0x0), fifo(0) {

  SetBuffSize(eventBuffSize);

  fEvt = new AliAnalysisCorrelationEvent[fifo];  //allocate pointer array of AliAnalysisKPEvents

  for(Int_t ii = 0; ii < fifo; ii++){ //Initialize particle table pointers to NULL

    (fEvt + ii)->fReconstructedFirst = NULL;

    (fEvt + ii)->fNumberCandidateFirst = 0;

    (fEvt + ii)->fReconstructedFirst = new AliReconstructedFirstC[maxFirstMult];


    (fEvt + ii)->fReconstructedSecond = NULL;

    (fEvt + ii)->fNumberCandidateSecond = 0;

    (fEvt + ii)->fReconstructedSecond = new AliReconstructedSecondC[maxSecondMult];


  }

}


//_____________________________________________________________________________
void AliAnalysisCorrelationEventCollection::FifoShift() { //Shift elements in FIFO by one and clear last element in FIFO 

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

void AliAnalysisCorrelationEventCollection::FifoClear() { //Shift elements in FIFO by one and clear last element in FIFO 

   (fEvt)->fNumberCandidateFirst=0;
   (fEvt)->fNumberCandidateSecond=0;

   for(Int_t j=0; j<3; j++) {

     (fEvt)->fPrimaryVertex[j] = 0.;

   }
}





