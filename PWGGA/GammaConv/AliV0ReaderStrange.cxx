#include <vector>
#include <TGeoGlobalMagField.h>

#include "AliV0ReaderStrange.h"
#include "AliKFParticle.h"
#include "AliV0ParticleStrange.h"
#include "AliAODv0.h"
#include "AliESDv0.h"
#include "AliV0.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliVEvent.h"
#include "AliPID.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliESDpid.h"
#include "AliESDtrackCuts.h"
#include "TRandom3.h"
#include "AliGenCocktailEventHeader.h"
#include "TList.h"
#include "AliKFConversionPhoton.h"
#include "AliAODConversionPhoton.h"
#include "AliConversionPhotonBase.h"
#include "TVector.h"
#include "AliKFVertex.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODHandler.h"
#include "AliPIDResponse.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TObjArray.h"
#include "AliVTrack.h"
#include "AliKFParticle.h"

class iostream;

using namespace std;

ClassImp(AliV0ReaderStrange)

//________________________________________________________________________
AliV0ReaderStrange::AliV0ReaderStrange(const char *name) : AliAnalysisTaskSE(name),
  fEventCuts(NULL),
  fV0Cuts(NULL),
  fVectorFoundGammas(0),
  fConversionGammas(NULL),
  fEventIsSelected(kFALSE)
{
  // Default constructor

  DefineInput(0, TChain::Class());
}

//________________________________________________________________________
AliV0ReaderStrange::~AliV0ReaderStrange()
{
  // default deconstructor

  if(fConversionGammas){
    fConversionGammas->Delete();// Clear Objects
    delete fConversionGammas;
    fConversionGammas=0x0;
  }
}


//________________________________________________________________________
void AliV0ReaderStrange::Init()
{
  // Initialize function to be called once before analysis
  if(fV0Cuts==NULL){
    if(fV0Cuts==NULL)AliError("No V0 Cut Selection initialized");
  }
  if(fEventCuts==NULL){
    if(fEventCuts==NULL)AliError("No Event Cut Selection initialized");
  }

  if(fConversionGammas != NULL){
    delete fConversionGammas;
    fConversionGammas=NULL;
  }

  if(fConversionGammas == NULL){
    fConversionGammas = new TClonesArray("AliKFParticle",1000);
//     if(kUseAODConversionPhoton){
//       fConversionGammas = new TClonesArray("AliAODConversionPhoton",100);}
//     else{
//       fConversionGammas = new TClonesArray("AliKFConversionPhoton",100);}
  }
  fConversionGammas->Delete();//Reset the TClonesArray
}

//________________________________________________________________________
void AliV0ReaderStrange::UserCreateOutputObjects()
{  
    fVectorFoundGammas.clear();
}  
  
//________________________________________________________________________
void AliV0ReaderStrange::UserExec(Option_t *option){

  AliESDEvent * esdEvent = dynamic_cast<AliESDEvent*>(fInputEvent);
  if(esdEvent) {
    if (!TGeoGlobalMagField::Instance()->GetField()) esdEvent->InitMagneticField();
  }

  // Check if correctly initialized
  if(!fConversionGammas)Init();

  // User Exec
  fEventIsSelected=ProcessEvent(fInputEvent,fMCEvent);
}

//________________________________________________________________________
Bool_t AliV0ReaderStrange::ProcessEvent(AliVEvent *inputEvent,AliMCEvent *mcEvent)
{
  //Reset the TClonesArray
  fConversionGammas->Delete();

  fInputEvent=inputEvent;
  fMCEvent=mcEvent;

  if(!fInputEvent){
    AliError("No Input event");
    return kFALSE;
  }
  if(!fEventCuts){AliError("No EventCuts");return kFALSE;}
  if(!fV0Cuts){AliError("No V0 Cuts");return kFALSE;}

  // Event Cuts
  if(!fEventCuts->EventIsSelected(fInputEvent,fMCEvent))return kFALSE;

  // Set Magnetic Field
  AliKFParticle::SetField(fInputEvent->GetMagneticField());

  if(fInputEvent->IsA()==AliESDEvent::Class()){
    ProcessESDV0s();
  }
  if(fInputEvent->IsA()==AliAODEvent::Class()){
    ProcessAODV0s();
  }
  
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliV0ReaderStrange::ProcessESDV0s()
{
  // Process ESD V0s for V0 reconstruction  
  
  AliV0ParticleStrange *fCurrentMotherLambdaCandidate=NULL;
  
  AliESDEvent *fESDEvent=dynamic_cast<AliESDEvent*>(fInputEvent);
    
  for(Int_t iV0=0; iV0<fESDEvent->GetNumberOfV0s(); ++iV0){
    AliESDv0 *fCurrentV0=(AliESDv0*)(fESDEvent->GetV0(iV0));
    if(!fCurrentV0){
      printf("Requested V0 does not exist");
      continue;
    }
    
    fCurrentMotherLambdaCandidate = ReconstructV0(fESDEvent, fCurrentV0, iV0);
    if(fCurrentMotherLambdaCandidate){
      new((*fConversionGammas)[fConversionGammas->GetEntriesFast()]) AliV0ParticleStrange(*fCurrentMotherLambdaCandidate);
      delete fCurrentMotherLambdaCandidate;
      fCurrentMotherLambdaCandidate=NULL;
    }
  } 
  
  return kTRUE;
}

///________________________________________________________________________
AliV0ParticleStrange *AliV0ReaderStrange::ReconstructV0(AliESDEvent *fESDEvent, AliESDv0 *fCurrentV0,Int_t currentV0Index)
{
  fV0Cuts->FillV0CutIndex(AliV0CutsStrange::kV0In);
  
  //checks if on the fly mode is set
  if(!fV0Cuts->SelectV0Finder(fCurrentV0->GetOnFlyStatus())){
    fV0Cuts->FillV0CutIndex(AliV0CutsStrange::kOnFly);
    return 0x0;
  }

//       Double_t lV0CosineOfPointingAngle = fCurrentV0->GetV0CosineOfPointingAngle();
//       GetV0CosineOfPointingAngle(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2]);
 
    if(fCurrentV0->GetV0CosineOfPointingAngle() < 0.98) return 0x0;
  
    AliVTrack * pos = (AliVTrack*)fESDEvent->GetTrack(fCurrentV0->GetPindex());
    AliVTrack * neg = (AliVTrack*)fESDEvent->GetTrack(fCurrentV0->GetNindex());
    const AliExternalTrackParam * paramPos = fCurrentV0->GetParamP();
    const AliExternalTrackParam * paramNeg = fCurrentV0->GetParamN();
    
    if(pos->GetSign() <0){//change tracks
      pos=neg;
      neg=fESDEvent->GetTrack(fCurrentV0->GetPindex());
      paramPos=paramNeg;
      paramNeg=fCurrentV0->GetParamP();
    }
        
  if(!pos || !neg ) {
    fV0Cuts->FillV0CutIndex(AliV0CutsStrange::kNoTracks);
    return 0x0;
//     return;
  }
    //cuts ------------    //remove like sign pairs
    if(pos->GetSign() == neg->GetSign()){
      fV0Cuts->FillV0CutIndex(AliV0CutsStrange::kSameSign);
      return 0x0;
    }
    
//         if(pos->Pt() < .16) return 0x0;
//     if(neg->Pt() < .16) return 0x0;
//     if(TMath::Abs(pos->Eta()) > 0.8) return 0x0;
//     if(TMath::Abs(neg->Eta()) > 0.8) return 0x0;
    
//     if( neg->GetKinkIndex(0) > 0 || pos->GetKinkIndex(0) > 0) return 0x0 ;        
//     if( !(pos->GetStatus() & AliESDtrack::kTPCrefit) || !(neg->GetStatus() & AliESDtrack::kTPCrefit) ) return 0x0;
  
  Bool_t isPiPlus = fV0Cuts->GetPIDpion(pos);
  Bool_t isPiMinus = fV0Cuts->GetPIDpion(neg);
  Bool_t isProton = fV0Cuts->GetPIDproton(pos);
  Bool_t isAntiProton = fV0Cuts->GetPIDproton(neg);

  
  if(!((isProton && isPiMinus) || (isAntiProton && isPiPlus))) return 0x0;
  
  // pi plus: 211, pi minus: -211, proton: 2212, anti-proton: -2212
    
  Int_t pdgCodePos = 2211; //proton
  Int_t pdgCodeNeg = -211; // pi minus
  if(isPiPlus){
    pdgCodePos = -211;
  }
  if(isAntiProton){
    pdgCodeNeg = -2212;
  }
  
//   AliVVertex *vertex = fESDEvent->GetPrimaryVertex();
//   fCurrentV0->DecayLength(vertex);
    
  const AliKFParticle fCurrentPositiveKFParticle(*(paramPos), pdgCodePos);
  const AliKFParticle fCurrentNegativeKFParticle(*(paramNeg), pdgCodeNeg);
  const Bool_t gamma = kFALSE;
  AliKFParticle *fCurrentMother = new AliKFParticle(fCurrentPositiveKFParticle, fCurrentNegativeKFParticle, gamma);
   
  AliV0ParticleStrange *fCurrentMotherV0 = new AliV0ParticleStrange(fCurrentMother);
  
   if(fMCEvent){

    Int_t labelp=TMath::Abs(fV0Cuts->GetTrack(fESDEvent,fCurrentMotherV0->GetTrackLabelPositive())->GetLabel());
    Int_t labeln=TMath::Abs(fV0Cuts->GetTrack(fESDEvent,fCurrentMotherV0->GetTrackLabelNegative())->GetLabel());

//     cout << "rec: " <<  currentTrackLabels[0] << "\t" << currentTrackLabels[1] << endl;
//     cout << "recProp: " <<  fCurrentMotherKF->GetTrackLabelPositive() << "\t" << fCurrentMotherKF->GetTrackLabelNegative() << endl;
//     cout << "MC: " <<  labeln << "\t" << labelp << endl;

    TParticle *fNegativeMCParticle = 0x0;
    if(labeln>-1) fNegativeMCParticle = fMCEvent->Particle(labeln);
    TParticle *fPositiveMCParticle = 0x0;
    if(labelp>-1) fPositiveMCParticle = fMCEvent->Particle(labelp);

    if(fPositiveMCParticle&&fNegativeMCParticle){
      fCurrentMotherV0->SetMCLabelPositive(labelp);
      fCurrentMotherV0->SetMCLabelNegative(labeln);
    }
  }
  
  fV0Cuts->FillV0CutIndex(AliV0CutsStrange::kV0Out);
  return fCurrentMotherV0;
}


///____________________________________________________________________________________________________________
Bool_t AliV0ReaderStrange::ProcessAODV0s()
{
  // Process ESD V0s for V0 reconstruction  
  
  AliV0ParticleStrange *fCurrentMotherLambdaCandidate=NULL;
  
  AliAODEvent *fAODEvent=dynamic_cast<AliAODEvent*>(fInputEvent);
    
  for(Int_t iV0=0; iV0<fAODEvent->GetNumberOfV0s(); ++iV0){
    AliAODv0 *fCurrentV0=(AliAODv0*)(fAODEvent->GetV0(iV0));
    if(!fCurrentV0){
      printf("Requested V0 does not exist");
      continue;
    }
    
    fCurrentMotherLambdaCandidate = ReconstructV0(fAODEvent, fCurrentV0, iV0);
    if(fCurrentMotherLambdaCandidate){
      new((*fConversionGammas)[fConversionGammas->GetEntriesFast()]) AliV0ParticleStrange(*fCurrentMotherLambdaCandidate);
      delete fCurrentMotherLambdaCandidate;
      fCurrentMotherLambdaCandidate=NULL;
    }
  } 
  
  return kTRUE;
}

///________________________________________________________________________
AliV0ParticleStrange *AliV0ReaderStrange::ReconstructV0(AliAODEvent *fAODEvent, AliAODv0 *fCurrentV0,Int_t currentV0Index)
{
  fV0Cuts->FillV0CutIndex(AliV0CutsStrange::kV0In);
  
  //checks if on the fly mode is set
  if(!fV0Cuts->SelectV0Finder(fCurrentV0->GetOnFlyStatus())){
    fV0Cuts->FillV0CutIndex(AliV0CutsStrange::kOnFly);
    return 0x0;
  }
  
//       Double_t lV0CosineOfPointingAngle = fCurrentV0->GetV0CosineOfPointingAngle();
//       GetV0CosineOfPointingAngle(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2]);
 
    AliAODVertex *vertex = fAODEvent->GetPrimaryVertex();
    if(fCurrentV0->CosPointingAngle(vertex) < 0.98) return 0x0;
    
    AliVTrack* pos = (AliVTrack*)fCurrentV0->GetDaughter(0);
    AliVTrack* neg = (AliVTrack*)fCurrentV0->GetDaughter(1);
    
    AliExternalTrackParam paramPos;
    paramPos.CopyFromVTrack(pos);
    AliExternalTrackParam paramNeg;
    paramPos.CopyFromVTrack(neg);
    
    if(pos->GetSign() <0){//change tracks
      pos=neg;
      neg=(AliVTrack*)fCurrentV0->GetDaughter(0);
      paramPos=paramNeg;
      paramNeg.CopyFromVTrack(pos);
    }
        
  if(!pos || !neg ) {
    fV0Cuts->FillV0CutIndex(AliV0CutsStrange::kNoTracks);
    return 0x0;
  }

  //cuts ------------    //remove like sign pairs
    if(pos->GetSign() == neg->GetSign()){
      fV0Cuts->FillV0CutIndex(AliV0CutsStrange::kSameSign);
      return 0x0;
    }
    
//         if(pos->Pt() < .16) return 0x0;
//     if(neg->Pt() < .16) return 0x0;
//     if(TMath::Abs(pos->Eta()) > 0.8) return 0x0;
//     if(TMath::Abs(neg->Eta()) > 0.8) return 0x0;
    
//     if( neg->GetKinkIndex(0) > 0 || pos->GetKinkIndex(0) > 0) return 0x0 ;        
//     if( !(pos->GetStatus() & AliESDtrack::kTPCrefit) || !(neg->GetStatus() & AliESDtrack::kTPCrefit) ) return 0x0;
  
  Bool_t isPiPlus = fV0Cuts->GetPIDpion(pos);
  Bool_t isPiMinus = fV0Cuts->GetPIDpion(neg);
  Bool_t isProton = fV0Cuts->GetPIDproton(pos);
  Bool_t isAntiProton = fV0Cuts->GetPIDproton(neg);

  
  if(!((isProton && isPiMinus) || (isAntiProton && isPiPlus))) return 0x0;
  
  // pi plus: 211, pi minus: -211, proton: 2212, anti-proton: -2212
    
  Int_t pdgCodePos = 2211; //proton
  Int_t pdgCodeNeg = -211; // pi minus
  if(isPiPlus){
    pdgCodePos = -211;
  }
  if(isAntiProton){
    pdgCodeNeg = -2212;
  }
  
  const AliKFParticle fCurrentPositiveKFParticle(paramPos, pdgCodePos);
  const AliKFParticle fCurrentNegativeKFParticle(paramNeg, pdgCodeNeg);
  const Bool_t gamma = kFALSE;
  AliKFParticle *fCurrentMother = new AliKFParticle(fCurrentPositiveKFParticle, fCurrentNegativeKFParticle, gamma);
   
  AliV0ParticleStrange *fCurrentMotherV0 = new AliV0ParticleStrange(fCurrentMother);
  
   if(fMCEvent){

    Int_t labelp=TMath::Abs(fV0Cuts->GetTrack(fAODEvent,fCurrentMotherV0->GetTrackLabelPositive())->GetLabel());
    Int_t labeln=TMath::Abs(fV0Cuts->GetTrack(fAODEvent,fCurrentMotherV0->GetTrackLabelNegative())->GetLabel());

//     cout << "rec: " <<  currentTrackLabels[0] << "\t" << currentTrackLabels[1] << endl;
//     cout << "recProp: " <<  fCurrentMotherKF->GetTrackLabelPositive() << "\t" << fCurrentMotherKF->GetTrackLabelNegative() << endl;
//     cout << "MC: " <<  labeln << "\t" << labelp << endl;

    TParticle *fNegativeMCParticle = 0x0;
    if(labeln>-1) fNegativeMCParticle = fMCEvent->Particle(labeln);
    TParticle *fPositiveMCParticle = 0x0;
    if(labelp>-1) fPositiveMCParticle = fMCEvent->Particle(labelp);

    if(fPositiveMCParticle&&fNegativeMCParticle){
      fCurrentMotherV0->SetMCLabelPositive(labelp);
      fCurrentMotherV0->SetMCLabelNegative(labeln);
    }
  }
  
  
  fV0Cuts->FillV0CutIndex(AliV0CutsStrange::kV0Out);
  return fCurrentMotherV0;
}
