/*
***********************************************************
    Variable definitions for event plane correction framework
    Contact: Jaap Onderwaater, j.onderwaater@gsi.de, jacobus.onderwaater@cern.ch
    Based on work of Ionut-Cristian Arsene
***********************************************************
*/

#include <iostream>
#include <fstream>


#include "AliReducedQnFillEvent.h"
//#include "AliQnCorrectionsVarManager.h"
#include "AliReducedVarManager.h"
#include "AliReducedEventInfo.h"
#include "AliReducedTrackInfo.h"
#include "AliReducedFMDInfo.h"
#include "AliHistogramManager.h"

#include "AliQnCorrectionsDataVector.h"
#include <AliQnCorrectionsConfiguration.h>
#include <AliQnCorrectionsManager.h>


#include <TChain.h>
#include <TH1D.h>
#include <TFile.h>


using std::cout;
using std::endl;
ClassImp(AliReducedQnFillEvent)

#ifdef ALIREDUCEDVARMANAGER_H
#define VAR AliReducedVarManager
#endif
//#ifdef ALIQNCORRECTIONS_VARMANAGER_H
//#define VAR AliQnCorrectionsVarManager
//#endif


AliReducedQnFillEvent::AliReducedQnFillEvent() :
  TNamed("AliReducedQnFillEvent","Fill functions"),
  fEvent(0x0),
  fEventPlaneManager(0x0),
  fEventPlaneHistos(0x0),
  fFillVZERO(kFALSE),
  fFillTPC(kFALSE),    
  fFillZDC(kFALSE),
  fFillTZERO(kFALSE),
  fFillFMD(kFALSE)
{
  //
  // Default constructor
  //

}

//_____________________________________________________________________________
AliReducedQnFillEvent::~AliReducedQnFillEvent()
{
  //
  // Destructor
  //
}



//__________________________________________________________________
void AliReducedQnFillEvent::SetDetectors() {
  //
  // determine which detectors are used (to call only the necessary detector fill functions)

  if(fEventPlaneManager->GetDetectorId(VAR::kTPC  )!=-1)    fFillTPC = kTRUE;
  if(fEventPlaneManager->GetDetectorId(VAR::kVZERO)!=-1)  fFillVZERO = kTRUE;
  if(fEventPlaneManager->GetDetectorId(VAR::kTZERO)!=-1)  fFillTZERO = kTRUE;
  if(fEventPlaneManager->GetDetectorId(VAR::kZDC  )!=-1)    fFillZDC = kTRUE;
  if(fEventPlaneManager->GetDetectorId(VAR::kFMD  )!=-1)    fFillFMD = kTRUE;


}


//__________________________________________________________________
void AliReducedQnFillEvent::Process(AliReducedEventInfo* event, Float_t* values) {

  //FillEventInfo(values);
  fEvent=event;
  FillDetectors(values);

}

//__________________________________________________________________
//void AliReducedQnFillEvent::FillEventInfo(Float_t* values) {
//  //
//  // fill event info
//  //
//
//
//  values[VAR::kRunNo]       = fEvent->GetRunNumber();
//  //values[VAR::kTriggerMask] = event->TriggerMask();
//  //values[VAR::kIsPhysicsSelection]  = (event->IsPhysicsSelection() ? 1.0 : 0.0);
//  //values[VAR::kNVtxTPCContributors] = event->VertexTPCContributors();
//  values[VAR::kVtxX]        = -999.;
//  values[VAR::kVtxY]        = -999.;
//  values[VAR::kVtxZ]        = -999.;
//  const AliVVertex *primVtx = fEvent->GetPrimaryVertex();
//  if (primVtx){
//    values[VAR::kVtxX]        = primVtx->GetX();
//    values[VAR::kVtxY]        = primVtx->GetY();
//    values[VAR::kVtxZ]        = primVtx->GetZ();
//    values[VAR::kNVtxContributors]    = primVtx->GetNContributors();
//  }
//
//  AliESDEvent* esdEvent = static_cast<AliESDEvent*>(fEvent);
//  AliCentrality* cent = esdEvent->GetCentrality();
//  if(cent){
//    values[VAR::kCentVZERO]   = cent->GetCentralityPercentile("V0M");
//    values[VAR::kCentSPD]     = cent->GetCentralityPercentile("CL1");
//    values[VAR::kCentTPC]     = cent->GetCentralityPercentile("TRK");
//    values[VAR::kCentQuality] = cent->GetQuality();
//  }
//    
//
//  AliVVZERO* vzero = fEvent->GetVZEROData();
//  values[VAR::kVZEROATotalMult]     = vzero->GetMTotV0A();
//  values[VAR::kVZEROCTotalMult]     = vzero->GetMTotV0C();
//  values[VAR::kVZEROTotalMult]      = values[VAR::kVZEROATotalMult]+values[VAR::kVZEROCTotalMult];
//
//  //fEventPlaneManager->EventCuts()->IsSelected(values);
//}



//__________________________________________________________________
//void AliReducedQnFillEvent::FillTrackInfo(AliReducedTrackInfo* particle, Float_t* values) {
//
//  Float_t dcaxy=0.0;
//  Float_t dcaz=0.0;
//  //particle->GetImpactParameters(dcaxy,dcaz);
//
//  values[VAR::kPx]        = particle->Px();
//  values[VAR::kPy]        = particle->Py();
//  values[VAR::kPz]        = particle->Pz();
//  values[VAR::kPt]        = particle->Pt();
//  values[VAR::kP]         = particle->P();
//  values[VAR::kPhi]       = particle->Phi();
//  values[VAR::kTheta]     = particle->Theta();
//  values[VAR::kEta]       = particle->Eta();
//  values[VAR::kCharge]    = particle->Charge();
//  values[VAR::kDcaXY]     = dcaxy;
//  values[VAR::kDcaZ]      = dcaz;
//
//  AliAODTrack* aodTrack=static_cast<AliAODTrack*>(particle);
//
//  //values[VAR::kITSncls]       = particle->GetNcls(0); 
//  values[VAR::kTPCncls]       = aodTrack->GetTPCNcls();
//  values[VAR::kTPCchi2]       = aodTrack->Chi2perNDF();
//  values[VAR::kTPCsignal]     = aodTrack->GetTPCsignal();
//  for(Int_t ibit=0; ibit<9; ibit++) values[VAR::kFilterBit+ibit]     = aodTrack->TestFilterBit(BIT(ibit));
//
//  //delete aodTrack;
//
//}


//_________________________________
void AliReducedQnFillEvent::FillDetectors(Float_t* values){

  if(fFillTPC)   FillTPC(values);
  if(fFillVZERO) FillVZERO();
  if(fFillZDC)   FillZDC();
  if(fFillTZERO) FillTZERO();
  if(fFillFMD)   FillFMD();

}


//_________________________________
void AliReducedQnFillEvent::FillTPC(Float_t* values){
  //
  // fill TPC info
  //
  Int_t nTrack=-1;
  AliQnCorrectionsConfiguration* QnConf =0x0;
  TClonesArray* QnConfList = fEventPlaneManager->GetQnConfigurations(VAR::kTPC);

  AliReducedTrackInfo* track = 0x0;
  TClonesArray* trackList = fEvent->GetTracks();
  TIter nextTrack(trackList);
  for(Int_t it=0; it<fEvent->NTracks(); ++it) {
    track = (AliReducedTrackInfo*)nextTrack();

    VAR::FillTrackInfo(track, values);

    fEventPlaneManager->AddDataVector(VAR::kTPC, track->Phi());


    //for(Int_t iconf=0; iconf<QnConfList->GetEntriesFast(); iconf++){
    //  QnConf = (AliQnCorrectionsConfiguration*) QnConfList->At(iconf);
    //  if(QnConf->PassCuts(values)){
    //      fEventPlaneHistos->FillHistClass("TrackQA_"+QnConf->QnConfigurationName(), values);
    //  }
    //}


  }



}


//_________________________________
void AliReducedQnFillEvent::FillVZERO(){
  //
  // fill VZERO info
  //

  Double_t weight=0.;
  const Double_t kX[8] = {0.92388, 0.38268, -0.38268, -0.92388, -0.92388, -0.38268, 0.38268, 0.92388};    // cosines of the angles of the VZERO sectors (8 per ring)
  const Double_t kY[8] = {0.38268, 0.92388, 0.92388, 0.38268, -0.38268, -0.92388, -0.92388, -0.38268};    // sines     -- " --
  AliQnCorrectionsConfiguration* QnConf = 0x0;

  for(Int_t ich=0; ich<64; ich++){
    weight=fEvent->MultChannelVZERO(ich);
    if(weight<0.01) weight=0.;

    fEventPlaneManager->AddDataVector(VAR::kVZERO, TMath::ATan2(kY[ich%8],kX[ich%8]), weight, ich);   // 1st ich is position in array, 2nd ich is channel id

  }
}



//_________________________________
void AliReducedQnFillEvent::FillTZERO(){
  //
  // fill ESD TZERO info
  //

  Double_t weight=0.0;
  const Double_t kX[24] = {/* Cside */ 0.905348,0.571718,0.0848977,-0.424671,-0.82045,-0.99639,-0.905348,-0.571718,-0.0848977,0.424671,0.82045,0.99639, /* Aside */ 0.99995,0.870982,0.508635,0.00999978,-0.491315,-0.860982,-0.99995,-0.870982,-0.508635,-0.0100001,0.491315,0.860982};
  const Double_t kY[24] = {/* Cside */ 0.424671,0.82045,0.99639,0.905348,0.571718,0.0848976,-0.424671,-0.82045,-0.99639,-0.905348,-0.571719,-0.0848975, /* Aside */ -0.00999983,0.491315,0.860982,0.99995,0.870982,0.508635,0.00999974,-0.491315,-0.860982,-0.99995,-0.870982,-0.508635};

  AliQnCorrectionsConfiguration* QnConf = 0x0;

  for(Int_t ich=0; ich<24; ich++){
    weight=fEvent->AmplitudeTZEROch(ich);
    if(weight<0.01) weight=0.;

    fEventPlaneManager->AddDataVector(VAR::kTZERO, TMath::ATan2(kY[ich%8],kX[ich%8]), weight, ich);   // 1st ich is position in array, 2nd ich is channel id


  }
}




//_________________________________
void AliReducedQnFillEvent::FillZDC(){
  //
  // fill ZDC info
  //

  const Double_t kX[10] = { /* Cside */ 0.0,  -1.75,  1.75, -1.75, 1.75, /* Aside */  0.0,  1.75, -1.75, 1.75, -1.75  };
  const Double_t kY[10] = { /* Cside */ 0.0,  -1.75, -1.75,  1.75, 1.75, /* Aside */  0.0, -1.75, -1.75, 1.75,  1.75  };


  AliQnCorrectionsConfiguration* QnConf = 0x0;

  Double_t weight=0.0;
  Float_t ZDCenergy[10];
  for(Int_t i=0; i<10; ++i)    ZDCenergy[i]  = fEvent->EnergyZDCnTree(i);

  for(Int_t ich=0; ich<10; ich++){
    if(weight<0.01) weight=0.;
    fEventPlaneManager->AddDataVector(VAR::kZDC, TMath::ATan2(kY[ich%8],kX[ich%8]), weight, ich);   // 1st ich is position in array, 2nd ich is channel id
  }

}




//_________________________________
void AliReducedQnFillEvent::FillFMD(){
  //
  // fill FMD info
  //
  Int_t nTrack=-1;

  AliReducedFMDInfo* fmd = 0x0;
  TClonesArray* fmdList = fEvent->GetFMD();
  TIter nextTrack(fmdList);
  for(Int_t it=0; it<fmdList->GetEntriesFast(); ++it) {
    fmd = (AliReducedFMDInfo*)nextTrack();
    if(!fmd) continue;
    //cout<<it<<"  "<<fmd->Multiplicity()<<"  "<<fmdList->GetEntriesFast()<<endl;
    //else cout<<fmd<<"  "<<it<<"  "<<fmd->Id()<<endl;
    fEventPlaneManager->AddDataVector(VAR::kFMD, fmd->Phi(), fmd->Multiplicity(), TMath::Abs(fmd->Id()));
  }
}
