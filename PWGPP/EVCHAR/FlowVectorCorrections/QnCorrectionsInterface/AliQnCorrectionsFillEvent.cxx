/*
 ***********************************************************
 Variable definitions for event plane correction framework
Contact: Jaap Onderwaater, j.onderwaater@gsi.de, jacobus.onderwaater@cern.ch
Based on work of Ionut-Cristian Arsene
 ***********************************************************
 */

#include <iostream>
#include <fstream>


#include "AliQnCorrectionsFillEvent.h"
#include "AliQnCorrectionsVarManager.h"

#include <AliQnCorrectionsDataVector.h>
#include <AliQnCorrectionsConfiguration.h>
#include <AliQnCorrectionsManager.h>
#include "AliQnCorrectionsHistos.h"

#include <TChain.h>
#include <TH1D.h>
#include <TFile.h>

#include <AliCFContainer.h>
#include <AliInputEventHandler.h>
#include <AliAnalysisManager.h>
#include <AliCentrality.h>
#include <AliMultSelection.h>
#include <AliAnalysisManager.h>
#include <AliAnalysisTaskSE.h>
#include <AliAODForwardMult.h>
#include <AliForwardUtil.h>

#include <AliVEvent.h>
#include <AliVZDC.h>
#include <AliVVZERO.h>
#include <AliVParticle.h>

#include <AliESDInputHandler.h>
#include <AliESDEvent.h>
#include <AliESDHeader.h>
#include <AliESDtrack.h>
#include <AliESDtrackCuts.h>
#include <AliESDFMD.h>

#include <AliAODInputHandler.h>
#include <AliAODEvent.h>
#include <AliAODHeader.h>
#include <AliAODTrack.h>


using std::cout;
using std::endl;
ClassImp(AliQnCorrectionsFillEvent)

#ifdef ALIQNCORRECTIONS_VARMANAGER_H
#define VAR AliQnCorrectionsVarManager
#endif


  AliQnCorrectionsFillEvent::AliQnCorrectionsFillEvent() :
    TNamed("AliQnCorrectionsFillEvent","Fill functions"),
    fEvent(0x0),
    fEventPlaneManager(0x0),
    fQAhistos(0x0),
    fUseTPCStandaloneTracks(kFALSE),
    fFillVZERO(kFALSE),
    fFillTPC(kFALSE),    
    fFillZDC(kFALSE),
    fFillTZERO(kFALSE),
    fFillFMD(kFALSE),
    fFillRawFMD(kFALSE),
    fFillSPD(kFALSE),
    fIsAOD(kFALSE),
    fIsESD(kFALSE)
{
  //
  // Default constructor
  //

}

//_____________________________________________________________________________
AliQnCorrectionsFillEvent::~AliQnCorrectionsFillEvent()
{
  //
  // Destructor
  //
}



//__________________________________________________________________
void AliQnCorrectionsFillEvent::SetDetectors() {
  //
  // determine which detectors are used (to call only the necessary detector fill functions)

  if(fEventPlaneManager->GetDetectorId(VAR::kTPC  )!=-1)    fFillTPC = kTRUE;
  if(fEventPlaneManager->GetDetectorId(VAR::kVZERO)!=-1)  fFillVZERO = kTRUE;
  if(fEventPlaneManager->GetDetectorId(VAR::kTZERO)!=-1)  fFillTZERO = kTRUE;
  if(fEventPlaneManager->GetDetectorId(VAR::kZDC  )!=-1)    fFillZDC = kTRUE;
  if(fEventPlaneManager->GetDetectorId(VAR::kFMD  )!=-1)    fFillFMD = kTRUE;
  if(fEventPlaneManager->GetDetectorId(VAR::kFMDraw)!=-1)fFillRawFMD = kTRUE;
  if(fEventPlaneManager->GetDetectorId(VAR::kSPD)!=-1)      fFillSPD = kTRUE;


}


//__________________________________________________________________
void AliQnCorrectionsFillEvent::Process(AliAnalysisTaskSE* task, AliVEvent* event, Float_t* values) {

  fEvent=event;

  TString aod = "AOD";
  TString esd = "ESD";

  fIsAOD = ( aod.EqualTo(fEvent->Whoami()) ? kTRUE : kFALSE );
  fIsESD = ( esd.EqualTo(fEvent->Whoami()) ? kTRUE : kFALSE );

  FillEventInfo(values);
  FillDetectors(task, values);

}

//__________________________________________________________________
void AliQnCorrectionsFillEvent::FillEventInfo(Float_t* values) {
  //
  // fill event info
  //


  values[VAR::kRunNo]       = fEvent->GetRunNumber();
  //values[VAR::kTriggerMask] = event->TriggerMask();
  //values[VAR::kIsPhysicsSelection]  = (event->IsPhysicsSelection() ? 1.0 : 0.0);
  //values[VAR::kNVtxTPCContributors] = event->VertexTPCContributors();
  values[VAR::kVtxX]        = -999.;
  values[VAR::kVtxY]        = -999.;
  values[VAR::kVtxZ]        = -999.;
  const AliVVertex *primVtx = fEvent->GetPrimaryVertex();
  if (primVtx){
    values[VAR::kVtxX]        = primVtx->GetX();
    values[VAR::kVtxY]        = primVtx->GetY();
    values[VAR::kVtxZ]        = primVtx->GetZ();
    values[VAR::kNVtxContributors]    = primVtx->GetNContributors();
  }

  AliMultSelection *MultSelection = (AliMultSelection * ) fEvent->FindListObject("MultSelection");
  if(MultSelection){
    values[VAR::kVZEROMultPercentile] = MultSelection->GetMultiplicityPercentile("V0M", kTRUE);
  }

  AliESDEvent* esdEvent = static_cast<AliESDEvent*>(fEvent);
  AliCentrality* cent = esdEvent->GetCentrality();
  if(cent){
    values[VAR::kCentVZERO]   = cent->GetCentralityPercentile("V0M");
    values[VAR::kCentSPD]     = cent->GetCentralityPercentile("CL1");
    values[VAR::kCentTPC]     = cent->GetCentralityPercentile("TRK");
    values[VAR::kCentQuality] = cent->GetQuality();
  }


  AliVVZERO* vzero = fEvent->GetVZEROData();
  values[VAR::kVZEROATotalMult]     = vzero->GetMTotV0A();
  values[VAR::kVZEROCTotalMult]     = vzero->GetMTotV0C();
  values[VAR::kVZEROTotalMult]      = values[VAR::kVZEROATotalMult]+values[VAR::kVZEROCTotalMult];

  AliMultiplicity* spdmult = (AliMultiplicity*) fEvent->GetMultiplicity();
  values[VAR::kSPDntracklets]      = spdmult->GetNumberOfTracklets();
  values[VAR::kSPDnSingleClusters] = spdmult->GetNumberOfSingleClusters();
  //fEventPlaneManager->EventCuts()->IsSelected(values);
}



//__________________________________________________________________
void AliQnCorrectionsFillEvent::FillTrackInfo(AliVParticle* particle, Float_t* values) {

  Float_t dcaxy=0.0;
  Float_t dcaz=0.0;
  //particle->GetImpactParameters(dcaxy,dcaz);

  values[VAR::kPx]        = particle->Px();
  values[VAR::kPy]        = particle->Py();
  values[VAR::kPz]        = particle->Pz();
  values[VAR::kPt]        = particle->Pt();
  values[VAR::kP]         = particle->P();
  values[VAR::kPhi]       = particle->Phi();
  values[VAR::kTheta]     = particle->Theta();
  values[VAR::kEta]       = particle->Eta();
  values[VAR::kCharge]    = particle->Charge();
  values[VAR::kDcaXY]     = dcaxy;
  values[VAR::kDcaZ]      = dcaz;

  AliAODTrack* aodTrack=static_cast<AliAODTrack*>(particle);

  //values[VAR::kITSncls]       = particle->GetNcls(0); 
  values[VAR::kTPCncls]       = aodTrack->GetTPCNcls();
  values[VAR::kTPCchi2]       = aodTrack->Chi2perNDF();
  values[VAR::kTPCsignal]     = aodTrack->GetTPCsignal();
  for(Int_t ibit=0; ibit<9; ibit++) values[VAR::kFilterBit+ibit]     = aodTrack->TestFilterBit(BIT(ibit));
  values[VAR::kFilterBitMask768]  = (aodTrack->TestFilterBit(BIT(8))||aodTrack->TestFilterBit(BIT(9)));

  //delete aodTrack;

}


//__________________________________________________________________
void AliQnCorrectionsFillEvent::FillTrackInfo(AliESDtrack* particle, Float_t* values) {

  Float_t dcaxy=0.0;
  Float_t dcaz=0.0;
  particle->GetImpactParameters(dcaxy,dcaz);

  values[VAR::kPx]        = particle->Px();
  values[VAR::kPy]        = particle->Py();
  values[VAR::kPz]        = particle->Pz();
  values[VAR::kPt]        = particle->Pt();
  values[VAR::kP]         = particle->P();
  values[VAR::kPhi]       = particle->Phi();
  values[VAR::kTheta]     = particle->Theta();
  values[VAR::kEta]       = particle->Eta();
  values[VAR::kCharge]    = particle->Charge();
  values[VAR::kDcaXY]     = dcaxy;
  values[VAR::kDcaZ]      = dcaz;

  //values[VAR::kITSncls]       = particle->GetNcls(0); 
  values[VAR::kTPCncls]       = particle->GetTPCNcls();
  values[VAR::kTPCnclsIter1]  = particle->GetTPCNclsIter1();
  values[VAR::kTPCchi2]       = values[VAR::kTPCncls]>0 ? particle->GetTPCchi2()/values[VAR::kTPCncls] : 0.0;
  values[VAR::kTPCchi2Iter1]  = values[VAR::kTPCnclsIter1]>0 ? particle->GetTPCchi2Iter1()/values[VAR::kTPCnclsIter1] : 0.0;
  values[VAR::kTPCsignal]     = particle->GetTPCsignal();



}

//_________________________________
void AliQnCorrectionsFillEvent::FillDetectors(AliAnalysisTaskSE* task, Float_t* values){

  if(fFillTPC)   FillTPC(values);
  if(fFillVZERO) FillVZERO();
  if(fFillZDC)   FillZDC();
  if(fFillTZERO) FillTZERO();
  if(fFillFMD)   FillFMD(task);
  if(fFillRawFMD)FillRawFMD(values);
  if(fFillSPD) FillSPDTracklets(values);

}


//_________________________________
void AliQnCorrectionsFillEvent::FillTPC(Float_t* values){
  //
  // fill TPC info
  //

  if(fIsAOD) FillAodTPC(values);
  if(fIsESD) FillEsdTPC(values);


}


//_________________________________
void AliQnCorrectionsFillEvent::FillAodTPC(Float_t* values){
  //
  // fill AOD TPC info
  //

  AliVParticle* vTrack;
  AliQnCorrectionsConfiguration* QnConf = 0x0;
  TClonesArray* QnConfList = fEventPlaneManager->GetQnConfigurations(fEventPlaneManager->GetDetectorId(VAR::kTPC));

  for (Int_t iTrack = 0; iTrack < fEvent->GetNumberOfTracks(); ++iTrack)
  {
    vTrack = fEvent->GetTrack(iTrack); //carefull do not modify it othwise  need to work with a copy 
    if (!vTrack) continue;

    FillTrackInfo(vTrack, values);
    fQAhistos->FillHistClass("TrackQA_NoCuts", values);

    fEventPlaneManager->AddDataVector(VAR::kTPC, vTrack->Phi());

    for(Int_t iconf=0; iconf<QnConfList->GetEntriesFast(); iconf++){
      QnConf = (AliQnCorrectionsConfiguration*) QnConfList->At(iconf);
      if(QnConf->PassCuts(values)){
        fQAhistos->FillHistClass("TrackQA_"+QnConf->QnConfigurationName(), values);
      }
    }

  }



}




//_________________________________
void AliQnCorrectionsFillEvent::FillEsdTPC(Float_t* values){
  //
  // fill ESD TPC info
  //

  AliESDtrack* esdTrack;
  AliQnCorrectionsConfiguration* QnConf = 0x0;
  TClonesArray* QnConfList = fEventPlaneManager->GetQnConfigurations(fEventPlaneManager->GetDetectorId(VAR::kTPC));

  const AliESDEvent* esdEvent = static_cast<AliESDEvent*>(fEvent);
  const AliESDEvent& esd = *esdEvent;

  for (Int_t iTrack = 0; iTrack < fEvent->GetNumberOfTracks(); ++iTrack)
  {
    AliESDtrack* track = 0x0;
    esdTrack = esd.GetTrack(iTrack); //carefull do not modify it othwise  need to work with a copy 
    if(fUseTPCStandaloneTracks) AliESDtrack* track = AliESDtrackCuts::GetTPCOnlyTrack(const_cast<AliESDEvent*>(&esd),esdTrack->GetID());
    else track = esdTrack;
    if (!track) continue;

    FillTrackInfo(track, values);
    fQAhistos->FillHistClass("TrackQA_NoCuts", values);

    fEventPlaneManager->AddDataVector(VAR::kTPC, track->Phi());

    for(Int_t iconf=0; iconf<QnConfList->GetEntriesFast(); iconf++){
      QnConf = (AliQnCorrectionsConfiguration*) QnConfList->At(iconf);
      if(QnConf->PassCuts(values)){
        fQAhistos->FillHistClass("TrackQA_"+QnConf->QnConfigurationName(), values);
      }
    }

    if(fUseTPCStandaloneTracks) delete track;
  }
}



//_________________________________________________________________________________
void AliQnCorrectionsFillEvent::FillSPDTracklets(Float_t* values) {
  //
  // fill SPD info
  //

  AliQnCorrectionsConfiguration* QnConf = 0x0;
  TClonesArray* QnConfList = fEventPlaneManager->GetQnConfigurations(fEventPlaneManager->GetDetectorId(VAR::kSPD));

  Int_t nTracklets = 0;

  AliMultiplicity* mult = (AliMultiplicity*) fEvent->GetMultiplicity();
  nTracklets = mult->GetNumberOfTracklets();
  for(Int_t iTracklet=0; iTracklet<nTracklets; ++iTracklet) {
    values[VAR::kSPDtrackletEta]    = mult->GetEta(iTracklet);
    values[VAR::kSPDtrackletPhi]    = mult->GetPhi(iTracklet);

    fEventPlaneManager->AddDataVector(VAR::kSPD, values[VAR::kSPDtrackletPhi]);

    for(Int_t iconf=0; iconf<QnConfList->GetEntriesFast(); iconf++){
      QnConf = (AliQnCorrectionsConfiguration*) QnConfList->At(iconf);
      if(QnConf->Cuts()){
      if(QnConf->PassCuts(values)){
        fQAhistos->FillHistClass("TrackletQA_"+QnConf->QnConfigurationName(), values);
      }}
      else{
        fQAhistos->FillHistClass("TrackletQA_"+QnConf->QnConfigurationName(), values);
      }
    }
  }

}





//_________________________________
void AliQnCorrectionsFillEvent::FillVZERO(){
  //
  // fill VZERO info
  //

  Double_t weight=0.;
  const Double_t kX[8] = {0.92388, 0.38268, -0.38268, -0.92388, -0.92388, -0.38268, 0.38268, 0.92388};    // cosines of the angles of the VZERO sectors (8 per ring)
  const Double_t kY[8] = {0.38268, 0.92388, 0.92388, 0.38268, -0.38268, -0.92388, -0.92388, -0.38268};    // sines     -- " --

  AliVVZERO* vzero = fEvent->GetVZEROData();

  for(Int_t ich=0; ich<64; ich++){
    weight=vzero->GetMultiplicity(ich);
    if(weight>0.01) {
      fEventPlaneManager->AddDataVector(VAR::kVZERO, TMath::ATan2(kY[ich%8],kX[ich%8]), weight, ich);   // 1st ich is position in array, 2nd ich is channel id
    }

  }
}



//_________________________________
void AliQnCorrectionsFillEvent::FillTZERO(){
  //
  // fill ESD TZERO info
  //

  Double_t weight=0.0;
  const Double_t kX[24] = {/* Cside */ 0.905348,0.571718,0.0848977,-0.424671,-0.82045,-0.99639,-0.905348,-0.571718,-0.0848977,0.424671,0.82045,0.99639, /* Aside */ 0.99995,0.870982,0.508635,0.00999978,-0.491315,-0.860982,-0.99995,-0.870982,-0.508635,-0.0100001,0.491315,0.860982};
  const Double_t kY[24] = {/* Cside */ 0.424671,0.82045,0.99639,0.905348,0.571718,0.0848976,-0.424671,-0.82045,-0.99639,-0.905348,-0.571719,-0.0848975, /* Aside */ -0.00999983,0.491315,0.860982,0.99995,0.870982,0.508635,0.00999974,-0.491315,-0.860982,-0.99995,-0.870982,-0.508635};

  const AliESDTZERO* tzero= ((AliESDEvent*)fEvent)->GetESDTZERO();


  for(Int_t ich=0; ich<24; ich++){
    weight=tzero->GetT0amplitude()[ich];
    if(weight>0.01) {
      fEventPlaneManager->AddDataVector(VAR::kTZERO, TMath::ATan2(kY[ich],kX[ich]), weight, ich);   // 1st ich is position in array, 2nd ich is channel id
    }

  }
}




//_________________________________
void AliQnCorrectionsFillEvent::FillZDC(){
  //
  // fill ZDC info
  //


  Double_t weight=0.0;
  const Double_t kX[10] = { /* Cside */ 0.0,  -1.75,  1.75, -1.75, 1.75, /* Aside */  0.0,  1.75, -1.75, 1.75, -1.75  };
  const Double_t kY[10] = { /* Cside */ 0.0,  -1.75, -1.75,  1.75, 1.75, /* Aside */  0.0, -1.75, -1.75, 1.75,  1.75  };


  AliVZDC* zdc = (AliVZDC*) fEvent->GetZDCData();

  Double_t ZDCenergy[10];
  for(Int_t i=0; i<5; ++i)    ZDCenergy[i]  = zdc->GetZNCTowerEnergy()[i];
  for(Int_t i=5; i<10; ++i)   ZDCenergy[i]  = zdc->GetZNATowerEnergy()[i-5];

  for(Int_t ich=1; ich<10; ich++){
    if(ich==5) continue;
    weight=ZDCenergy[ich];
    if(weight>100.) {
      fEventPlaneManager->AddDataVector(VAR::kZDC, TMath::ATan2(kY[ich],kX[ich]), weight, ich);   // 1st ich is position in array, 2nd ich is channel id
    }
  }
}


//_________________________________________________________________________________
void AliQnCorrectionsFillEvent::FillFMD(AliAnalysisTaskSE* task)
{
  //
  // fill ESD FMD info
  //

  Float_t m,eta,phi;

  AliAODEvent* aodEvent = AliForwardUtil::GetAODEvent(task);
  if (!aodEvent) {cout<<"didn't get AOD"<<endl; return;}


  TObject* obj = aodEvent->FindListObject("Forward");  
  if (!obj) return;

  AliAODForwardMult* aodForward = static_cast<AliAODForwardMult*>(obj);

  // //if (!aodForward->CheckEvent(mask,ipZmin,ipZmax,cMin,cMax)) return 0;

  //Double_t ret = 0;
  const TH2D& d2Ndetadphi = aodForward->GetHistogram();

  //TH2D* fFMDhist=(TH2D*)d2Ndetadphi.Clone("fmdmap");


  Float_t FMDtotalmult=0.0;
  Int_t nEta = d2Ndetadphi.GetXaxis()->GetNbins();
  Int_t nPhi = d2Ndetadphi.GetYaxis()->GetNbins();


  // Loop over eta 
  Int_t nFMD=-1;
  for (Int_t iEta = 1; iEta <= nEta; iEta++) {
    Int_t valid = d2Ndetadphi.GetBinContent(iEta, 0);
    if (!valid) continue; // No data expected for this eta 

    eta = d2Ndetadphi.GetXaxis()->GetBinCenter(iEta);
    // Loop over phi 
    for (Int_t iPhi = 1; iPhi <= nPhi; iPhi++) {
      phi = d2Ndetadphi.GetYaxis()->GetBinCenter(iPhi);
      m     =  d2Ndetadphi.GetBinContent(iEta, iPhi);
      if(m<0.01) continue;
      nFMD++;

      fEventPlaneManager->AddDataVector(VAR::kFMD, phi, m, iEta*nPhi+iPhi);   // 1st ich is position in array, 2nd ich is channel id


    }
  }
}



//_________________________________
void AliQnCorrectionsFillEvent::FillRawFMD(Float_t* values)
{
  //
  // fill Raw FMD info
  //
  Bool_t isESD = (fEvent->IsA()==AliESDEvent::Class());
  if(!isESD) return;

  AliESDEvent* esdEvent = static_cast<AliESDEvent*>(fEvent);

  AliESDFMD* esdFmd = esdEvent->GetFMDData();

  Int_t id=-1;
  Int_t maxDet=3;
  Int_t maxRing=2;
  Int_t maxSector;
  Int_t maxStrip;
  Float_t m=0.0;
  Double_t phi,eta;
  Char_t ring;

  for(UShort_t det = 1; det <= maxDet; ++det) {
    (det == 1 ? maxRing=1 : maxRing=2);
    for(UShort_t ir = 0; ir < maxRing; ++ir) {
      ring = (ir == 0 ? 'I' : 'O');
      (ir == 0 ? maxSector=20 : maxSector=40);
      (ir == 0 ? maxStrip=512 : maxStrip=256);
      for(UShort_t sec = 0; sec < maxSector; ++sec) {
        phi  =  esdFmd->Phi(det, ring, sec, 0)/180.*TMath::Pi();
        for(UShort_t str = 0; str < maxStrip; ++str) {
          ++id;
          eta  =  esdFmd->Eta(det, ring, sec, str);
          m    =  esdFmd->Multiplicity(det, ring, sec, str);
          if(m ==  AliESDFMD::kInvalidMult) m=0;
          values[AliQnCorrectionsVarManager::kFMDEta] = eta;
          fEventPlaneManager->AddDataVector(VAR::kFMDraw, phi, m, id);   // 1st ich is position in array, 2nd ich is channel id
        }  // end loop over strips
      }  // end loop over sectors      
    }  // end loop over rings
  } // end loop over detectors
}


