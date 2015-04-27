/*
***********************************************************
    Variable definitions for event plane correction framework
    Contact: Jaap Onderwaater, j.onderwaater@gsi.de, jacobus.onderwaater@cern.ch
    Based on work of Ionut-Cristian Arsene
***********************************************************
*/

#include <iostream>
#include <fstream>


#include "AliEventPlaneDetector.h"
#include "AliEventPlaneConfiguration.h"
#include "AliEventPlaneManager.h"
#include "AliEventPlaneVarManager.h"
#include "AliEventPlaneHistos.h"
#include <AliCentrality.h>
#include <AliVEvent.h>
#include <AliESDEvent.h>
#include <AliESDHeader.h>
#include <AliESDtrack.h>
#include <AliESDtrackCuts.h>
#include "AliEventPlaneESDVarManager.h"

ClassImp(AliEventPlaneESDVarManager)

#ifdef ALIEVENTPLANEVARMANAGER_H
#define VAR AliEventPlaneVarManager
#endif

//using namespace AliEventPlaneVarManager;
//void AliEventPlaneESDVarManager::FillEventInfo(TObject* event, Float_t* values);
//void AliEventPlaneESDVarManager::FillTrackInfo(TObject* particle, Float_t* values);

AliEventPlaneESDVarManager::AliEventPlaneESDVarManager() :
  TNamed("AliEventPlaneESDVarManager","Fill functions")
{
  //
  // Default constructor
  //
}

//_____________________________________________________________________________
AliEventPlaneESDVarManager::~AliEventPlaneESDVarManager()
{
  //
  // Destructor
  //
}


//__________________________________________________________________
void AliEventPlaneESDVarManager::FillEventInfo(AliVEvent* event, Float_t* values) {
  //
  // fill event wise info
  //
  values[VAR::kRunNo]       = event->GetRunNumber();
  //values[VAR::kTriggerMask] = event->TriggerMask();
  //values[VAR::kIsPhysicsSelection]  = (event->IsPhysicsSelection() ? 1.0 : 0.0);
  //values[VAR::kNVtxTPCContributors] = event->VertexTPCContributors();
  values[VAR::kVtxX]        = -999.;
  values[VAR::kVtxY]        = -999.;
  values[VAR::kVtxZ]        = -999.;
  const AliVVertex *primVtx = event->GetPrimaryVertex();
  if (primVtx){
    values[VAR::kVtxX]        = primVtx->GetX();
    values[VAR::kVtxY]        = primVtx->GetY();
    values[VAR::kVtxZ]        = primVtx->GetZ();
    values[VAR::kNVtxContributors]    = primVtx->GetNContributors();
  }

  AliESDEvent* esdEvent = static_cast<AliESDEvent*>(event);
  AliCentrality* cent = esdEvent->GetCentrality();
  if(cent){
    values[VAR::kCentVZERO]   = cent->GetCentralityPercentile("V0M");
    values[VAR::kCentSPD]     = cent->GetCentralityPercentile("CL1");
    values[VAR::kCentTPC]     = cent->GetCentralityPercentile("TRK");
    values[VAR::kCentQuality] = cent->GetQuality();
  }
    
}






//__________________________________________________________________
void AliEventPlaneESDVarManager::FillTrackInfo(AliESDtrack* particle, Float_t* values) {

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

  values[VAR::kITSncls]       = particle->GetNcls(0); 
  values[VAR::kTPCncls]       = particle->GetTPCNcls();
  values[VAR::kTPCnclsIter1]  = particle->GetTPCNclsIter1();
  values[VAR::kTPCchi2]       = particle->GetTPCchi2()/values[VAR::kTPCncls];
  values[VAR::kTPCchi2Iter1]  = particle->GetTPCchi2Iter1()/values[VAR::kTPCnclsIter1];
  values[VAR::kTPCsignal]     = particle->GetTPCsignal();
  ////values[VAR::kNclsTPCiter1]  = particle->GetTPCNclsIter1(); // TODO: get rid of the plain numbers
  //values[VAR::kNFclsTPC]      = particle->GetTPCNclsF();
  //values[VAR::kNFclsTPCr]     = particle->GetTPCClusterInfo(2,1);
  //values[VAR::kNFclsTPCrFrac] = particle->GetTPCClusterInfo(2);
  //values[VAR::kTPCsignalNfrac]= tpcNcls>0?tpcSignalN/tpcNcls:0;
  //values[VAR::kNclsTRD]       = particle->GetNcls(2); // TODO: get rid of the plain numbers
  //values[VAR::kTRDntracklets] = particle->GetTRDntracklets(); // TODO: GetTRDtracklets/GetTRDntracklets?
  //values[VAR::kTRDpidQuality] = particle->GetTRDpidQuality();
  //values[VAR::kTrackStatus]   = (Double_t)particle->GetStatus();
  //if (tpcNcls>0) values[AliDielectronVarContainer::kTPCchi2Cl] = particle->GetTPCchi2() / tpcNcls;

}



//_________________________________
void AliEventPlaneESDVarManager::FillTPC(AliEventPlaneManager* EPmanager, AliVEvent* event, Float_t* values){

  Int_t nTrack=-1;
  AliEventPlaneConfiguration* EPconf = 0x0;
  TClonesArray* epConfList = EPmanager->GetEventPlaneConfigurations(AliEventPlaneManager::kTPC);
  TClonesArray& detector = *(EPmanager->GetReducedDetector(AliEventPlaneManager::kTPC));

  const AliESDEvent* esdEvent = static_cast<AliESDEvent*>(event);
  const AliESDEvent& esd = *esdEvent;

  for (Int_t iTrack = 0; iTrack < esd.GetNumberOfTracks(); ++iTrack)
  {
    AliESDtrack* esdTrack = esd.GetTrack(iTrack); //carefull do not modify it othwise  need to work with a copy 
    AliESDtrack* track = AliESDtrackCuts::GetTPCOnlyTrack(const_cast<AliESDEvent*>(&esd),esdTrack->GetID());
    if (!track) continue;

    nTrack++;
    
    AliEventPlaneDetector *reducedDetector=new(detector[nTrack]) AliEventPlaneDetector();
    FillTrackInfo(track, values);

    reducedDetector->SetPhi(track->Phi());
    reducedDetector->SetX(TMath::Cos(track->Phi()));
    reducedDetector->SetY(TMath::Sin(track->Phi()));
    reducedDetector->SetWeight(1.);



    Bool_t once = kTRUE;
    Bool_t trackUsed = kFALSE;

    for(Int_t iconf=0; iconf<epConfList->GetEntriesFast(); iconf++){
      EPconf = (AliEventPlaneConfiguration*) epConfList->At(iconf);
      if(EPconf->IsTrackSelected(values)){
          reducedDetector->SetEventPlaneDetector( EPconf->LocalIndex() );
          AliEventPlaneHistos::Instance()->FillHistClass("Tracks_"+EPconf->EventPlaneDetectorName(), values);
          AliEventPlaneHistos::Instance()->FillHistClass("TrackQA_"+EPconf->EventPlaneDetectorName(), values);
          //if(fRunLightWeight) {if(EPconf->CalibrationStep()==0) AliEventPlaneHistos::FillHistClass("TrackQA_"+EPconf->EventPlaneDetectorName(), values);}
          //else AliEventPlaneHistos::FillHistClass("TrackQA_"+EPconf->EventPlaneDetectorName(), values, AliReducedVarManager::GetUsedVars());
          trackUsed = kTRUE;
      }
    }

    if(!trackUsed) nTrack--;

  }



}




//_________________________________
void AliEventPlaneESDVarManager::FillVZERO(AliEventPlaneManager* EPmanager,AliVEvent* event){

  Double_t weight=0.;
  const Double_t kX[8] = {0.92388, 0.38268, -0.38268, -0.92388, -0.92388, -0.38268, 0.38268, 0.92388};    // cosines of the angles of the VZERO sectors (8 per ring)
  const Double_t kY[8] = {0.38268, 0.92388, 0.92388, 0.38268, -0.38268, -0.92388, -0.92388, -0.38268};    // sines     -- " --
  AliEventPlaneConfiguration* EPconf = 0x0;
  AliVVZERO* vzero = event->GetVZEROData();

  for(Int_t ich=0; ich<64; ich++){
    weight=vzero->GetMultiplicity(ich);
    if(weight<EPmanager->VZEROminMult()) weight=0.;

    TClonesArray& detector = *(EPmanager->GetReducedDetector(AliEventPlaneManager::kVZERO));
    AliEventPlaneDetector *reducedDetector=new(detector[ich]) AliEventPlaneDetector();
    // copy vzero data and set respective coordinates
    reducedDetector->SetX(kX[ich%8]);
    reducedDetector->SetY(kY[ich%8]);
    reducedDetector->SetWeight(weight);
    reducedDetector->SetId(ich);
    reducedDetector->SetBin(0);

    //// set event plane subdetectors
    TClonesArray* epConfList = EPmanager->GetEventPlaneConfigurations(AliEventPlaneManager::kVZERO);
    for(Int_t iconf=0; iconf<epConfList->GetEntriesFast(); iconf++){
    EPconf = (AliEventPlaneConfiguration*) epConfList->At(iconf);
     if(!EPconf) continue;

    //TIter nextEPconf(GetEventPlaneConfigurations(AliEventPlaneManager::kVZERO));
    //while((EPconf=static_cast<AliEventPlaneConfiguration*>(nextEPconf()))) {
    //  if(!EPconf) continue;
      if(EPconf->UseChannel(ich)) reducedDetector->SetEventPlaneDetector(EPconf->LocalIndex());
    }
  }
}



//_________________________________
void AliEventPlaneESDVarManager::FillTZERO(AliEventPlaneManager* EPmanager,AliVEvent* event){

  Double_t weight=0.0;
  const Double_t kX[24] = {/* Cside */ 0.905348,0.571718,0.0848977,-0.424671,-0.82045,-0.99639,-0.905348,-0.571718,-0.0848977,0.424671,0.82045,0.99639, /* Aside */ 0.99995,0.870982,0.508635,0.00999978,-0.491315,-0.860982,-0.99995,-0.870982,-0.508635,-0.0100001,0.491315,0.860982};
  const Double_t kY[24] = {/* Cside */ 0.424671,0.82045,0.99639,0.905348,0.571718,0.0848976,-0.424671,-0.82045,-0.99639,-0.905348,-0.571719,-0.0848975, /* Aside */ -0.00999983,0.491315,0.860982,0.99995,0.870982,0.508635,0.00999974,-0.491315,-0.860982,-0.99995,-0.870982,-0.508635};

  AliEventPlaneConfiguration* EPconf = 0x0;

  const AliESDTZERO* tzero= ((AliESDEvent*)event)->GetESDTZERO();

  for(Int_t ich=0; ich<24; ich++){
    weight=tzero->GetT0amplitude()[ich];
    if(weight<EPmanager->TZEROminMult()) weight=0.;

    TClonesArray& detector = *(EPmanager->GetReducedDetector(AliEventPlaneManager::kTZERO));
    AliEventPlaneDetector *reducedDetector=new(detector[ich]) AliEventPlaneDetector();
    // copy tzero data and set respective coordinates
    reducedDetector->SetX(kX[ich]);
    reducedDetector->SetY(kY[ich]);
    reducedDetector->SetWeight(weight);
    reducedDetector->SetId(ich);
    reducedDetector->SetBin(0);



    // set event plane subdetectors
    TIter nextEPconf(EPmanager->GetEventPlaneConfigurations(AliEventPlaneManager::kTZERO));
    while((EPconf=static_cast<AliEventPlaneConfiguration*>(nextEPconf()))) {
      if(!EPconf) continue;
      if(EPconf->UseChannel(ich)) reducedDetector->SetEventPlaneDetector(EPconf->LocalIndex());
    }
  }
}



