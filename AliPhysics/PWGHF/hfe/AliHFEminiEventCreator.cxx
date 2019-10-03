/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

//------------------------------------------------------------------------------------------------
//  This is a derived version of  AliHFEreducedEventCreatorESD.cxx/h by M.Fasel <M.Fasel@gsi.de>
//  The contents are much less than  AliHFEreducedEventCreatorESD class.
//  Contributors: Nirbhay K. Behera, Jiyeon Kwon, Jonghan Park
//  At the moment, keep information only for HFE analysis meant for IP fit method.
//  Works both for ESD and AOD. Collision system: pp, pA and AA.
//------------------------------------------------------------------------------------------------


#include <TBits.h>
#include <TFile.h>
#include <TParticle.h>
#include <TArrayI.h>
#include <TString.h>
#include <TTree.h>
#include <TH1D.h>

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODMCHeader.h"
#include "AliAnalysisUtils.h"
#include "AliCentrality.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliAODHeader.h"
#include "AliMCEventHandler.h"
#include "AliMultiplicity.h"
#include "AliMultSelection.h"
#include "AliPIDResponse.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVCluster.h"
#include "AliVVZERO.h"
#include "AliVZDC.h"
#include "AliMagF.h"
#include <TGeoGlobalMagField.h>
#include "TTreeStream.h"
#include "AliEventplane.h"
#include "AliLog.h"

#include "AliHFEcuts.h"
#include "AliHFEextraCuts.h"
#include "AliHFEmcQA.h"
#include "AliHFEpidTPC.h"
#include "AliHFEsignalCuts.h"
#include "AliHFEV0taginfo.h"
#include "AliHFEminiEvent.h"
#include "AliHFEminiTrack.h"


#include "AliHFEminiEventCreator.h"

using std::cout;
using std::endl;

ClassImp(AliHFEminiEventCreator)

AliHFEminiEventCreator::AliHFEminiEventCreator():
  AliAnalysisTaskSE(),
  fEventNumber(0),
  fNclustersTPC(70),
  fNclustersTPCPID(0),
  fNclustersITS(2),
  fTPCChiSquare(4.),
  fRatioTPCNcluster(.6),
  fPtMin(0.5),
  fPtMax(100.),
  fEta(0.8),
  fVz(10.),
  fDCAxy(1.),
  fDCAz(2.),
  fProdVz(0.5),
  fSPDResolution(0.25),
  fNsigmaTPClow(-1.),
  fNsigmaTPChigh(3.),
  fNsigmaTOF(3.),
  fITSLayerStatus(kFALSE),
  fHasMCdata(kFALSE),
  fAODanalysis(kFALSE),
  fRemoveFirstEvent(kFALSE),
  fSelectSignalOnly(kFALSE),
  fCollisionSystem("pp"),
  fList(NULL),
  fHFEtree(NULL),
  fNevents(NULL),
  fNtracks(NULL),
  fESD(NULL),
  fAOD(NULL),
  fVevent(NULL),
  fMCEvent(NULL),
  fStack(NULL),
  fPIDResponse(NULL),
  fAODMCHeader(NULL),
  fAODArrayMCInfo(NULL),
  fInputHandler(NULL),
  fMCEventHandler(NULL),
  fAnalysisUtils(NULL),
  fHFEevent(NULL),
  fTrackCuts(NULL),
  fExtraCuts(NULL),
  fSignalCuts(NULL),
  fTPCpid(NULL),
  fV0Tagger(NULL)

{
  //
  // Default constructor
  //
}

AliHFEminiEventCreator::AliHFEminiEventCreator(const char *name):
  AliAnalysisTaskSE(name),
  fEventNumber(0),
  fNclustersTPC(70),
  fNclustersTPCPID(0),
  fNclustersITS(2),
  fTPCChiSquare(4.),
  fRatioTPCNcluster(.6),
  fPtMin(0.5),
  fPtMax(100.),
  fEta(0.8),
  fVz(10.),
  fDCAxy(1.),
  fDCAz(2.),
  fProdVz(0.5),
  fSPDResolution(0.25),
  fNsigmaTPClow(-1.),
  fNsigmaTPChigh(3.),
  fNsigmaTOF(3.),
  fITSLayerStatus(kFALSE),
  fHasMCdata(kFALSE),
  fAODanalysis(kFALSE),
  fRemoveFirstEvent(kFALSE),
  fSelectSignalOnly(kFALSE),
  fCollisionSystem("pp"),
  fList(NULL),
  fHFEtree(NULL),
  fNevents(NULL),
  fNtracks(NULL),
  fESD(NULL),
  fAOD(NULL),
  fVevent(NULL),
  fMCEvent(NULL),
  fStack(NULL),
  fPIDResponse(NULL),
  fAODMCHeader(NULL),
  fAODArrayMCInfo(NULL),
  fInputHandler(NULL),
  fMCEventHandler(NULL),
  fAnalysisUtils(NULL),
  fHFEevent(NULL),
  fTrackCuts(NULL),
  fExtraCuts(NULL),
  fSignalCuts(NULL),
  fTPCpid(NULL),
  fV0Tagger(NULL)
{
  //
  // Default constructor
  //
  fTPCpid = new AliHFEpidTPC("QAtpcPID");
  fAnalysisUtils = new AliAnalysisUtils;
  fV0Tagger = new AliHFEV0taginfo("Tagger");
  DefineOutput(1, TList::Class());
  
}

AliHFEminiEventCreator::~AliHFEminiEventCreator(){

  // Default destructor
  
  if(fList) delete fList;
  if(fTPCpid) delete fTPCpid;
  if(fAnalysisUtils) delete fAnalysisUtils;
  if(fV0Tagger) delete fV0Tagger;
  if(fHFEevent) delete fHFEevent;
  if(fSignalCuts) delete fSignalCuts;
  if(fTrackCuts) delete fTrackCuts;
  
}

void AliHFEminiEventCreator::UserCreateOutputObjects(){
  // Create tree, signal cuts and track cuts
  fList = new TList();
  fList->SetOwner(kTRUE);
  
  fInputHandler = dynamic_cast<AliVEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  
  fPIDResponse = dynamic_cast<AliInputEventHandler *>(fInputHandler)->GetPIDResponse();
  if(!fPIDResponse){
    AliError("No PID response task found !!");
  }

  //if(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()) fHasMCdata = kTRUE;
  //else fHasMCdata = kFALSE;

  //Set the analysis type: AOD or ESD
  if(!TString(fInputHandler->IsA()->GetName()).CompareTo("AliAODInputHandler")) fAODanalysis = kTRUE;
  else fAODanalysis = kFALSE;
  
  fAnalysisUtils = new AliAnalysisUtils();
  
  fSignalCuts = new AliHFEsignalCuts("HFEsignalCuts", "HFE MC Signal definition");
 
  fTrackCuts = new AliHFEcuts("fTrackCuts", "Basic HFE track cuts");
  fTrackCuts->CreateStandardCuts();
  if(fAODanalysis) fTrackCuts->SetAOD();
  fTrackCuts->SetPIDResponse(fPIDResponse);
  fTrackCuts->SetMaxChi2perClusterTPC( fTPCChiSquare );   // switch off this cut (for the moment);
  fTrackCuts->SetMinNClustersTPC( fNclustersTPC );
  fTrackCuts->SetMinNClustersTPCPID( fNclustersTPCPID );
  fTrackCuts->SetMinRatioTPCclusters( fRatioTPCNcluster );
  fTrackCuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);
  fTrackCuts->SetMinNClustersITS( fNclustersITS );
  fTrackCuts->SetCutITSpixel(AliHFEextraCuts::kBoth);
  fTrackCuts->SetCheckITSLayerStatus( fITSLayerStatus );
  fTrackCuts->SetEtaRange( -fEta, fEta );
  fTrackCuts->SetPtRange( fPtMin, fPtMax); // (min, max) = (0.1, 100.)
  fTrackCuts->SetRejectKinkDaughters();
  fTrackCuts->SetAcceptKinkMothers();
  fTrackCuts->SetUseCorrelationVertex();
  fTrackCuts->SetSPDVtxResolutionCut();
  fTrackCuts->SetMaxImpactParam(fDCAxy, fDCAz);
  fTrackCuts->SetUseMixedVertex(kTRUE);
  fTrackCuts->SetVertexRange( fVz );
  fTrackCuts->Initialize();

  fExtraCuts = new AliHFEextraCuts("hfeExtraCuts","HFE Extra Cuts");

  fHFEevent = new AliHFEminiEvent;
  //OpenFile(1);
  fHFEtree = new TTree("HFEtree", "HFE event tree");
  fHFEtree->Branch("HFEevent", "AliHFEminiEvent", fHFEevent, 256000, 0);
  fHFEtree->SetAutoSave(100000000);
  fList->Add(fHFEtree);
  
  fNevents = new TH1D("Nevents","Number of events", 4, 0., 4.);
  fList->Add(fNevents);
  fNtracks = new TH1D("Ntracks","Number of tracks", 7, 0., 7.);
  fList->Add(fNtracks);
  
  PostData(1, fList);
  
}

void AliHFEminiEventCreator::UserExec(Option_t *){
  //
  // User Exec: Fill debug Tree

  fNevents->Fill(1);

  if(!fInputHandler)
    fInputHandler = dynamic_cast<AliVEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  
  if(!fInputHandler){
    AliError("No InputHandler");
    return;
  }

 
  //Event setups--------------------------
  if(fHasMCdata){
    SetUpMCEvent();
  }
  
  if(fAODanalysis){//AOD events----
    SetUpAODAnalysis();
  }
  else{//ESD events-------
    SetUpESDAnalysis();
  }
  
  if(fRemoveFirstEvent){
    if(fAnalysisUtils->IsFirstEventInChunk(fVevent)) return;
  }
  
  AliDebug(1, "Event Selected");
  
 
  // Reject pile up------------------------------------- 
  if( RemovePileupEvents(fVevent) ) return;
  AliDebug(1, "No pile up-> Ready for Event analysis");

  
  //Event analysis-------->
  if(fAODanalysis) ExecAODEvent();
  else ExecESDEvent();

  AliDebug(1, "Event analysis finished--> next event !!");
 
}

//-------------------------------------------------------------------
//====================================================================
void AliHFEminiEventCreator::SetUpMCEvent(){
  
  fMCEventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if(!fMCEventHandler){
    AliError("No MC Event Handler available");
    return;
  }
  
  fMCEvent = fMCEventHandler->MCEvent();
  if(!fMCEvent){
    AliError("No MC Event available");
    return;
  }
  
  fSignalCuts->SetMCEvent(fMCEvent);
  fTrackCuts->SetMCEvent(fMCEvent);
  
  fStack = fMCEvent->Stack();
  if(!fStack){
    AliError("No MC Event Stack available");
    return;
  }
  
  
  
}
//-----------------------------------------------------------------------



//====================================================================
void AliHFEminiEventCreator::SetUpESDAnalysis(){
  
  fESD = dynamic_cast<AliESDEvent*>(fInputHandler->GetEvent());
  if(!fESD){
    AliError("No ESD Event available");
    return;
  }

  fVevent = dynamic_cast<AliVEvent*>(fInputHandler->GetEvent());
  if (!fVevent) {
    printf("ERROR: fVEvent not available\n");
    return;
  }

   fTrackCuts->SetRecEvent(fESD);

   if(!fTrackCuts->CheckEventCuts("fCutsEvRec", fESD)){
     AliDebug(1, "Event rejected by the event cuts\n");
     return;
   }
   
   if(fV0Tagger){
     fV0Tagger->Reset();
     fV0Tagger->TagV0Tracks(fESD);
   }
   

}

//-------------------------------------------------------------------

//====================================================================
void AliHFEminiEventCreator::SetUpAODAnalysis(){
  
  fAOD = dynamic_cast<AliAODEvent*>(fInputHandler->GetEvent());
  if(!fAOD){
    AliError("No AOD Event available");
    return;
  }

  fVevent = dynamic_cast<AliVEvent*>(fInputHandler->GetEvent());
  if (!fVevent) {
    AliError("ERROR: fVEvent not available");
    return;
  }

  if( fHasMCdata ){
    
    fAODMCHeader = dynamic_cast<AliAODMCHeader *>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
    if(!fAODMCHeader){ 
      AliError("No AliAODMCHeader");
      return;
    }
    
    fAODArrayMCInfo = dynamic_cast<TClonesArray *>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
    if(!fAODArrayMCInfo){ 
      AliError("No AOD MC particles");
      return;
    }
    
    if( fAODArrayMCInfo ){
      fSignalCuts->SetMCAODInfo(fAODArrayMCInfo);
      fTrackCuts->SetMCEvent(fAOD);
    }
    
  }
  
  fTrackCuts->SetRecEvent(fAOD);

  if(!fTrackCuts->CheckEventCuts("fCutsEvRec", fAOD)){
    AliDebug(1, "Event rejected by the event cuts\n");
    return;
  }

  
  if(fV0Tagger){
    fV0Tagger->Reset();
    fV0Tagger->TagV0Tracks(fAOD);
  }
  
}
//-------------------------------------------------------------------

//===================================================================
void AliHFEminiEventCreator::ExecAODEvent(){

 // Make Mini Event 
  fHFEevent->~AliHFEminiEvent(); 
  new(fHFEevent) AliHFEminiEvent();

  if(!fExtraCuts){
    fExtraCuts = new AliHFEextraCuts("hfeExtraCuts","HFE Extra Cuts");
  }
  fExtraCuts->SetRecEventInfo(fAOD);

  //Check the vertex cut ----
  if( !ProperVertex(fAOD) ) return;
  
  if(fCollisionSystem == "pp"){ //pp events
    GetSPDnVZEROMultiplicity(fAOD);
  }
  else{//only pPb or PbPb events ----
    GetEventCentrality(fAOD);
  }
  
  //Get the Nch multiplicity--
  if( fHasMCdata ){
    Int_t nMCparticles = GetMCPrimaryNch();
    if(nMCparticles != 0) fHFEevent->SetPrimaryNchMC(nMCparticles);
  }
  
  
  AliAODTrack copyTrack;
  AliAODMCParticle *mctrack(NULL);
  AliAODTrack *aodTrack = 0x0;  //track->aodTrack
  
  for(Int_t itrack = 0; itrack < fAOD->GetNumberOfTracks(); itrack++){
    
    fNtracks->Fill(1); 
    aodTrack = dynamic_cast<AliAODTrack *>(fAOD->GetTrack(itrack));
    if(!aodTrack) continue;

    if(fHasMCdata){
      Int_t tempLabel = aodTrack->GetLabel();
      if(tempLabel >= 0){
      AliAODMCParticle *aodMCparticle = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(tempLabel));
      if(!aodMCparticle) continue; 
      if( TMath::Abs(aodMCparticle->PdgCode()) != 11) continue;// only electrons
      }
    }

    if (!aodTrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue;
    // Cut track (Only basic track cuts)
    if(!fTrackCuts->CheckParticleCuts(AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepRecKineITSTPC, aodTrack)) continue;
    fNtracks->Fill(2);

    //ITS layer cut-----------
    //if( ! (aodTrack->HasPointOnITSLayer(0) && aodTrack->HasPointOnITSLayer(1)) ) continue; //straight forward way--
    if(!fTrackCuts->CheckParticleCuts(AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepHFEcutsITS, aodTrack)) continue;
    fNtracks->Fill(3);
    
    //HFE DCA
    Float_t dcaxy = -999., dcaz = -999.;
    fExtraCuts->GetImpactParameters((AliVTrack *)aodTrack, dcaxy, dcaz);
    if(TMath::Abs(dcaxy) > fDCAxy) continue;
    if(TMath::Abs(dcaz)  > fDCAz) continue;

    if(!fHasMCdata){ //PID only for data

      //TOF nsigma cut---
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTOF(aodTrack, AliPID::kElectron)) > fNsigmaTOF) continue;
      fNtracks->Fill(4);
      
      if(fTPCpid->HasEtaCorrection()) fTPCpid->ApplyEtaCorrection(&copyTrack, AliHFEpidObject::kAODanalysis);
      copyTrack.~AliAODTrack();
      new(&copyTrack) AliAODTrack(*aodTrack);
      //TPC nsigma cut--
      if(fPIDResponse->NumberOfSigmasTPC(&copyTrack, AliPID::kElectron) < fNsigmaTPClow || fPIDResponse->NumberOfSigmasTPC(&copyTrack, AliPID::kElectron) > fNsigmaTPChigh ) continue;
      
    }//PID
    
    fNtracks->Fill(5);
    
    //Store track info in the minitree---
    AliHFEminiTrack hfetrack;
    hfetrack.SetSignedPt(aodTrack->Pt(), aodTrack->Charge() > 0);
    
    Double_t hfeImpactParam(-999.), hfeImpactParamResol(-999.);
    fExtraCuts->GetHFEImpactParameters((AliVTrack *)aodTrack, hfeImpactParam, hfeImpactParamResol);
    hfetrack.SetHFEImpactParam(hfeImpactParam, hfeImpactParamResol);
    
    fNtracks->Fill(6);
    
    if(fHasMCdata){
      //Fill Monte-Carlo Information
      Int_t label = TMath::Abs(aodTrack->GetLabel());
      if(label < fAODArrayMCInfo->GetEntriesFast())
	mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(label));
      
      if(mctrack){
        AliDebug(2, "Associated MC particle found");
	//        if(fTrackCuts->CheckParticleCuts(static_cast<UInt_t>(AliHFEcuts::kStepMCGenerated), mctrack)) hfetrack.SetMCSignal();
	
        // Kinematics
        hfetrack.SetMCSignedPt(mctrack->Pt(),mctrack->Charge() > 0.);
        hfetrack.SetMCPDG(mctrack->PdgCode());
	
        // Get Mother PDG code of the particle
        Int_t motherlabel = TMath::Abs(mctrack->GetMother());
        if(motherlabel >= 0 && motherlabel < fAODArrayMCInfo->GetEntriesFast()){
          AliAODMCParticle *mother = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(motherlabel));
          if(mother) 
	    hfetrack.SetMCMotherPdg(mother->PdgCode());   
        }
        
        // derive source
        hfetrack.SetMCSource(static_cast<Int_t>(fSignalCuts->GetSignalSource(aodTrack))); 
        Double_t mcmpt = -1;
        Int_t mcelectronSource = fSignalCuts->GetMCQAObject()->GetElecSource(mctrack,kTRUE, mcmpt);
        hfetrack.SetMCElectronSource(mcelectronSource);
        hfetrack.SetMCElectronSourcePt(mcmpt);

	fNtracks->Fill(7);
	
      } else {
        AliDebug(2, "Associated MC particle not found");
      }
      //cout <<" HFE task done " << endl; 
    }
    
    
    fHFEevent->AddTrack(&hfetrack);
  }
  
  fEventNumber++;
  
  fHFEtree->Fill();
  fNevents->Fill(2);
  
  PostData(1, fList);
  
}

//-------------------------------------------------------------------

//===================================================================
void AliHFEminiEventCreator::ExecESDEvent(){

  // Make Mini Event 
  fHFEevent->~AliHFEminiEvent(); 
  new(fHFEevent)AliHFEminiEvent();

   if(!fExtraCuts){
    fExtraCuts = new AliHFEextraCuts("hfeExtraCuts","HFE Extra Cuts");
  }
  fExtraCuts->SetRecEventInfo(fESD);

  //Check the vertex cut ----
  if( !ProperVertex(fESD) ) return;
  
  if(fCollisionSystem == "pp"){ //pp events
    GetSPDnVZEROMultiplicity(fESD);
  }
  else{//only pPb or PbPb events ----
    GetEventCentrality(fESD);
  }
  
  //Get the Nch multiplicity--
  if( fHasMCdata ){
    Int_t nMCparticles = GetMCPrimaryNch();
    if(nMCparticles != 0) fHFEevent->SetPrimaryNchMC(nMCparticles);
  }
  
  
  AliESDtrack copyTrack;
  AliESDtrack *track = 0x0;
  AliMCParticle *mctrack(NULL);
  
  for(Int_t itrack = 0; itrack < fESD->GetNumberOfTracks(); itrack++){
    
    fNtracks->Fill(1); 
    track = dynamic_cast<AliESDtrack *>(fESD->GetTrack(itrack));
    if(!track) continue;
    
    if(fHasMCdata){
      Int_t tempLabel = track->GetLabel();
      if(tempLabel >= 0){
	AliMCParticle *esdMCparticle = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(tempLabel));
	if(!esdMCparticle) continue;
	if( TMath::Abs(esdMCparticle->PdgCode()) != 11) continue;// store only electrons
      }
    }

    // Cut track (Only basic track cuts)
    if(!fTrackCuts->CheckParticleCuts(AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepRecKineITSTPC, track)) continue;
    fNtracks->Fill(2);

    //ITS layer cut-----------
    //if( ! (track->HasPointOnITSLayer(0) && track->HasPointOnITSLayer(1)) ) continue; //straight forward way--
    if(!fTrackCuts->CheckParticleCuts(AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepHFEcutsITS, track)) continue;
    fNtracks->Fill(3);

    
    // HFE DCA
    Float_t dcaxy = -999., dcaz = -999.;
    fExtraCuts->GetImpactParameters((AliVTrack *)track, dcaxy, dcaz);
    if(TMath::Abs(dcaxy) > fDCAxy) continue;
    if(TMath::Abs(dcaz)  > fDCAz) continue;
    
    if(!fHasMCdata){//PID only for data----
      
      //TOF nsigma cut---
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kElectron)) > fNsigmaTOF) continue;
      fNtracks->Fill(4);
      
      copyTrack.~AliESDtrack();
      new(&copyTrack) AliESDtrack(*track);
      if(fTPCpid->HasEtaCorrection()) fTPCpid->ApplyEtaCorrection(&copyTrack, AliHFEpidObject::kESDanalysis);
      //TPC nsigma cut--
      if(fPIDResponse->NumberOfSigmasTPC(&copyTrack, AliPID::kElectron) < fNsigmaTPClow || fPIDResponse->NumberOfSigmasTPC(&copyTrack, AliPID::kElectron) > fNsigmaTPChigh ) continue;
      
    }//PID
    
    fNtracks->Fill(5);
    
    // Kinematics
    AliHFEminiTrack hfetrack;
    hfetrack.SetSignedPt(track->Pt(), track->Charge() > 0);
    
    Double_t hfeImpactParam(-999.), hfeImpactParamResol(-999.);
    fExtraCuts->GetHFEImpactParameters((AliVTrack *)track, hfeImpactParam, hfeImpactParamResol);
    hfetrack.SetHFEImpactParam(hfeImpactParam, hfeImpactParamResol);
    
    fNtracks->Fill(6);
    
    if(fHasMCdata){
      // Fill Monte-Carlo Information
      Int_t label = TMath::Abs(track->GetLabel());
      if(label < fMCEvent->GetNumberOfTracks())
	mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(label));
      
      if(mctrack){
        AliDebug(2, "Associated MC particle found");
	//        if(fTrackCuts->CheckParticleCuts(static_cast<UInt_t>(AliHFEcuts::kStepMCGenerated), mctrack)) hfetrack.SetMCSignal();
        // Kinematics
        hfetrack.SetMCSignedPt(mctrack->Pt(),mctrack->Charge() > 0.);
        hfetrack.SetMCPDG(mctrack->PdgCode());
	
        // Get Mother PDG code of the particle
        Int_t motherlabel = TMath::Abs(mctrack->GetMother());
        if(motherlabel >= 0 && motherlabel < fMCEvent->GetNumberOfTracks()){
          AliMCParticle *mother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(motherlabel));
          if(mother) 
	    hfetrack.SetMCMotherPdg(mother->PdgCode());   
        }
        
        // derive source
        hfetrack.SetMCSource(static_cast<Int_t>(fSignalCuts->GetSignalSource(track))); 
        Double_t mcmpt = -1;
        Int_t mcelectronSource = fSignalCuts->GetMCQAObject()->GetElecSource(mctrack,kTRUE, mcmpt);
        hfetrack.SetMCElectronSource(mcelectronSource);
        hfetrack.SetMCElectronSourcePt(mcmpt);
	
	fNtracks->Fill(7);
	
      } else {
        AliDebug(2, "Associated MC particle not found");
      }
    }
    
    
    fHFEevent->AddTrack(&hfetrack);
  }
  
  fEventNumber++;
  
  fHFEtree->Fill();
  fNevents->Fill(2);
  
  PostData(1, fList);
 
}
//-------------------------------------------------------------------
//===================================================================
Bool_t AliHFEminiEventCreator::RemovePileupEvents(AliVEvent *event){

  if(!fAnalysisUtils) return kTRUE;
  
  fAnalysisUtils->SetUseMVPlpSelection(kTRUE);
  fAnalysisUtils->SetUseOutOfBunchPileUp(kTRUE);
  
  if(fAnalysisUtils->IsPileUpMV(event)){
    AliDebug(1, "Pileup Event with multiple vertex\n");
    return kTRUE;
  }
  
  if(fAnalysisUtils->IsOutOfBunchPileUp(event)){
    AliDebug(1, "Pileup Event with out of Bunch\n");
    return kTRUE;
  }
  if(fAnalysisUtils->IsSPDClusterVsTrackletBG(event)){
    AliDebug(1, "SPD cluster vs Trklt BG cut\n");
    return kTRUE;
  }

  return kFALSE; //no pile up

}

//-------------------------------------------------------------------
//===================================================================
Bool_t AliHFEminiEventCreator::ProperVertex(AliVEvent *event){

//Get Primary Vertex---
  const AliVVertex *vertex = event->GetPrimaryVertex();
  if(!vertex){
    AliError("Primary vertex not found");
    return kFALSE;
  }
  
  Double_t vtx[3];
  vertex->GetXYZ(vtx);
  if( !vtx[2]) return kFALSE;
  if(TMath::Abs(vtx[2]) > fVz) return kFALSE;
  Double_t tempVZ;
  tempVZ = vtx[2];
  
  // Get Primary Vertex from SPD
  Double_t vtxSPD[3];
  Double_t vcovSPD[6];
  const AliVVertex *vertexSPD = event->GetPrimaryVertexSPD();
  if(!vertexSPD){
    AliError("Primary vertex SPD not found");
    return kFALSE;
  }
  else{
    vertexSPD->GetXYZ(vtxSPD);
    if(!(vertexSPD->GetNContributors())) return kFALSE;
    if(TMath::Abs(tempVZ-vtxSPD[2]) > fProdVz ) return kFALSE;
    vertex->GetCovarianceMatrix(vcovSPD);
    if(TMath::Sqrt(vcovSPD[5]) > fSPDResolution ) return kFALSE;
  }

  fHFEevent->SetVZ(vertex->GetZ()); //fill the Vz position.

  return kTRUE;

}
  
//===================================================================
//-------------------------------------------------------------------
void AliHFEminiEventCreator::GetSPDnVZEROMultiplicity(AliVEvent *event){
  
  //Get VZERO multiplicity--
  if(fAODanalysis ){
    AliAODVZERO *vzeroinfoAOD = dynamic_cast<AliAODVZERO *>( dynamic_cast<AliAODEvent *>(event)->GetVZEROData());
    if(vzeroinfoAOD) fHFEevent->SetV0Multiplicity(vzeroinfoAOD->GetMTotV0A(), vzeroinfoAOD->GetMTotV0C());
  }
  else {
    AliVVZERO *vzeroinfo = event->GetVZEROData();
    if(vzeroinfo) fHFEevent->SetV0Multiplicity(vzeroinfo->GetMTotV0A(), vzeroinfo->GetMTotV0C());
  }
  
  //Set SPD multiplicity
  Int_t rawTracklets = 0;
  Int_t nTracklets = 0;
  if( fAODanalysis ){ //J. Wangner ways
    AliAODTracklets *tracklets = ((AliAODEvent*)event)->GetTracklets();
    nTracklets = tracklets->GetNumberOfTracklets();
    for(Int_t imult = 0; imult < nTracklets; imult++){
      Double_t theta = tracklets->GetTheta(imult);
      Double_t eta = -TMath::Log(TMath::Tan(theta/2.0));
      if( TMath::Abs(eta) > 1. ) continue;
      rawTracklets += 1;
    }
    
  }
  else{
    nTracklets = ((AliESDEvent*)fESD)->GetMultiplicity()->GetNumberOfTracklets();	
    for(Int_t imult = 0; imult < nTracklets; imult++){
      Double_t eta = ((AliESDEvent*)fESD)->GetMultiplicity()->GetEta(imult);
      if( TMath::Abs(eta) > 1. ) continue;
      rawTracklets += 1;
    }
  }

  if(rawTracklets > 0) fHFEevent->SetSPDMultiplicity(rawTracklets); // skip zero SPD tracklet events 
  
}
//-------------------------------------------------------------------

//===================================================================
void AliHFEminiEventCreator::GetEventCentrality(AliVEvent *event){

  Double_t CentPercentile = 0.;
  AliCentrality *fCentrality = NULL;

  AliMultSelection *fMultSelection = fMultSelection = (AliMultSelection *) event->FindListObject("MultSelection");
 
  if(!fMultSelection ){
    AliWarning("AliMultSelection object not found!");
  }
  else CentPercentile = fMultSelection->GetMultiplicityPercentile("V0M", false); 
  
  //cout << "Centrality = " << CentPercentile << endl;
  
  Double_t centRange[12] = { 0.,5.,10., 20., 30., 40., 50., 60.,70.,80., 90., 100.00001};
  Int_t CentBin = -1;
  
  for(Int_t ibin = 0; ibin < 11; ibin++){
    if( CentPercentile >= centRange[ibin] && CentPercentile < centRange[ibin+1] ){
      CentBin = ibin;
      break;
    }//if
  }
  
  if(  CentBin == -1 ) CentBin = 11;
  
  fHFEevent->SetCentralityBin(CentBin);
  
}
//-------------------------------------------------------------------


//===================================================================
Int_t AliHFEminiEventCreator::GetMCPrimaryNch(){

  //Loop on MC tracks only to calculate PhyPrimary tacks only
  Int_t nMCparticles = 0;
  
  if( fAODanalysis ){ 
    for (Int_t igen = 0; igen < fAODArrayMCInfo->GetEntriesFast(); igen++){
      AliAODMCParticle *mctrack=(AliAODMCParticle*)fAODArrayMCInfo->UncheckedAt(igen);
      if(!mctrack) continue;
      if( !mctrack->IsPhysicalPrimary() ) continue;
      if( mctrack->Charge() == 0 ) continue;
      if( TMath::Abs(mctrack->Eta()) > 1.0 ) continue;
      nMCparticles += 1; 
    }
  }
  else{
    for(Int_t iMC = 0; iMC < fStack->GetNprimary(); ++iMC){
      AliVParticle* particle = fMCEvent->GetTrack(iMC);
      if(!particle) continue;
      if(!fStack->IsPhysicalPrimary(iMC))  continue;
      if( particle->Charge() == 0 ) continue;
      if( TMath::Abs(particle->Eta()) > 1. ) continue;
      nMCparticles += 1;
    }     
  }//else

  return nMCparticles;
  
}
//-------------------------------------------------------------------  


//===================================================================
Bool_t AliHFEminiEventCreator::IsTOFmismatch(const AliVTrack *const track, const AliPIDResponse *const pid) const {
  //
  // Is TOF mismatch
  //
  Double_t probs[AliPID::kSPECIESC];
  AliPIDResponse::EDetPidStatus status = pid->ComputeTOFProbability(track, AliPID::kSPECIESC, probs);
  return status == AliPIDResponse::kDetMismatch;
}


//===================================================================
void AliHFEminiEventCreator::Terminate(Option_t *){
  //
  // Terminate
  //
  AliInfo("terminating...\n");
  
}

