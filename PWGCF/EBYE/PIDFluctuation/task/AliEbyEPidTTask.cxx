/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: ALICE Offline.                                                 *
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


//=========================================================================//
//             AliEbyE Analysis for Particle Ratio Fluctuation             //
//                   Deepika Rathee  | Satyajit Jena                       //
//                   drathee@cern.ch | sjena@cern.ch                       //
//                 Dealing with Wide pT Window Modified to ESDs            //
//=========================================================================//

#include "TChain.h"
#include "TList.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliVEvent.h"
#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliAODHeader.h"
#include "AliAODpidUtil.h"
#include "AliHelperPID.h"

#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliESDtrackCuts.h"
#include "AliESDInputHandler.h"
#include "AliESDpid.h"
#include "AliAODInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliCentrality.h"
#include "AliTracker.h"


using std::endl;
using std::cout;

#include "AliEbyEPidTTask.h"

ClassImp(AliEbyEPidTTask)

//-----------------------------------------------------------------------
AliEbyEPidTTask::AliEbyEPidTTask( const char *name ) : 
AliAnalysisTaskSE( name ), 
  fInputEventHandler(NULL),
  fESD(NULL),
  fAOD(NULL),
  fMCEvent(NULL),
  fStack(NULL), 
  fAODHandler(NULL),
  fESDHandler(NULL),
  fMCStack(NULL),
  fArrayMC(NULL),
  fESDtrackCuts(NULL),
  
  fThnList(NULL), 
  fAODtrackCutBit(128),
  fHelperPID(0x0),
  fEventCounter(NULL), 
  fPidCont(NULL),
  
  fVxMax(3.), 
  fVyMax(3.), 
  fVzMax(10.), 
  fPtMin(0.2),   
  fPtMax(3.), 
  fEtaMin(-1.), 
  fEtaMax(1.),  
  fDcaXy(10.),
  fDcaZ(10.),  
  
  fIsMC(kFALSE),
  fIsAOD(kFALSE),
  fDebug(kFALSE),
  fIsQa(kFALSE),
  
  fRunNumber(0),
  fNumberOfTracks(12000),
  fNumberOfTracksM(12000),
  fNTracks(0) {
   
  DefineOutput(1, TList::Class()); //! Connect Outpput....
  DefineOutput(2, TTree::Class()); //! Connect Outpput....
}

AliEbyEPidTTask::~AliEbyEPidTTask() {
  //!   Cleaning up
  if (fThnList)   delete fThnList;
  if (fHelperPID) delete fHelperPID;
  if (fPidCont)    delete fPidCont;
}

//---------------------------------------------------------------------------------
void AliEbyEPidTTask::UserCreateOutputObjects() {
  fThnList = new TList();
  fThnList->SetOwner(kTRUE);

 if (!fIsAOD) {
    if(!fESDtrackCuts)
      fESDtrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    else 
      Printf(" >>>>  User Track Cuts <<<< ");
    fESDtrackCuts->Print();
    
    Printf(" >>> DCAxy in TC [%8.4f:%8.4f]", 
	   fESDtrackCuts->GetMinDCAToVertexXY(), fESDtrackCuts->GetMaxDCAToVertexXY());
    Printf(" >>> DCAz in TC  [%8.4f:%8.4f]", 
	   fESDtrackCuts->GetMinDCAToVertexZ(), fESDtrackCuts->GetMaxDCAToVertexZ());
	       
    Float_t r1,r2;
    fESDtrackCuts->GetPtRange(r1,r2);
    Printf(" >>> Pt in TC  [%10.4f:%10.4f]",r1,r2);

    fESDtrackCuts->GetRapRange(r1,r2);
    Printf(" >>> Rap in TC [%10.4f:%10.4f]",r1,r2);

    fESDtrackCuts->GetEtaRange(r1,r2);
    Printf(" >>> Eta in TC [%10.4f:%10.4f]",r1,r2);
  }     



  fEventCounter = new TH1D("fEventCounter","EventCounter", 100, 0.5,100.5);
  fThnList->Add(fEventCounter);
  
  if (fHelperPID) {
    fThnList->Add(new TList);
    TList *list =  static_cast<TList*>(fThnList->Last());
    list->SetName("HelperPID");
    list->SetOwner(kTRUE);
    TList *ll = (TList*)fHelperPID->GetOutputList();
    for (Int_t ikey = 0; ikey < ll->GetEntries(); ikey++) {
      list->Add(ll->At(ikey));
    }
  }

  TDirectory *owd = gDirectory;
  OpenFile(1);
  fPidCont = new TTree("Event","fPidCont");
  owd->cd();

  fPidCont->Branch("cent",fCentrality,"fCentrality[6]/F");
  fPidCont->Branch("fRunNumber", &fRunNumber,  "fRunNumber/I");
  fPidCont->Branch("Trigger", fTrigMask,  "fTrigMask[5]/I");
  fPidCont->Branch("vertex",&fVtx,"fVtx[3]/F");
  fPidCont->Branch("fNumberOfTracks", &fNumberOfTracks,"fNumberOfTracks/I");
  fPidCont->Branch("fTrackPt",   fTrackPt,        "fTrackPt[fNumberOfTracks]/F");
  fPidCont->Branch("fTrackPhi",  fTrackPhi,       "fTrackPhi[fNumberOfTracks]/F");
  fPidCont->Branch("fTrackEta",  fTrackEta,       "fTrackEta[fNumberOfTracks]/F");
  fPidCont->Branch("fTrackDxy",  fTrackDxy,       "fTrackDxy[fNumberOfTracks]/I");
  fPidCont->Branch("fTrackDz",   fTrackDxy,       "fTrackDz[fNumberOfTracks]/F");
  fPidCont->Branch("fTrackPid",  fTrackPid,       "fTrackPid[fNumberOfTracks]/I");
  if(fIsMC) {
    fPidCont->Branch("fNumberOfTracksM", &fNumberOfTracksM, "fNumberOfTracksM/I");
    fPidCont->Branch("fTrackPtM",        fTrackPtM,        "fTrackPtM[fNumberOfTracksM]/F");
    fPidCont->Branch("fTrackPhiM",       fTrackPhiM,       "fTrackPhiM[fNumberOfTracksM]/F");
    fPidCont->Branch("fTrackEtaM",       fTrackEtaM,       "fTrackEtaM[fNumberOfTracksM]/F");
    fPidCont->Branch("fTrackPidM",       fTrackPidM,       "fTrackPidM[fNumberOfTracksM]/I");
  }

  PostData(1, fThnList);
  PostData(2, fPidCont);  
}

//----------------------------------------------------------------------------------
void AliEbyEPidTTask::UserExec( Option_t * ){
  fEventCounter->Fill(1);
 if (SetupEvent() < 0) {
   fEventCounter->Fill(2);
   PostData(1, fThnList); 
   PostData(2, fPidCont);
   return;
 }
 fEventCounter->Fill(3);
 fNTracks  = (fESD) ? fESD->GetNumberOfTracks() : fAOD->GetNumberOfTracks();  
 if (fIsMC && fIsAOD) {
   fArrayMC = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
   if (!fArrayMC)
     AliFatal("No array of MC particles found !!!"); 
 }
 fEventCounter->Fill(4);
 ExecEvents();
 fEventCounter->Fill(5);
 PostData(1, fThnList); 
 PostData(2, fPidCont);

}

//----------------------------------------------------------------------------------
void AliEbyEPidTTask::ExecEvents(){
  fEventCounter->Fill(6);
  Int_t iTracks = 0; 
  for (Int_t idxTrack = 0; idxTrack < fNTracks; ++idxTrack) {
    AliVTrack *track = (fESD) ? static_cast<AliVTrack*>(fESD->GetTrack(idxTrack)) : static_cast<AliVTrack*>(fAOD->GetTrack(idxTrack)); 
    
    if(!AcceptTrackL(track)) continue;
    
    Float_t dca[2], cov[3]; // 
    if (fESD)
      (dynamic_cast<AliESDtrack*>(track))->GetImpactParameters(dca, cov);
    else  {
      Double_t dcaa[2] = {-999,-999};
      Double_t cova[3] = {-999,-999,-999};
      AliAODTrack* clone =dynamic_cast<AliAODTrack*>(track->Clone("trk_clone"));
      Bool_t propagate = clone->PropagateToDCA(fAOD->GetPrimaryVertex(),fAOD->GetMagneticField(),100.,dcaa,cova);
      delete clone;  
      if (!propagate) dca[0] = -999;
      else dca[0] = Float_t(dcaa[0]);
    }
    
    if ( TMath::Abs(dca[0]) > fDcaXy ) continue;
    if ( TMath::Abs(dca[1]) > fDcaZ )  continue;

    Int_t icharge = track->Charge() < 0 ? -1 : 1;
    Int_t a = fHelperPID->GetParticleSpecies(track,kTRUE);
    Int_t b = -999;
    if (a == 0 )      b = 1;
    else if (a == 1 ) b = 2;
    else if (a == 2 ) b = 3;
    else              b = 4;    
    fTrackPt[iTracks]  = (Float_t)track->Pt();
    fTrackPhi[iTracks] = (Float_t)track->Phi();
    fTrackEta[iTracks] = (Float_t)track->Eta();
    fTrackDxy[iTracks] = dca[0];
    fTrackDz[iTracks]  = dca[1];
    fTrackPid[iTracks] = icharge*b;
    iTracks++;
  }
  fNumberOfTracks = iTracks;
  fEventCounter->Fill(7);
  //---- - -- - - - - -   -  -- - - - ---- - - - ---
  if (fIsMC) {
    Int_t mTracks = 0;
    fEventCounter->Fill(8);
    if (fIsAOD) {
      for (Int_t idxMC = 0; idxMC < fArrayMC->GetEntries(); idxMC++) {
	AliAODMCParticle *particle = static_cast<AliAODMCParticle*>(fArrayMC->At(idxMC));
	if (!particle) 
	  continue;
	if (!AcceptTrackLMC((AliVParticle*)particle, idxMC)) continue;
	Int_t icharge = (particle->PdgCode() < 0) ? -1 : 1;
	Int_t iPid = -999;  
	if      (TMath::Abs(particle->PdgCode()) ==  211) iPid = 1; // pion
	else if (TMath::Abs(particle->PdgCode()) ==  321) iPid = 2; // kaon
	else if (TMath::Abs(particle->PdgCode()) == 2212) iPid = 3; // proton
	else    iPid = 4;

	fTrackPtM[mTracks]     = (Float_t)particle->Pt();
	fTrackPhiM[mTracks]    = (Float_t)particle->Phi();
	fTrackEtaM[mTracks]    = (Float_t)particle->Eta();
	fTrackPidM[mTracks] = icharge*iPid;
      
	mTracks++;
      }
      fEventCounter->Fill(9);
      fNumberOfTracksM = mTracks;
    } else if (fESD) {
      fEventCounter->Fill(10);
      for (Int_t idxMC = 0; idxMC < fStack->GetNprimary(); ++idxMC) {
	AliVParticle* particle = fMCEvent->GetTrack(idxMC);
	if (!particle) 
	  continue;
	if (!AcceptTrackLMC(particle, idxMC)) continue;
	   Int_t icharge = (particle->PdgCode() < 0) ? -1 : 1;
	  
	   Int_t iPid = -999;  
	   if      (TMath::Abs(particle->PdgCode()) ==  211) iPid = 1; // pion
	   else if (TMath::Abs(particle->PdgCode()) ==  321) iPid = 2; // kaon
	   else if (TMath::Abs(particle->PdgCode()) == 2212) iPid = 3; // proton
	   else  iPid = 4;
	   
	   fTrackPtM[mTracks]     = (Float_t)particle->Pt();
	   fTrackPhiM[mTracks]    = (Float_t)particle->Phi();
	   fTrackEtaM[mTracks]    = (Float_t)particle->Eta();
	   fTrackPidM[mTracks] = icharge*iPid;
	   
	   mTracks++;
      }
      fEventCounter->Fill(11);
      fNumberOfTracksM = mTracks;
    }
  }
  fEventCounter->Fill(12);
  fPidCont->Fill();  
}

//________________________________________________________________________
Int_t AliEbyEPidTTask::SetupEvent() {
  ResetCurrentEvent();
  fEventCounter->Fill(13);
  if (!fIsAOD && (SetupESD() < 0)) {
    fEventCounter->Fill(14);
    AliError("ESD Event failed");
    return -1;
  }
  fEventCounter->Fill(15);
  if (fIsAOD && (SetupAOD() < 0)) {
    fEventCounter->Fill(16);
    AliError("AOD Event failed");
    return -1;
  }
  fEventCounter->Fill(17);
  if ((fIsMC && !fIsAOD) && (SetupMC() < 0)) {
    fEventCounter->Fill(18);
      AliError("MC Event failed");
      return -1;
  }
  
  fEventCounter->Fill(19);
  if(SetupEventCR(fESDHandler, fAODHandler, fMCEvent) < 0) {
    fEventCounter->Fill(20);
    AliError("MC Event failed");
    return -1;
  }
  fEventCounter->Fill(21);
  return RejectedEvent() ? -2 : 0;
}

//________________________________________________________________________
void AliEbyEPidTTask::ResetCurrentEvent() { 
  fEventCounter->Fill(22);
  fESD = NULL; 
  fAOD = NULL;
  if (fIsMC && !fIsAOD) 
    fMCEvent = NULL; 
  else if(fIsMC && fIsAOD) 
    fArrayMC = NULL; 
  return;
}


//________________________________________________________________________
Int_t AliEbyEPidTTask::SetupMC() {
  fEventCounter->Fill(23);
  AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*> 
    (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  
  if (!mcH) {
    AliError("MC event handler not available");
    return -1;
  }

  fEventCounter->Fill(24);
  fMCEvent = mcH->MCEvent();
  if (!fMCEvent) {
    AliError("MC event not available");
    return -1;
  }
  fEventCounter->Fill(25);
  AliHeader* header = fMCEvent->Header();
  if (!header) {
    AliError("MC header not available");
    return -1;
  }
  fEventCounter->Fill(26);
  fMCStack = fMCEvent->Stack(); 
  if (!fMCStack) {
    AliError("MC stack not available");
    return -1;
  }
  fEventCounter->Fill(27);
  if (!header->GenEventHeader()) {
    AliError("Could not retrieve genHeader from header");
    return -1;
  }
  fEventCounter->Fill(28);
  if (!fMCEvent->GetPrimaryVertex()){
    AliError("Could not get MC vertex");
    return -1;
  }
  fEventCounter->Fill(29);
  return 0;
}


//________________________________________________________________________
Int_t AliEbyEPidTTask::SetupAOD() {
  fAODHandler= dynamic_cast<AliAODInputHandler*> 
    (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  fEventCounter->Fill(30);
  if (!fAODHandler) {
    AliError("Could not get AOD input handler");
    return -1;
  } 
  fEventCounter->Fill(31);
  fAOD = fAODHandler->GetEvent();
  if (!fAOD) {
    AliError("Could not get AOD event");
    return -1;
  }
  fEventCounter->Fill(32);
  if (!fAODHandler->GetPIDResponse()) {
    AliError("Could not get PID response");
    return -1;
  } 
  fEventCounter->Fill(33);
  if (!fAOD->GetPrimaryVertex()) {
    AliError("Could not get primary vertex");
    return -1;
  }
  fEventCounter->Fill(34);
  if (!((AliVAODHeader*)fAOD->GetHeader())->GetCentralityP()) {
    AliError("Could not get centrality");
    return -1;
  }
  fEventCounter->Fill(35);
  return 0;
}


//________________________________________________________________________
Int_t AliEbyEPidTTask::SetupESD() {
  fEventCounter->Fill(36);
  fESDHandler= dynamic_cast<AliESDInputHandler*> 
    (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  fEventCounter->Fill(37);
  if (!fESDHandler) {
    AliError("Could not get ESD input handler");
    return -1;
  } 
  fEventCounter->Fill(38);
  fESD = fESDHandler->GetEvent();
  if (!fESD) {
    AliError("Could not get ESD event");
    return -1;
  }
  fEventCounter->Fill(39);
  if (!fESDHandler->GetPIDResponse()) {
    AliError("Could not get PID response");
    return -1;
  } 
  fEventCounter->Fill(40);
  if (!fESD->GetPrimaryVertexTracks()) {
    AliError("Could not get vertex from tracks");
    return -1;
  }
  fEventCounter->Fill(41);
  if (!fESD->GetCentrality()) {
    AliError("Could not get centrality");
    return -1;
  }
  fEventCounter->Fill(42);
  return 0;
}


//___________________________________________________________
Bool_t AliEbyEPidTTask::AcceptTrackL(AliVTrack *track) const {
 if (!track) 
   return kFALSE; 
 if (track->Charge() == 0) 
   return kFALSE; 
  
  
 if (fIsAOD) {  // AOD
 AliAODTrack * trackAOD = dynamic_cast<AliAODTrack*>(track);
 if (!trackAOD) {
   AliError("Pointer to dynamic_cast<AliAODTrack*>(track) = ZERO");
   return kFALSE; 
 }
 if (!trackAOD->TestFilterBit(fAODtrackCutBit))
   return kFALSE;
 } else {      // ESDs
   if(!fESDtrackCuts->AcceptTrack(dynamic_cast<AliESDtrack*>(track)))  return kFALSE;
 }

 if(track->Pt() < fPtMin || track->Pt() > fPtMax )  return kFALSE; 
 if (TMath::Abs(track->Eta()) > fEtaMax) return kFALSE; 

 return kTRUE;
}


//___________________________________________________________
Bool_t AliEbyEPidTTask::AcceptTrackLMC(AliVParticle *particle, Int_t idxMC) const {
  if(!particle) return kFALSE;

  if (particle->Charge() == 0.0) 
    return kFALSE;
  
  if (fIsAOD) {
    if(!(static_cast<AliAODMCParticle*>(particle))->IsPhysicalPrimary()) return kFALSE;
  } else {
    if(!fStack->IsPhysicalPrimary(idxMC)) return kFALSE;
  }
 
  if (particle->Pt() < fPtMin || particle->Pt() > fPtMax) return kFALSE;
  if (TMath::Abs(particle->Eta()) > fEtaMax) return kFALSE;

  return kTRUE;
}



//________________________________________________________________________

Int_t AliEbyEPidTTask::SetupEventCR(AliESDInputHandler *esdHandler, 
						   AliAODInputHandler *aodHandler, 
						   AliMCEvent *mcEvent) {
  fEventCounter->Fill(43);  
if(esdHandler){
    fInputEventHandler = static_cast<AliInputEventHandler*>(esdHandler);
    fESD               = dynamic_cast<AliESDEvent*>(fInputEventHandler->GetEvent());
    if (!fESD) {
      AliError("ESD event handler not available");
      return -1;
    }
  }
 else if(aodHandler){
    fInputEventHandler = static_cast<AliInputEventHandler*>(aodHandler);
    fAOD               = dynamic_cast<AliAODEvent*>(fInputEventHandler->GetEvent());
    if (!fAOD) {
      AliError("AOD event handler not available");
      return -1;
    }
  }
 fEventCounter->Fill(45);
  //  fPIDResponse = fInputEventHandler->GetPIDResponse();
  fMCEvent     = mcEvent;
  if (fMCEvent)
    fStack     = fMCEvent->Stack();
  
  AliCentrality *centrality = NULL;
  fEventCounter->Fill(46);
  if(esdHandler)
    centrality = fESD->GetCentrality();
  else if(aodHandler)
    centrality = ((AliVAODHeader*)fAOD->GetHeader())->GetCentralityP();

  fEventCounter->Fill(47);
  if (!centrality) {
    AliError("Centrality not available");
    return -1;
  }
  fEventCounter->Fill(48);
  fCentrality[0] = centrality->GetCentralityPercentile("V0M");
  fCentrality[1] = centrality->GetCentralityPercentile("V0A");
  fCentrality[2] = centrality->GetCentralityPercentile("V0C");
  fCentrality[3] = centrality->GetCentralityPercentile("CL1");
  fCentrality[4] = centrality->GetCentralityPercentile("ZNA");
  fCentrality[5] = centrality->GetCentralityPercentile("ZNC");
    
  return 0;
}
//________________________________________________________________________
Bool_t AliEbyEPidTTask::TriggeredEvents() {
  fEventCounter->Fill(49);
  Bool_t *aTriggerFired = new Bool_t[5];
  for (Int_t ii = 0; ii < 5; ++ii)
    aTriggerFired[ii] = 0;
  if ((fInputEventHandler->IsEventSelected() & AliVEvent::kMB))          aTriggerFired[0] = 1;
  if ((fInputEventHandler->IsEventSelected() & AliVEvent::kCentral))     aTriggerFired[1] = 1;
  if ((fInputEventHandler->IsEventSelected() & AliVEvent::kSemiCentral)) aTriggerFired[2] = 1;
  if ((fInputEventHandler->IsEventSelected() & AliVEvent::kEMCEJE))      aTriggerFired[3] = 1;
  if ((fInputEventHandler->IsEventSelected() & AliVEvent::kEMCEGA))      aTriggerFired[4] = 1;
  for (Int_t ii=0; ii<5; ++ii) {
    fTrigMask[ii] = aTriggerFired[ii];
  }

 Bool_t isTriggered = kFALSE;

  for (Int_t ii=0; ii<5; ++ii) {
    if(aTriggerFired[ii]) {
      isTriggered = kTRUE;
    }
  }
  
  delete[] aTriggerFired;
  return isTriggered;
}

//________________________________________________________________________
Bool_t AliEbyEPidTTask::RejectedEvent() {
  fEventCounter->Fill(50);
  if (!TriggeredEvents()) return kFALSE;
  const AliESDVertex* vtxESD = NULL;
  const AliAODVertex* vtxAOD = NULL;
  if (fESD){
    vtxESD = fESD->GetPrimaryVertexTracks();
    if (!vtxESD) return kFALSE;

    if(TMath::Abs(vtxESD->GetX()) > fVxMax) return kFALSE;
    if(TMath::Abs(vtxESD->GetY()) > fVyMax) return kFALSE;
    if(TMath::Abs(vtxESD->GetZ()) > fVzMax) return kFALSE;
    fVtx[0] = vtxESD->GetX();
    fVtx[1] = vtxESD->GetY();
    fVtx[2] = vtxESD->GetZ();
  }
  else if (fAOD){
    vtxAOD = fAOD->GetPrimaryVertex();
    if (!vtxAOD) return kFALSE;
    
    if(TMath::Abs(vtxAOD->GetX()) > fVxMax)  return kFALSE;
    if(TMath::Abs(vtxAOD->GetY()) > fVyMax)  return kFALSE;
    if(TMath::Abs(vtxAOD->GetZ()) > fVzMax)  return kFALSE;
    fVtx[0] = vtxAOD->GetX();
    fVtx[1] = vtxAOD->GetY();
    fVtx[2] = vtxAOD->GetZ();
  } else 
    return kFALSE;
  
  return kTRUE;
}

//___________________________________________________________
void AliEbyEPidTTask::Terminate( Option_t * ){
  Info("AliEbyEPidTTask"," Task Successfully finished");
}
