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

/* $Id$ */
 
#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TChain.h>
#include <TFile.h>
#include <TList.h>

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisCuts.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"

#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliESD.h"
#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliAODVZERO.h"
#include "AliTOFHeader.h"
#include "AliAODTracklets.h"
#include "AliAODCaloCells.h"
#include "AliAODCaloTrigger.h"
#include "AliAODMCParticle.h"
#include "AliVEvent.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliInputEventHandler.h"
#include "AliMultiInputEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliTrackSelectionFactory.h"
#include "AliVTrackSelection.h"
#include "AliLog.h"
#include "AliAODDimuon.h"


ClassImp(AliAnalysisTaskSE)

////////////////////////////////////////////////////////////////////////
AliVHeader*      AliAnalysisTaskSE::fgAODHeader         = NULL;
AliTOFHeader*    AliAnalysisTaskSE::fgTOFHeader         = NULL;
AliAODVZERO*     AliAnalysisTaskSE::fgAODVZERO          = NULL;
TClonesArray*    AliAnalysisTaskSE::fgAODTracks         = NULL;
TClonesArray*    AliAnalysisTaskSE::fgAODVertices       = NULL;
TClonesArray*    AliAnalysisTaskSE::fgAODV0s            = NULL;
TClonesArray*    AliAnalysisTaskSE::fgAODPMDClusters    = NULL;
TClonesArray*    AliAnalysisTaskSE::fgAODJets           = NULL;
TClonesArray*    AliAnalysisTaskSE::fgAODFMDClusters    = NULL;
TClonesArray*    AliAnalysisTaskSE::fgAODCaloClusters   = NULL;
TClonesArray*    AliAnalysisTaskSE::fgAODMCParticles    = NULL;
AliAODTracklets* AliAnalysisTaskSE::fgAODTracklets      = NULL;
AliAODCaloCells* AliAnalysisTaskSE::fgAODEmcalCells     = NULL;
AliAODCaloCells* AliAnalysisTaskSE::fgAODPhosCells      = NULL;
AliAODCaloTrigger* AliAnalysisTaskSE::fgAODEMCALTrigger = NULL;
AliAODCaloTrigger* AliAnalysisTaskSE::fgAODPHOSTrigger  = NULL;
TClonesArray*    AliAnalysisTaskSE::fgAODDimuons        = NULL;
TClonesArray*    AliAnalysisTaskSE::fgAODHmpidRings     = NULL;

AliAnalysisTaskSE::AliAnalysisTaskSE():
    AliAnalysisTask(),
    fDebug(0),
    fEntry(0),
    fInputEvent(0x0),
    fESDfriend(0x0),
    fInputHandler(0x0),
    fOutputAOD(0x0),
    fMCEvent(0x0),
    fTreeA(0x0),
    fCurrentRunNumber(-1),
    fHistosQA(0x0),
    fOfflineTriggerMask(0),
    fMultiInputHandler(0),
    fMCEventHandler(0),
    fTrackSelectionFactory(0),
    fTrackSelection(0)
{
  // Default constructor
}

AliAnalysisTaskSE::AliAnalysisTaskSE(const char* name):
    AliAnalysisTask(name, "AnalysisTaskSE"),
    fDebug(0),
    fEntry(0),
    fInputEvent(0x0),
    fESDfriend(0x0),
    fInputHandler(0x0),
    fOutputAOD(0x0),
    fMCEvent(0x0),
    fTreeA(0x0),
    fCurrentRunNumber(-1),
    fHistosQA(0x0),
    fOfflineTriggerMask(0),
    fMultiInputHandler(0),
    fMCEventHandler(0),
    fTrackSelectionFactory(0),
    fTrackSelection(0)
{
  // Default constructor
    DefineInput (0, TChain::Class());
    DefineOutput(0,  TTree::Class());
}

AliAnalysisTaskSE::AliAnalysisTaskSE(const AliAnalysisTaskSE& obj):
    AliAnalysisTask(obj),
    fDebug(0),
    fEntry(0),
    fInputEvent(0x0),
    fESDfriend(0x0),
    fInputHandler(0x0),
    fOutputAOD(0x0),
    fMCEvent(0x0),
    fTreeA(0x0),
    fCurrentRunNumber(-1),
    fHistosQA(0x0),
    fOfflineTriggerMask(0),
    fMultiInputHandler(obj.fMultiInputHandler),
    fMCEventHandler(obj.fMCEventHandler),
    fTrackSelectionFactory(obj.fTrackSelectionFactory),
    fTrackSelection(obj.fTrackSelection)
{
// Copy constructor
    fDebug            = obj.fDebug;
    fEntry            = obj.fEntry;
    fInputEvent       = obj.fInputEvent;
    fESDfriend        = obj.fESDfriend;
    fInputHandler     = obj.fInputHandler;
    fOutputAOD        = obj.fOutputAOD;
    fMCEvent          = obj.fMCEvent;
    fTreeA            = obj.fTreeA;    
    fCurrentRunNumber = obj.fCurrentRunNumber;
    fHistosQA         = obj.fHistosQA;

}


AliAnalysisTaskSE& AliAnalysisTaskSE::operator=(const AliAnalysisTaskSE& other)
{
// Assignment
  if(&other == this) return *this;
  AliAnalysisTask::operator=(other);

    AliAnalysisTask::operator=(other);
    fDebug            = other.fDebug;
    fEntry            = other.fEntry;
    fInputEvent       = other.fInputEvent;
    fESDfriend        = other.fESDfriend;
    fInputHandler     = other.fInputHandler;
    fOutputAOD        = other.fOutputAOD;
    fMCEvent          = other.fMCEvent;
    fTreeA            = other.fTreeA;    
    fCurrentRunNumber = other.fCurrentRunNumber;
    fHistosQA         = other.fHistosQA;
    fOfflineTriggerMask = other.fOfflineTriggerMask;
    fMultiInputHandler  = other.fMultiInputHandler;
    fMCEventHandler     = other.fMCEventHandler;
    fTrackSelectionFactory = other.fTrackSelectionFactory;
    fTrackSelection = other.fTrackSelection;
    return *this;
}

//______________________________________________________________________________
void AliAnalysisTaskSE::ConnectInputData(Option_t* /*option*/)
{
// Connect the input data
    if (fDebug > 1) printf("AnalysisTaskSE::ConnectInputData() \n");

   // Connect input handlers (multi input handler is handled)
    ConnectMultiHandler();
    
    if (fInputHandler && fInputHandler->GetTree()) {
	if (fInputHandler->GetTree()->GetBranch("ESDfriend."))
	    fESDfriend = ((AliESDInputHandler*)fInputHandler)->GetESDfriend();

	fInputEvent = fInputHandler->GetEvent();
    } else if( fMCEvent ) {
         AliWarning("No Input Event Handler connected, only MC Truth Event Handler") ; 
    } else {
         AliError("No Input Event Handler connected") ; 
         return ; 
    }
    // Disconnect multi handler
    DisconnectMultiHandler();
}

void AliAnalysisTaskSE::CreateOutputObjects()
{
// Create the output container
//
//  Default AOD
    if (fDebug > 1) printf("AnalysisTaskSE::CreateOutPutData() \n");

    AliAODHandler* handler = dynamic_cast<AliAODHandler*> 
         ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    
    Bool_t merging = kFALSE;
    AliAODInputHandler* aodIH = static_cast<AliAODInputHandler*>((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
    if (aodIH) {
	if (aodIH->GetMergeEvents()) merging = kTRUE;
    }

    // Check if AOD replication has been required
    if (handler) {
	fOutputAOD   = handler->GetAOD();
	fTreeA = handler->GetTree();
	if (fOutputAOD && !(handler->IsStandard())) {
	    if ((handler->NeedsHeaderReplication() || merging) && !(fgAODHeader)) 
		{
		 if (fDebug > 1) AliInfo("Replicating header");
		 fgAODHeader = new AliAODHeader;
		 handler->AddBranch("AliAODHeader", &fgAODHeader);
		}
            if ((handler->NeedsTOFHeaderReplication() || merging) && !(fgTOFHeader))
                {
                 if (fDebug > 1) AliInfo("Replicating TOFheader");
                 fgTOFHeader = new AliTOFHeader;
                 handler->AddBranch("AliTOFHeader", &fgTOFHeader);
                }
            if ((handler->NeedsVZEROReplication() || merging) && !(fgAODVZERO))
                {
                 if (fDebug > 1) AliInfo("Replicating VZERO");
                 fgAODVZERO = new AliAODVZERO;
                 handler->AddBranch("AliAODVZERO", &fgAODVZERO);
                }

	    if ((handler->NeedsTracksBranchReplication() || merging) && !(fgAODTracks))      
	    {   
		if (fDebug > 1) AliInfo("Replicating track branch\n");
		fgAODTracks = new TClonesArray("AliAODTrack",500);
		fgAODTracks->SetName("tracks");
		handler->AddBranch("TClonesArray", &fgAODTracks);
	    }    
	    if ((handler->NeedsVerticesBranchReplication() || merging) && !(fgAODVertices))
	    {
		if (fDebug > 1) AliInfo("Replicating vertices branch\n");
		fgAODVertices = new TClonesArray("AliAODVertex",500);
		fgAODVertices->SetName("vertices");
		handler->AddBranch("TClonesArray", &fgAODVertices);
	    }	
	    if ((handler->NeedsV0sBranchReplication()) && !(fgAODV0s))	  
	    {   
		if (fDebug > 1) AliInfo("Replicating V0s branch\n");
		fgAODV0s = new TClonesArray("AliAODv0",500);
		fgAODV0s->SetName("v0s");
		handler->AddBranch("TClonesArray", &fgAODV0s);
	    }
	    if ((handler->NeedsTrackletsBranchReplication()) && !(fgAODTracklets))	  
	    {   
		if (fDebug > 1) AliInfo("Replicating Tracklets branch\n");
		fgAODTracklets = new AliAODTracklets("tracklets","tracklets");
		handler->AddBranch("AliAODTracklets", &fgAODTracklets);
	    }
	    if ((handler->NeedsPMDClustersBranchReplication()) && !(fgAODPMDClusters))	  
	    {   
		if (fDebug > 1) AliInfo("Replicating PMDClusters branch\n");
		fgAODPMDClusters = new TClonesArray("AliAODPmdCluster",500);
		fgAODPMDClusters->SetName("pmdClusters");
		handler->AddBranch("TClonesArray", &fgAODPMDClusters);
	    }
	    if ((handler->NeedsJetsBranchReplication() || merging) && !(fgAODJets))	  
	    {   
		if (fDebug > 1) AliInfo("Replicating Jets branch\n");
		fgAODJets = new TClonesArray("AliAODJet",500);
		fgAODJets->SetName("jets");
		handler->AddBranch("TClonesArray", &fgAODJets);
	    }
	    if ((handler->NeedsFMDClustersBranchReplication()) && !(fgAODFMDClusters))	  
	    {   
		AliInfo("Replicating FMDClusters branch\n");
		fgAODFMDClusters = new TClonesArray("AliAODFmdCluster",500);
		fgAODFMDClusters->SetName("fmdClusters");
		handler->AddBranch("TClonesArray", &fgAODFMDClusters);
	    }
	    if ((handler->NeedsCaloClustersBranchReplication() || merging) && !(fgAODCaloClusters))	  
	    {   
		if (fDebug > 1) AliInfo("Replicating CaloClusters branch\n");
		fgAODCaloClusters = new TClonesArray("AliAODCaloCluster",500);
		fgAODCaloClusters->SetName("caloClusters");
		handler->AddBranch("TClonesArray", &fgAODCaloClusters);

		fgAODEmcalCells = new AliAODCaloCells("emcalCells","emcalCells",AliVCaloCells::kEMCALCell);
		handler->AddBranch("AliAODCaloCells", &fgAODEmcalCells);

		fgAODPhosCells = new AliAODCaloCells("phosCells","phosCells",AliVCaloCells::kPHOSCell);
		handler->AddBranch("AliAODCaloCells", &fgAODPhosCells);
		}
	    if ((handler->NeedsCaloTriggerBranchReplication() || merging) && !(fgAODEMCALTrigger))	  
	    {   
		if (fDebug > 1) AliInfo("Replicating EMCAL Calo Trigger branches\n");
		fgAODEMCALTrigger = new AliAODCaloTrigger("emcalTrigger","emcalTrigger");
		handler->AddBranch("AliAODCaloTrigger", &fgAODEMCALTrigger);
		}
		if ((handler->NeedsCaloTriggerBranchReplication() || merging) && !(fgAODPHOSTrigger))	  
		{   
		if (fDebug > 1) AliInfo("Replicating PHOS Calo Trigger branches\n");
		fgAODPHOSTrigger = new AliAODCaloTrigger("phosTrigger","phosTrigger");
		handler->AddBranch("AliAODCaloTrigger", &fgAODPHOSTrigger);
	    }
	    if ((handler->NeedsMCParticlesBranchReplication() || merging) && !(fgAODMCParticles))	  
	    {   
		if (fDebug > 1) AliInfo("Replicating MCParticles branch\n");
		fgAODMCParticles = new TClonesArray("AliAODMCParticle",500);
		fgAODMCParticles->SetName("mcparticles");
		handler->AddBranch("TClonesArray", &fgAODMCParticles);
	    }
            if ((handler->NeedsDimuonsBranchReplication() || merging) && !(fgAODDimuons))      
	    {   
		if (fDebug > 1) AliInfo("Replicating dimuon branch\n");
		fgAODDimuons = new TClonesArray("AliAODDimuon",0);
		fgAODDimuons->SetName("dimuons");
		handler->AddBranch("TClonesArray", &fgAODDimuons);
	    }    
	    if ((handler->NeedsHMPIDBranchReplication() || merging) && !(fgAODHmpidRings))      
	    {   
		if (fDebug > 1) AliInfo("Replicating HMPID branch\n");
		fgAODHmpidRings = new TClonesArray("AliAODHMPIDrings",0);
		fgAODHmpidRings->SetName("hmpidRings");
		handler->AddBranch("TClonesArray", &fgAODHmpidRings);
	    }
                

	    // cache the pointerd in the AODEvent
	    fOutputAOD->GetStdContent();
	}
    }

    if(fTrackSelectionFactory && !fTrackSelection)
      fTrackSelection = fTrackSelectionFactory->CreateTrackCuts(
          AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->IsA() == AliAODInputHandler::Class() ?
              AliTrackSelectionFactory::kAOD :
              AliTrackSelectionFactory::kESD
              );
    ConnectMultiHandler();
    UserCreateOutputObjects();
    DisconnectMultiHandler();
}

void AliAnalysisTaskSE::Exec(Option_t* option)
{
//
// Exec analysis of one event
    
    ConnectMultiHandler();
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (mgr->GetDebugLevel() > 1) {
       if (!mgr->GetTopTasks()->FindObject(this))
          printf("    -> Executing sub-task %s\n", GetName());
    }   
    if ( fDebug >= 10)
      printf("Task is active %5d\n", IsActive());
    
    if (fDebug > 1) AliInfo("AliAnalysisTaskSE::Exec() \n");
//
    AliAODHandler* handler = dynamic_cast<AliAODHandler*>(mgr->GetOutputEventHandler());

    AliAODInputHandler* aodH = dynamic_cast<AliAODInputHandler*>(fInputHandler);
//
// Was event selected ? If no event selection mechanism, the event SHOULD be selected (AG)
    UInt_t isSelected = AliVEvent::kAny;
    if( fInputHandler && (fInputHandler->GetEventSelection() || aodH)) {
      // Get the actual offline trigger mask for the event and AND it with the
      // requested mask. If no mask requested select by default the event.
      if (fOfflineTriggerMask)
	isSelected = fOfflineTriggerMask & fInputHandler->IsEventSelected();
    }
//  Functionality below moved in the filter tasks (AG)
//    if (handler) handler->SetFillAOD(isSelected);

    if( fInputHandler ) {
      fEntry = fInputHandler->GetReadEntry();
      fESDfriend = ((AliESDInputHandler*)fInputHandler)->GetESDfriend();
    }
    

    // Notify the change of run number
    if (InputEvent() && (InputEvent()->GetRunNumber() != fCurrentRunNumber)) {
	fCurrentRunNumber = InputEvent()->GetRunNumber();
	NotifyRun();
    } else if( fMCEvent )
      fEntry = fMCEvent->Header()->GetEvent(); 

    if ( !((Entry()-1)%100) && fDebug > 0) 
      AliInfo(Form("%s ----> Processing event # %lld", CurrentFileName(), Entry()));
    
    if (aodH) fMCEvent = aodH->MCEvent();

    if (handler && aodH) {
	Bool_t merging = aodH->GetMergeEvents();
	// Do not analyze merged events if last embedded file has less events than normal event, 
	// skip analysis after last embeded event 
	if(merging){
	  if(aodH->GetReadEntry() + aodH->GetMergeOffset() >= aodH->GetTreeToMerge()->GetEntriesFast()){
	    // Do I need to add the lines before the return?
	    // Added protection in case the derived task is not an AOD producer.
	    AliAnalysisDataSlot *out0 = GetOutputSlot(0);
	    if (out0 && out0->IsConnected()) PostData(0, fTreeA);    
	    DisconnectMultiHandler();
	    return;
	  }
	}
      
	AliAODEvent* aod = dynamic_cast<AliAODEvent*>(InputEvent());

	if (aod && !(handler->IsStandard()) && !(handler->AODIsReplicated())) {
	    if ((handler->NeedsHeaderReplication() || merging) && (fgAODHeader))
	    {
	      // copy the contents by assigment
	      *fgAODHeader =  *(aod->GetHeader());
	    }
            if ((handler->NeedsTOFHeaderReplication() || merging) && (fgTOFHeader))
            {
              if (aod->GetTOFHeader()) *fgTOFHeader =  *(aod->GetTOFHeader());
            }
            if ((handler->NeedsVZEROReplication() || merging) && (fgAODVZERO) && aod->GetVZEROData())
            {
              *fgAODVZERO = *(aod->GetVZEROData());
            }

	    if ((handler->NeedsTracksBranchReplication() || (merging &&  aodH->GetMergeTracks())) && (fgAODTracks))
	    {
		TClonesArray* tracks = aod->GetTracks();
		new (fgAODTracks) TClonesArray(*tracks);
	    }
	    if ((handler->NeedsVerticesBranchReplication() || merging) && (fgAODVertices))
	    {
		TClonesArray* vertices = aod->GetVertices();
		new (fgAODVertices) TClonesArray(*vertices);
	    }
	    if ((handler->NeedsV0sBranchReplication()) && (fgAODV0s))
	    {
		TClonesArray* v0s = aod->GetV0s();
		new (fgAODV0s) TClonesArray(*v0s);
	    }
	    if ((handler->NeedsTrackletsBranchReplication()) && (fgAODTracklets))
	    {
	      *fgAODTracklets = *aod->GetTracklets();
	    }
	    if ((handler->NeedsPMDClustersBranchReplication()) && (fgAODPMDClusters))
	    {
		TClonesArray* pmdClusters = aod->GetPmdClusters();
		new (fgAODPMDClusters) TClonesArray(*pmdClusters);
	    }
	    if ((handler->NeedsJetsBranchReplication() || (merging &&aodH->GetMergeTracks())) && (fgAODJets))
	    {
		TClonesArray* jets = aod->GetJets();
		new (fgAODJets) TClonesArray(*jets);
	    }
	    if ((handler->NeedsFMDClustersBranchReplication()) && (fgAODFMDClusters))
	    {
		TClonesArray* fmdClusters = aod->GetFmdClusters();
		new (fgAODFMDClusters) TClonesArray(*fmdClusters);
	    }
	    if ((handler->NeedsCaloClustersBranchReplication() || 
		 (merging && (aodH->GetMergeEMCALClusters() || aodH->GetMergePHOSClusters()))) 
		&& (fgAODCaloClusters))
	    {
		TClonesArray* caloClusters = aod->GetCaloClusters();
		new (fgAODCaloClusters) TClonesArray(*caloClusters);
	    }

	    if ((handler->NeedsMCParticlesBranchReplication() || merging) && (fgAODMCParticles))
	    {
		TClonesArray* mcParticles = (TClonesArray*) (aod->FindListObject("mcparticles"));
		if( mcParticles )
		  new (fgAODMCParticles) TClonesArray(*mcParticles);
	    }
	    
	    if ((handler->NeedsDimuonsBranchReplication() || (merging && aodH->GetMergeTracks())) && (fgAODDimuons))
	    {
	        fgAODDimuons->Clear();
		TClonesArray& dimuons = *fgAODDimuons;
		TClonesArray& tracksnew = *fgAODTracks;
		
                Int_t nMuonTrack[100]; 
                for(Int_t imuon = 0; imuon < 100; imuon++) nMuonTrack[imuon] = 0;
                Int_t nMuons=0;
		for(Int_t ii=0; ii < fgAODTracks->GetEntries(); ii++){
		    AliAODTrack *track = (AliAODTrack*) fgAODTracks->At(ii);
		    if(track->IsMuonTrack()) {
			nMuonTrack[nMuons]= ii;
			nMuons++;
		    }  
		}
                Int_t jDimuons=0;
		if(nMuons >= 2){
		    for(Int_t i = 0; i < nMuons; i++){
			Int_t index0 = nMuonTrack[i];
			for(Int_t j = i+1; j < nMuons; j++){
			    Int_t index1 = nMuonTrack[j];
			    tracksnew.At(index0)->ResetBit(kIsReferenced);
			    tracksnew.At(index0)->SetUniqueID(0); 
			    tracksnew.At(index1)->ResetBit(kIsReferenced);
			    tracksnew.At(index1)->SetUniqueID(0);
			    new(dimuons[jDimuons++]) AliAODDimuon(tracksnew.At(index0),tracksnew.At(index1));
			}
		    }    
		}
	    }
            if ((handler->NeedsHMPIDBranchReplication()) && (fgAODHmpidRings))
	    {
		TClonesArray* hmpidRings = aod->GetHMPIDrings();
		new (fgAODHmpidRings) TClonesArray(*hmpidRings);
	    }
            
            
            
	    // Additional merging if needed
	    if (merging) {
	      Int_t nc;

		// mcParticles
		TClonesArray* mcparticles = (TClonesArray*) ((aodH->GetEventToMerge())->FindListObject("mcparticles"));
		if( mcparticles ){
		  Int_t npart = mcparticles->GetEntries();
		  nc = fgAODMCParticles->GetEntries();
		  
		  for (Int_t i = 0; i < npart; i++) {
		    AliAODMCParticle* particle = (AliAODMCParticle*) mcparticles->At(i);
		    new((*fgAODMCParticles)[nc++]) AliAODMCParticle(*particle);
		  }
		}

		// tracks
		TClonesArray* tracks = aodH->GetEventToMerge()->GetTracks();
		if(tracks && aodH->GetMergeTracks()){
		  Int_t ntr = tracks->GetEntries();
		  nc  = fgAODTracks->GetEntries();	
		  for (Int_t i = 0; i < ntr; i++) {
		    AliAODTrack*    track = (AliAODTrack*) tracks->At(i);
		    AliAODTrack* newtrack = new((*fgAODTracks)[nc++]) AliAODTrack(*track);
		    newtrack->SetLabel(newtrack->GetLabel() + fgAODMCParticles->GetEntries());
		  }
		  
		  for (Int_t i = 0; i < nc; i++) 
		    {
		      AliAODTrack* track = (AliAODTrack*) fgAODTracks->At(i);
		      track->ResetBit(kIsReferenced);
		      track->SetUniqueID(0);
		    }
		}
		
		// clusters
		TClonesArray* clusters = aodH->GetEventToMerge()->GetCaloClusters();
		if( clusters  && (aodH->GetMergeEMCALClusters() || aodH->GetMergePHOSClusters())) {
		  Int_t ncl  = clusters->GetEntries();
		  nc         =  fgAODCaloClusters->GetEntries();
		  for (Int_t i = 0; i < ncl; i++) {
		    AliAODCaloCluster*    cluster = (AliAODCaloCluster*) clusters->At(i);
		    if(cluster->IsEMCAL() && !aodH->GetMergeEMCALClusters() ) continue;
		    if(cluster->IsPHOS()  && !aodH->GetMergePHOSClusters()  ) continue;   
		    new((*fgAODCaloClusters)[nc++]) AliAODCaloCluster(*cluster);
		  }
		}

		// EMCAL cells
		//*fgAODEmcalCells =  *(aod->GetEMCALCells()); // This will be valid after 10.Mar.2011.
		if(aodH->GetMergeEMCALCells()) 
		  {
		    AliAODCaloCells* copycells = aod->GetEMCALCells();
		    fgAODEmcalCells->CreateContainer(copycells->GetNumberOfCells());
		    nc  = copycells->GetNumberOfCells();
		    
		    while( nc-- ){ fgAODEmcalCells->SetCell(nc,copycells->GetCellNumber(nc),copycells->GetAmplitude(nc),
							    copycells->GetTime(nc),copycells->GetMCLabel(nc),copycells->GetEFraction(nc)); }
		    
		    AliAODCaloCells* cellsA = aodH->GetEventToMerge()->GetEMCALCells();
		    if( cellsA )
		      {
			Int_t ncells  = cellsA->GetNumberOfCells();
			nc = fgAODEmcalCells->GetNumberOfCells();
			
			for (Int_t i  = 0; i < ncells; i++) 
			  {
			    Int_t cn  = cellsA->GetCellNumber(i);
			    Int_t pos = fgAODEmcalCells->GetCellPosition(cn);
			    
			    if (pos >= 0) 
			      {
				Double_t amp = cellsA->GetAmplitude(i) + fgAODEmcalCells->GetAmplitude(pos);
				
				//Check if it is MC, depending on that assing the mc lable, time and e fraction
//				Double_t time    = 0;
				Int_t    mclabel =-1;
			        Double_t efrac   = 0;
				if(cellsA->GetMCLabel(i) >= 0 && fgAODEmcalCells->GetMCLabel(pos) < 0)
				  {
				    mclabel = cellsA->GetMCLabel(i) ;
//				    time    = fgAODEmcalCells->GetTime(pos) ; // Time from data
				    if(amp > 0) efrac = cellsA->GetAmplitude(i) / amp;
				  }
				else if(fgAODEmcalCells->GetMCLabel(pos) >= 0 &&  cellsA->GetMCLabel(i) < 0)
				  {
				    mclabel = fgAODEmcalCells->GetMCLabel(pos) ;
//				    time    = cellsA->GetTime(i) ; // Time from data
				    if(amp > 0) efrac = fgAODEmcalCells->GetAmplitude(pos) / amp;
				  }
				else 
				  { // take all from input
				    mclabel = cellsA->GetMCLabel(i) ;
//				    time    = cellsA->GetTime(i) ; 
				    if(amp > 0) efrac = cellsA->GetAmplitude(i) / amp;  
				  }
				
				fgAODEmcalCells->SetCell(pos, cn, amp,cellsA->GetTime(i),mclabel,efrac);
				
			      } else 
			      {
				AliAODCaloCells* copycells1 = new AliAODCaloCells(*fgAODEmcalCells);
				fgAODEmcalCells->CreateContainer(nc+1);
				Int_t nn = copycells1->GetNumberOfCells();
				
				while( nn-- ){ fgAODEmcalCells->SetCell(nn,copycells1->GetCellNumber(nn),copycells1->GetAmplitude(nn),
									copycells1->GetTime(nn),copycells1->GetMCLabel(nn),copycells1->GetEFraction(nn)); }
				
				fgAODEmcalCells->SetCell(nc++,cn,cellsA->GetAmplitude(i),cellsA->GetTime(i), cellsA->GetMCLabel(i),1.);
				
				delete copycells1;
			      }
			  }
			fgAODEmcalCells->Sort();
		      }
		  } // merge emcal cells
		
		
		// PHOS cells
		//*fgAODPhosCells =  *(aod->GetPHOSCells()); // This will be valid after 10.Mar.2011.
		if(aodH->GetMergePHOSCells()) 
		  {
		    AliAODCaloCells* copycells = aod->GetPHOSCells();
		    fgAODPhosCells->CreateContainer(copycells->GetNumberOfCells());
		    nc  = copycells->GetNumberOfCells();
		    
		    while( nc-- ){ fgAODPhosCells->SetCell(nc,copycells->GetCellNumber(nc),copycells->GetAmplitude(nc),
							   copycells->GetTime(nc),copycells->GetMCLabel(nc),copycells->GetEFraction(nc)); }
		    
		    AliAODCaloCells* cellsP = aodH->GetEventToMerge()->GetPHOSCells();
		    if( cellsP )
		      {
			Int_t ncellsP  = cellsP->GetNumberOfCells();
			nc = fgAODPhosCells->GetNumberOfCells();
			
			for (Int_t i  = 0; i < ncellsP; i++) 
			  {
			    Int_t cn  = cellsP->GetCellNumber(i);
			    Int_t pos = fgAODPhosCells->GetCellPosition(cn);
			    
			    if (pos >= 0) 
			      {
				Double_t amp = cellsP->GetAmplitude(i) + fgAODPhosCells->GetAmplitude(pos);
				
				//Check if it is MC, depending on that assing the mc lable, time and e fraction
//				Double_t time    = 0;
				Int_t    mclabel =-1;
				Double_t    efrac   = 0;
				if(cellsP->GetMCLabel(i) >= 0 && fgAODPhosCells->GetMCLabel(pos) < 0)
				  {
				    mclabel = cellsP->GetMCLabel(i) ;
//				    time    = fgAODPhosCells->GetTime(pos) ; // Time from data
				    if(amp > 0) efrac = cellsP->GetAmplitude(i) / amp;
				  }
				else if(fgAODPhosCells->GetMCLabel(pos) >= 0 &&  cellsP->GetMCLabel(i) < 0)
				  {
				    mclabel = fgAODPhosCells->GetMCLabel(pos) ;
//				    time    = cellsP->GetTime(i) ; // Time from data
				    if(amp > 0) efrac = fgAODPhosCells->GetAmplitude(pos) / amp;
				  }
				else 
				  { // take all from input
				    mclabel = cellsP->GetMCLabel(i) ;
//				    time    = cellsP->GetTime(i) ; 
				    if(amp > 0) efrac = cellsP->GetAmplitude(i) / amp;  
				  }
				
				fgAODPhosCells->SetCell(pos, cn, amp,cellsP->GetTime(i),mclabel,efrac);                
				
			      } else 
			      {
				AliAODCaloCells* copycells1 = new AliAODCaloCells(*fgAODPhosCells);
				fgAODPhosCells->CreateContainer(nc+1);
				Int_t nn = copycells1->GetNumberOfCells();
				
				while( nn-- ){ fgAODPhosCells->SetCell(nn,copycells1->GetCellNumber(nn),copycells1->GetAmplitude(nn), 
								       copycells1->GetTime(nn),copycells1->GetMCLabel(nn),copycells1->GetEFraction(nn)); }
				
				fgAODPhosCells->SetCell(nc++,cn,cellsP->GetAmplitude(i),cellsP->GetTime(i), cellsP->GetMCLabel(i),1.);
				
				delete copycells1;
			      }
			  }
			fgAODPhosCells->Sort();
		      }
		  } // Merge PHOS Cells
		
		if (aodH->GetMergeEMCALTrigger() && aod->GetCaloTrigger("EMCAL")) 
		{
			Int_t   tsEMCAL[48][64], px, py, ts;
			Float_t foEMCAL[48][64], am;
			
			for (Int_t i = 0; i < 48; i++) for (Int_t j = 0; j < 64; j++) 
			{
				tsEMCAL[i][j] = 0;
				foEMCAL[i][j] = 0.;
			}
      
      AliAODCaloTrigger& trg0 = *(aod->GetCaloTrigger("EMCAL"));
      trg0.Reset();
      while (trg0.Next())
      {
        trg0.GetPosition(px, py);
        
        if (px > -1 && py > -1) 
        {
          trg0.GetL1TimeSum(ts);
          if (ts > -1) tsEMCAL[px][py] += ts;
          
          trg0.GetAmplitude(am);
          if (am > -1) foEMCAL[px][py] += am;
        }
      }
      
      AliAODCaloTrigger& trg1 = *((aodH->GetEventToMerge())->GetCaloTrigger("EMCAL"));
      
      trg1.Reset();
      while (trg1.Next())
      {
        trg1.GetPosition(px, py);
        
        if (px > -1 && py > -1) 
        {
          trg1.GetL1TimeSum(ts);
          if (ts > -1) tsEMCAL[px][py] += ts;
          
          trg1.GetAmplitude(am);
          if (am > -1) foEMCAL[px][py] += am;
        }
      }
      
      int nEntries = 0;
      for (Int_t i = 0; i < 48; i++) 
        for (Int_t j = 0; j < 64; j++) 
          if (tsEMCAL[i][j] || foEMCAL[i][j]) nEntries++;
      
      fgAODEMCALTrigger->Allocate(nEntries);
      Int_t timesL0[10]; for (int i = 0; i < 10; i++) timesL0[i] = -1;
      
      for (Int_t i = 0; i < 48; i++) 
        for (Int_t j = 0; j < 64; j++) 
          if (tsEMCAL[i][j] || foEMCAL[i][j]) 
            fgAODEMCALTrigger->Add(i, j, foEMCAL[i][j], -1., timesL0, 0, tsEMCAL[i][j], 0);
    }
		
		if (aodH->GetMergePHOSTrigger()) 
		{
			// To be implemented by PHOS
		}
	    } // merging
	    
	    handler->SetAODIsReplicated();
	}
    }


// Call the user analysis    
	if (isSelected) UserExec(option);
    
// Added protection in case the derived task is not an AOD producer.
    AliAnalysisDataSlot *out0 = GetOutputSlot(0);
    if (out0 && out0->IsConnected()) PostData(0, fTreeA);    

    DisconnectMultiHandler();
}

const char* AliAnalysisTaskSE::CurrentFileName()
{
// Returns the current file name    
    if( fInputHandler )
      return fInputHandler->GetTree()->GetCurrentFile()->GetName();
    else if( fMCEvent )
      return ((AliMCEventHandler*) ((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler()))->TreeK()->GetCurrentFile()->GetName();
    else return "";
}

void AliAnalysisTaskSE::AddAODBranch(const char* cname, void* addobj, const char *fname)
{
    // Add a new branch to the aod tree
    AliAODHandler* handler = dynamic_cast<AliAODHandler*> 
	((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    if (handler) {
	handler->AddBranch(cname, addobj, fname);
    }
}

Bool_t AliAnalysisTaskSE::IsStandardAOD() const
{
// Check if the output AOD handler is configured for standard or delta AOD.
// Users should first check that AODEvent() returns non-null.
    AliAODHandler* handler = dynamic_cast<AliAODHandler*> 
         ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    if (!handler) {
       Error("IsStandardAOD", "No AOD handler. Please use AODEvent() to check this first");
       return kTRUE;
    }
    return handler->IsStandard();   
}

Bool_t AliAnalysisTaskSE::Notify()
{
    return (UserNotify());
}

const AliEventTag *AliAnalysisTaskSE::EventTag() const
{
// Returns tag for the current event, if any. The return value should always be checked by the user.
   if (!fInputHandler) {
      Error("EventTag", "Input handler not yet available. Call this in UserExec");
      return NULL;
   }
   return fInputHandler->GetEventTag();
}

void AliAnalysisTaskSE::LoadBranches() const
{
// Load all branches declared in fBranchNames data member of the parent class.
// Should be called in UserExec.
  if (!fInputHandler) {
     Error("LoadBranches", "Input handler not available yet. Call this in UserExec");
     return;
  }
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (mgr->GetAutoBranchLoading()) return;
  TString taskbranches;
  GetBranches(fInputHandler->GetDataType(), taskbranches);
  if (taskbranches.IsNull()) return;
  TObjArray *arr = taskbranches.Tokenize(",");
  TIter next(arr);
  TObject *obj;
  while ((obj=next())) mgr->LoadBranch(obj->GetName());
  delete arr;
}


//_________________________________________________________________________________________________
void AliAnalysisTaskSE::ConnectMultiHandler()
{
   //
   // Connect MultiHandler
   //
   fInputHandler = (AliInputEventHandler *)((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
   fMultiInputHandler = dynamic_cast<AliMultiInputEventHandler *>(fInputHandler);
   if (fMultiInputHandler) {
      fInputHandler = dynamic_cast<AliInputEventHandler *>(fMultiInputHandler->GetFirstInputEventHandler());
      fMCEventHandler = dynamic_cast<AliInputEventHandler *>(fMultiInputHandler->GetFirstMCEventHandler());
   } else { 
      fMCEventHandler = dynamic_cast<AliInputEventHandler *>((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());
   }
   if (fMCEventHandler) fMCEvent = fMCEventHandler->MCEvent();
}

//_________________________________________________________________________________________________
void AliAnalysisTaskSE::DisconnectMultiHandler()
{
   //
   // Disconnect MultiHandler
   //
   if (fMultiInputHandler) fInputHandler = fMultiInputHandler;
}

/**
 * Perform track selection in case a virtual track selection is provided.
 * @return List of accepted tracks in case a track selection is provided, null otherwise
 */
TObjArray *AliAnalysisTaskSE::GetAcceptedTracks(){
  TObjArray *result = NULL;
  if(fTrackSelection) result = fTrackSelection->GetAcceptedTracks(InputEvent());
  return result;
}
