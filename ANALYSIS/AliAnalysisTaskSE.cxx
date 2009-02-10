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
#include "AliAnalysisDataSlot.h"
#include "AliESDEvent.h"
#include "AliESD.h"
#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliAODTracklets.h"
#include "AliVEvent.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliLog.h"


ClassImp(AliAnalysisTaskSE)

////////////////////////////////////////////////////////////////////////
AliAODHeader*    AliAnalysisTaskSE::fgAODHeader        = NULL;
TClonesArray*    AliAnalysisTaskSE::fgAODTracks        = NULL;
TClonesArray*    AliAnalysisTaskSE::fgAODVertices      = NULL;
TClonesArray*    AliAnalysisTaskSE::fgAODV0s           = NULL;
TClonesArray*    AliAnalysisTaskSE::fgAODPMDClusters   = NULL;
TClonesArray*    AliAnalysisTaskSE::fgAODJets          = NULL;
TClonesArray*    AliAnalysisTaskSE::fgAODFMDClusters   = NULL;
TClonesArray*    AliAnalysisTaskSE::fgAODCaloClusters  = NULL;
AliAODTracklets* AliAnalysisTaskSE::fgAODTracklets     = NULL;


AliAnalysisTaskSE::AliAnalysisTaskSE():
    AliAnalysisTask(),
    fDebug(0),
    fEntry(0),
    fInputEvent(0x0),
    fInputHandler(0x0),
    fOutputAOD(0x0),
    fMCEvent(0x0),
    fTreeA(0x0)
{
  // Default constructor
}

AliAnalysisTaskSE::AliAnalysisTaskSE(const char* name):
    AliAnalysisTask(name, "AnalysisTaskSE"),
    fDebug(0),
    fEntry(0),
    fInputEvent(0x0),
    fInputHandler(0x0),
    fOutputAOD(0x0),
    fMCEvent(0x0),
    fTreeA(0x0)
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
    fInputHandler(0x0),
    fOutputAOD(0x0),
    fMCEvent(0x0),
    fTreeA(0x0)
{
// Copy constructor
    fDebug        = obj.fDebug;
    fEntry        = obj.fEntry;
    fInputEvent   = obj.fInputEvent;
    fInputHandler = obj.fInputHandler;
    fOutputAOD    = obj.fOutputAOD;
    fMCEvent      = obj.fMCEvent;
    fTreeA        = obj.fTreeA;    
    printf("Constructor (3) \n");
}


AliAnalysisTaskSE& AliAnalysisTaskSE::operator=(const AliAnalysisTaskSE& other)
{
// Assignment
    AliAnalysisTask::operator=(other);
    fDebug        = other.fDebug;
    fEntry        = other.fEntry;
    fInputEvent   = other.fInputEvent;
    fInputHandler = other.fInputHandler;
    fOutputAOD    = other.fOutputAOD;
    fMCEvent      = other.fMCEvent;
    fTreeA        = other.fTreeA;    
    return *this;
}


void AliAnalysisTaskSE::ConnectInputData(Option_t* /*option*/)
{
// Connect the input data
    if (fDebug > 1) printf("AnalysisTaskSE::ConnectInputData() \n");
//
//  ESD
//
    fInputHandler = (AliInputEventHandler*) 
         ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
//
//  Monte Carlo
//
    AliMCEventHandler*    mcH = 0;
    mcH = (AliMCEventHandler*) ((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());
    if (mcH) fMCEvent = mcH->MCEvent();
    
    
    if (fInputHandler) {
         fInputEvent = fInputHandler->GetEvent();
    } else if( fMCEvent ) {
         AliWarning("No Input Event Handler connected, only MC Truth Event Handler") ; 
    } else {
         AliError("No Input Event Handler connected") ; 
         return ; 
    }
}

void AliAnalysisTaskSE::CreateOutputObjects()
{
// Create the output container
//
//  Default AOD
    if (fDebug > 1) printf("AnalysisTaskSE::CreateOutPutData() \n");

    AliAODHandler* handler = (AliAODHandler*) 
         ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    
    // Check if AOD replication has been required
    
    if (handler) {
	fOutputAOD   = handler->GetAOD();
	fTreeA = handler->GetTree();
	if (!(handler->IsStandard())) {
	    if ((handler->NeedsHeaderReplication()) && !(fgAODHeader)) 
		{
		 if (fDebug > 1) AliInfo("Replicating header");
		 fgAODHeader = new AliAODHeader;
		 handler->AddBranch("AliAODHeader", &fgAODHeader);
		}
	    if ((handler->NeedsTracksBranchReplication()) && !(fgAODTracks))      
	    {   
		if (fDebug > 1) AliInfo("Replicating track branch\n");
		fgAODTracks = new TClonesArray("AliAODTrack",500);
		fgAODTracks->SetName("tracks");
		handler->AddBranch("TClonesArray", &fgAODTracks);
	    }    
	    if ((handler->NeedsVerticesBranchReplication()) && !(fgAODVertices))
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
	    if ((handler->NeedsJetsBranchReplication()) && !(fgAODJets))	  
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
	    if ((handler->NeedsCaloClustersBranchReplication()) && !(fgAODCaloClusters))	  
	    {   
		if (fDebug > 1) AliInfo("Replicating CaloClusters branch\n");
		fgAODCaloClusters = new TClonesArray("AliAODCaloCluster",500);
		fgAODCaloClusters->SetName("caloClusters");
		handler->AddBranch("TClonesArray", &fgAODCaloClusters);
	    }
	}
    } else {
	AliWarning("No AOD Event Handler connected.") ; 
    }
    UserCreateOutputObjects();
}

void AliAnalysisTaskSE::Exec(Option_t* option)
{
//
// Exec analysis of one event
    if (fDebug > 1) AliInfo("AliAnalysisTaskSE::Exec() \n");
    if( fInputHandler ) 
       fEntry = fInputHandler->GetReadEntry();
    else if( fMCEvent )
       fEntry = fMCEvent->Header()->GetEvent(); 
    if ( !((Entry()-1)%100) && fDebug > 0) 
         AliInfo(Form("%s ----> Processing event # %lld", CurrentFileName(), Entry()));

    AliAODHandler* handler = (AliAODHandler*) 
	((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    AliAODInputHandler* aodH = dynamic_cast<AliAODInputHandler*>(fInputHandler);

    if (handler && aodH) {
	if (!(handler->IsStandard()) && !(handler->AODIsReplicated())) {
	    if ((handler->NeedsHeaderReplication()) && (fgAODHeader))
	    {
	      // copy the contents by assigment
	      *fgAODHeader =  *(dynamic_cast<AliAODHeader*>(InputEvent()->GetHeader()));
	    }
	    if ((handler->NeedsTracksBranchReplication()) && (fgAODTracks))
	    {
		TClonesArray* tracks = (dynamic_cast<AliAODEvent*>(InputEvent()))->GetTracks();
		new (fgAODTracks) TClonesArray(*tracks);
	    }
	    if ((handler->NeedsVerticesBranchReplication()) && (fgAODVertices))
	    {
		TClonesArray* vertices = (dynamic_cast<AliAODEvent*>(InputEvent()))->GetVertices();
		new (fgAODVertices) TClonesArray(*vertices);
	    }
	    if ((handler->NeedsV0sBranchReplication()) && (fgAODV0s))
	    {
		TClonesArray* V0s = (dynamic_cast<AliAODEvent*>(InputEvent()))->GetV0s();
		new (fgAODV0s) TClonesArray(*V0s);
	    }
	    if ((handler->NeedsTrackletsBranchReplication()) && (fgAODTracklets))
	    {
	      *fgAODTracklets = *(dynamic_cast<AliAODEvent*>(InputEvent()))->GetTracklets();
	    }
	    if ((handler->NeedsPMDClustersBranchReplication()) && (fgAODPMDClusters))
	    {
		TClonesArray* PMDClusters = (dynamic_cast<AliAODEvent*>(InputEvent()))->GetPmdClusters();
		new (fgAODPMDClusters) TClonesArray(*PMDClusters);
	    }
	    if ((handler->NeedsJetsBranchReplication()) && (fgAODJets))
	    {
		TClonesArray* Jets = (dynamic_cast<AliAODEvent*>(InputEvent()))->GetJets();
		new (fgAODJets) TClonesArray(*Jets);
	    }
	    if ((handler->NeedsFMDClustersBranchReplication()) && (fgAODFMDClusters))
	    {
		TClonesArray* FMDClusters = (dynamic_cast<AliAODEvent*>(InputEvent()))->GetFmdClusters();
		new (fgAODFMDClusters) TClonesArray(*FMDClusters);
	    }
	    if ((handler->NeedsCaloClustersBranchReplication()) && (fgAODCaloClusters))
	    {
		TClonesArray* CaloClusters = (dynamic_cast<AliAODEvent*>(InputEvent()))->GetCaloClusters();
		new (fgAODCaloClusters) TClonesArray(*CaloClusters);
	    }
	    //
	    handler->SetAODIsReplicated();
	    
	}
    }

// Call the user analysis    
    UserExec(option);
    // Added protection in case the derived task is not an AOD producer.
    AliAnalysisDataSlot *out0 = GetOutputSlot(0);
    if (out0 && out0->IsConnected()) PostData(0, fTreeA);    
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

void AliAnalysisTaskSE::AddAODBranch(const char* cname, void* addobj)
{
    // Add a new branch to the aod tree
    AliAODHandler* handler = (AliAODHandler*) 
	((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    if (handler) {
	handler->AddBranch(cname, addobj);
    }
}
