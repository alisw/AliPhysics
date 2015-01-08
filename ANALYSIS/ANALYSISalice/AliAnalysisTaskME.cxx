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

#include "AliAnalysisTaskME.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataSlot.h"
#include "AliAODEvent.h"
#include "AliVEvent.h"
#include "AliAODHandler.h"
#include "AliMultiEventInputHandler.h"
#include "AliLog.h"


ClassImp(AliAnalysisTaskME)

////////////////////////////////////////////////////////////////////////

AliAnalysisTaskME::AliAnalysisTaskME():
    AliAnalysisTask(),
    fDebug(0),
    fEntry(0),
    fFreshBufferOnly(kFALSE),
    fInputHandler(0x0),
    fOutputAOD(0x0),
    fTreeA(0x0),
    fOfflineTriggerMask(0)
{
  // Default constructor
}

AliAnalysisTaskME::AliAnalysisTaskME(const char* name):
    AliAnalysisTask(name, "AnalysisTaskME"),
    fDebug(0),
    fEntry(0),
    fFreshBufferOnly(kFALSE),
    fInputHandler(0x0),
    fOutputAOD(0x0),
    fTreeA(0x0),
    fOfflineTriggerMask(0)
{
  // Default constructor
    DefineInput (0, TChain::Class());
    DefineOutput(0,  TTree::Class());
}

AliAnalysisTaskME::AliAnalysisTaskME(const AliAnalysisTaskME& obj):
    AliAnalysisTask(obj),
    fDebug(0),
    fEntry(0),
    fFreshBufferOnly(kFALSE),
    fInputHandler(0x0),
    fOutputAOD(0x0),
    fTreeA(0x0),
    fOfflineTriggerMask(0)
{
// Copy constructor
    fDebug        = obj.fDebug;
    fEntry        = obj.fEntry;
    fInputHandler = obj.fInputHandler;
    fOutputAOD    = obj.fOutputAOD;
    fTreeA        = obj.fTreeA; 
    fOfflineTriggerMask = obj.fOfflineTriggerMask;
}


AliAnalysisTaskME& AliAnalysisTaskME::operator=(const AliAnalysisTaskME& other)
{
// Assignment
    if (this != &other) {
	AliAnalysisTask::operator=(other);
	fDebug           = other.fDebug;
	fEntry           = other.fEntry;
	fFreshBufferOnly = other.fFreshBufferOnly;
	fInputHandler    = other.fInputHandler;
	fOutputAOD       = other.fOutputAOD;
	fTreeA           = other.fTreeA;    
	fOfflineTriggerMask = other.fOfflineTriggerMask;
    }
    return *this;
}


void AliAnalysisTaskME::ConnectInputData(Option_t* /*option*/)
{
// Connect the input data
    if (fDebug > 1) printf("AnalysisTaskME::ConnectInputData() \n");
//
//  Multi AOD
//
    fInputHandler = dynamic_cast<AliMultiEventInputHandler*> 
	((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
    if (fInputHandler == 0) {
	AliFatal("Event Handler has to be MultiEventInputHandler !");
    } else {
	// Check that we have an event pool
	if (!fInputHandler->GetEventPool()) {
	    fInputHandler->SetEventPool(AliAnalysisManager::GetAnalysisManager()->GetEventPool());
	    if (!fInputHandler->GetEventPool()) 
		AliFatal("MultiEventInputHandler has no EventPool connected !");
	}
    }
}

void AliAnalysisTaskME::CreateOutputObjects()
{
// Create the output container
//
//  Default AOD
    if (fDebug > 1) printf("AnalysisTaskME::CreateOutPutData() \n");

    AliAODHandler* handler = (AliAODHandler*) 
	((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    
    if (handler) {
	fOutputAOD   = handler->GetAOD();
	fTreeA = handler->GetTree();
    } else {
	AliWarning("No AOD Event Handler connected.") ; 
    }
    UserCreateOutputObjects();
}

void AliAnalysisTaskME::Exec(Option_t* option)
{
//
// Exec analysis of one event

    if (fDebug > 1) AliInfo("AliAnalysisTaskME::Exec() \n");
    if( fInputHandler ) 
       fEntry = fInputHandler->GetReadEntry();
    if ( !((Entry()-1)%100) && fDebug > 0) 
         AliInfo(Form("%s ----> Processing event # %lld", CurrentFileName(), Entry()));

    AliAODHandler* outputHandler = (AliAODHandler*) 
	((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());         
//
// Was event selected ? If no event selection mechanism, the event SHOULD be selected (AG)
    UInt_t isSelected = AliVEvent::kAny;
    if(fInputHandler && fInputHandler->GetEventSelection()) {
      // Get the actual offline trigger mask for the event and AND it with the
      // requested mask. If no mask requested select by default the event.
      if (fOfflineTriggerMask)
	isSelected = fOfflineTriggerMask & fInputHandler->IsEventSelected();
    }
    
    if (!isSelected) { 
	if (fDebug > 1) AliInfo("Event rejected \n");
	fInputHandler->EventSkipped();
	return;
    }
// Call the user analysis    
    
    if (fInputHandler && fInputHandler->IsBufferReady()) {
	if ((fFreshBufferOnly && fInputHandler->IsFreshBuffer()) || !fFreshBufferOnly)
	{
	    if (outputHandler) outputHandler->SetFillAOD(kTRUE);
	    UserExec(option);
	    // Added protection in case the derived task is not an AOD producer.
	    AliAnalysisDataSlot *out0 = GetOutputSlot(0);
	    if (out0 && out0->IsConnected()) PostData(0, fTreeA);
	} else {
	    if (outputHandler) outputHandler->SetFillAOD(kFALSE);
	}
    } else {
	AliInfo(Form("Waiting for buffer to be ready !\n"));
    }
}

const char* AliAnalysisTaskME::CurrentFileName()
{
// Returns the current file name    
    if(fInputHandler )
	return fInputHandler->GetTree()->GetCurrentFile()->GetName();
    else return "";
}

void AliAnalysisTaskME::AddAODBranch(const char* cname, void* addobj, const char *fname)
{
    // Add a new branch to the aod tree
    AliAODHandler* handler = (AliAODHandler*) 
	((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    if (handler) {
	handler->AddBranch(cname, addobj, fname);
    }
}

AliVEvent*  AliAnalysisTaskME::GetEvent(Int_t iev)
{
    // Get an event from the input handler
    return (fInputHandler->GetEvent(iev));
}

