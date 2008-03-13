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
#include "AliESDEvent.h"
#include "AliESD.h"
#include "AliAODEvent.h"
#include "AliVEvent.h"
#include "AliAODHandler.h"
#include "AliMCEventHandler.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliLog.h"


ClassImp(AliAnalysisTaskSE)

////////////////////////////////////////////////////////////////////////

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
    
    if (handler) {
	fOutputAOD   = handler->GetAOD();
	fTreeA = handler->GetTree();
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
         
// Call the user analysis    
    UserExec(option);
    PostData(0, fTreeA);
    
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

