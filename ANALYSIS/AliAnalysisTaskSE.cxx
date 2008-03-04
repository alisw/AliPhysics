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


ClassImp(AliAnalysisTaskSE)

////////////////////////////////////////////////////////////////////////

AliAnalysisTaskSE::AliAnalysisTaskSE():
    AliAnalysisTask(),
    fDebug(0),
    fInputEvent(0x0),
    fOutputAOD(0x0),
    fMCEvent(0x0),
    fTreeA(0x0)
{
  // Default constructor
}

AliAnalysisTaskSE::AliAnalysisTaskSE(const char* name):
    AliAnalysisTask(name, "AnalysisTaskSE"),
    fDebug(0),
    fInputEvent(0x0),
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
    fInputEvent(0x0),
    fOutputAOD(0x0),
    fMCEvent(0x0),
    fTreeA(0x0)
{
// Copy constructor
    fDebug      = obj.fDebug;
    fInputEvent = obj.fInputEvent;
    fOutputAOD  = obj.fOutputAOD;
    fMCEvent    = obj.fMCEvent;
    fTreeA      = obj.fTreeA;    
}


AliAnalysisTaskSE& AliAnalysisTaskSE::operator=(const AliAnalysisTaskSE& other)
{
// Assignment
    AliAnalysisTask::operator=(other);
    fDebug      = other.fDebug;
    fInputEvent = other.fInputEvent;
    fOutputAOD  = other.fOutputAOD;
    fMCEvent    = other.fMCEvent;
    fTreeA      = other.fTreeA;    
    return *this;
}


void AliAnalysisTaskSE::ConnectInputData(Option_t* /*option*/)
{
// Connect the input data
    if (fDebug > 1) printf("AnalysisTaskSE::ConnectInputData() \n");
//
//  ESD
//
    AliInputEventHandler* handler = (AliInputEventHandler*) 
	((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
    if (handler) fInputEvent = handler->GetEvent();
//
//  Monte Carlo
//
    AliMCEventHandler*    mcH = 0;
    mcH = (AliMCEventHandler*) ((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());
    if (mcH) fMCEvent = mcH->MCEvent();
}

void AliAnalysisTaskSE::CreateOutputObjects()
{
// Create the output container
//
//  Default AOD
    if (fDebug > 1) printf("AnalysisTaskSE::CreateOutPutData() \n");

    AliAODHandler* handler = (AliAODHandler*) 
	((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    
    fOutputAOD   = handler->GetAOD();
    fTreeA = handler->GetTree();
    UserCreateOutputObjects();
}

void AliAnalysisTaskSE::Exec(Option_t* option)
{
//
// Exec analysis of one event
    PostData(0, fTreeA);
    UserExec(option);
}



