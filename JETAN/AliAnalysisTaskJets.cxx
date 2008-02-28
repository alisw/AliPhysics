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

#include "AliAnalysisTaskJets.h"
#include "AliAnalysisManager.h"
#include "AliJetFinder.h"
#include "AliJetHistos.h"
#include "AliESDEvent.h"
#include "AliESD.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliMCEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"


ClassImp(AliAnalysisTaskJets)

////////////////////////////////////////////////////////////////////////

AliAnalysisTaskJets::AliAnalysisTaskJets():
    AliAnalysisTaskSE(),
    fDebug(0),
    fJetFinder(0x0),
    fTree(0x0),
    fHistos(0x0),
    fListOfHistos(0x0)
{
  // Default constructor
}

AliAnalysisTaskJets::AliAnalysisTaskJets(const char* name):
    AliAnalysisTaskSE(name),
    fDebug(0),
    fJetFinder(0x0),
    fTree(0x0),
    fHistos(0x0),
    fListOfHistos(0x0)
{
  // Default constructor
    DefineOutput(1, TList::Class());
}

void AliAnalysisTaskJets::UserCreateOutputObjects()
{
// Create the output container
//
    if (fDebug > 1) printf("AnalysisTaskJets::CreateOutPutData() \n");
//  Connec default AOD to jet finder

    fJetFinder->ConnectAOD(AODEvent());
//
//  Histograms
    OpenFile(1);
    fListOfHistos = new TList();
    fHistos       = new AliJetHistos();
    fHistos->AddHistosToList(fListOfHistos);
}

void AliAnalysisTaskJets::Init()
{
    // Initialization
    if (fDebug > 1) printf("AnalysisTaskJets::Init() \n");

    // Call configuration file
    gROOT->LoadMacro("ConfigJetAnalysis.C");
    fJetFinder = (AliJetFinder*) gInterpreter->ProcessLine("ConfigJetAnalysis()");
    // Initialise Jet Analysis
    fJetFinder->Init();
    // Write header information to local file
    fJetFinder->WriteHeaders();
}

    



void AliAnalysisTaskJets::UserExec(Option_t */*option*/)
{
// Execute analysis for current event
//
    
    if (MCEvent()) {
	AliStack* stack = MCEvent()->Stack();
    }
    fJetFinder->GetReader()->SetInputEvent(InputEvent(), AODEvent(), MCEvent());
    fJetFinder->ProcessEvent();

    // Fill control histos
    fHistos->FillHistos(AODEvent()->GetJets());
    // Post the data
    PostData(1, fListOfHistos);
}

void AliAnalysisTaskJets::Terminate(Option_t */*option*/)
{
// Terminate analysis
//
    if (fDebug > 1) printf("AnalysisJets: Terminate() \n");
//    if (fJetFinder) fJetFinder->FinishRun();
}

