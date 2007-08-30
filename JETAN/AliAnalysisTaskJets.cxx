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
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"


ClassImp(AliAnalysisTaskJets)

////////////////////////////////////////////////////////////////////////

AliAnalysisTaskJets::AliAnalysisTaskJets():
    fDebug(0),
    fJetFinder(0x0),
    fChain(0x0),
    fESD(0x0),
    fAOD(0x0),
    fTreeA(0x0),
    fHistos(0x0),
    fListOfHistos(0x0)
{
  // Default constructor
}

AliAnalysisTaskJets::AliAnalysisTaskJets(const char* name):
    AliAnalysisTask(name, "AnalysisTaskJets"),
    fDebug(0),
    fJetFinder(0x0),
    fChain(0x0),
    fESD(0x0),
    fAOD(0x0),
    fTreeA(0x0),
    fHistos(0x0),
    fListOfHistos(0x0)
{
  // Default constructor
    DefineInput (0, TChain::Class());
    DefineOutput(0, TTree::Class());
    DefineOutput(1, TList::Class());
}

void AliAnalysisTaskJets::CreateOutputObjects()
{
// Create the output container
//
//  Default AOD
    if (fDebug > 1) printf("AnalysisTaskJets::CreateOutPutData() \n");
    AliAODHandler* handler = (AliAODHandler*) ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    
    fAOD   = handler->GetAOD();
    fTreeA = handler->GetTree();
    fJetFinder->ConnectAOD(fAOD);
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

void AliAnalysisTaskJets::ConnectInputData(Option_t */*option*/)
{
// Connect the input data
    if (fDebug > 1) printf("AnalysisTaskJets::ConnectInputData() \n");
    fChain = (TChain*)GetInputData(0);
    fESD = new AliESDEvent();
    fESD->ReadFromTree(fChain);

    AliMCEventHandler*    mcTruth = (AliMCEventHandler*) 
	((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());

    fJetFinder->GetReader()->SetInputEvent(fESD, fAOD, mcTruth);
}

void AliAnalysisTaskJets::Exec(Option_t */*option*/)
{
// Execute analysis for current event
//
    AliMCEventHandler*    mctruth = (AliMCEventHandler*) 
	((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());
    if (mctruth) {
	AliStack* stack = mctruth->Stack();
	printf("AliAnalysisTaskJets: Number of tracks on stack %5d\n", stack->GetNtrack());
    }
    
    AliESD* old = fESD->GetAliESDOld();
    if (old) {
	fESD->CopyFromOldESD();
    }
    
    Long64_t ientry = fChain->GetReadEntry();
    if (fDebug > 1) printf("Analysing event # %5d\n", (Int_t) ientry);
    fJetFinder->ProcessEvent(ientry);

    // Fill control histos
    fHistos->FillHistos(fAOD->GetJets());
    
    // Post the data
    PostData(0, fTreeA);
    PostData(1, fListOfHistos);
}

void AliAnalysisTaskJets::Terminate(Option_t */*option*/)
{
// Terminate analysis
//
    if (fDebug > 1) printf("AnalysisJets: Terminate() \n");
    // if (fJetFinder) fJetFinder->FinishRun();
}

