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

/* $Id: AliAnalysisTaskJets.cxx 24457 2008-03-13 10:02:22Z morsch $ */
 
#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TChain.h>
#include <TFile.h>
#include <TList.h>

#include "AliAnalysisManager.h"
#include "AliCdfJetFinder.h"
#include "AliESDEvent.h"
#include "AliESD.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliMCEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"

#include "AliAnalysisTaskJetsCDF.h"

ClassImp(AliAnalysisTaskJetsCDF)

//______________________________________________________________________________
AliAnalysisTaskJetsCDF::AliAnalysisTaskJetsCDF():
    AliAnalysisTaskSE(),
    fJetFinder(0x0),
    fListOfHistos(0x0)
{
  // Default constructor
}

//______________________________________________________________________________
AliAnalysisTaskJetsCDF::AliAnalysisTaskJetsCDF(const char* name):
    AliAnalysisTaskSE(name),
    fJetFinder(0x0),
    fListOfHistos(0x0)
{
  // Default constructor
    DefineOutput(1, TList::Class());
}

//______________________________________________________________________________
void AliAnalysisTaskJetsCDF::UserCreateOutputObjects()
{
// Create the output container
//
    if (fDebug > 1) printf("AnalysisTaskJetsCDF::CreateOutPutData() \n");
//  Connect default AOD to jet finder

    fJetFinder->ConnectAOD(AODEvent());
//
//  Histograms
    OpenFile(1);
    fListOfHistos = new TList();
    if (!fJetFinder) {
       Error("UserCreateOutputObjects", "No jet finder connected. Aborting.");
       return;
    }
    fJetFinder->CreateOutputObjects(fListOfHistos);   
}

//______________________________________________________________________________
void AliAnalysisTaskJetsCDF::LocalInit()
{
    // Initialization
    if (fDebug > 1) printf("AnalysisTaskJets::Init() \n");

    // Call configuration file
    gROOT->LoadMacro("ConfigJetAnalysisCDF.C");
    fJetFinder = (AliCdfJetFinder*) gInterpreter->ProcessLine("ConfigJetAnalysis()");
    // Initialise Jet Analysis
    fJetFinder->Init();
    // Write header information to local file
    fJetFinder->WriteHeaders();
}

//______________________________________________________________________________
void AliAnalysisTaskJetsCDF::UserExec(Option_t */*option*/)
{
// Execute analysis for current event
//
    
    fJetFinder->GetReader()->SetInputEvent(InputEvent(), AODEvent(), MCEvent());
    fJetFinder->ProcessEvent();

    // Post the data
    PostData(1, fListOfHistos);
}

//______________________________________________________________________________
void AliAnalysisTaskJetsCDF::Terminate(Option_t */*option*/)
{
// Terminate analysis
//
    if (fDebug > 1) printf("AliAnalysisJetsCDF: Terminate() \n");
//    if (fJetFinder) fJetFinder->FinishRun();
}

