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
 
#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>

#include "AliAnalysisTaskJets.h"
#include "AliJetFinder.h"
#include "AliESD.h"

ClassImp(AliAnalysisTaskJets)

////////////////////////////////////////////////////////////////////////

AliAnalysisTaskJets::AliAnalysisTaskJets():
    fDebug(0),
    fJetFinder(0x0),
    fChain(0x0),
    fESD(0x0),
    fTreeJ(0x0)
{
  // Default constructor
}

AliAnalysisTaskJets::AliAnalysisTaskJets(const char* name):
    AliAnalysisTask(name, "AnalysisTaskJets"),
    fDebug(0),
    fJetFinder(0x0),
    fChain(0x0),
    fESD(0x0),
    fTreeJ(0x0)
{
  // Default constructor
    DefineInput (0, TChain::Class());
    DefineOutput(0, TTree::Class());
}

void AliAnalysisTaskJets::CreateOutputObjects()
{
// Create the output container
    OpenFile(0);
    fTreeJ = fJetFinder->MakeTreeJ("TreeJ");
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
//
    if (fDebug > 1) printf("AnalysisTaskJets::ConnectInputData() \n");
    fChain = (TChain*)GetInputData(0);

    char ** address = (char **)GetBranchAddress(0, "ESD");
    if (address)     {
	
// Branch has been already connected
	fESD = (AliESD*)(*address);
    }
    else     {
// First task taking the branch enables it
	fESD = new AliESD();
	SetBranchAddress(0, "ESD", &fESD);
    }
    
    fJetFinder->ConnectTree(fChain, fESD);
}

void AliAnalysisTaskJets::Exec(Option_t */*option*/)
{
// Execute analysis for current event
//
    Long64_t ientry = fChain->GetReadEntry();
    if (fDebug > 1) printf("Analysing event # %5d\n", (Int_t) ientry);
    fJetFinder->ProcessEvent(ientry);
    PostData(0, fTreeJ);
}

void AliAnalysisTaskJets::Terminate(Option_t */*option*/)
{
// Terminate analysis
//
    if (fDebug > 1) printf("AnalysisJets: Terminate() \n");
    // if (fJetFinder) fJetFinder->FinishRun();
}

