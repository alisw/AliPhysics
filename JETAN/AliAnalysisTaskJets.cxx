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
#include "AliJetHeader.h"
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
  fConfigFile("ConfigJetAnalysis.C"),  
  fNonStdBranch(""),  
  fJetFinder(0x0),
  fHistos(0x0),
  fListOfHistos(0x0)
{
  // Default constructor
}

AliAnalysisTaskJets::AliAnalysisTaskJets(const char* name):
    AliAnalysisTaskSE(name),
    fConfigFile("ConfigJetAnalysis.C"),  
    fNonStdBranch(""),  
    fJetFinder(0x0),
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

    

    if(fNonStdBranch.Length()==0){
      //  Connec default AOD to jet finder
      fJetFinder->ConnectAOD(AODEvent());
    }
    else{
      // Create a new branch for jets...
      // how is this is reset cleared in the UserExec....
      // Can this be handled by the framework?
      TClonesArray *tca = new TClonesArray("AliAODJet", 0);
      tca->SetName(fNonStdBranch);
      AddAODBranch("TClonesArray",&tca);
      fJetFinder->ConnectAODNonStd(AODEvent(),fNonStdBranch.Data());
    }
    //  Histograms
    OpenFile(1);
    fListOfHistos = new TList();
    fHistos       = new AliJetHistos();
    fHistos->AddHistosToList(fListOfHistos);
    
    // Add the JetFinderInforamtion to the Outputlist
    AliJetHeader *fH = fJetFinder->GetHeader();
    // Compose a characteristic output name
    // with the name of the output branch
    if(fH){
      if(fNonStdBranch.Length()==0){
	fH->SetName("AliJetHeader_jets");
      }
      else{
	fH->SetName(Form("AliJetHeader_%s",fNonStdBranch.Data()));
      }
    }
    OutputTree()->GetUserInfo()->Add(fH);
}

void AliAnalysisTaskJets::Init()
{
    // Initialization
    if (fDebug > 1) printf("AnalysisTaskJets::Init() \n");

    // Call configuration file
    gROOT->LoadMacro(fConfigFile);
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


  // Fill control histos
  TClonesArray* jarray = 0;
  if(fNonStdBranch.Length()==0){
    jarray = AODEvent()->GetJets();
  }
  else{
    jarray =  dynamic_cast<TClonesArray*>(AODEvent()->FindListObject(fNonStdBranch.Data()));
    jarray->Delete();    // this is our responsibility, clear before filling again
  }

  fJetFinder->GetReader()->SetInputEvent(InputEvent(), AODEvent(), MCEvent());
  fJetFinder->ProcessEvent();

  fHistos->FillHistos(jarray);
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

