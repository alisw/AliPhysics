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

//----------------------------------------------------------------
// Analysis task for interfacing the jet reader with the analysis framework
//
// Author: magali.estienne@subatech.in2p3.fr 
//         alexandre.shabetai@cern.ch
//----------------------------------------------------------------

#include <Riostream.h> 
#include <TROOT.h>
#include <TInterpreter.h>
#include <TTree.h>

#include "AliAnalysisTaskJetsReader.h"
#include "AliAnalysisManager.h"
#include "AliJetReader.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliJetCalTrk.h"

ClassImp(AliAnalysisTaskJetsReader)

////////////////////////////////////////////////////////////////////////

AliAnalysisTaskJetsReader::AliAnalysisTaskJetsReader():
  AliAnalysisTaskSE(),
  fConfigFile("ConfigJetReaderAnalysis.C"),
  fJetReader(0x0),
  fReadAODFromOutput(0),
  fReaderEvent(0x0),
  fExchangeTree(0x0)
{
  // Default constructor
}

//----------------------------------------------------------------
AliAnalysisTaskJetsReader::AliAnalysisTaskJetsReader(const char* name):
  AliAnalysisTaskSE(name),
  fConfigFile("ConfigJetReaderAnalysis.C"),
  fJetReader(0x0),
  fReadAODFromOutput(0),
  fReaderEvent(0x0),
  fExchangeTree(0x0)
{
  // Default constructor
  DefineOutput(1, TTree::Class());

}

//----------------------------------------------------------------
AliAnalysisTaskJetsReader::~AliAnalysisTaskJetsReader()
{
  // Destructor

  if (fExchangeTree)
   delete fExchangeTree;

}

//----------------------------------------------------------------
void AliAnalysisTaskJetsReader::UserCreateOutputObjects()
{
  // Create the TTree to be exchanged between the reader and finder
  //
   
  if (fDebug > 1) printf("AnalysisTaskJetsReader::CreateOutPutData() \n");

  fExchangeTree = new TTree("jets_ExchangeContainer","ExchangeTree");

  fJetReader->InitTasks();

  fReaderEvent = fJetReader->GetCalTrkEvent();

  fExchangeTree->Branch("AliJetCalTrkEvent", &fReaderEvent);

  PostData(1, fExchangeTree);

}

//----------------------------------------------------------------
void AliAnalysisTaskJetsReader::Init()
{
  // Initialization
  if (fDebug > 1) printf("AnalysisTaskJets::Init() \n");
 
  // Call configuration file
  if (fConfigFile.Length()) {
    gROOT->LoadMacro(fConfigFile);
    fJetReader = (AliJetReader*) gInterpreter->ProcessLine("ConfigJetReaderAnalysis()");
  }

}

//----------------------------------------------------------------
void AliAnalysisTaskJetsReader::UserExec(Option_t */*option*/)
{

  // Execute analysis for current event
  //

  // Clear current CalTrkEvent 
  fReaderEvent->Clear();

  fExchangeTree->Reset();

  // Give InputEvent to the reader
  if (dynamic_cast<AliAODEvent*>(InputEvent()) !=  0 && !fReadAODFromOutput) {
    // AOD is input event..........................................V                                       
    fJetReader->SetInputEvent(InputEvent(), InputEvent(), MCEvent());
  } else {
    // AOD is read from output ....................................V     
    fJetReader->SetInputEvent(InputEvent(), AODEvent(), MCEvent());
  }

  // Process current event
  fJetReader->ProcessEvent();

  // Fill object to be exchanged between reader and finder tasks
  fExchangeTree->Fill();

  // Post the data
  PostData(1, fExchangeTree);

  return;

}

//----------------------------------------------------------------
void AliAnalysisTaskJetsReader::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  if (fDebug > 1) printf("AnalysisJets: Terminate() \n");

}

