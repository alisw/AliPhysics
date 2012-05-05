// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE Project            * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  Sedat Altinpinar <Sedat.Altinpinar@cern.ch>           *
//*                  Hege Erdal       <hege.erdal@gmail.com>               *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/// @file   AliAnalysisTaskDxHFEParticleSelection.cxx
/// @author Sedat Altinpinar, Hege Erdal, Matthias Richter
/// @date   2012-03-19
/// @brief  AnalysisTask electron selection for D0 - HFE correlation
///

#include "AliAnalysisTaskDxHFEParticleSelection.h"
#include "AliDxHFEParticleSelection.h"
#include "AliAnalysisManager.h"
#include "AliLog.h"
#include "AliESDInputHandler.h"
#include "TChain.h"
#include "TSystem.h"
#include "TFile.h"

/// ROOT macro for the implementation of ROOT specific class methods
ClassImp(AliAnalysisTaskDxHFEParticleSelection)

AliAnalysisTaskDxHFEParticleSelection::AliAnalysisTaskDxHFEParticleSelection(const char* opt)
  : AliAnalysisTaskSE("AliAnalysisTaskDxHFEParticleSelection")
  , fOutput(0)
  , fOption(opt)
  , fSelector(NULL)
{
  // constructor
  //
  //

  DefineSlots();
}

int AliAnalysisTaskDxHFEParticleSelection::DefineSlots()
{
  // define the data slots
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  return 0;
}

AliAnalysisTaskDxHFEParticleSelection::~AliAnalysisTaskDxHFEParticleSelection()
{
  // destructor
  //
  //

  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor

  if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutput;
    fOutput = 0;
  }

  if (fSelector) {
    fSelector->Clear();
    delete fSelector;
    fSelector=NULL;
  }
}

void AliAnalysisTaskDxHFEParticleSelection::UserCreateOutputObjects()
{
  // create result objects and add to output list

  fOutput = new TList;
  fOutput->SetOwner();
  
  // all tasks must post data once for all outputs
  PostData(1, fOutput);
}

void AliAnalysisTaskDxHFEParticleSelection::UserExec(Option_t* /*option*/)
{
  // process the event

  // TODO: implement correct input, this is likely not to be the
  // ESD
  TObject* pInput=InputEvent();
  if (!pInput) {
    AliError("failed to get input");
    return;
  }
  AliVEvent *pEvent = dynamic_cast<AliVEvent*>(pInput);
  if(!pEvent){
    AliError(Form("input of wrong class type %s, expecting AliVEvent", pInput->ClassName()));
    return;
  }

  PostData(1, fOutput);
}

void AliAnalysisTaskDxHFEParticleSelection::FinishTaskOutput()
{
  // end of the processing
}

void AliAnalysisTaskDxHFEParticleSelection::Terminate(Option_t *)
{
  // last action on the client
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {
    AliFatal("failed to get output container");
    return;
  }

}
