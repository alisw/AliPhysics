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

/// @file   AliAnalysisTaskDxHFECorrelation.cxx
/// @author Sedat Altinpinar, Hege Erdal, Matthias Richter
/// @date   2012-03-19
/// @brief  AnalysisTask D0 - HFE correlation
///

#include "AliAnalysisTaskDxHFECorrelation.h"
#include "AliDxHFECorrelation.h"
#include "AliDxHFEParticleSelectionD0.h"
#include "AliDxHFEParticleSelectionEl.h"
#include "AliAnalysisManager.h"
#include "AliLog.h"
#include "AliESDInputHandler.h"
#include "TChain.h"
#include "TSystem.h"
#include "TFile.h"
#include <memory>

/// ROOT macro for the implementation of ROOT specific class methods
ClassImp(AliAnalysisTaskDxHFECorrelation)

AliAnalysisTaskDxHFECorrelation::AliAnalysisTaskDxHFECorrelation(const char* opt)
  : AliAnalysisTaskSE("AliAnalysisTaskDxHFECorrelation")
  , fOutput(0)
  , fOption(opt)
  , fCorrelation(NULL)
  , fD0s(NULL)
  , fElectrons(NULL)
{
  // constructor
  //
  //

  DefineSlots();
}

int AliAnalysisTaskDxHFECorrelation::DefineSlots()
{
  // define the data slots
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  return 0;
}

AliAnalysisTaskDxHFECorrelation::~AliAnalysisTaskDxHFECorrelation()
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
  if (fD0s) delete fD0s;
  fD0s=NULL;
  if (fElectrons) delete fElectrons;
  fElectrons=NULL;
  if (fCorrelation) delete fCorrelation;
  fCorrelation=NULL;
}

void AliAnalysisTaskDxHFECorrelation::UserCreateOutputObjects()
{
  // create result objects and add to output list

  std::auto_ptr<TList> Output(new TList);
  std::auto_ptr<AliDxHFEParticleSelection> D0s(new AliDxHFEParticleSelectionD0);
  std::auto_ptr<AliDxHFEParticleSelection> Electrons(new AliDxHFEParticleSelectionEl);
  std::auto_ptr<AliDxHFECorrelation> Correlation(new AliDxHFECorrelation);

  if (!Output     .get() ||
      !D0s        .get() ||
      !Electrons  .get() ||
      !Correlation.get()) {
    AliFatal("allocation of worker classes failed");
    return;
  }

  fOutput      = Output     .release();
  fD0s         = D0s        .release();
  fElectrons   = Electrons  .release();
  fCorrelation = Correlation.release();

  fOutput->SetOwner();

  // all tasks must post data once for all outputs
  PostData(1, fOutput);
}

void AliAnalysisTaskDxHFECorrelation::UserExec(Option_t* /*option*/)
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

  std::auto_ptr<TObjArray> pSelectedD0s(fD0s->Select(pEvent));
  std::auto_ptr<TObjArray> pSelectedElectrons(fElectrons->Select(pEvent));
  int iResult=fCorrelation->Fill(pSelectedD0s.get(), pSelectedElectrons.get());
  if (iResult<0) {
    AliError(Form("%s processing failed with error %d", fCorrelation->GetName(), iResult));
  }

  PostData(1, fOutput);
}

void AliAnalysisTaskDxHFECorrelation::FinishTaskOutput()
{
  // end of the processing
}

void AliAnalysisTaskDxHFECorrelation::Terminate(Option_t *)
{
  // last action on the client
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {
    AliFatal("failed to get output container");
    return;
  }

}
