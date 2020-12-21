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

/* $Id: AliAnalysisTaskHistogram.cxx 46301 2011-01-06 14:25:27Z agheata $ */

#include "AliAnalysisTaskHistogram.h"
#include "TChain.h"
#include "AliAnalysisManager.h"

ClassImp(AliAnalysisTaskHistogram)

AliAnalysisTaskHistogram::AliAnalysisTaskHistogram() 
   :AliAnalysisTaskSE()
{
  // constructor
}

AliAnalysisTaskHistogram::AliAnalysisTaskHistogram(const char *name)
   :AliAnalysisTaskSE(name)
{
  /// constructor

  DefineInput(0, TChain::Class());
  DefineOutput(1, TH1::Class());
}

AliAnalysisTaskHistogram::~AliAnalysisTaskHistogram()
{
   /// destructor

}

void AliAnalysisTaskHistogram::CreateOutputObjects()
{
  fPtHist = new TH1F("fPtHist", "fPtHist", 100, 0, 20);
  PostData(1, fPtHist);
}


void AliAnalysisTaskHistogram::UserExec(Option_t *) 
{
  Int_t nTracks = fInputEvent->GetNumberOfTracks();
  
  for (int i=0; i<nTracks; i++)
    fPtHist->Fill(fInputEvent->GetTrack(i)->Pt());
}

AliAnalysisTaskHistogram *AliAnalysisTaskHistogram::AddTask(TString suffix)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    return nullptr;
  }
  // get the input event handler, again via a static method.
  // this handler is part of the managing system and feeds events
  // to your task
  if (!mgr->GetInputEventHandler())
  {
    return nullptr;
  }
  // by default, a file is open for writing. here, we get the filename
  TString fileName = "AO2D.root";
  if (!suffix.IsNull())
    fileName += ":" + suffix; // create a subfolder in the file
  // now we create an instance of your task
  AliAnalysisTaskHistogram *task = new AliAnalysisTaskHistogram((TString("AliAnalysisTaskHistogram") + suffix).Data());
  if (!task)
    return nullptr;
  // add your task to the manager
  mgr->AddTask(task);
  // your task needs input: here we connect the manager to your task
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  // same for the output
  mgr->ConnectOutput(task, 1, mgr->CreateContainer("HIST", TH1::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName()));
  return task;
}
