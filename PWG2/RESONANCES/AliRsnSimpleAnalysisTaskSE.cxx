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

//----------------------------------------------------------------------------------
//  Class AliRsnSimpleAnalysisTaskSE
// ------------------------
// Reader for conversion of ESD output into the internal format
// used for resonance study.
// ---
// original author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
// ---
// adapted for Analysis Framework
// by    : R. Vernet                          (email: renaud.vernet@cern.ch)
//----------------------------------------------------------------------------------

#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TROOT.h>
#include <TObjArray.h>

#include "AliLog.h"
#include "AliVEvent.h"
#include "AliMCEventHandler.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"

#include "AliMCEvent.h"

#include "AliRsnPairDef.h"
#include "AliRsnEvent.h"
#include "AliRsnReader.h"
#include "AliRsnPID.h"
#include "AliRsnSimpleFunction.h"
#include "AliRsnSimpleAnalyzer.h"
#include "AliRsnSimpleAnalysisTaskSE.h"

ClassImp(AliRsnSimpleAnalysisTaskSE)

//_____________________________________________________________________________
AliRsnSimpleAnalysisTaskSE::AliRsnSimpleAnalysisTaskSE() :
    AliAnalysisTaskSE(),
    fReader(0x0),
    fPID(0x0),
    fAnalyzer(0x0),
    fRsnEvent(0x0),
    fHistograms(0x0)
{
//
// Default constructor (not recommended)
//
}

//_____________________________________________________________________________
AliRsnSimpleAnalysisTaskSE::AliRsnSimpleAnalysisTaskSE(const char *name) :
    AliAnalysisTaskSE(name),
    fReader(0x0),
    fPID(0x0),
    fAnalyzer(0x0),
    fRsnEvent(0x0),
    fHistograms(0x0)
{
//
// Working constructor (recommended)
//

  DefineOutput(1, TList::Class());
}

//_____________________________________________________________________________
void AliRsnSimpleAnalysisTaskSE::UserCreateOutputObjects()
{
//
// Create the output container
//

  // check for presence of NECESSARY data-members
  if (!fReader)
  {
    AliFatal("Event reader not initialized. Impossible to continue. Aborting with fatal error.");
    return;
  }
  if (!fPID)
  {
    AliFatal("PID manager not initialized. Impossible to continue. Aborting with fatal error.");
    return;
  }
  if (!fAnalyzer)
  {
    AliFatal("Analysis manager not initialized. Impossible to continue. Aborting with fatal error.");
    return;
  }

  // output histogram list
  fHistograms = new TList;

  // initialize analyzer
  fAnalyzer->Init();

  // store all histograms in the functions into the list
  TObjArray *array = fAnalyzer->GetSingle();
  AliRsnSimpleFunction *fcn;
  TH1D *h1D;
  TH2D *h2D;
  if (array)
  {
    TObjArrayIter iter(array);
    while ((fcn = (AliRsnSimpleFunction*)iter.Next()))
    {
      h1D = fcn->GetHistogram1D();
      h2D = fcn->GetHistogram2D();
      if (h1D) fHistograms->AddLast(h1D);
      if (h2D) fHistograms->AddLast(h2D);

    }
  }
  else
  {
    AliWarning("No single-event functions in analyzer");
  }
  array = fAnalyzer->GetMix();
  if (array)
  {
    TObjArrayIter iter(array);
    while ((fcn = (AliRsnSimpleFunction*)iter.Next()))
    {
      h1D = fcn->GetHistogram1D();
      h2D = fcn->GetHistogram2D();
      if (h1D) fHistograms->AddLast(h1D);
      if (h2D) fHistograms->AddLast(h2D);
    }
  }
  else
  {
    AliWarning("No mixing functions in analyzer");
  }
}

//_____________________________________________________________________________
void AliRsnSimpleAnalysisTaskSE::UserExec(Option_t */*option*/)
{
//
// Loops on input container to store data of all tracks.
// Uses the AliRsnReader methods to save them in the output.
//

  // clear previous event
  if (!fRsnEvent)
  {
    fRsnEvent = new AliRsnEvent;
    fRsnEvent->Init();
  }
  else fRsnEvent->Clear();

  // read event
  if (!fReader->Fill(fRsnEvent, fInputEvent, fMCEvent)) AliWarning("Failed reading");
  AliInfo(Form("Collected %d tracks", fRsnEvent->GetMultiplicity()));

  // identify event if the class is available
  if (fPID) fPID->Process(fRsnEvent);

  // process event with analyzer
  fAnalyzer->Process(fRsnEvent);

  // post histograms in slot #1
  PostData(1, fHistograms);
}

//_____________________________________________________________________________
Bool_t AliRsnSimpleAnalysisTaskSE::Configure(const char *configFile)
{
//
// Configure this object using an external macro which creates
// all required objects and stores them in the appropriate data members
// Returns kFALSE if not all the required objects were created.
//

  gROOT->LoadMacro(configFile);
  gROOT->ProcessLine("RsnConfig()");

  fReader = (AliRsnReader*)gDirectory->Get("RsnReader");
  fPID = (AliRsnPID*)gDirectory->Get("RsnPID");
  fAnalyzer = (AliRsnSimpleAnalyzer*)gDirectory->Get("RsnSimpleAnalyzer");

  // check for presence of NECESSARY data-members
  if (!fReader)
  {
    AliError("Event reader not initialized. Impossible to continue. Aborting with fatal error.");
    return kFALSE;
  }
  if (!fPID)
  {
    AliError("PID manager not initialized. Impossible to continue. Aborting with fatal error.");
    return kFALSE;
  }
  if (!fAnalyzer)
  {
    AliError("Analysis manager not initialized. Impossible to continue. Aborting with fatal error.");
    return kFALSE;
  }

  return kTRUE;
}

//_____________________________________________________________________________
void AliRsnSimpleAnalysisTaskSE::PrintSettings()
{
//
// Print analysis settings
//

  AliInfo("==== ANALYSIS TASK INFO =======================================================");

  // reader
  AliInfo(Form("Reader address, name: %x %s", fReader, fReader->GetName()));
  switch (fReader->GetSource())
  {
    case AliRsnReader::kESD:
      AliInfo("Reader source: kESD");
      break;
    case AliRsnReader::kESDTPC:
      AliInfo("Reader source: kESDTPC");
      break;
    case AliRsnReader::kAOD:
      AliInfo("Reader source: kAOD");
      break;
    case AliRsnReader::kMC:
      AliInfo("Reader source: kMC");
      break;
    default:
      AliInfo("Reader source not properly set");
  }
  AliInfo(Form("Reader->CheckSplit  = %s", (fReader->CheckSplit() ? "true" : "false")));
  AliInfo(Form("Reader->RejectFakes = %s", (fReader->RejectFakes() ? "true" : "false")));

  // PID
  Int_t i;
  AliRsnPID::EType type;
  AliInfo(Form("PID address, name: %x", fPID, fPID->GetName()));
  for (i = 0; i < AliRsnPID::kSpecies; i++)
  {
    type = (AliRsnPID::EType)i;
    AliInfo(Form("Prior probability [%d] = %f (%s)", i, fPID->GetPriorProbability(type), AliRsnPID::ParticleName(type)));
  }
  AliInfo(Form("PID momentum    threshold = %f", fPID->GetMaxPt()));
  AliInfo(Form("PID probability threshold = %f", fPID->GetMinProb()));

  // analyzer
  AliRsnSimpleFunction *fcn;
  AliInfo(Form("Analyzer address, name: %x", fAnalyzer, fAnalyzer->GetName()));
  TObjArrayIter iter1(fAnalyzer->GetSingle());
  while ((fcn = (AliRsnSimpleFunction*)iter1.Next()))
  {
    AliInfo(Form("Single-event function: %s [%s]", fcn->GetName(), fcn->ClassName()));
  }
  if (fAnalyzer->GetMix())
  {
    TObjArrayIter iter2(fAnalyzer->GetMix());
    while ((fcn = (AliRsnSimpleFunction*)iter2.Next()))
    {
      AliInfo(Form("Mix function: %s [%s]", fcn->GetName(), fcn->ClassName()));
    }
  }

  AliInfo("===============================================================================");
}
