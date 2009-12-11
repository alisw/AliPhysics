/* $Id: AliPhysicsSelection.cxx 35782 2009-10-22 11:54:31Z jgrosseo $ */

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

//-------------------------------------------------------------------------
//                      Implementation of   Class AliPhysicsSelection
// This class selects collision candidates from data runs, applying selection cuts on triggers 
// and background rejection based on the content of the ESD
//
// Usage:
//
// Create the object and initialize it with the correct run number:
//   fPhysicsSelection = new AliPhysicsSelection;
//   fPhysicsSelection->Initialize(104160);
//
// To check if an event is a collision candidate, use:
//   fPhysicsSelection->IsCollisionCandidate(fESD)
//
// After processing save the resulting histograms to a file with (a folder physics_selection 
//   will be created that contains the histograms):
//   fPhysicsSelection->SaveHistograms("physics_selection")
//
// To print statistics after processing use:
//   fPhysicsSelection->Print();
//
//   Origin: Jan Fiete Grosse-Oetringhaus, CERN
//-------------------------------------------------------------------------

#include <Riostream.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TIterator.h>
#include <TDirectory.h>

#include <AliPhysicsSelection.h>

#include <AliTriggerAnalysis.h>
#include <AliLog.h>

#include <AliESDEvent.h>

ClassImp(AliPhysicsSelection)

AliPhysicsSelection::AliPhysicsSelection() :
  fTriggerAnalysis(0),
  fHistStatistics(0),
  fHistBunchCrossing(0)
{
  // constructor
  
  AliLog::SetClassDebugLevel("AliPhysicsSelection", AliLog::kWarning);
}
    
AliPhysicsSelection::~AliPhysicsSelection()
{
  // destructor
  
  if (fTriggerAnalysis)
  {
    delete fTriggerAnalysis;
    fTriggerAnalysis = 0;
  }

  if (fHistStatistics)
  {
    delete fHistStatistics;
    fHistStatistics = 0;
  }
  
  if (fHistBunchCrossing)
  {
    delete fHistBunchCrossing;
    fHistBunchCrossing = 0;
  }
}
    
Bool_t AliPhysicsSelection::IsCollisionCandidate(const AliESDEvent* aEsd)
{
  // checks if the given event is a collision candidate
  
  if (!fTriggerAnalysis)
    AliFatal("Not initialized!");
    
  const AliESDHeader* esdHeader = aEsd->GetHeader();
  if (!esdHeader)
  {
    AliError("ESD Header could not be retrieved");
    return kFALSE;
  }
  
  // check event type (should be PHYSICS = 7)
  if (esdHeader->GetEventType() != 7)
    return kFALSE;  
  
  fHistStatistics->Fill(1);
  
  fTriggerAnalysis->FillTriggerClasses(aEsd);
    
  for (Int_t i=0; i < fRequTrigClasses.GetEntries(); i++)
  {
    const char* triggerClass = ((TObjString*) fRequTrigClasses.At(i))->String();
    if (!aEsd->IsTriggerClassFired(triggerClass))
    {
      AliDebug(AliLog::kDebug, Form("Rejecting event because trigger class %s is not present", triggerClass));
      return kFALSE;
    }
  }      
  
  for (Int_t i=0; i < fRejTrigClasses.GetEntries(); i++)
  {
    const char* triggerClass = ((TObjString*) fRejTrigClasses.At(i))->String();
    if (aEsd->IsTriggerClassFired(triggerClass))
    {
      AliDebug(AliLog::kDebug, Form("Rejecting event because trigger class %s is present", triggerClass));
      return kFALSE;
    }
  }
  
  fTriggerAnalysis->FillHistograms(aEsd);
  
  fHistStatistics->Fill(2);
    
  Bool_t fastOR = fTriggerAnalysis->IsOfflineTriggerFired(aEsd, AliTriggerAnalysis::kSPDGFO);
  Bool_t v0BB = fTriggerAnalysis->IsOfflineTriggerFired(aEsd, AliTriggerAnalysis::kV0A) || fTriggerAnalysis->IsOfflineTriggerFired(aEsd, AliTriggerAnalysis::kV0C);
  Bool_t v0BG = fTriggerAnalysis->IsOfflineTriggerFired(aEsd, AliTriggerAnalysis::kV0ABG) || fTriggerAnalysis->IsOfflineTriggerFired(aEsd, AliTriggerAnalysis::kV0CBG);
 
  if (fastOR)
    fHistStatistics->Fill(3);
  if (v0BB)
    fHistStatistics->Fill(4);
  if (v0BG)
    fHistStatistics->Fill(5);
    
  if (fastOR || v0BB)
    fHistStatistics->Fill(6);
    
  if (!fastOR && !v0BB)
  {
    AliDebug(AliLog::kDebug, "Rejecting event because neither FO nor V0 has triggered");
    return kFALSE;
  }
  
  if (v0BG)
  {
    AliDebug(AliLog::kDebug, "Rejecting event because of V0 BG flag");
    return kFALSE;
  }
      
  fHistStatistics->Fill(7);
  
  // TODO additional background identification
  
  fHistStatistics->Fill(9);
  
  fHistBunchCrossing->Fill(aEsd->GetBunchCrossNumber());
  
  AliDebug(AliLog::kDebug, "Accepted event");
  
  return kTRUE;
}
    
Bool_t AliPhysicsSelection::Initialize(UInt_t runNumber)
{
  // initializes the object for the given run
  // TODO having the run number here and parameters hardcoded is clearly temporary, a way needs to be found to have a CDB-like configuration also for analysis
  
  AliInfo(Form("Initializing for run %d", runNumber));
  
  fRequTrigClasses.Clear();
  fRejTrigClasses.Clear();
  
  fRequTrigClasses.Add(new TObjString("CINT1B-ABCE-NOPF-ALL"));
  
  if (!fTriggerAnalysis)
  {
    fTriggerAnalysis = new AliTriggerAnalysis;
    fTriggerAnalysis->EnableHistograms();
  }
    
  fTriggerAnalysis->SetSPDGFOThreshhold(1);
  fTriggerAnalysis->SetV0TimeOffset(0);
  
  if (runNumber == 104321)
    fTriggerAnalysis->SetV0TimeOffset(7.5);
  
  if (fHistStatistics)
  {
    fHistStatistics->Reset();
  }
  else
  {
    fHistStatistics = new TH1F("fHistStatistics", "fHistStatistics;;event count", 10, 0.5, 10.5);
    
    Int_t n = 1;
    fHistStatistics->GetXaxis()->SetBinLabel(n++, "Total");
    fHistStatistics->GetXaxis()->SetBinLabel(n++, "Correct trigger class(es)");
    fHistStatistics->GetXaxis()->SetBinLabel(n++, "FO");
    fHistStatistics->GetXaxis()->SetBinLabel(n++, "V0 BB");
    fHistStatistics->GetXaxis()->SetBinLabel(n++, "V0 BG");
    fHistStatistics->GetXaxis()->SetBinLabel(n++, "FO | V0 BB");
    fHistStatistics->GetXaxis()->SetBinLabel(n++, "(FO | V0 BB) & !V0 BG");
    fHistStatistics->GetXaxis()->SetBinLabel(n++, "Background identification");
    fHistStatistics->GetXaxis()->SetBinLabel(n++, "Accepted");
  }
  
  if (fHistBunchCrossing)
  {
    fHistBunchCrossing->Reset();
  }
  else
    fHistBunchCrossing = new TH1F("fHistBunchCrossing", "fHistBunchCrossing;bunch crossing number;accepted events", 4000, -0.5, 3999.5);
    
  return kTRUE;
}

void AliPhysicsSelection::Print(Option_t* /* option */) const
{
  // print the configuration
  
  AliInfo("Configuration:");
  
  TString str("Required trigger classes: ");
  for (Int_t i=0; i < fRequTrigClasses.GetEntries(); i++)
    str += ((TObjString*) fRequTrigClasses.At(i))->String() + " ";
  AliInfo(str);
  
  str = "Rejected trigger classes: ";
  for (Int_t i=0; i < fRejTrigClasses.GetEntries(); i++)
    str += ((TObjString*) fRejTrigClasses.At(i))->String() + " ";
  AliInfo(str);
  
  AliInfo(Form("Requiring %d FO chips (offline) or V0A or V0C and no V0 BG flag", fTriggerAnalysis->GetSPDGFOThreshhold()));
  
  if (fTriggerAnalysis->GetV0TimeOffset() > 0)
    AliInfo(Form("V0 time offset active: %.2f ns", fTriggerAnalysis->GetV0TimeOffset()));
    
  AliInfo("");
  
  AliInfo("Selection statistics:");
    
  AliInfo("Total available events:");
  fTriggerAnalysis->PrintTriggerClasses();
  
  AliInfo(Form("Total events: %d", (Int_t) fHistStatistics->GetBinContent(1)));
  AliInfo(Form("Correct trigger class: %d", (Int_t) fHistStatistics->GetBinContent(2)));
  AliInfo(Form("Selected collision candidates: %d", (Int_t) fHistStatistics->GetBinContent(9)));
}

Long64_t AliPhysicsSelection::Merge(TCollection* list)
{
  // Merge a list of AliMultiplicityCorrection objects with this (needed for
  // PROOF).
  // Returns the number of merged objects (including this).

  if (!list)
    return 0;

  if (list->IsEmpty())
    return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj;

  // collections of all histograms
  const Int_t nHists = 9;
  TList collections[nHists];

  Int_t count = 0;
  while ((obj = iter->Next())) {

    AliPhysicsSelection* entry = dynamic_cast<AliPhysicsSelection*> (obj);
    if (entry == 0) 
      continue;

    Int_t n = 0;
    collections[n++].Add(entry->fTriggerAnalysis);
    collections[n++].Add(entry->fHistStatistics);
    collections[n++].Add(entry->fHistBunchCrossing);

    count++;
  }

  Int_t n = 0;
  fTriggerAnalysis->Merge(&collections[n++]);
  fHistStatistics->Merge(&collections[n++]);
  fHistBunchCrossing->Merge(&collections[n++]);
  
  delete iter;

  return count+1;
}

void AliPhysicsSelection::SaveHistograms(const char* folder) const
{
  // write histograms to current directory
  
  if (!fHistStatistics)
    return;
    
  if (folder)
  {
    gDirectory->mkdir(folder);
    gDirectory->cd(folder);
  }
  
  fHistStatistics->Write();
  fHistBunchCrossing->Write();
  
  gDirectory->mkdir("trigger_histograms");
  gDirectory->cd("trigger_histograms");
  
  fTriggerAnalysis->SaveHistograms();
  
  gDirectory->cd("..");
  
  if (folder)
    gDirectory->cd("..");
}
