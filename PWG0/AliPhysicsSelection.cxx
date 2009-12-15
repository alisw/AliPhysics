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
#include <TObjArray.h>

#include <AliPhysicsSelection.h>

#include <AliTriggerAnalysis.h>
#include <AliLog.h>

#include <AliESDEvent.h>

ClassImp(AliPhysicsSelection)

AliPhysicsSelection::AliPhysicsSelection() :
  AliAnalysisCuts("AliPhysicsSelection", "AliPhysicsSelection"),
  fCurrentRun(-1),
  fCollTrigClasses(),
  fBGTrigClasses(),
  fTriggerAnalysis(),
  fBackgroundIdentification(0),
  fHistStatistics(0),
  fHistBunchCrossing(0)
{
  // constructor
  
  fCollTrigClasses.SetOwner(1);
  fBGTrigClasses.SetOwner(1);
  fTriggerAnalysis.SetOwner(1);
  
  AliLog::SetClassDebugLevel("AliPhysicsSelection", AliLog::kWarning);
}
    
AliPhysicsSelection::~AliPhysicsSelection()
{
  // destructor
  
  fCollTrigClasses.Delete();
  fBGTrigClasses.Delete();
  fTriggerAnalysis.Delete();

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

Bool_t AliPhysicsSelection::CheckTriggerClass(const AliESDEvent* aEsd, const char* trigger) const
{
  // checks if the given trigger class(es) are found for the current event
  // format of trigger: +TRIGGER1 -TRIGGER2
  //   requires TRIGGER1 and rejects TRIGGER2
  
  TString str(trigger);
  TObjArray* tokens = str.Tokenize(" ");
  
  for (Int_t i=0; i < tokens->GetEntries(); i++)
  {
    TString str2(((TObjString*) tokens->At(i))->String());
    
    if (str2[0] != '+' && str2[0] != '-')
      AliFatal(Form("Invalid trigger syntax: %s", trigger));
      
    Bool_t flag = (str2[0] == '+');
    
    str2.Remove(0, 1);
    
    if (flag && !aEsd->IsTriggerClassFired(str2))
    {
      AliDebug(AliLog::kDebug, Form("Rejecting event because trigger class %s is not present", str2.Data()));
      delete tokens;
      return kFALSE;
    }
    if (!flag && aEsd->IsTriggerClassFired(str2))
    {
      AliDebug(AliLog::kDebug, Form("Rejecting event because trigger class %s is present", str2.Data()));
      delete tokens;
      return kFALSE;
    }
  }
  
  delete tokens;
  return kTRUE;
}
    
Bool_t AliPhysicsSelection::IsCollisionCandidate(const AliESDEvent* aEsd)
{
  // checks if the given event is a collision candidate
  
  const AliESDHeader* esdHeader = aEsd->GetHeader();
  if (!esdHeader)
  {
    AliError("ESD Header could not be retrieved");
    return kFALSE;
  }
  
  // check event type (should be PHYSICS = 7)
  if (esdHeader->GetEventType() != 7)
    return kFALSE;  
    
  if (fCurrentRun != aEsd->GetRunNumber())
    if (!Initialize(aEsd->GetRunNumber()))
      AliFatal(Form("Could not initialize for run %d", aEsd->GetRunNumber()));
    
  Bool_t accept = kFALSE;
    
  Int_t count = fCollTrigClasses.GetEntries() + fBGTrigClasses.GetEntries();
  for (Int_t i=0; i < count; i++)
  {
    const char* triggerClass = 0;
    if (i < fCollTrigClasses.GetEntries())
      triggerClass = ((TObjString*) fCollTrigClasses.At(i))->String();
    else
      triggerClass = ((TObjString*) fBGTrigClasses.At(i - fCollTrigClasses.GetEntries()))->String();
  
    AliDebug(AliLog::kDebug, Form("Processing trigger class %s", triggerClass));
  
    AliTriggerAnalysis* triggerAnalysis = static_cast<AliTriggerAnalysis*> (fTriggerAnalysis.At(i));
  
    triggerAnalysis->FillTriggerClasses(aEsd);
    
    if (CheckTriggerClass(aEsd, triggerClass))
    {
      triggerAnalysis->FillHistograms(aEsd);
  
      fHistStatistics->Fill(1, i);
    
      Int_t fastOR = triggerAnalysis->SPDFiredChips(aEsd, 0);
      Bool_t v0A = triggerAnalysis->IsOfflineTriggerFired(aEsd, AliTriggerAnalysis::kV0A);
      Bool_t v0C = triggerAnalysis->IsOfflineTriggerFired(aEsd, AliTriggerAnalysis::kV0C);
      Bool_t v0ABG = triggerAnalysis->IsOfflineTriggerFired(aEsd, AliTriggerAnalysis::kV0ABG);
      Bool_t v0CBG = triggerAnalysis->IsOfflineTriggerFired(aEsd, AliTriggerAnalysis::kV0CBG);
      Bool_t v0BG = v0ABG || v0CBG;
    
      if (fastOR > 0)
        fHistStatistics->Fill(2, i);
      if (fastOR > 1)
        fHistStatistics->Fill(3, i);
        
      if (v0A)
        fHistStatistics->Fill(4, i);
      if (v0C)
        fHistStatistics->Fill(5, i);
      if (v0ABG)
        fHistStatistics->Fill(6, i);
      if (v0CBG)
        fHistStatistics->Fill(7, i);
        
      if ((fastOR > 0 || v0A || v0C) && !v0BG)
        fHistStatistics->Fill(8, i);
    
      if (fastOR > 0 && (v0A || v0C) && !v0BG)
        fHistStatistics->Fill(9, i);
  
      if (v0A && v0C && !v0BG)
        fHistStatistics->Fill(10, i);
        
      if (fastOR > 1 || (fastOR > 0 && (v0A || v0C)) || (v0A && v0C))
      {
        if (!v0BG)
        {
          fHistStatistics->Fill(11, i);
      
          if (fBackgroundIdentification && !fBackgroundIdentification->IsSelected(const_cast<AliESDEvent*> (aEsd)))
          {
            AliDebug(AliLog::kDebug, "Rejecting event because of background identification");
            fHistStatistics->Fill(12, i);
          }
          else
          {
            AliDebug(AliLog::kDebug, "Accepted event for histograms");
            
            fHistStatistics->Fill(13, i);
            fHistBunchCrossing->Fill(aEsd->GetBunchCrossNumber(), i);
            if (i < fCollTrigClasses.GetEntries())
              accept = kTRUE;
          }
        }
        else
          AliDebug(AliLog::kDebug, "Rejecting event because of V0 BG flag");
      }
      else
        AliDebug(AliLog::kDebug, "Rejecting event because trigger condition is not fulfilled");
    }
  }
 
  if (accept)
    AliDebug(AliLog::kDebug, "Accepted event as collision candidate");
  
  return accept;
}
    
Bool_t AliPhysicsSelection::Initialize(UInt_t runNumber)
{
  // initializes the object for the given run
  // TODO having the run number here and parameters hardcoded is clearly temporary, a way needs to be found to have a CDB-like configuration also for analysis
  
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  if (fCurrentRun != -1)
    AliFatal("Processing several runs is not supported, yet");
  
  AliInfo(Form("Initializing for run %d", runNumber));
  fCurrentRun = runNumber;
  
  fTriggerAnalysis.Delete();
  fCollTrigClasses.Delete();
  fBGTrigClasses.Delete();
  
  fCollTrigClasses.Add(new TObjString("+CINT1B-ABCE-NOPF-ALL"));
  fBGTrigClasses.Add(new TObjString("+CINT1A-ABCE-NOPF-ALL"));
  fBGTrigClasses.Add(new TObjString("+CINT1C-ABCE-NOPF-ALL"));
  fBGTrigClasses.Add(new TObjString("+CINT1-E-NOPF-ALL"));
  
  Int_t count = fCollTrigClasses.GetEntries() + fBGTrigClasses.GetEntries();
  
  for (Int_t i=0; i<count; i++)
  {
    AliTriggerAnalysis* triggerAnalysis = new AliTriggerAnalysis;
    triggerAnalysis->EnableHistograms();
    triggerAnalysis->SetSPDGFOThreshhold(1);
    triggerAnalysis->SetV0TimeOffset(0);
    
    switch (runNumber)
    {
      case 104316:
      case 104320:
      case 104321: //OK
      case 104439:
        triggerAnalysis->SetV0TimeOffset(7.5);
        break;
    }
  
    fTriggerAnalysis.Add(triggerAnalysis);
  }
      
  if (fHistStatistics)
    delete fHistStatistics;

  fHistStatistics = new TH2F("fHistStatistics", "fHistStatistics;;", 13, 0.5, 13.5, count, -0.5, -0.5 + count);
    
  Int_t n = 1;
  fHistStatistics->GetXaxis()->SetBinLabel(n++, "Correct trigger class(es)");
  fHistStatistics->GetXaxis()->SetBinLabel(n++, "FO >= 1");
  fHistStatistics->GetXaxis()->SetBinLabel(n++, "FO >= 2");
  fHistStatistics->GetXaxis()->SetBinLabel(n++, "V0A");
  fHistStatistics->GetXaxis()->SetBinLabel(n++, "V0C");
  fHistStatistics->GetXaxis()->SetBinLabel(n++, "V0A BG");
  fHistStatistics->GetXaxis()->SetBinLabel(n++, "V0C BG");
  fHistStatistics->GetXaxis()->SetBinLabel(n++, "(FO >= 1 | V0A | VOC) & !V0 BG");
  fHistStatistics->GetXaxis()->SetBinLabel(n++, "FO >= 1 & (V0A | VOC) & !V0 BG");
  fHistStatistics->GetXaxis()->SetBinLabel(n++, "V0A & VOC & !V0 BG");
  fHistStatistics->GetXaxis()->SetBinLabel(n++, "(FO >= 2 | (FO >= 1 & (V0A | VOC)) | (V0A & VOC)) & !V0 BG");
  fHistStatistics->GetXaxis()->SetBinLabel(n++, "Background identification");
  fHistStatistics->GetXaxis()->SetBinLabel(n++, "Accepted");
  
  if (fHistBunchCrossing)
    delete fHistBunchCrossing;

  fHistBunchCrossing = new TH2F("fHistBunchCrossing", "fHistBunchCrossing;bunch crossing number;", 4000, -0.5, 3999.5,  count, -0.5, -0.5 + count);
    
  n = 1;
  for (Int_t i=0; i < fCollTrigClasses.GetEntries(); i++)
  {
    fHistStatistics->GetYaxis()->SetBinLabel(n, ((TObjString*) fCollTrigClasses.At(i))->String());
    fHistBunchCrossing->GetYaxis()->SetBinLabel(n, ((TObjString*) fCollTrigClasses.At(i))->String());
    n++;
  }
  for (Int_t i=0; i < fBGTrigClasses.GetEntries(); i++)
  {
    fHistStatistics->GetYaxis()->SetBinLabel(n, ((TObjString*) fBGTrigClasses.At(i))->String());
    fHistBunchCrossing->GetYaxis()->SetBinLabel(n, ((TObjString*) fBGTrigClasses.At(i))->String());
    n++;
  }
    
  TH1::AddDirectory(oldStatus);
  
  return kTRUE;
}

void AliPhysicsSelection::Print(Option_t* /* option */) const
{
  // print the configuration
  
  Printf("Configuration initialized for run %d:", fCurrentRun);
  
  Printf("Collision trigger classes:");
  for (Int_t i=0; i < fCollTrigClasses.GetEntries(); i++)
    Printf("%s", ((TObjString*) fCollTrigClasses.At(i))->String().Data());
  
  Printf("Background trigger classes:");
  for (Int_t i=0; i < fBGTrigClasses.GetEntries(); i++)
    Printf("%s", ((TObjString*) fBGTrigClasses.At(i))->String().Data());

  AliTriggerAnalysis* triggerAnalysis = dynamic_cast<AliTriggerAnalysis*> (fTriggerAnalysis.At(0));
  
  if (triggerAnalysis)
  {
    if (triggerAnalysis->GetV0TimeOffset() > 0)
      Printf("V0 time offset active: %.2f ns", triggerAnalysis->GetV0TimeOffset());
    
    Printf("\nTotal available events:");
    
    triggerAnalysis->PrintTriggerClasses();
  }
  
  if (fHistStatistics)
  {
    Printf("\nSelection statistics for first collision trigger:");
    
    Printf("Total events with correct trigger class: %d", (Int_t) fHistStatistics->GetBinContent(1, 1));
    Printf("Selected collision candidates: %d", (Int_t) fHistStatistics->GetBinContent(13, 1));
  }
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
    collections[n++].Add(&(entry->fTriggerAnalysis));
    collections[n++].Add(entry->fHistStatistics);
    collections[n++].Add(entry->fHistBunchCrossing);
    if (entry->fBackgroundIdentification)
      collections[n++].Add(entry->fBackgroundIdentification);

    count++;
  }

  Int_t n = 0;
  fTriggerAnalysis.Merge(&collections[n++]);
  fHistStatistics->Merge(&collections[n++]);
  fHistBunchCrossing->Merge(&collections[n++]);
  if (fBackgroundIdentification)
    fBackgroundIdentification->Merge(&collections[n++]);
  
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
  
  Int_t count = fCollTrigClasses.GetEntries() + fBGTrigClasses.GetEntries();
  for (Int_t i=0; i < count; i++)
  {
    TString triggerClass = "trigger_histograms_";
    if (i < fCollTrigClasses.GetEntries())
      triggerClass += ((TObjString*) fCollTrigClasses.At(i))->String();
    else
      triggerClass += ((TObjString*) fBGTrigClasses.At(i - fCollTrigClasses.GetEntries()))->String();
  
    gDirectory->mkdir(triggerClass);
    gDirectory->cd(triggerClass);
  
    static_cast<AliTriggerAnalysis*> (fTriggerAnalysis.At(i))->SaveHistograms();
    
    gDirectory->cd("..");
  }
 
  if (fBackgroundIdentification)
  {
    gDirectory->mkdir("background_identification");
    gDirectory->cd("background_identification");
      
    fBackgroundIdentification->GetOutput()->Write();
      
    gDirectory->cd("..");
  }
  
  if (folder)
    gDirectory->cd("..");
}
