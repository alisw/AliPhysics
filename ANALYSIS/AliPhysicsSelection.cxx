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
// Create the object:
//   fPhysicsSelection = new AliPhysicsSelection;
//
// For MC data, call
//   fPhysicsSelection->SetAnalyzeMC()
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
// Usually the class selects the trigger scheme by itself depending on the run number.
// Nevertheless, you can do that manually by calling AddCollisionTriggerClass() and AddBGTriggerClass()
// Example:
// To define the class CINT1B-ABCE-NOPF-ALL as collision trigger (those will be accepted as  
// collision candidates when they pass the selection):
//   AddCollisionTriggerClass("+CINT1B-ABCE-NOPF-ALL #769 #3119");
// To select on bunch crossing IDs in addition, use:
//   AddCollisionTriggerClass("+CINT1B-ABCE-NOPF-ALL #769 #3119");
// To define the class CINT1A-ABCE-NOPF-ALL as a background trigger (those will only be counted
// for the control histograms):
//   AddBGTriggerClass("+CINT1A-ABCE-NOPF-ALL");
// You can also specify more than one trigger class in a string or you can require that some are *not*
// present. The following line would require CSMBA-ABCE-NOPF-ALL, but CSMBB-ABCE-NOPF-ALL is not allowed
// to be present:
//   AddBGTriggerClass("+CSMBA-ABCE-NOPF-ALL -CSMBB-ABCE-NOPF-ALL");
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
  fMC(kFALSE),
  fCollTrigClasses(),
  fBGTrigClasses(),
  fTriggerAnalysis(),
  fBackgroundIdentification(0),
  fHistStatistics(0),
  fHistBunchCrossing(0),
  fSkipTriggerClassSelection(0),
  fUsingCustomClasses(0)
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
  
  Bool_t foundBCRequirement = kFALSE;
  Bool_t foundCorrectBC = kFALSE;
  
  TString str(trigger);
  TObjArray* tokens = str.Tokenize(" ");
  
  for (Int_t i=0; i < tokens->GetEntries(); i++)
  {
    TString str2(((TObjString*) tokens->At(i))->String());
    
    if (str2[0] == '+' || str2[0] == '-')
    {
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
    else if (str2[0] == '#')
    {
      foundBCRequirement = kTRUE;
    
      str2.Remove(0, 1);
      
      Int_t bcNumber = str2.Atoi();
      AliDebug(AliLog::kDebug, Form("Checking for bunch crossing number %d", bcNumber));
      
      if (aEsd->GetBunchCrossNumber() == bcNumber)
      {
        foundCorrectBC = kTRUE;
        AliDebug(AliLog::kDebug, Form("Found correct bunch crossing %d", bcNumber));
      }
    }
    else
      AliFatal(Form("Invalid trigger syntax: %s", trigger));
  }
  
  delete tokens;
  
  if (foundBCRequirement && !foundCorrectBC)
    return kFALSE;
  
  return kTRUE;
}
    
Bool_t AliPhysicsSelection::IsCollisionCandidate(const AliESDEvent* aEsd)
{
  // checks if the given event is a collision candidate
  
  if (fCurrentRun != aEsd->GetRunNumber())
    if (!Initialize(aEsd->GetRunNumber()))
      AliFatal(Form("Could not initialize for run %d", aEsd->GetRunNumber()));
    
  const AliESDHeader* esdHeader = aEsd->GetHeader();
  if (!esdHeader)
  {
    AliError("ESD Header could not be retrieved");
    return kFALSE;
  }
  
  // check event type; should be PHYSICS = 7 for data and 0 for MC
  if (!fMC)
  {
    if (esdHeader->GetEventType() != 7)
      return kFALSE;
  }
  else
  {
    if (esdHeader->GetEventType() != 0)
      AliFatal(Form("Invalid event type for MC: %d", esdHeader->GetEventType()));
  }
  
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
      
      // hardware trigger (should only remove events for MC)
      // replay CINT1B hardware trigger
      // TODO this has to depend on the actual hardware trigger (and that depends on the run...)
      Int_t fastORHW = triggerAnalysis->SPDFiredChips(aEsd, 1); // SPD number of chips from trigger bits (!)
      
      AliTriggerAnalysis::V0Decision v0ADecision = triggerAnalysis->V0Trigger(aEsd, AliTriggerAnalysis::kASide);
      AliTriggerAnalysis::V0Decision v0CDecision = triggerAnalysis->V0Trigger(aEsd, AliTriggerAnalysis::kCSide);
      
      Bool_t v0A = (v0ADecision == AliTriggerAnalysis::kV0BB);
      Bool_t v0C = (v0CDecision == AliTriggerAnalysis::kV0BB);
        
      if (fastORHW == 0 && !v0A && !v0C)
      {
        AliDebug(AliLog::kDebug, "Rejecting event because hardware trigger is not fired");
        continue;
      }
      
      fHistStatistics->Fill(2, i);
    
      // offline trigger
      Int_t fastOROffline = triggerAnalysis->SPDFiredChips(aEsd, 0); // SPD number of chips from clusters (!)
      Bool_t v0ABG = (v0ADecision == AliTriggerAnalysis::kV0BG);
      Bool_t v0CBG = (v0CDecision == AliTriggerAnalysis::kV0BG);
      Bool_t v0BG = v0ABG || v0CBG;
    
      if (fastOROffline > 0)
        fHistStatistics->Fill(3, i);
      if (fastOROffline > 1)
        fHistStatistics->Fill(4, i);
        
      if (v0A)
        fHistStatistics->Fill(5, i);
      if (v0C)
        fHistStatistics->Fill(6, i);
      if (v0ABG)
        fHistStatistics->Fill(7, i);
      if (v0CBG)
        fHistStatistics->Fill(8, i);
        
      if (v0ADecision == AliTriggerAnalysis::kV0Fake)
        fHistStatistics->Fill(9, i);
      if (v0CDecision == AliTriggerAnalysis::kV0Fake)
        fHistStatistics->Fill(10, i);
      
      if (fastOROffline > 1 && !v0BG)
        fHistStatistics->Fill(11, i);
        
      if ((fastOROffline > 0 || v0A || v0C) && !v0BG)
        fHistStatistics->Fill(12, i);
    
      if (fastOROffline > 0 && (v0A || v0C) && !v0BG)
        fHistStatistics->Fill(13, i);
  
      if (v0A && v0C && !v0BG)
        fHistStatistics->Fill(14, i);
      
      if (fastOROffline > 1 || (fastOROffline > 0 && (v0A || v0C)) || (v0A && v0C))
      {
        if (!v0BG)
        {
          fHistStatistics->Fill(15, i);
      
          if (fBackgroundIdentification && !fBackgroundIdentification->IsSelected(const_cast<AliESDEvent*> (aEsd)))
          {
            AliDebug(AliLog::kDebug, "Rejecting event because of background identification");
            fHistStatistics->Fill(16, i);
          }
          else
          {
            AliDebug(AliLog::kDebug, "Accepted event for histograms");
            
            fHistStatistics->Fill(17, i);
            fHistBunchCrossing->Fill(aEsd->GetBunchCrossNumber(), i);
            if (i < fCollTrigClasses.GetEntries() || fSkipTriggerClassSelection)
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

Int_t AliPhysicsSelection::GetTriggerScheme(UInt_t runNumber)
{
  // returns the current trigger scheme (classes that are accepted/rejected)
  
  if (fMC)
    return 0;
    
  // TODO dependent on run number
  
  switch (runNumber)
  {
    // CSMBB triggers
    case 104044:
    case 105054:
    case 105057:
      return 2;
  }
  
  // default: CINT1 suite
  return 1;
}  
    
Bool_t AliPhysicsSelection::Initialize(UInt_t runNumber)
{
  // initializes the object for the given run
  // TODO having the run number here and parameters hardcoded is clearly temporary, a way needs to be found to have a CDB-like configuration also for analysis
  
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  Int_t triggerScheme = GetTriggerScheme(runNumber);
  
  if (!fUsingCustomClasses && fCurrentRun != -1 && triggerScheme != GetTriggerScheme(fCurrentRun))
    AliFatal("Processing several runs with different trigger schemes is not supported");
  
  AliInfo(Form("Initializing for run %d", runNumber));
  
  // initialize first time?
  if (fCurrentRun == -1)
  {
    if (fUsingCustomClasses) {
      AliInfo("Using user-provided trigger classes");
    } else {
      switch (triggerScheme)
      {
      case 0:
        fCollTrigClasses.Add(new TObjString(""));
        break;
        
      case 1:
        fCollTrigClasses.Add(new TObjString("+CINT1B-ABCE-NOPF-ALL"));
        fBGTrigClasses.Add(new TObjString("+CINT1A-ABCE-NOPF-ALL"));
        fBGTrigClasses.Add(new TObjString("+CINT1C-ABCE-NOPF-ALL"));
        fBGTrigClasses.Add(new TObjString("+CINT1-E-NOPF-ALL"));
        break;
        
      case 2:
        fCollTrigClasses.Add(new TObjString("+CSMBB-ABCE-NOPF-ALL"));
        fBGTrigClasses.Add(new TObjString("+CSMBA-ABCE-NOPF-ALL -CSMBB-ABCE-NOPF-ALL"));
        fBGTrigClasses.Add(new TObjString("+CSMBC-ABCE-NOPF-ALL -CSMBB-ABCE-NOPF-ALL"));
        break;
        
      default:
        AliFatal(Form("Unsupported trigger scheme %d", triggerScheme));
      }
    }
    
    Int_t count = fCollTrigClasses.GetEntries() + fBGTrigClasses.GetEntries();
    
    for (Int_t i=0; i<count; i++)
    {
      AliTriggerAnalysis* triggerAnalysis = new AliTriggerAnalysis;
      triggerAnalysis->SetAnalyzeMC(fMC);
      triggerAnalysis->EnableHistograms();
      triggerAnalysis->SetSPDGFOThreshhold(1);
      fTriggerAnalysis.Add(triggerAnalysis);
    }
      
    if (fHistStatistics)
      delete fHistStatistics;
  
    fHistStatistics = new TH2F("fHistStatistics", "fHistStatistics;;", 17, 0.5, 17.5, count, -0.5, -0.5 + count);
      
    Int_t n = 1;
    fHistStatistics->GetXaxis()->SetBinLabel(n++, "Trigger class");
    fHistStatistics->GetXaxis()->SetBinLabel(n++, "Hardware trigger");
    fHistStatistics->GetXaxis()->SetBinLabel(n++, "FO >= 1");
    fHistStatistics->GetXaxis()->SetBinLabel(n++, "FO >= 2");
    fHistStatistics->GetXaxis()->SetBinLabel(n++, "V0A");
    fHistStatistics->GetXaxis()->SetBinLabel(n++, "V0C");
    fHistStatistics->GetXaxis()->SetBinLabel(n++, "V0A BG");
    fHistStatistics->GetXaxis()->SetBinLabel(n++, "V0C BG");
    fHistStatistics->GetXaxis()->SetBinLabel(n++, "V0A Fake");
    fHistStatistics->GetXaxis()->SetBinLabel(n++, "V0C Fake");
    fHistStatistics->GetXaxis()->SetBinLabel(n++, "FO >= 2 &!V0 BG");
    fHistStatistics->GetXaxis()->SetBinLabel(n++, "(FO >= 1 | V0A | V0C) & !V0 BG");
    fHistStatistics->GetXaxis()->SetBinLabel(n++, "FO >= 1 & (V0A | V0C) & !V0 BG");
    fHistStatistics->GetXaxis()->SetBinLabel(n++, "V0A & V0C & !V0 BG");
    fHistStatistics->GetXaxis()->SetBinLabel(n++, "(FO >= 2 | (FO >= 1 & (V0A | V0C)) | (V0A & V0C)) & !V0 BG");
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
  }
    
  Int_t count = fCollTrigClasses.GetEntries() + fBGTrigClasses.GetEntries();
  for (Int_t i=0; i<count; i++)
  {
    AliTriggerAnalysis* triggerAnalysis = static_cast<AliTriggerAnalysis*> (fTriggerAnalysis.At(i));
  
    switch (runNumber)
    {
      case 104315:
      case 104316:
      case 104320:
      case 104321:
      case 104439:
        triggerAnalysis->SetV0TimeOffset(7.5);
        break;
      default:
        triggerAnalysis->SetV0TimeOffset(0);
    }
  }
    
  fCurrentRun = runNumber;
  
  TH1::AddDirectory(oldStatus);
  
  return kTRUE;
}

void AliPhysicsSelection::Print(Option_t* /* option */) const
{
  // print the configuration
  
  Printf("Configuration initialized for run %d (MC: %d):", fCurrentRun, fMC);
  
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
  
  if (fHistStatistics && fCollTrigClasses.GetEntries() > 0)
  {
    Printf("\nSelection statistics for first collision trigger (%s):", ((TObjString*) fCollTrigClasses.First())->String().Data());
    
    Printf("Total events with correct trigger class: %d", (Int_t) fHistStatistics->GetBinContent(1, 1));
    Printf("Selected collision candidates: %d", (Int_t) fHistStatistics->GetBinContent(fHistStatistics->GetXaxis()->FindBin("Accepted"), 1));
  }
  
  if (fHistBunchCrossing)
  {
    Printf("\nBunch crossing statistics:");
    
    for (Int_t i=1; i<=fHistBunchCrossing->GetNbinsY(); i++)
    {
      TString str;
      str.Form("Trigger %s has accepted events in the bunch crossings: ", fHistBunchCrossing->GetYaxis()->GetBinLabel(i));
      
      for (Int_t j=1; j<=fHistBunchCrossing->GetNbinsX(); j++)
        if (fHistBunchCrossing->GetBinContent(j, i) > 0)
          str += Form("%d, ", (Int_t) fHistBunchCrossing->GetXaxis()->GetBinCenter(j));
       
      Printf("%s", str.Data());
    }
    
    for (Int_t j=1; j<=fHistBunchCrossing->GetNbinsX(); j++)
    {
      Int_t count = 0;
      for (Int_t i=1; i<=fHistBunchCrossing->GetNbinsY(); i++)
      {
        if (fHistBunchCrossing->GetBinContent(j, i) > 0)
          count++;
      }
      if (count > 1)
        Printf("WARNING: Bunch crossing %d has more than one trigger class active. Check BPTX functioning for this run!", (Int_t) fHistBunchCrossing->GetXaxis()->GetBinCenter(j));
    }
  }

  if (fUsingCustomClasses)        
    Printf("WARNING: Using custom trigger classes!");
  if (fSkipTriggerClassSelection) 
    Printf("WARNING: Skipping trigger class selection!");
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
      
    collections[0].Add(&(entry->fTriggerAnalysis));
    if (entry->fHistStatistics)
      collections[1].Add(entry->fHistStatistics);
    if (entry->fHistBunchCrossing)
      collections[2].Add(entry->fHistBunchCrossing);
    if (entry->fBackgroundIdentification)
      collections[3].Add(entry->fBackgroundIdentification);

    count++;
  }

  fTriggerAnalysis.Merge(&collections[0]);
  if (fHistStatistics)
    fHistStatistics->Merge(&collections[1]);
  if (fHistBunchCrossing)
    fHistBunchCrossing->Merge(&collections[2]);
  if (fBackgroundIdentification)
    fBackgroundIdentification->Merge(&collections[3]);
  
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
