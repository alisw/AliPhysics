/* $Id: AliPhysicsSelection.h 35782 2009-10-22 11:54:31Z jgrosseo $ */

#ifndef ALIPHYSICSSELECTION_H
#define ALIPHYSICSSELECTION_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                      Implementation of   Class AliPhysicsSelection
// 
// This class selects collision candidates from data runs, applying selection cuts on triggers 
// and background rejection based on the content of the ESD
//
//   Origin: Jan Fiete Grosse-Oetringhaus, CERN
//-------------------------------------------------------------------------

#include <AliAnalysisCuts.h>
#include <TList.h>
#include "TObjString.h"

class AliESDEvent;
class TH2F;
class TCollection;
class AliTriggerAnalysis;

class AliPhysicsSelection : public AliAnalysisCuts
{
  public:
    AliPhysicsSelection();
    virtual ~AliPhysicsSelection();
    
    // AliAnalysisCuts interface
    virtual Bool_t IsSelected(TObject* obj) { return IsCollisionCandidate((const AliESDEvent*) obj); }
    virtual Bool_t IsSelected(TList*) { return kFALSE; }
    
    Bool_t IsCollisionCandidate(const AliESDEvent* aEsd);
    Bool_t Initialize(UInt_t runNumber);
    
    void SetAnalyzeMC(Bool_t flag = kTRUE) { fMC = flag; }
    void SetSkipTriggerClassSelection(Bool_t flag = kTRUE) { fSkipTriggerClassSelection = flag; }
   
    void AddBackgroundIdentification(AliAnalysisCuts* background) { fBackgroundIdentification = background; }
    
    virtual void Print(Option_t* option = "") const;
    virtual Long64_t Merge(TCollection* list);
    void SaveHistograms(const char* folder = 0) const;
    
    const TList* GetCollisionTriggerClasses() const { return &fCollTrigClasses; }
    const TList* GetBGTriggerClasses()        const { return &fBGTrigClasses; }
    void AddCollisionTriggerClass(const char* className){ fCollTrigClasses.Add(new TObjString(className)); fUsingCustomClasses = kTRUE; }
    void AddBGTriggerClass(const char* className)       { fBGTrigClasses.Add(new TObjString(className));  fUsingCustomClasses = kTRUE; }
  
    AliTriggerAnalysis* GetTriggerAnalysis() { return (fTriggerAnalysis.GetEntries() > 0) ? (AliTriggerAnalysis*) fTriggerAnalysis.At(0) : 0; }    
    
    const TH2F* GetStatisticsHistogram() const { return fHistStatistics; }
    const TH2F* GetBunchCrossingHistogram() const { return fHistBunchCrossing; }
    
  protected:
    Bool_t CheckTriggerClass(const AliESDEvent* aEsd, const char* trigger) const;
    Int_t GetTriggerScheme(UInt_t runNumber);
    
    Int_t fCurrentRun;      // run number for which the object is initialized
    Bool_t fMC;             // flag if MC is analyzed
    TList fCollTrigClasses; // trigger class identifying collision candidates
    TList fBGTrigClasses;   // trigger classes identifying background events
    
    TList fTriggerAnalysis; // list of offline trigger objects (several are needed to keep the control histograms separate per trigger class)
  
    AliAnalysisCuts* fBackgroundIdentification; // class that performs additional background identification
    
    TH2F* fHistStatistics;      // how many events are cut away why
    TH2F* fHistBunchCrossing;   // histograms of accepted bunch crossing numbers
    
    Bool_t fSkipTriggerClassSelection;  // flag that determines if the trigger classs selection is skipped
    Bool_t fUsingCustomClasses;         // flag that is set if costum trigger classes are defined

    ClassDef(AliPhysicsSelection, 2)
    
  private:
    AliPhysicsSelection(const AliPhysicsSelection&);
    AliPhysicsSelection& operator=(const AliPhysicsSelection&);
};

#endif
