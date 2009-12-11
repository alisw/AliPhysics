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

#include <TObject.h>
#include <TList.h>

class AliESDEvent;
class TH1;
class TH2;
class TCollection;
class AliTriggerAnalysis;

class AliPhysicsSelection : public TObject
{
  public:
    AliPhysicsSelection();
    virtual ~AliPhysicsSelection();
    
    Bool_t IsCollisionCandidate(const AliESDEvent* aEsd);
    Bool_t Initialize(UInt_t runNumber);
    
    virtual void Print(Option_t* option = "") const;
    virtual Long64_t Merge(TCollection* list);
    void SaveHistograms(const char* folder = 0) const;
    
    TList* GetRequiredTriggerClasses() { return &fRequTrigClasses; }
    TList* GetRejectedTriggerClasses() { return &fRejTrigClasses; }
    
  protected:
    AliTriggerAnalysis* fTriggerAnalysis; // offline trigger object
  
    TList fRequTrigClasses; // list of trigger class required for collision candidates
    TList fRejTrigClasses;  // list of trigger classes not allowed for collision candidates
    
    TH1* fHistStatistics;      // how many events are cut away why
    TH1* fHistBunchCrossing;   // histograms of accepted bunch crossing numbers
    
    ClassDef(AliPhysicsSelection, 1)
    
  private:
    AliPhysicsSelection(const AliPhysicsSelection&);
    AliPhysicsSelection& operator=(const AliPhysicsSelection&);
};

#endif
