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
//           Michele Floris, CERN
//-------------------------------------------------------------------------

#include <AliAnalysisCuts.h>
#include <TList.h>
#include "TObjString.h"

#define VERBOSE_STAT

class AliESDEvent;
class TH2F;
class TH1F;
class TCollection;
class AliTriggerAnalysis;
class AliAnalysisTaskSE;

class AliPhysicsSelection : public AliAnalysisCuts
{
public:

  enum {kStatTriggerClass=1,kStatHWTrig,kStatV0ABG,kStatV0CBG,kStatMB1,kStatMB1Prime,kStatFMD,kStatFO1,kStatFO2,kStatV0A,kStatV0C,kStatSSD1,kStatFO1AndV0,kStatV0,kStatAny2Hits,kStatOffline,kStatBG,kStatAccepted};

#ifdef VERBOSE_STAT
  enum {kStatRowBG=5,kStatRowAcc,kStatRowBGFrac,kStatRowAccFrac,kStatRowErrGoodFrac,kStatRowGoodFrac,kStatRowErrGood,kStatRowGood};
#else
  enum {kStatRowBG=5,kStatRowAcc,kStatRowGood};
#endif

  enum {kStatIdxAll=0,kStatIdxBin0=1};


  typedef Bool_t (*Bin0Callback_t)(const AliESDEvent *);

  AliPhysicsSelection();
  virtual ~AliPhysicsSelection();
    
  // AliAnalysisCuts interface
  virtual Bool_t IsSelected(TObject* obj) { return IsCollisionCandidate((const AliESDEvent*) obj); }
  virtual Bool_t IsSelected(TList*) { return kFALSE; }
    
  Bool_t IsCollisionCandidate(const AliESDEvent* aEsd);
  Bool_t Initialize(Int_t runNumber);
    
  void SetAnalyzeMC(Bool_t flag = kTRUE) { fMC = flag; }
  void SetSkipTriggerClassSelection(Bool_t flag = kTRUE) { fSkipTriggerClassSelection = flag; }
  void SetSkipV0(Bool_t flag=kTRUE) { fSkipV0 = flag;}
   
  void AddBackgroundIdentification(AliAnalysisCuts* background) { fBackgroundIdentification = background; }
    
  virtual void Print(Option_t* option = "") const;
  virtual Long64_t Merge(TCollection* list);
  void SaveHistograms(const char* folder = 0) const;
    
  const TList* GetCollisionTriggerClasses() const { return &fCollTrigClasses; }
  const TList* GetBGTriggerClasses()        const { return &fBGTrigClasses; }
  void AddCollisionTriggerClass(const char* className){ fCollTrigClasses.Add(new TObjString(className)); fUsingCustomClasses = kTRUE; }
  void AddBGTriggerClass(const char* className)       { fBGTrigClasses.Add(new TObjString(className));  fUsingCustomClasses = kTRUE; }
  
  AliTriggerAnalysis* GetTriggerAnalysis(Int_t i = 0) { return (fTriggerAnalysis.GetEntries() > 0) ? (AliTriggerAnalysis*) fTriggerAnalysis.At(i) : 0; }    
    
  const TH2F* GetStatisticsHistogram(Int_t idx=kStatIdxAll) const { return fHistStatistics[idx]; }
  const TH2F* GetBunchCrossingHistogram() const { return fHistBunchCrossing; }
    
  void SetBIFactors(Int_t run);
  
  void SetUseBXNumbers(Bool_t flag = kTRUE) {fUseBXNumbers = flag;}
  void SetComputeBG   (Bool_t flag = kTRUE) {fComputeBG    = flag; if(flag) fUseBXNumbers = flag;}
  void SetUseMuonTriggers(Bool_t flag = kTRUE) {fUseMuonTriggers = flag;}
  void SetBin0Callback( const char * cb) {fBin0CallBack = cb;} 
  void SetBin0CallbackViaPointer( Bin0Callback_t cb) {fBin0CallBackPointer = cb;}// WARNING: THIS SHOULD NOT BE USED, WILL BE REMOVED SOON
  

protected:
  Bool_t CheckTriggerClass(const AliESDEvent* aEsd, const char* trigger) const;
  Int_t GetTriggerScheme(UInt_t runNumber) const;
  const char * GetBXIDs(UInt_t runNumber, const char * trigger ) ;
  const char * GetFillingScheme(UInt_t runNumber) ;
  Int_t GetRatioBBBE(Int_t runNumber);
  TH2F * BookHistStatistics(const char * tag) ;

    
  Int_t fCurrentRun;      // run number for which the object is initialized
  Bool_t fMC;             // flag if MC is analyzed
  TList fCollTrigClasses; // trigger class identifying collision candidates
  TList fBGTrigClasses;   // trigger classes identifying background events
    
  TList fTriggerAnalysis; // list of offline trigger objects (several are needed to keep the control histograms separate per trigger class)
  
  AliAnalysisCuts* fBackgroundIdentification; // class that performs additional background identification
    
  TH2F* fHistStatistics[2];      // how many events are cut away why {all,bin 0}
  TH2F* fHistBunchCrossing;   // histograms of accepted bunch crossing numbers
  TH1F* fHistTriggerPattern;  // Pattern of the individual detectors in the MB1 trigger. Can reveal inconsistencies/inefficiencies in the trigger 
    
  Bool_t fSkipTriggerClassSelection;  // flag that determines if the trigger classs selection is skipped
  Bool_t fUsingCustomClasses;         // flag that is set if costum trigger classes are defined
  Bool_t fSkipV0;                     // ignore information from v0

  Float_t fBIFactorA;                 // ratio of interacting over non interacting bunch intensities for beam 1
  Float_t fBIFactorC;                 // ratio of interacting over non interacting bunch intensities for beam 2

  Int_t fRatioBEEE; // ratio between the number of BX in the Beam-Empty and the Empty-Empty. Depends on the filling and on the trigger scheme

  Bool_t fComputeBG; // Switch on computation of background and filling of relevant stat table entries. If you enable this you can only process one run at a time (the relative bunch intensity used to compute this chages from run to run)
  Bool_t fUseBXNumbers;	// Explicitely select "good" bunch crossing numbers (exclude pilot, afterpulses and fakes). If you anable this you can only process  runs within the same filling scheme.
  Bool_t fUseMuonTriggers;	// if true, also use the muon triggers
  TString fFillingScheme; // stores the filling scheme of the current run.

  TString fBin0CallBack; // callback used to determine if an event is in the bin0 (name of the task where the callback is implemented);
  Bin0Callback_t fBin0CallBackPointer; //! don't stream this. TO BE REMOVED SOON

  ClassDef(AliPhysicsSelection, 8)
    
    private:
  AliPhysicsSelection(const AliPhysicsSelection&);
  AliPhysicsSelection& operator=(const AliPhysicsSelection&);
};

#endif
