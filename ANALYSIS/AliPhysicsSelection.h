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
#include "AliVEvent.h"
#include "AliLog.h"
#include "AliAnalysisManager.h"

//#define VERBOSE_STAT

class AliESDEvent;
class TH2F;
class TH1F;
class TCollection;
class AliTriggerAnalysis;
class AliAnalysisTaskSE;
class AliOADBPhysicsSelection ;
class AliOADBFillingScheme    ;
class AliOADBTriggerAnalysis;
class TPRegexp;

class AliPhysicsSelection : public AliAnalysisCuts
{
public:

public:
  enum {kStatTriggerClass=1,kStatHWTrig,kStatV0ABG,kStatV0CBG,kStatMB1,kStatMB1Prime,kStatLaserCut,kStatFO1,kStatFO2,kStatFO2L1,kStatV0A,kStatV0C, kStatT0BB, kStatT0BG, kStatT0PileUp, kStatZDCA,kStatZDCC,kStatZDCAC,kStatZDCTime,kStatV0,kStatOffline,kStatBG, kStatAcceptedPileUp, kStatAccepted};

#ifdef VERBOSE_STAT
  enum {kStatRowAllB=0, kStatRowAllAC, kStatRowAllE, kStatRowBG,kStatRowAcc,kStatRowBGFrac,kStatRowAccFrac,kStatRowErrGoodFrac,kStatRowGoodFrac,kStatRowErrGood,kStatRowGood}; // offset wrt fBGStatOffset
#else
  enum {kStatRowAllB=0, kStatRowAllAC, kStatRowAllE, kStatRowBG, kStatRowAcc,kStatRowGood}; // offset wrt fBGStatOffset
#endif

  enum {kStatIdxAll=0,kStatIdxBin0=1};
  enum ETriggerLogic { kCINT1 = 0, kCMBS2A, kCMBS2C, kCMBAC, kCMBACS2, kHighMultL1 };

  typedef Bool_t (*Bin0Callback_t)(const AliESDEvent *);

  AliPhysicsSelection();
  virtual ~AliPhysicsSelection();
    
  // AliAnalysisCuts interface
  virtual UInt_t GetSelectionMask(const TObject* obj) { return IsCollisionCandidate((const AliESDEvent*) obj); }
  virtual Bool_t IsSelected(TList*) { return kFALSE; }
  virtual Bool_t IsSelected(TObject*)  {return kFALSE;}
    
  Int_t  GetCurrentRun() const {return fCurrentRun;}
  UInt_t IsCollisionCandidate(const AliESDEvent* aEsd);
  Bool_t Initialize(const AliESDEvent* aEsd);
  Bool_t Initialize(Int_t runNumber);
    
  void SetAnalyzeMC(Bool_t flag = kTRUE) { fMC = flag; }
  void SetSkipTriggerClassSelection(Bool_t flag = kTRUE) { fSkipTriggerClassSelection = flag; }
  void SetSkipV0(Bool_t flag=kTRUE) { AliError ("This method is deprecated! Please change the OADB configuration instead") ; fSkipV0 = flag;}
   
  void AddBackgroundIdentification(AliAnalysisCuts*) { AliFatal ("This method is deprecated. Please use the kSPDClsVsTrkBG tag in the OADB configuration instead"); }
    
  virtual void Print(const Option_t* option = "") const;
  virtual Long64_t Merge(TCollection* list);
  void SaveHistograms(const char* folder = 0);
    
  const TList* GetCollisionTriggerClasses() const { return &fCollTrigClasses; }
  const TList* GetBGTriggerClasses()        const { return &fBGTrigClasses; }
  void AddCollisionTriggerClass(const char* className);
  void AddBGTriggerClass(const char* className)       ;
  void SetCustomOADBObjects(AliOADBPhysicsSelection * oadbPS, AliOADBFillingScheme * oadbFS, AliOADBTriggerAnalysis * oadbTA = 0) { fPSOADB = oadbPS; fFillOADB = oadbFS; fTriggerOADB = oadbTA; fUsingCustomClasses = kTRUE;}
  const AliOADBPhysicsSelection * GetOADBPhysicsSelection() const {return fPSOADB;  }
  const AliOADBFillingScheme    * GetOADBFillingScheme()    const {return fFillOADB;}
  const AliOADBTriggerAnalysis  * GetOADBTriggerAnalysis()  const {return fTriggerOADB;}

  Int_t GetBGStatOffset() const { return fBGStatOffset; }


  AliTriggerAnalysis* GetTriggerAnalysis(Int_t i = 0) { return (fTriggerAnalysis.GetEntries() > 0) ? (AliTriggerAnalysis*) fTriggerAnalysis.At(i) : 0; }    
    
  const TH2F* GetStatisticsHistogram(Int_t idx=kStatIdxAll) const { return fHistStatistics[idx]; }
  const TH2F* GetBunchCrossingHistogram() const { return fHistBunchCrossing; }
  virtual TObject *GetStatistics(const Option_t *option) const;
    
  void SetBIFactors(const AliESDEvent * aESD);
  
  void SetUseBXNumbers(Bool_t flag = kTRUE) {fUseBXNumbers = flag;}
  void SetComputeBG   (UInt_t flag = AliVEvent::kMB) {fComputeBG    = flag; if(flag) fUseBXNumbers = flag;}
  void SetUseMuonTriggers(Bool_t) { Printf("ERROR: Use of AliPhysicsSelection::SetUseMuonTriggers is deprecated. Use AliVEvent::kMUON selection instead."); }
  void SetBin0Callback( const char * cb) {fBin0CallBack = cb;} 
  void SetBin0CallbackViaPointer( Bin0Callback_t cb) {fBin0CallBackPointer = cb;}// WARNING: THIS SHOULD NOT BE USED, WILL BE REMOVED SOON
  
  static const char * GetOADBFileName() {   static TString filename; filename.Form("%s/COMMON/PHYSICSSELECTION/data/physicsSelection.root", AliAnalysisManager::GetOADBPath()); return filename.Data();};

protected:
  UInt_t CheckTriggerClass(const AliESDEvent* aEsd, const char* trigger, Int_t& triggerLogic) const;
  Bool_t EvaluateTriggerLogic(const AliESDEvent* aEsd, AliTriggerAnalysis* triggerAnalysis, const char* triggerLogic, Bool_t offline);
  TH2F * BookHistStatistics(const char * tag) ;
  Int_t GetStatRow(const char * triggerBXClass, UInt_t offlineTriggerType, UInt_t ** rowIDs) const;
  const char * GetTriggerString(TObjString * obj);


  Int_t fCurrentRun;      // run number for which the object is initialized
  Bool_t fMC;             // flag if MC is analyzed
  TList fCollTrigClasses; // trigger class identifying collision candidates
  TList fBGTrigClasses;   // trigger classes identifying background events
    
  TList fTriggerAnalysis; // list of offline trigger objects (several are needed to keep the control histograms separate per trigger class)
      
  TH2F* fHistStatistics[2];      // how many events are cut away why {all,bin 0}
  //  TH2F* fHistStatisticsTokens;   // how many events are cut away why (new version, bins for all tokens. Only used tokens are filled)
  TH2F* fHistBunchCrossing;   // histograms of accepted bunch crossing numbers
  TH1F* fHistTriggerPattern;  // Pattern of the individual detectors in the MB1 trigger. Can reveal inconsistencies/inefficiencies in the trigger 
    
  Bool_t fSkipTriggerClassSelection;  // flag that determines if the trigger classs selection is skipped
  Bool_t fUsingCustomClasses;         // flag that is set if costum trigger classes are defined
  Bool_t fSkipV0;                     // ignore information from v0

  Float_t fBIFactorA;                 // ratio of interacting over non interacting bunch intensities for beam 1
  Float_t fBIFactorC;                 // ratio of interacting over non interacting bunch intensities for beam 2
  Float_t fBIFactorAC;                // ratio of interacting over non interacting bunch intensities for beam 1&2: 
                                      // for some trigger classes we do not distinguish betweem 1 and 2


  UInt_t fComputeBG; // Switch on computation of background and filling of relevant stat table entries. If you enable this you can only process one run at a time (the relative bunch intensity used to compute this chages from run to run); This is a mask of the trigger types to be used in the selection, as defined in AliVEvent
  Int_t fBGStatOffset; // Offset of the BG statistics computed at the end of the processing. It depends on how many trigger classes you have selected
  Bool_t fUseBXNumbers;	// Explicitely select "good" bunch crossing numbers (exclude pilot, afterpulses and fakes). If you anable this you can only process  runs within the same filling scheme.
  Bool_t fUseMuonTriggers;	// if true, also use the muon triggers
  TString fFillingScheme; // stores the filling scheme of the current run.

  TString fBin0CallBack; // callback used to determine if an event is in the bin0 (name of the task where the callback is implemented);
  Bin0Callback_t fBin0CallBackPointer; //! don't stream this. TO BE REMOVED SOON

  Bool_t fIsPP; // True if processing pp run, false if heavy ion

  AliOADBPhysicsSelection * fPSOADB; // Physics selection OADB object
  AliOADBFillingScheme    * fFillOADB; // Filling scheme OADB object
  AliOADBTriggerAnalysis  * fTriggerOADB; // Trigger analysis OADB object

  TPRegexp* fRegexp; //! regular expression for trigger tokens
  TList* fCashedTokens; //! trigger token lookup list

  ClassDef(AliPhysicsSelection, 15)
    
    private:
  AliPhysicsSelection(const AliPhysicsSelection&);
  AliPhysicsSelection& operator=(const AliPhysicsSelection&);
};

#endif
