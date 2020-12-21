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
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliTriggerAnalysis.h"
// In case of ROOT6 it is necessary to stay for the moment with the v5
// version of TFormula as the v6 version produces a large amount of
// warnings at runtime.
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,3,0)
#include <v5/TFormula.h>
typedef ROOT::v5::TFormula R5TFormula;
#else
#include <TFormula.h>
typedef TFormula R5TFormula;
#endif

class AliVEvent;
class TH2F;
class TH1F;
class TCollection;
class AliAnalysisTaskSE;
class AliOADBPhysicsSelection;
class AliOADBFillingScheme;
class AliOADBTriggerAnalysis;
class TPRegexp;
class StringToRegexp;

typedef std::pair<R5TFormula, std::vector<AliTriggerAnalysis::Trigger>> FormulaAndBits;
typedef std::map<std::string, FormulaAndBits> StringToFormula;

class AliPhysicsSelection : public AliAnalysisCuts{
public:
  // These enums are deprecated
  enum {kStatRowAllB=0, kStatRowAllAC, kStatRowAllE, kStatRowBG, kStatRowAcc,kStatRowGood};
  
  typedef Bool_t (*Bin0Callback_t)(const AliESDEvent *);

  AliPhysicsSelection();
  AliPhysicsSelection(const char *name);
  virtual ~AliPhysicsSelection();
    
  virtual UInt_t GetSelectionMask(const TObject* obj) { return IsCollisionCandidate((const AliVEvent*) obj); }
  virtual Bool_t IsSelected(TList*)   { return kFALSE; }
  virtual Bool_t IsSelected(TObject*) { return kFALSE; }
    
  Int_t  GetCurrentRun() const {return fCurrentRun;}
  UInt_t IsCollisionCandidate(const AliVEvent* event);
  Bool_t Initialize(const AliVEvent* event);
  Bool_t Initialize(Int_t runNumber);
  
  void FillStatistics();
  
  virtual void Print(const Option_t* option = "") const;
  virtual Long64_t Merge(TCollection* list);
  void SaveHistograms(const char* folder = 0);
  void ApplyPileupCuts(Bool_t val) { fPileupCutsEnabled = val; }
  
  const TList* GetCollisionTriggerClasses() const { return &fCollTrigClasses; }
  const TList* GetBGTriggerClasses()        const { return &fBGTrigClasses; }
  const AliOADBPhysicsSelection * GetOADBPhysicsSelection() const {return fPSOADB;  }
  const AliOADBFillingScheme    * GetOADBFillingScheme()    const {return fFillOADB;}
  const AliOADBTriggerAnalysis  * GetOADBTriggerAnalysis()  const {return fTriggerOADB;}

  AliTriggerAnalysis* GetTriggerAnalysis(Int_t i = 0) { return (fTriggerAnalysis.GetEntries() > 0) ? (AliTriggerAnalysis*) fTriggerAnalysis.At(i) : 0; }    
  
  void SetAnalyzeMC(Bool_t flag = kTRUE) { fMC = flag; }
  void SetUseBXNumbers(Bool_t flag = kTRUE) {fUseBXNumbers = flag;}
  void SetCustomOADBObjects(AliOADBPhysicsSelection * oadbPS, AliOADBFillingScheme * oadbFS, AliOADBTriggerAnalysis * oadbTA = 0) { fPSOADB = oadbPS; fFillOADB = oadbFS; fTriggerOADB = oadbTA; fUsingCustomClasses = kTRUE;}
  
  virtual TObject *GetStatistics(const Option_t *option) const { return fHistList.FindObject("fHistStat"); }
  void SetBin0Callback( const char * cb) { AliError("This method is deprecated"); } 
  void SetBin0CallbackViaPointer( Bin0Callback_t cb) { AliError("This method is deprecated"); }
  void SetSkipTriggerClassSelection(Bool_t flag = kTRUE) { AliError("This method is deprecated"); }
  
  static const char * GetOADBFileName() { static TString filename; filename.Form("%s/COMMON/PHYSICSSELECTION/data/physicsSelection.root", AliAnalysisManager::GetOADBPath()); return filename.Data();};

  void SetPassName(const TString passName) { fPassName = passName; }
  void DetectPassName();
  void ReadOCDB(Bool_t val) { fReadOCDB=val; }
  Bool_t IsMC() const { return fMC; }
protected:
  UInt_t CheckTriggerClass(const AliVEvent* event, const char* trigger, Int_t& triggerLogic) const;
  Bool_t EvaluateTriggerLogic(const AliVEvent* event, AliTriggerAnalysis* triggerAnalysis, const char* triggerLogic, Bool_t offline);
  const char * GetTriggerString(TObjString * obj);

  TString fPassName;          // pass name for current run
  Int_t fCurrentRun;          // run number for which the object is initialized
  Bool_t fMC;                 // flag if MC is analyzed
  Bool_t fPileupCutsEnabled;  // flag to enable/disable cuts sensitive to in/out-of-bunch pileup
  Bool_t fIsPP;               // True if processing pp run, false if heavy ion
  Bool_t fReadOCDB;           // Flag to read thresholds from OCDB
  Bool_t fUseBXNumbers;       // Explicitly select "good" bunch crossing numbers
  Bool_t fUsingCustomClasses; // flag that is set if custom trigger classes are defined
  TList fCollTrigClasses;     // trigger class identifying collision candidates
  TList fBGTrigClasses;       // trigger classes identifying background events
  TList fTriggerAnalysis;     // list of AliTriggerAnalysis objects (several are needed to keep the control histograms separate per trigger class)
  TList fHistList;            // list of output histos
  TH2F* fHistStat;            //!
  
  AliOADBPhysicsSelection* fPSOADB;      // Physics selection OADB object
  AliOADBFillingScheme*    fFillOADB;    // Filling scheme OADB object
  AliOADBTriggerAnalysis*  fTriggerOADB; // Trigger analysis OADB object

  StringToFormula *fTriggerToFormula; //! Map trigger strings to TFormulas
  FormulaAndBits& FindForumla(const char* triggerLogic); //! Returns pair of TFormula and trigger bits

  StringToRegexp* fTriggerToRegexp; //!
  TPRegexp& FindRegexp(const std::string& triggers) const;

  ClassDef(AliPhysicsSelection, 24)
private:
  AliPhysicsSelection(const AliPhysicsSelection&);
  AliPhysicsSelection& operator=(const AliPhysicsSelection&);
};

#endif
