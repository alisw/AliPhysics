#ifndef ALIANALYSISTASKMULTIDIELECTRON_H
#define ALIANALYSISTASKMULTIDIELECTRON_H
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//#####################################################
//#                                                   # 
//#        Basic Analysis task for Dielectron         #
//#          single event analysis                    #
//#                                                   #
//#  by WooJin J. Park, GSI / W.J.Park@gsi.de         #
//#     Ionut C. Arsene, GSI / I.C.Arsene@gsi.de      #
//#     Magnus Mager, CERN / Magnus.Mager@cern.ch     #
//#     Jens Wiechula, Uni HD / Jens.Wiechula@cern.ch #
//#                                                   #
//#####################################################

#include "TList.h"

#include "AliAnalysisTaskSE.h"

// #include "AliDielectronPID.h"

class AliDielectron;
class TH1D;
class AliAnalysisCuts;
class AliTriggerAnalysis;

class AliAnalysisTaskMultiDielectron : public AliAnalysisTaskSE {
  
public:
  AliAnalysisTaskMultiDielectron();
  AliAnalysisTaskMultiDielectron(const char *name);
  virtual ~AliAnalysisTaskMultiDielectron();

  enum ETriggerLogig {kAny, kExact};

  virtual void UserExec(Option_t *option);
  virtual void UserCreateOutputObjects();
  virtual void FinishTaskOutput();
  //temporary
//   virtual void NotifyRun(){AliDielectronPID::SetCorrVal((Double_t)fCurrentRunNumber);}
  
  void UsePhysicsSelection(Bool_t phy=kTRUE) {fSelectPhysics=phy;}
  void SetTriggerMask(ULong64_t mask) {fTriggerMask=mask;}
  UInt_t GetTriggerMask() const { return fTriggerMask; }
  void SetExcludeTriggerMask(ULong64_t mask) {fExcludeTriggerMask=mask;}
  UInt_t GetExcludeTriggerMask() const { return fExcludeTriggerMask; }
  void SetTriggerLogic(ETriggerLogig log) {fTriggerLogic=log;}
  ETriggerLogig GetTriggerLogic() const {return fTriggerLogic;}
  void SetFiredTriggerName(const char* select, Bool_t exclude=kFALSE) {fFiredTrigger=select; fFiredExclude=exclude;}

  void SetEventFilter(AliAnalysisCuts * const filter) {fEventFilter=filter;}
  void SetTriggerOnV0AND(Bool_t v0and=kTRUE)    { fTriggerOnV0AND=v0and;    }
  void SetRejectPileup(Bool_t pileup=kTRUE)     { fRejectPileup=pileup;     }
  void AddDielectron(AliDielectron * const die) { fListDielectron.Add(die); }
  void SetBeamEnergy(Double_t beamEbyHand=-1.)  { fBeamEnergy=beamEbyHand;  }
  
protected:
  enum {kAllEvents=0, kSelectedEvents, kV0andEvents, kFilteredEvents, kPileupEvents, kNbinsEvent};
  TList fListDielectron;             // List of dielectron framework instances
  TList fListHistos;                 //! List of histogram manager lists in the framework classes
  TList fListCF;                     //! List with CF Managers

  Bool_t fSelectPhysics;             // Whether to use physics selection
  UInt_t fTriggerMask;               // Event trigger mask
  UInt_t fExcludeTriggerMask;        // Triggers to exclude from the analysis
  Bool_t fTriggerOnV0AND;            // if to trigger on V0and
  TString fFiredTrigger;             // online trigger class name
  Bool_t fFiredExclude;              // cut logic: select(default) or exclude
  Bool_t fRejectPileup;              // pileup rejection wanted
  Double_t fBeamEnergy;              // beam energy in GeV (set by hand)

  ETriggerLogig fTriggerLogic;       // trigger logic: any or all bits need to be matching
  
  AliTriggerAnalysis *fTriggerAnalysis; //! trigger analysis class

  AliAnalysisCuts *fEventFilter;     // event filter
  
  TH1D *fEventStat;                  //! Histogram with event statistics
  
  AliAnalysisTaskMultiDielectron(const AliAnalysisTaskMultiDielectron &c);
  AliAnalysisTaskMultiDielectron& operator= (const AliAnalysisTaskMultiDielectron &c);
  
  ClassDef(AliAnalysisTaskMultiDielectron, 3); //Analysis Task handling multiple instances of AliDielectron
};
#endif
