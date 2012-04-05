#ifndef ALIANALYSISTASKDIELECTRONSE_H
#define ALIANALYSISTASKDIELECTRONSE_H
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

#include "AliAnalysisTaskSE.h"

class AliDielectron;
class AliTriggerAnalysis;
class TH1D;

class AliAnalysisTaskDielectronSE : public AliAnalysisTaskSE {
  
public:
  AliAnalysisTaskDielectronSE();
  AliAnalysisTaskDielectronSE(const char *name);
  virtual ~AliAnalysisTaskDielectronSE(){;}

  virtual void  UserExec(Option_t *option);
  virtual void  UserCreateOutputObjects();
  
  void UsePhysicsSelection(Bool_t phy=kTRUE) {fSelectPhysics=phy;}
  void SetTriggerMask(UInt_t mask) {fTriggerMask=mask;}
  UInt_t GetTriggerMask() const { return fTriggerMask; }

  void SetEventFilter(AliAnalysisCuts * const filter) {fEventFilter=filter;}
  void SetTriggerOnV0AND(Bool_t v0and=kTRUE)    { fTriggerOnV0AND=v0and;    }
  void SetRejectPileup(Bool_t pileup=kTRUE)     { fRejectPileup=pileup;     }
  
  void SetDielectron(AliDielectron * const die) { fDielectron = die; }
  
private:
  enum {kAllEvents=0, kSelectedEvents, kV0andEvents, kFilteredEvents, kPileupEvents, kNbinsEvent};
  
  AliDielectron *fDielectron;             // Dielectron framework object

  Bool_t fSelectPhysics;             // Whether to use physics selection
  UInt_t fTriggerMask;               // Event trigger mask
  Bool_t fTriggerOnV0AND;            // if to trigger on V0and
  Bool_t fRejectPileup;              // pileup rejection wanted

  AliTriggerAnalysis *fTriggerAnalysis; //! trigger analysis class

  AliAnalysisCuts *fEventFilter;     // event filter

  TH1D *fEventStat;                  //! Histogram with event statistics
  
  AliAnalysisTaskDielectronSE(const AliAnalysisTaskDielectronSE &c);
  AliAnalysisTaskDielectronSE& operator= (const AliAnalysisTaskDielectronSE &c);
  
  ClassDef(AliAnalysisTaskDielectronSE, 1);
};
#endif
