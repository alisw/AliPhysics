#ifndef ALIANALYSISTASKDIELECTRONME_H
#define ALIANALYSISTASKDIELECTRONME_H
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

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

#include "AliAnalysisTaskME.h"

#include "AliDielectronPID.h"

class AliDielectron;
class TH1D;

class AliAnalysisTaskDielectronME : public AliAnalysisTaskME {
  
public:
  AliAnalysisTaskDielectronME();
  AliAnalysisTaskDielectronME(const char *name);
  virtual ~AliAnalysisTaskDielectronME(){  }

  virtual void UserExec(Option_t *option);
  virtual void UserCreateOutputObjects();
  virtual void FinishTaskOutput();
  //temporary
  //virtual void NotifyRun(){AliDielectronPID::SetCorrVal((Double_t)fCurrentRunNumber);}
  virtual void NotifyRun(){AliDielectronPID::SetCorrVal((Double_t)GetEvent(0)->GetRunNumber());}
  
  void UsePhysicsSelection(Bool_t phy=kTRUE) {fSelectPhysics=phy;}
  void SetTriggerMask(UInt_t mask) {fTriggerMask=mask;}
  UInt_t GetTriggerMask() const { return fTriggerMask; }
  void SetPoolDepth(Int_t depth=2){fPoolDepth=depth;}
  
  void AddDielectron(AliDielectron * const die) { fListDielectron.Add(die); }
  
private:
  
  TList fListDielectron;             // List of dielectron framework instances
  TList fListHistos;                 //! List of histogram manager lists in the framework classes
  TList fListCF;                     //! List with CF Managers

  Int_t fPoolDepth;                  // Pool depth for event mixing
  Bool_t fSelectPhysics;             // Whether to use physics selection
  UInt_t fTriggerMask;               // Event trigger mask

  TH1D *fEventStat;                  //! Histogram with event statistics
  
  AliAnalysisTaskDielectronME(const AliAnalysisTaskDielectronME &c);
  AliAnalysisTaskDielectronME& operator= (const AliAnalysisTaskDielectronME &c);
  
  ClassDef(AliAnalysisTaskDielectronME, 1); //Analysis Task handling multiple instances of AliDielectron
};
#endif
