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

class AliDielectron;

class AliAnalysisTaskMultiDielectron : public AliAnalysisTaskSE {
  
public:
  AliAnalysisTaskMultiDielectron();
  AliAnalysisTaskMultiDielectron(const char *name);
  virtual ~AliAnalysisTaskMultiDielectron(){  }

  virtual void UserExec(Option_t *option);
  virtual void UserCreateOutputObjects();
  virtual void FinishTaskOutput();
  
  void UsePhysicsSelection(Bool_t phy=kTRUE) {fSelectPhysics=phy;}
  
  void AddDielectron(AliDielectron * const die) { fListDielectron.Add(die); }
  
private:
  
  TList fListDielectron;             // List of dielectron framework instances
  TList fListHistos;                 //! List of histogram manager lists in the framework classes
  TList fListCF;                     //! List with CF Managers

  Bool_t fSelectPhysics;             // Whether to use physics selection
  
  AliAnalysisTaskMultiDielectron(const AliAnalysisTaskMultiDielectron &c);
  AliAnalysisTaskMultiDielectron& operator= (const AliAnalysisTaskMultiDielectron &c);
  
  ClassDef(AliAnalysisTaskMultiDielectron, 1); //Analysis Task handling multiple instances of AliDielectron
};
#endif
