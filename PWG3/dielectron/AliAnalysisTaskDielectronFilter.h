#ifndef ALIANALYSISTASKDIELECTRONFILTER_H
#define ALIANALYSISTASKDIELECTRONFILTER_H
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//#####################################################
//#                                                   # 
//#        Dielectron even filter task                #
//#                                                   #
//#                                                   #
//#  by WooJin J. Park, GSI / W.J.Park@gsi.de         #
//#     Ionut C. Arsene, GSI / I.C.Arsene@gsi.de      #
//#     Magnus Mager, CERN / Magnus.Mager@cern.ch     #
//#     Jens Wiechula, Uni HD / Jens.Wiechula@cern.ch #
//#                                                   #
//#####################################################
/*
Filter Event based on cuts provided in the AliDielectron class.

Write an AOD file containing events with Dielectron candidates.
Add a sattelite AOD with the array of candidates.
*/



#include "AliAnalysisTaskSE.h"

class AliDielectron;

class AliAnalysisTaskDielectronFilter : public AliAnalysisTaskSE {
  
public:
  AliAnalysisTaskDielectronFilter();
  AliAnalysisTaskDielectronFilter(const char *name);
  virtual ~AliAnalysisTaskDielectronFilter(){}

  virtual void UserExec(Option_t *option);
  virtual void Init();
  virtual void LocalInit() {Init();}

  void UsePhysicsSelection(Bool_t phy=kTRUE) {fSelectPhysics=phy;}
  
  void SetDielectron(AliDielectron * const die) { fDielectron = die; }
  
private:
  
  AliDielectron *fDielectron;             // J/psi framework object

  Bool_t fSelectPhysics;                  // Whether to use physics selection
  
  AliAnalysisTaskDielectronFilter(const AliAnalysisTaskDielectronFilter &c);
  AliAnalysisTaskDielectronFilter& operator= (const AliAnalysisTaskDielectronFilter &c);
  
  ClassDef(AliAnalysisTaskDielectronFilter, 1);
};
#endif
