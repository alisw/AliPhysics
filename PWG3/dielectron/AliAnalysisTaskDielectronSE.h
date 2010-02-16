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

class AliAnalysisTaskDielectronSE : public AliAnalysisTaskSE {
  
public:
  AliAnalysisTaskDielectronSE();
  AliAnalysisTaskDielectronSE(const char *name);
  virtual ~AliAnalysisTaskDielectronSE(){;}

  virtual void  UserExec(Option_t *option);
  virtual void  UserCreateOutputObjects();
  
  void SetDielectron(AliDielectron * const die) { fDielectron = die; }
  
private:
  
  AliDielectron *fDielectron;             // Dielectron framework object

  AliAnalysisTaskDielectronSE(const AliAnalysisTaskDielectronSE &c);
  AliAnalysisTaskDielectronSE& operator= (const AliAnalysisTaskDielectronSE &c);
  
  ClassDef(AliAnalysisTaskDielectronSE, 1);
};
#endif
