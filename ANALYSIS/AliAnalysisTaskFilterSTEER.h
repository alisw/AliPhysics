#ifndef ALIANALYSISTASKFILTERSTEER_H
#define ALIANALYSISTASKFILTERSTEER_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

//*************************************************************************
// Class AliAnalysisTaskFilterSTEER
// 
//*************************************************************************

#include "AliAnalysisTaskFilter.h"

class AliAnalysisTaskFilterSTEER : public AliAnalysisTaskFilter
{
 public:

  AliAnalysisTaskFilterSTEER();
  AliAnalysisTaskFilterSTEER(const char *name, Double_t ptCut, Double_t fractionCut);
  virtual ~AliAnalysisTaskFilterSTEER();

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual Bool_t UserSelectESDfriendForCurrentEvent();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

 private:

  AliAnalysisTaskFilterSTEER(const AliAnalysisTaskFilterSTEER &);
  AliAnalysisTaskFilterSTEER& operator=(const AliAnalysisTaskFilterSTEER&);

  Double_t      fPtCut;           // pt cut
  Double_t      fKeepFraction;    // fraction of tracks to keep
  //
  AliESDEvent  *fESDInput;        // ESD input object
  AliESDfriend *fESDfriendInput;  // ESD input friend object
  ClassDef(AliAnalysisTaskFilterSTEER,1); 
};

#endif

