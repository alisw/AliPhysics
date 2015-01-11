#ifndef ALIANALYSISTASKFILTERSTEER_H
#define ALIANALYSISTASKFILTERSTEER_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

//*************************************************************************
// Class AliAnalysisTaskFilterSteer
// Implementing the filtering of the friend in order to reduce them to 1%
// 
//*************************************************************************

#include "AliAnalysisTaskFilter.h"

class AliAnalysisTaskFilterSteer : public AliAnalysisTaskFilter
{
 public:

  AliAnalysisTaskFilterSteer();
  AliAnalysisTaskFilterSteer(const char *name);
  virtual ~AliAnalysisTaskFilterSteer();

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual Bool_t UserSelectESDfriendForCurrentEvent();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
  void SetFraction(Double_t fraction) {fFraction = fraction;}
  Double_t GetFraction() const {return fFraction;}

 private:

  AliAnalysisTaskFilterSteer(const AliAnalysisTaskFilterSteer &);
  AliAnalysisTaskFilterSteer& operator=(const AliAnalysisTaskFilterSteer&);

  Double_t      fFraction;    // fraction of events for which to keep the friends
  //
  AliESDEvent  *fESDInput;        // ESD input object
  AliESDfriend *fESDfriendInput;  // ESD input friend object
  ClassDef(AliAnalysisTaskFilterSteer,1); 
};

#endif

