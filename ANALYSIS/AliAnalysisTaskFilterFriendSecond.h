#ifndef ALIANALYSISTASKFILTERFRIENDSECOND_H
#define ALIANALYSISTASKFILTERFRIENDSECOND_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

//*************************************************************************
// Class AliAnalysisTaskFilterFriendSecond
// Test Task 
//*************************************************************************

#include <TString.h>
#include "AliAnalysisTaskFilter.h"

class AliAnalysisTaskFilterFriendSecond : public AliAnalysisTaskFilter
{
 public:

  AliAnalysisTaskFilterFriendSecond();
  AliAnalysisTaskFilterFriendSecond(const char *name);
  virtual ~AliAnalysisTaskFilterFriendSecond();


  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual Bool_t UserSelectESDfriendForCurrentEvent();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

 private:

  AliAnalysisTaskFilterFriendSecond(const AliAnalysisTaskFilterFriendSecond &);
  AliAnalysisTaskFilterFriendSecond& operator=(const AliAnalysisTaskFilterFriendSecond&);

  AliESDEvent  *fESDInput;        // ESD input object
  AliESDfriend *fESDfriendInput;  // ESD input friend object
  ClassDef(AliAnalysisTaskFilterFriendSecond,1); 
};

#endif

