#ifndef ALIANALYSISTASKFILTERFRIEND_H
#define ALIANALYSISTASKFILTERFRIEND_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

//*************************************************************************
// Class AliAnalysisTaskFilterFriend
// Test Task 
//*************************************************************************

#include "AliAnalysisTaskFilter.h"

class AliAnalysisTaskFilterFriend : public AliAnalysisTaskFilter
{
 public:

  AliAnalysisTaskFilterFriend();
  AliAnalysisTaskFilterFriend(const char *name);
  virtual ~AliAnalysisTaskFilterFriend();

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual Bool_t UserSelectESDfriendForCurrentEvent();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

 private:

  AliAnalysisTaskFilterFriend(const AliAnalysisTaskFilterFriend &);
  AliAnalysisTaskFilterFriend& operator=(const AliAnalysisTaskFilterFriend&);

  AliESDEvent  *fESDInput;        // ESD input object
  AliESDfriend *fESDfriendInput;  // ESD input friend object
  ClassDef(AliAnalysisTaskFilterFriend,1); 
};

#endif

