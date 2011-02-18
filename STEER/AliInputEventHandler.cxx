/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//-------------------------------------------------------------------------
//     Event handler for event input 
//     Author: Andreas Morsch, CERN
//-------------------------------------------------------------------------


#include "AliInputEventHandler.h"
#include "AliVEvent.h"
#include "AliVCuts.h"
#include "AliLog.h"


ClassImp(AliInputEventHandler)

//______________________________________________________________________________
AliInputEventHandler::AliInputEventHandler() :
    AliVEventHandler(),
    fTree(0),
    fBranches(""),
    fBranchesOn(""),
    fNewEvent(kTRUE),
    fEventCuts(0),
    fIsSelectedResult(0),
    fMixingHandler(0),
    fParentHandler(0)
{
  // default constructor
}

//______________________________________________________________________________
AliInputEventHandler::~AliInputEventHandler() 
{
// destructor
}

//______________________________________________________________________________
AliInputEventHandler::AliInputEventHandler(const char* name, const char* title):
  AliVEventHandler(name, title),
  fTree(0),
  fBranches(""),
  fBranchesOn(""),
  fNewEvent(kTRUE),
  fEventCuts(0),
  fIsSelectedResult(0),
  fMixingHandler(0),
  fParentHandler(0)
{
// Named constructor.
}

void AliInputEventHandler::SwitchOffBranches() const {
  //
  // Switch of branches on user request
    TObjArray * tokens = fBranches.Tokenize(" ");
    Int_t ntok = tokens->GetEntries();
    for (Int_t i = 0; i < ntok; i++)  {
	TString str = ((TObjString*) tokens->At(i))->GetString();
	if (str.Length() == 0)
	    continue;
	fTree->SetBranchStatus(Form("%s%s%s","*", str.Data(), "*"), 0);
	AliDebug(1,Form("Branch %s switched off", str.Data()));
    }
  delete tokens;
}

void AliInputEventHandler::SwitchOnBranches() const {
  //
  // Switch of branches on user request
  TObjArray * tokens = fBranchesOn.Tokenize(" ");
  Int_t ntok = tokens->GetEntries();

  for (Int_t i = 0; i < ntok; i++)  {
      TString str = ((TObjString*) tokens->At(i))->GetString();
      if (str.Length() == 0)
	  continue;
      fTree->SetBranchStatus(Form("%s%s%s","*", str.Data(), "*"), 1);
      AliDebug(1,Form("Branch %s switched on", str.Data()));
  }
  delete tokens;
}

TObject *AliInputEventHandler::GetStatistics(Option_t *) const
{
// Returns the statistics object(s) (TH2F histogram) produced by the physics
// selection. Implementations both for ESD and AOD input handlers.
  return NULL;
}
   
Long64_t AliInputEventHandler::GetReadEntry() const 
{
  // Get the current entry.
  return fTree->GetReadEntry();
}
