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
//     Event handler for multiple VEvent input.
//     This class handles multiple inputs for event mixing. 
//     Author: Andreas Morsch, CERN
//-------------------------------------------------------------------------

#include "AliMultiEventInputHandler.h"
#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliVEventPool.h"
#include "AliVCuts.h"
#include "AliLog.h"
#include <TObjArray.h>
#include <TTree.h>
#include <TList.h>
#include <TEntryList.h>


ClassImp(AliMultiEventInputHandler)

AliMultiEventInputHandler::AliMultiEventInputHandler() :
    AliInputEventHandler(),
    fBufferSize(0),
    fFormat(1),
    fNBuffered(0),
    fIndex(0),
    fCurrentBin(0),
    fCurrentEvt(0),
    fInit(0),
    fEventPool(0),
    fEventBuffer(0),
    fEventSkipped(0)
{
  // Default constructor
}

//______________________________________________________________________________
AliMultiEventInputHandler::AliMultiEventInputHandler(Int_t size, Int_t format) :
    AliInputEventHandler(),
    fBufferSize(size),
    fFormat(format),
    fNBuffered(0),
    fIndex(0),
    fCurrentBin(0),
    fCurrentEvt(0),
    fInit(0),
    fEventPool(0),
    fEventBuffer(0),
    fEventSkipped(0)
{
  // constructor
}

//______________________________________________________________________________
AliMultiEventInputHandler::AliMultiEventInputHandler(const char* name, const char* title, Int_t size, Int_t format):
    AliInputEventHandler(name, title),
    fBufferSize(size),
    fFormat(format),
    fNBuffered(0),
    fIndex(0),
    fCurrentBin(0),
    fCurrentEvt(0),
    fInit(0),
    fEventPool(0),
    fEventBuffer(0),
    fEventSkipped(0)
{
    // Constructor

}

//______________________________________________________________________________
AliMultiEventInputHandler::~AliMultiEventInputHandler() 
{
// Destructor
}

Bool_t AliMultiEventInputHandler::Init(TTree* tree, Option_t* /*opt*/)
{
    // Initialisation necessary for each new tree
    if (!fEventBuffer) {
	fEventBuffer = new AliVEvent*[fBufferSize];
	
	for (Int_t i = 0; i < fBufferSize; i++) 
	    if (fFormat == 1) {
		fEventBuffer[i] = new AliAODEvent();
	    } else if (fFormat == 0) {
		fEventBuffer[i] = new AliESDEvent();
	    } else{
		AliWarning(Form("Unknown Format %5d", fFormat));
	    }
    }
    

    fTree = tree;
    fInit = 1;
    
    if (!fTree) return kFALSE;
    for (Int_t i = 0; i < fBufferSize; i++) 
	fEventBuffer[i]->Clear();
    fIndex     = 0;
    fNBuffered = 1;
    return kTRUE;
}


Bool_t AliMultiEventInputHandler::Notify(const char */*path*/)
{
    // Connect to new tree

    TList* connectedList = (TList*) (fTree->GetUserInfo()->FindObject("AODObjectsConnectedToTree"));   
    if (connectedList && !fInit) {
	fEventBuffer[0]->ReadFromTree(fTree, "reconnect");
    } else {
	if (fInit) fEventBuffer[0]->ReadFromTree(fTree, "");
    }
    
    fCurrentEvt = 0;
    fInit = 0;
    
    return (kTRUE);
}

Bool_t AliMultiEventInputHandler::BeginEvent(Long64_t /*entry*/)
{
    // Actions before analysis of each event 
    //
    // Reset the number of events buffered for this bin to 0
    
    if (fCurrentBin != (fEventPool->BinNumber())) {
	fCurrentBin = fEventPool->BinNumber();
	fNBuffered = 0;
    }
  //
  // Event selection
  // 
    if (fFormat == 0) {
      fIsSelectedResult = 0;
      if (fEventCuts && !IsUserCallSelectionMask())
	fIsSelectedResult = 
	  fEventCuts->GetSelectionMask((AliESDEvent*)fEventBuffer[fIndex]); 
    }
    
    return kTRUE;
}

Bool_t AliMultiEventInputHandler::FinishEvent()
{
    // 
    // Connect the next event in the buffer to the tree
    if (!fEventSkipped) fIndex++;
    fIndex %= fBufferSize;
    AliInfo(Form("Connecting buffer entry %5d", fIndex));
    fEventBuffer[fIndex]->Clear();
    fCurrentEvt++;
    if (fEventBuffer[fIndex]->GetList() && fCurrentEvt > (fBufferSize - 1))
	fEventBuffer[fIndex]->GetList()->Delete();

    fEventBuffer[fIndex]->ReadFromTree(fTree, "reconnect");

    fNBuffered++;
    if (fNBuffered > fBufferSize) fNBuffered = fBufferSize;
    
    Int_t nmax = fTree->GetEntries();
    if (fTree->GetEntryList()) {
	nmax = (fTree->GetEntryList()->GetN());
    } else {
	if (fTree->GetTree()) nmax = fTree->GetTree()->GetEntries();
    }
    
    if (fCurrentEvt == nmax)
    {
	for (Int_t i = 0; i < fBufferSize; i++) {
	    fEventBuffer[i]->Clear();
	}
    }
    
    return (kTRUE);
}

AliVEvent* AliMultiEventInputHandler::GetEvent(Int_t iev) const
{
    // Get event number iev from buffer
    if ((iev < 0) || (iev >= fBufferSize))
    {
	AliWarning(Form("Event number out of range: %10d", iev));
	return 0;
    }
	
    iev = fIndex - (fBufferSize - 1 - iev);
    if (iev < 0) iev += fBufferSize;
    AliInfo(Form("Event index in buffer is %5d", iev));
    return (fEventBuffer[iev]);
}

