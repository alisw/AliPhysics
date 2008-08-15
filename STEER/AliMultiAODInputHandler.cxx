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
//     Event handler for AOD input.
//     This class handles multiple inputs for event mixing. 
//     Author: Andreas Morsch, CERN
//-------------------------------------------------------------------------

#include "AliMultiAODInputHandler.h"
#include "AliAODEvent.h"
#include "AliVEventPool.h"
#include "AliLog.h"
#include <TObjArray.h>
#include <TTree.h>


ClassImp(AliMultiAODInputHandler)

AliMultiAODInputHandler::AliMultiAODInputHandler() :
    AliInputEventHandler(),
    fBufferSize(0),
    fNBuffered(0),
    fIndex(0),
    fCurrentBin(0),
    fTree(0),
    fEventPool(0),
    fEventBuffer(0)
{
  // Default constructor
}

//______________________________________________________________________________
AliMultiAODInputHandler::AliMultiAODInputHandler(Int_t size) :
    AliInputEventHandler(),
    fBufferSize(size),
    fNBuffered(0),
    fIndex(0),
    fCurrentBin(0),
    fTree(0),
    fEventPool(0),
    fEventBuffer(new AliAODEvent*[size])
{
  // Default constructor
    for (Int_t i = 0; i < size; i++) 
	fEventBuffer[i] = new AliAODEvent();
}

//______________________________________________________________________________
AliMultiAODInputHandler::AliMultiAODInputHandler(const char* name, const char* title, Int_t size):
    AliInputEventHandler(name, title),
    fBufferSize(size),
    fNBuffered(0),
    fIndex(0),
    fCurrentBin(0),
    fTree(0),
    fEventPool(0),
    fEventBuffer(new AliAODEvent*[size])
{
    // Constructor
    for (Int_t i = 0; i < size; i++) 
	fEventBuffer[i] = new AliAODEvent();
}

//______________________________________________________________________________
AliMultiAODInputHandler::~AliMultiAODInputHandler() 
{
// Destructor
}

Bool_t AliMultiAODInputHandler::Init(TTree* tree, Option_t* /*opt*/)
{
    // Initialisation necessary for each new tree
    fTree = tree;
    if (!fTree) return kFALSE;
    // Get pointer to AOD event
    fEventBuffer[0]->ReadFromTree(fTree);
    fIndex     = 0;
    fNBuffered = 1;
    return kTRUE;
}

Bool_t AliMultiAODInputHandler::BeginEvent(Long64_t /*entry*/)
{
    // Actions before analysis of each event 
    //
    // Reset the number of events buffered for this bin to 0
    if (fCurrentBin != (fEventPool->BinNumber())) {
	fCurrentBin = fEventPool->BinNumber();
	fNBuffered = 0;
    }
    return kTRUE;
}

Bool_t AliMultiAODInputHandler::FinishEvent()
{
    // 
    // Connect the next event in the buffer to the tree
    fIndex++;
    
    fIndex %= fBufferSize;
    AliInfo(Form("Connecting buffer entry %5d", fIndex));

    fEventBuffer[fIndex]->ReadFromTree(fTree, "reconnect");

    fNBuffered++;
    if (fNBuffered > fBufferSize) fNBuffered = fBufferSize;

    return (kTRUE);
}


AliAODEvent* AliMultiAODInputHandler::GetEvent(Int_t iev) const
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

