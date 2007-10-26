/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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
//---------------------------------------------------------------------------------
//                          Class AliMCEventHandler
// This class gives access to MC truth during the analysis.
// Monte Carlo truth is containe in the kinematics tree (produced particles) and 
// the tree of reference hits.
//      
// Origin: Andreas Morsch, CERN, andreas.morsch@cern.ch 
//---------------------------------------------------------------------------------



#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliTrackReference.h"
#include "AliHeader.h"
#include "AliStack.h"
#include "AliLog.h"

#include <TTree.h>
#include <TFile.h>
#include <TParticle.h>
#include <TString.h>
#include <TClonesArray.h>
#include <TDirectoryFile.h>

ClassImp(AliMCEventHandler)

AliMCEventHandler::AliMCEventHandler() :
    AliVEventHandler(),
    fMCEvent(new AliMCEvent()),
    fFileE(0),
    fFileK(0),
    fFileTR(0),
    fTreeE(0),
    fTreeK(0),
    fTreeTR(0),
    fDirK(0),
    fDirTR(0),
    fNEvent(-1),
    fEvent(-1),
    fPathName(new TString("./")),
    fExtension(""),
    fFileNumber(0),
    fEventsPerFile(0)
{
    // Default constructor
}

AliMCEventHandler::AliMCEventHandler(const char* name, const char* title) :
    AliVEventHandler(name, title),
    fMCEvent(new AliMCEvent()),
    fFileE(0),
    fFileK(0),
    fFileTR(0),
    fTreeE(0),
    fTreeK(0),
    fTreeTR(0),
    fDirK(0),
    fDirTR(0),
    fNEvent(-1),
    fEvent(-1),
    fPathName(new TString("./")),
    fExtension(""),
    fFileNumber(0),
    fEventsPerFile(0)
{
    // Constructor
}
AliMCEventHandler::~AliMCEventHandler()
{ 
    // Destructor
    delete fMCEvent;
    delete fFileE;
    delete fFileK;
    delete fFileTR;
}

Bool_t AliMCEventHandler::InitIO(Option_t* opt) 
{ 
    // Initialize input
    //
    if (!(strcmp(opt, "proof")) || !(strcmp(opt, "local"))) return kTRUE;
    //

    fFileE = TFile::Open(Form("%sgalice.root", fPathName->Data()));
    if (!fFileE) AliFatal(Form("AliMCEventHandler:galice.root not found in directory %s ! \n", fPathName->Data()));

    //
    // Tree E
    fFileE->GetObject("TE", fTreeE);
    // Connect Tree E to the MCEvent
    fMCEvent->ConnectTreeE(fTreeE);
    fNEvent = fTreeE->GetEntries();
    //
    // Tree K
    fFileK = TFile::Open(Form("%sKinematics%s.root", fPathName->Data(), fExtension));
    if (!fFileK) AliFatal(Form("AliMCEventHandler:Kinematics.root not found in directory %s ! \n", fPathName));
    fEventsPerFile = fFileK->GetNkeys() - fFileK->GetNProcessIDs();
    //
    // Tree TR
    fFileTR = TFile::Open(Form("%sTrackRefs%s.root", fPathName->Data(), fExtension));
    if (!fFileTR) AliWarning(Form("AliMCEventHandler:TrackRefs.root not found in directory %s ! \n", fPathName->Data()));
    //
    // Reset the event number
    fEvent      = -1;
    fFileNumber =  0;
    
    AliInfo(Form("AliMCEventHandler:Number of events in this directory %5d \n", fNEvent));
    return kTRUE;
}

Bool_t AliMCEventHandler::GetEvent(Int_t iev)
{
    // Load the event number iev
    //
    // Calculate the file number
    Int_t inew  = iev / fEventsPerFile;
    if (inew != fFileNumber) {
	fFileNumber = inew;
	if (!OpenFile(fFileNumber)){
	    return kFALSE;
	}
    }
    // Folder name
    char folder[20];
    sprintf(folder, "Event%d", iev);
    // TreeE
    fTreeE->GetEntry(iev);
    // Tree K
    fFileK->GetObject(folder, fDirK);
    if (!fDirK) {
	AliWarning(Form("AliMCEventHandler: Event #%5d not found\n", iev));
	return kFALSE;
    }
    fDirK ->GetObject("TreeK", fTreeK);
    // Connect TreeK to MCEvent
    fMCEvent->ConnectTreeK(fTreeK);
    //Tree TR 
    if (fFileTR) {
	// Check which format has been read
	fFileTR->GetObject(folder, fDirTR);
	fDirTR->GetObject("TreeTR", fTreeTR);
	//
	// Connect TR to MCEvent
	fMCEvent->ConnectTreeTR(fTreeTR);
    }
    //
    return kTRUE;
}

Bool_t AliMCEventHandler::OpenFile(Int_t i)
{
    // Open file i
    Bool_t ok = kTRUE;
    if (i > 0) {
	fExtension = Form("%d", i);
    } else {
	fExtension = "";
    }
    
    
    delete fFileK;
    fFileK = TFile::Open(Form("%sKinematics%s.root", fPathName->Data(), fExtension));
    if (!fFileK) {
	AliFatal(Form("AliMCEventHandler:Kinematics%s.root not found in directory %s ! \n", fExtension, fPathName->Data()));
	ok = kFALSE;
    }
    
    delete fFileTR;
    fFileTR = TFile::Open(Form("%sTrackRefs%s.root", fPathName->Data(), fExtension));
    if (!fFileTR) {
	AliWarning(Form("AliMCEventHandler:TrackRefs%s.root not found in directory %s ! \n", fExtension, fPathName->Data()));
	ok = kFALSE;
    }
    
    return ok;
}

Bool_t AliMCEventHandler::BeginEvent()
{ 
    // Read the next event
    fEvent++;
    if (fEvent >= fNEvent) {
	AliWarning(Form("AliMCEventHandler: Event number out of range %5d\n", fEvent));
	return kFALSE;
    }
    return GetEvent(fEvent);
}

Int_t AliMCEventHandler::GetParticleAndTR(Int_t i, TParticle*& particle, TClonesArray*& trefs)
{
    // Retrieve entry i
    return (fMCEvent->GetParticleAndTR(i, particle, trefs));
}

void AliMCEventHandler::DrawCheck(Int_t i, Int_t search)
{
    // Retrieve entry i and draw momentum vector and hits
    fMCEvent->DrawCheck(i, search);
}

Bool_t AliMCEventHandler::Notify(const char *path)
{
  // Notify about directory change
  // The directory is taken from the 'path' argument
  // Reconnect trees
    TString fileName(path);
    if(fileName.Contains("AliESDs.root")){
	fileName.ReplaceAll("AliESDs.root", "");
    }
    else if(fileName.Contains("galice.root")){
	// for running with galice and kinematics alone...
	fileName.ReplaceAll("galice.root", "");
    }
    
    *fPathName = fileName;
    printf("AliMCEventHandler::Notify() Path: %s\n", fPathName->Data());
    
    ResetIO();
    InitIO("");
    return kTRUE;
}
    
void AliMCEventHandler::ResetIO()
{
//  Clear header and stack
    fMCEvent->Clean();
    
// Delete Tree E    
    delete fTreeE; fTreeE = 0;
    
// Reset files
    if (fFileE)  delete fFileE;
    if (fFileK)  delete fFileK;
    if (fFileTR) delete fFileTR;
}

			    
Bool_t AliMCEventHandler::FinishEvent()
{
    // Clean-up after each event
    delete fDirTR;  fDirTR = 0;
    delete fDirK;   fDirK  = 0;    
    fMCEvent->FinishEvent();
    return kTRUE;
}

Bool_t AliMCEventHandler::Terminate()
{ 
    // Dummy 
    return kTRUE;
}

Bool_t AliMCEventHandler::TerminateIO()
{ 
    // Dummy
    return kTRUE;
}
    

void AliMCEventHandler::SetInputPath(char* fname)
{
    // Set the input path name
    delete fPathName;
    fPathName = new TString(fname);
}
