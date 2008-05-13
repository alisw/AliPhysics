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
 
//-------------------------------------------------------------------------
//  Class AliRsnSelectorRL
// ------------------------
// Reader for conversion of ESD output into the internal format
// used for resonance study.
// ---
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
// ---
// adapted for TSelector compliance
// by    : R. Vernet                 (email: renaud.vernet@cern.ch)
//-------------------------------------------------------------------------

#include <Riostream.h>

#include <TFile.h>
#include <TChain.h>
#include <TParticle.h>
#include <TRandom.h>
#include <TObjString.h>
#include <TObjectTable.h>
#include <TOrdCollection.h>
#include "AliRun.h"
#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliStack.h"
#include "AliHeader.h"
#include "AliESDtrack.h"
#include "AliRunLoader.h"
#include "AliRsnDaughter.h"
#include "AliRsnEvent.h"
#include "AliRsnSelectorRL.h"

ClassImp(AliRsnSelectorRL)

//--------------------------------------------------------------------------------------------------------
AliRsnSelectorRL::AliRsnSelectorRL(TTree*) :
  AliSelectorRL(),
  AliRsnReader(),
  fOutputPath(0),
  fIsRunLoaderOpen(0),
  fRsnEventTree(0),
  fRsnEvent(0),
  fRsnEventBranch(0)
{
//
// Constructor.
// Initializes all working parameters to default values:
// - ESD particle identification
// - rejection of non-ITS-refitted tracks
// - maximum distance allowed from primary vertex = 3.0 cm (beam pipe)
//
}

AliRsnSelectorRL::~AliRsnSelectorRL() {
  Clear();
}

//--------------------------------------------------------------------------------------------------------
AliRsnSelectorRL::AliRsnSelectorRL(const AliRsnSelectorRL& obj) :
  AliSelectorRL(), // not implemented a copy constructor for AliRsnSelectorRL
  AliRsnReader(obj),
  fOutputPath(obj.fOutputPath),
  fIsRunLoaderOpen(obj.fIsRunLoaderOpen),
  fRsnEventTree(obj.fRsnEventTree),
  fRsnEvent(obj.fRsnEvent),
  fRsnEventBranch(obj.fRsnEventBranch)
{
//
// Copy constructor
//
}
//--------------------------------------------------------------------------------------------------------
AliRsnSelectorRL& AliRsnSelectorRL::operator=(const AliRsnSelectorRL& obj) 
{
//
// Assignment operator
// works in the same way as the copy constructor
//
	if (this!=&obj) {
		fOutputPath = obj.fOutputPath;
		fIsRunLoaderOpen = obj.fIsRunLoaderOpen;
		fRsnEventTree = obj.fRsnEventTree;
		fRsnEvent = obj.fRsnEvent;
		fRsnEventBranch = obj.fRsnEventBranch;
	}
	
	return *this;
}
//--------------------------------------------------------------------------------------------------------
void AliRsnSelectorRL::Clear(Option_t * /*option*/)
{
  //
  // Does nothing.
  //
}
//--------------------------------------------------------------------------------------------------------

//--------------------------------------------------------
//             The following is Selector stuff
//--------------------------------------------------------

//--------------------------------------------------------------------------------------------------------
void AliRsnSelectorRL::Begin(TTree *) const
{
//
// Implementation of BEGIN method
//
	Info("Begin", "");
	TString option = GetOption();
}
//--------------------------------------------------------------------------------------------------------
void AliRsnSelectorRL::SlaveBegin(TTree * tree)
{
//
// Implementation of secondary BEGIN which 
// is called by separate process managers
//
	Info("SlaveBegin", "");
	Init(tree);
	TString option = GetOption();
}
//--------------------------------------------------------------------------------------------------------
void AliRsnSelectorRL::Init(TTree *tree)
{
//
// Initializer
// Connects the selector to a TTree and links the branch
// which is used to make translation ESD --> Rsn
//
	Info("Init","");
	
	if (!tree) return;
	fTree = tree;
	if (fDebugFlag) {
		Info("Init", "fTree=%p   fTree->GetCurrentFile()=%p", fTree, fTree->GetCurrentFile());
	}
	fTree->SetBranchAddress("ESD", &fESD);
	fRsnEventTree = new TTree("selection", "AliRsnEvents");
	TTree::SetBranchStyle(1);
	fRsnEvent = new AliRsnEvent;
	fRsnEventBranch = fRsnEventTree->Branch("RsnEvents", "AliRsnEvent", &fRsnEvent, 32000, 1);
	fRsnEventBranch->SetAutoDelete(kFALSE);
}
//--------------------------------------------------------------------------------------------------------
Bool_t AliRsnSelectorRL::Process(Long64_t entry)
{
//
// Main core of the Selector processing
// Reads the ESD input and creates a TTree or AliRsnEvent objects
//
	if (fDebugFlag) Info("Process", "Processing event %d", entry);
	if (!AliSelectorRL::Process(entry)) return kFALSE;
	
	AliStack* stack = GetStack();
	if (!stack) {
		Warning("Process", "NULL stack: cannot get kinematics info");
	}
	
	/*
	AliESDEvent *esdEvent = new AliESDEvent(fESD);
	AliRsnEvent *event = ReadESDEvent(esdEvent);
	AddMCInfo(event, stack);
	AddEvent(fRsnTree, event);
	*/
	return kTRUE;
}
//--------------------------------------------------------------------------------------------------------
void AliRsnSelectorRL::SlaveTerminate()
{
//
// SlaveTerminate
// Partial termination method
//
// test have been performed only on AliEn, please contact the author in case of merging problem under CAF/PROOF.
//
//
	Info("SlaveTerminate", "");
	
	// Add the histograms to the output on each slave server
	fOutput->Add(fRsnEventTree);
}
//--------------------------------------------------------------------------------------------------------
void AliRsnSelectorRL::Terminate()
{
//
// Global termination method
//
	Info("Terminate","");
	// fRsnEventTree = dynamic_cast<TTree*>(fOutput->FindObject("fRsnEventTree"));
	
	//AliSelector::Terminate();
	cout << fOutputPath << endl;
	Info("Terminate", Form("Saving in: %s", fOutputPath->Data()));
	TFile* file = TFile::Open(fOutputPath->Data(), "RECREATE");
	fRsnEventTree->Write();
	file->Close();
	
	delete fRsnEventTree;
	delete file;
	delete fOutputPath;
}
