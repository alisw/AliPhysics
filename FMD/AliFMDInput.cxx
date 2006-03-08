/**************************************************************************
 * Copyright(c) 2004, ALICE Experiment at CERN, All rights reserved. *
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

//___________________________________________________________________
//
// The classes defined here, are utility classes for reading in data
// for the FMD.  They are  put in a seperate library to not polute the
// normal libraries.  The classes are intended to be used as base
// classes for customized class that do some sort of analysis on the
// various types of data produced by the FMD. 
//
// Latest changes by Christian Holm Christensen
//
#include "AliFMDInput.h"	// ALIFMDHIT_H
#include "AliLog.h"		// ALILOG_H
#include "AliLoader.h"          // ALILOADER_H
#include "AliRunLoader.h"       // ALIRUNLOADER_H
#include "AliRun.h"             // ALIRUN_H
#include "AliStack.h"           // ALISTACK_H
#include "AliFMD.h"             // ALIFMD_H
#include "AliFMDHit.h"		// ALIFMDHIT_H
#include "AliFMDDigit.h"	// ALIFMDDigit_H
#include "AliFMDMultStrip.h"	// ALIFMDMultStrip_H
#include "AliFMDMultRegion.h"	// ALIFMDMultRegion_H
#include <TTree.h>              // ROOT_TTree
#include <TParticle.h>          // ROOT_TParticle
#include <TString.h>            // ROOT_TString
#include <TDatabasePDG.h>       // ROOT_TDatabasePDG
#include <TMath.h>              // ROOT_TMath
#include <TGeoManager.h>        // ROOT_TGeoManager 
#include <Riostream.h>		// ROOT_Riostream

//____________________________________________________________________
ClassImp(AliFMDInput)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif


//____________________________________________________________________
AliFMDInput::AliFMDInput()
  : fGAliceFile(""), 
    fLoader(0),
    fRun(0), 
    fStack(0),
    fFMDLoader(0), 
    fFMD(0),
    fTreeE(0),
    fTreeH(0),
    fTreeD(0),
    fTreeS(0),
    fTreeR(0), 
    fArrayE(0),
    fArrayH(0),
    fArrayD(0),
    fArrayS(0), 
    fArrayN(0), 
    fArrayP(0), 
    fTreeMask(0), 
    fIsInit(kFALSE)
{
  // Constructor of an FMD input object.  Specify what data to read in
  // using the AddLoad member function.   Sub-classes should at a
  // minimum overload the member function Event.   A full job can be
  // executed using the member function Run. 
}

  

//____________________________________________________________________
AliFMDInput::AliFMDInput(const char* gAliceFile)
  : fGAliceFile(gAliceFile),
    fLoader(0),
    fRun(0), 
    fStack(0),
    fFMDLoader(0), 
    fFMD(0),
    fTreeE(0),
    fTreeH(0),
    fTreeD(0),
    fTreeS(0),
    fTreeR(0), 
    fArrayE(0),
    fArrayH(0),
    fArrayD(0),
    fArrayS(0), 
    fArrayN(0), 
    fArrayP(0), 
    fTreeMask(0), 
    fIsInit(kFALSE)
{
  // Constructor of an FMD input object.  Specify what data to read in
  // using the AddLoad member function.   Sub-classes should at a
  // minimum overload the member function Event.   A full job can be
  // executed using the member function Run. 
}

//____________________________________________________________________
Int_t
AliFMDInput::NEvents() const 
{
  // Get number of events
  if (fTreeE) return fTreeE->GetEntries();
  return -1;
}

//____________________________________________________________________
Bool_t
AliFMDInput::Init()
{
  // Initialize the object.  Get the needed loaders, and such. 

  // Check if we have been initialized
  if (fIsInit) { 
    AliWarning("Already initialized");
    return fIsInit;
  }
  if (fGAliceFile.IsNull()) fGAliceFile = "galice.root";
  // Get the loader
  fLoader = AliRunLoader::Open(fGAliceFile.Data(), "Alice", "read");
  if (!fLoader) {
    AliError(Form("Coulnd't read the file %s", fGAliceFile.Data()));
    return kFALSE;
  }
  
  // Get the run 
  if  (fLoader->LoadgAlice()) return kFALSE;
  fRun = fLoader->GetAliRun();
  
  // Get the FMD 
  fFMD = static_cast<AliFMD*>(fRun->GetDetector("FMD"));
  if (!fFMD) {
    AliError("Failed to get detector FMD from loader");
    return kFALSE;
  }
  
  // Get the FMD loader
  fFMDLoader = fLoader->GetLoader("FMDLoader");
  if (!fFMDLoader) {
    AliError("Failed to get detector FMD loader from loader");
    return kFALSE;
  }
  if (fLoader->LoadHeader()) { 
    AliError("Failed to get event header information from loader");
    return kFALSE;
  }
  fTreeE = fLoader->TreeE();

  // Optionally, get the geometry 
  if (TESTBIT(fTreeMask, kGeometry)) {
    TString fname(fRun->GetGeometryFileName());
    if (fname.IsNull()) {
      Warning("Init", "No file name for the geometry from AliRun");
      fname = gSystem->DirName(fGAliceFile);
      fname.Append("/geometry.root");
    }
    fGeoManager = TGeoManager::Import(fname.Data());
    if (!fGeoManager) {
      Fatal("Init", "No geometry manager found");
      return kFALSE;
    }
  }
  
  fIsInit = kTRUE;
  return fIsInit;
}

//____________________________________________________________________
Bool_t
AliFMDInput::Begin(Int_t event)
{
  // Called at the begining of each event.  Per default, it gets the
  // data trees and gets pointers to the output arrays.   Users can
  // overload this, but should call this member function in the
  // overloaded member function of the derived class. 

  // Check if we have been initialized
  if (!fIsInit) { 
    AliError("Not initialized");
    return fIsInit;
  }
  // Get the event 
  if (fLoader->GetEvent(event)) return kFALSE;
  AliInfo(Form("Now in event %d/%d", event, NEvents()));

  // Possibly load global kinematics information 
  if (TESTBIT(fTreeMask, kKinematics)) {
    AliInfo("Getting kinematics");
    if (fLoader->LoadKinematics()) return kFALSE;
    fStack = fLoader->Stack();
  }
  // Possibly load FMD Hit information 
  if (TESTBIT(fTreeMask, kHits)) {
    AliInfo("Getting FMD hits");
    if (fFMDLoader->LoadHits()) return kFALSE;
    fTreeH = fFMDLoader->TreeH();
    if (!fArrayH) fArrayH = fFMD->Hits(); 
  }
  // Possibly load FMD Digit information 
  if (TESTBIT(fTreeMask, kDigits)) {
    AliInfo("Getting FMD digits");
    if (fFMDLoader->LoadDigits()) return kFALSE;
    fTreeD = fFMDLoader->TreeD();
    if (!fArrayD) fArrayD = fFMD->Digits();
  }
  // Possibly load FMD Sdigit information 
  if (TESTBIT(fTreeMask, kSDigits)) {
    AliInfo("Getting FMD summable digits");
    if (fFMDLoader->LoadSDigits()) return kFALSE;
    fTreeS = fFMDLoader->TreeS();
    if (!fArrayS) fArrayS = fFMD->SDigits();
  }
  // Possibly load FMD RecPoints information 
  if (TESTBIT(fTreeMask, kRecPoints)) {
    AliInfo("Getting FMD reconstructed points");
    if (fFMDLoader->LoadRecPoints()) return kFALSE;
    fTreeR = fFMDLoader->TreeR();
    if (!fArrayN) fArrayN = new TClonesArray("AliFMDMultStrip");
    // if (!fArrayP) fArrayP = new TClonesArray("AliFMDMultRegion");
    fTreeR->SetBranchAddress("FMD",  &fArrayN);
    // fTreeR->SetBranchAddress("FMDPoisson", &fArrayP);
  }
  return kTRUE;
}

//____________________________________________________________________
Bool_t
AliFMDInput::End()
{
  // Called at the end of each event.  Per default, it unloads the
  // data trees and resets the pointers to the output arrays.   Users
  // can overload this, but should call this member function in the
  // overloaded member function of the derived class. 

  // Check if we have been initialized
  if (!fIsInit) { 
    AliError("Not initialized");
    return fIsInit;
  }
  // Possibly unload global kinematics information 
  if (TESTBIT(fTreeMask, kKinematics)) {
    fLoader->UnloadKinematics();
    // fTreeK = 0;
    fStack = 0;
  }
  // Possibly unload FMD Hit information 
  if (TESTBIT(fTreeMask, kHits)) {
    fFMDLoader->UnloadHits();
    fTreeH = 0;
  }
  // Possibly unload FMD Digit information 
  if (TESTBIT(fTreeMask, kDigits)) {
    fFMDLoader->UnloadDigits();
    fTreeD = 0;
  }
  // Possibly unload FMD Sdigit information 
  if (TESTBIT(fTreeMask, kSDigits)) {
    fFMDLoader->UnloadSDigits();
    fTreeS = 0;
  }
  // Possibly unload FMD RecPoints information 
  if (TESTBIT(fTreeMask, kRecPoints)) {
    fFMDLoader->UnloadRecPoints();
    fTreeR = 0;
  }
  AliInfo("Now out event");
  return kTRUE;
}

//____________________________________________________________________
Bool_t
AliFMDInput::Run()
{
  // Run over all events and files references in galice.root 

  Bool_t retval;
  if (!(retval = Init())) return retval;

  Int_t nEvents = NEvents();
  for (Int_t event = 0; event < nEvents; event++) {
    if (!(retval = Begin(event))) break;
    if (!(retval = Event())) break;
    if (!(retval = End())) break;
  }
  if (!retval) return retval;
  retval = Finish();
  return retval;
}

//====================================================================
ClassImp(AliFMDInputHits)
#if 0
  ;
#endif

Bool_t
AliFMDInputHits::Event()
{
  // Read the hit tree, and pass each hit to the member function
  // ProcessHit.  Optionally, if the user said `AddLoad(kKinematics)'
  // the track corresponding to the hit will also be passed to the
  // ProcessHit member function of the derived class.
  if (!fTreeH) {
    AliError("No hit tree defined");
    return kFALSE;
  }
  Int_t nTracks = fTreeH->GetEntries();
  for (Int_t i = 0; i < nTracks; i++) {
    Int_t hitRead  = fTreeH->GetEntry(i);
    if (hitRead <= 0) continue;
    if (!fArrayH) {
      AliError("No hit array defined");
      return kFALSE;
    }
    Int_t nHit = fArrayH->GetEntries();
    if (nHit <= 0) continue;
    for (Int_t j = 0; j < nHit; j++) {
      AliFMDHit* hit = static_cast<AliFMDHit*>(fArrayH->At(j));
      if (!hit) continue;
      TParticle* track = 0;
      if (TESTBIT(fTreeMask, kKinematics) && fStack) {
	Int_t trackno = hit->Track();
	track = fStack->Particle(trackno);
      }
      if (!ProcessHit(hit, track)) return kFALSE;
    }    
  }
  return kTRUE;
}

//====================================================================
ClassImp(AliFMDInputDigits)
#if 0
  ;
#endif

Bool_t
AliFMDInputDigits::Event()
{
  // Read the digit tree, and pass each digit to the member function
  // ProcessDigit.
  Int_t nEv = fTreeD->GetEntries();
  for (Int_t i = 0; i < nEv; i++) {
    Int_t digitRead  = fTreeD->GetEntry(i);
    if (digitRead <= 0) continue;
    Int_t nDigit = fArrayD->GetEntries();
    if (nDigit <= 0) continue;
    for (Int_t j = 0; j < nDigit; j++) {
      AliFMDDigit* digit = static_cast<AliFMDDigit*>(fArrayD->At(j));
      if (!digit) continue;
      if (!ProcessDigit(digit)) return kFALSE;
    }    
  }
  return kTRUE;
}

//====================================================================
ClassImp(AliFMDInputSDigits)
#if 0
  ;
#endif

Bool_t
AliFMDInputSDigits::Event()
{
  // Read the summable digit tree, and pass each sumable digit to the
  // member function ProcessSdigit.
  Int_t nEv = fTreeD->GetEntries();
  for (Int_t i = 0; i < nEv; i++) {
    Int_t sdigitRead  = fTreeS->GetEntry(i);
    if (sdigitRead <= 0) continue;
    Int_t nSdigit = fArrayS->GetEntries();
    if (nSdigit <= 0) continue;
    for (Int_t j = 0; j < nSdigit; j++) {
      AliFMDSDigit* sdigit = static_cast<AliFMDSDigit*>(fArrayS->At(j));
      if (!sdigit) continue;
      if (!ProcessSDigit(sdigit)) return kFALSE;
    }    
  }
  return kTRUE;
}

//====================================================================
ClassImp(AliFMDInputRecPoints)
#if 0
  ;
#endif

Bool_t
AliFMDInputRecPoints::Event()
{
  // Read the reconstrcted points tree, and pass each reconstruction
  // object to either ProcessStrip (for AliFMDMultStrip objects), or 
  // ProcessRegion (for AliFMDMultRegion objects).
  Int_t nEv = fTreeR->GetEntries();
  for (Int_t i = 0; i < nEv; i++) {
    Int_t recRead  = fTreeR->GetEntry(i);
    if (recRead <= 0) continue;
    Int_t nRecStrip = fArrayN->GetEntries();
    for (Int_t j = 0; j < nRecStrip; j++) {
      AliFMDMultStrip* strip = static_cast<AliFMDMultStrip*>(fArrayN->At(j));
      if (!strip) continue;
      if (!ProcessStrip(strip)) return kFALSE;
    }    
    Int_t nRecRegion = fArrayP->GetEntries();
    for (Int_t j = 0; j < nRecRegion; j++) {
      AliFMDMultRegion* region =static_cast<AliFMDMultRegion*>(fArrayP->At(j));
      if (!region) continue;
      if (!ProcessRegion(region)) return kFALSE;
    }    
  }
  return kTRUE;
}

    

//____________________________________________________________________
//
// EOF
//
