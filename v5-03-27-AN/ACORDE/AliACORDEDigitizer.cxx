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

///_________________________________________________________________________
///
/// This class constructs Digits out of Hits
///
///

// --- Standard library ---

// --- ROOT system ---
#include <TMath.h>
#include <TTree.h>
#include <TRandom.h>

// --- AliRoot header files ---
#include "AliRun.h"
#include "AliACORDE.h"
#include "AliACORDEhit.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliDigitizationInput.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliACORDECalibData.h"
#include "AliACORDEConstants.h"

#include "AliACORDEdigit.h"
#include "AliACORDEDigitizer.h"

ClassImp(AliACORDEDigitizer)

AliACORDEDigitizer::AliACORDEDigitizer()
  :AliDigitizer(),
   fCalibData(GetCalibData()),
   fNdigits(0),
   fDigits(0)
  
{
  // default constructor
}

AliACORDEDigitizer::AliACORDEDigitizer(AliDigitizationInput* digInput)
  :AliDigitizer(digInput),
   fCalibData(GetCalibData()),
   fNdigits(0),
   fDigits(0)

{
  // constructor

}

AliACORDEDigitizer::~AliACORDEDigitizer()
{
  // destructor
  
  if (fDigits) {
    fDigits->Delete();
    delete fDigits;
    fDigits=0;
  }
}


Bool_t AliACORDEDigitizer::Init()
{
  // Initialises the digitizer
  
  // Initialises the Digit array
  fDigits = new TClonesArray ("AliACORDEdigit", 1000);
  
  return kTRUE;
}

void AliACORDEDigitizer::Digitize(Option_t* /*option*/)
{

  // Creates digits from hits

  // 1.- create and initialize temporal variables
  // 2.- get loaders to access hots and digits
  // 3.- loop over all input. 
  //     3.1 for each input loop over all track
  //     3.2 for each track loop over all hits
  //     3.3 for each hit, check 
  //         if energy above threshold
  //         if survives efficiency
  //         then store the plastic, the time and track in temp arrays
  //         if a hit already survived for this plastic take the hit
  //         with the lowest time
  // 4.- loop over temporal array
  //     if both plastic have a surviving hit and the time
  //     difference is below the limit, add a new digit
  // 


  // 1.- temporal variables
  Float_t emin = AliACORDEConstants::Instance()->HitEnergyThreshold();
  Float_t td = AliACORDEConstants::Instance()->MaxHitTimeDifference();
  Int_t modules[60]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Int_t moduls[60]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; 
  Int_t mods;
  Float_t PlasticTimes[2][60]; 
  Int_t PlasticTracks[2][60];
  for (Int_t i=0;i<60;i++) {
    PlasticTimes[0][i]=-1.0;
    PlasticTimes[1][i]=-1.0;
    PlasticTracks[0][i]=-1;
    PlasticTracks[1][i]=-1;
  }

  // 2.- get loaders
  AliRunLoader* outRunLoader =
    AliRunLoader::GetRunLoader(fDigInput->GetOutputFolderName());
  if (!outRunLoader) {
    Error("Exec", "Can not get output Run Loader");
    return;}
  
  AliLoader* outLoader = outRunLoader->GetLoader("ACORDELoader");
  if (!outLoader) {
    Error("Exec", "Can not get output ACORDE Loader");
    return;}
  
  outLoader->LoadDigits("update");
  if (!outLoader->TreeD()) outLoader->MakeTree("D");
  outLoader->MakeDigitsContainer();
  TTree* treeD  = outLoader->TreeD();
  Int_t bufsize = 16000;
  treeD->Branch("ACORDEdigit", &fDigits, bufsize);
  
  // 3. loop over inputs
  for (Int_t iInput = 0; iInput < fDigInput->GetNinputs(); iInput++) {
    AliRunLoader* runLoader =
      AliRunLoader::GetRunLoader(fDigInput->GetInputFolderName(iInput));
    AliLoader* loader = runLoader->GetLoader("ACORDELoader");
    if (!loader) {
      Error("Exec", "Can not get ACORDE Loader for input %d", iInput);
      continue;}
    
    if (!runLoader->GetAliRun()) runLoader->LoadgAlice();
    
    AliACORDE* acorde = (AliACORDE*) runLoader->GetAliRun()->GetDetector("ACORDE");
    if (!acorde) {
      Error("Exec", "No ACORDE detector for input %d", iInput);
      continue;}
    
    loader->LoadHits();
    TTree* treeH = loader->TreeH();
    if (!treeH) {
      Error("Exec", "Cannot get TreeH for input %d", iInput);
      continue; }
    
    TClonesArray* hits = acorde->Hits();
    
    // 3.1 loop over all tracks
    Int_t nTracks = (Int_t) treeH->GetEntries();
    for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
      acorde->ResetHits();
      treeH->GetEvent(iTrack);
      Int_t nHits = hits->GetEntriesFast();
      // 3.2 loop over hits
      for (Int_t iHit = 0; iHit < nHits; iHit++) {
        AliACORDEhit* hit = (AliACORDEhit *)hits->UncheckedAt(iHit);
	// 3.3 select hit
	// get hit info
	Float_t eloss_mev = hit->Eloss()*1000.0;
	Int_t module = hit->GetModule();
        modules[module]=1;
	Int_t plastic = hit->GetPlastic();
	Float_t time_ns = hit->GetTime()*1e9;
	Float_t eff = TMath::Sqrt(fCalibData->GetEfficiency(module));
	// if enough energy and efficiency
	if( eloss_mev > emin && gRandom->Uniform() < eff ) {
	  // if first hit or earlier track
	  if ((PlasticTimes[plastic-1][module-1] == -1.0) ||
	      (PlasticTimes[plastic-1][module-1] > time_ns) ) {
	    PlasticTimes[plastic-1][module-1]= time_ns;
	    PlasticTracks[plastic-1][module-1]= hit->GetTrack();
	  } 
	}
      } // end of hit   loop
    } // end of track loop
    for(Int_t i=0;i<60;i++){moduls[i]=modules[i];}
    
    loader->UnloadHits();
    
  }  // end of input loop

  // 4.- loop over temporal arrays to add hits
  Int_t tracks[3]={-1,-1,-1};
  for (Int_t i=0; i<60; i++) {
    // if both modules have a hit
    // if time diff small enough
    Float_t diff = TMath::Abs(PlasticTimes[0][i]-PlasticTimes[1][i]);
    if (diff < td) {
      tracks[0] = PlasticTracks[0][i];
      if (PlasticTracks[0][i] != PlasticTracks[1][i]) 
	tracks[1] = PlasticTracks[1][i];
      if(moduls[i]==1) {
	mods = i;
	//	Float_t module_time = (PlasticTimes[0][i] > PlasticTimes[1][i] ? 
	//			       PlasticTimes[0][i] : PlasticTimes[1][i]);
	//	AddDigit(tracks, mods, module_time);
	AddDigit(tracks, mods, 0);
      }
    }
  }
  treeD->Fill();
  outLoader->WriteDigits("OVERWRITE");
  outLoader->UnloadDigits();
  ResetDigit();
}

//____________________________________________________________________________

void AliACORDEDigitizer::AddDigit(Int_t* track, Int_t module, Float_t time)
{
  
  // Adds Digit
  
  TClonesArray &ldigits = *fDigits;
  new(ldigits[fNdigits++]) AliACORDEdigit(track,module,time);
}

void AliACORDEDigitizer::AddDigit(Int_t* modul,Float_t time)
{
	// MRC Adds Digit
  TClonesArray &ldigits = *fDigits;
  new(ldigits[fNdigits++]) AliACORDEdigit(modul,time);
  

}
void AliACORDEDigitizer::ResetDigit()
{
//
// Clears Digits
//
  fNdigits = 0;
  if (fDigits) fDigits->Delete();
}


AliACORDECalibData* AliACORDEDigitizer::GetCalibData() const

{
  AliCDBManager *man = AliCDBManager::Instance();

  AliCDBEntry *entry=0;

  entry = man->Get("ACORDE/Calib/Data");

  if(!entry){
    AliWarning("Load of calibration data from default storage failed!");
    AliWarning("Calibration data will be loaded from local storage ($ALICE_ROOT)");
    Int_t runNumber = man->GetRun();
    entry = man->GetStorage("local://$ALICE_ROOT/OCDB")
      ->Get("ACORDE/Calib/Data",runNumber);

  }

  // Retrieval of data in directory ACORDE/Calib/Data:


  AliACORDECalibData *calibdata = 0;

  if (entry) calibdata = (AliACORDECalibData*) entry->GetObject();
  if (!calibdata)  AliError("No calibration data from calibration database !");


  return calibdata;

}

