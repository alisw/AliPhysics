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

#include <Riostream.h>
#include <TTree.h>
#include <TClonesArray.h>

#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliACORDETrigger.h"
#include "AliACORDEConstants.h"

//______________________________________________________________________
ClassImp(AliACORDETrigger)

AliACORDETrigger::AliACORDETrigger()
  :AliTriggerDetector(),
   fSingleMuon(0),
   fMultiMuon(0)
{

  cout << " ================>>>>>>>>>>>>>>>>> AliACORDETrigger" << endl;
   SetName("ACORDE");
   CreateInputs();
   for (Int_t i=0; i<60; i++) fModuleFired[i]=kFALSE;
}

void AliACORDETrigger::CreateInputs()
{
  // Do not create inputs again!!
  if( fInputs.GetEntriesFast() > 0 ) return;

  // two acorde triggers, single muon and multicoincidence
  fInputs.AddLast( new 
		   AliTriggerInput( "0ASL", 
				    "ACORDE", 0 ) );
  fInputs.AddLast( new 
		   AliTriggerInput( "0AMU",
				    "ACORDE", 0 ) );
}

void AliACORDETrigger::Trigger()
{
  
  // 1.- Get loaders and pointers
  // 2.- Loop over all entries
  //     set temporal variables to default values
  //     start the loop
  // 3.- Loop over all digits in an entrie
  //     Fill temporal arrays
  //     Find module with lowest time
  // 4.- Loop over temporal arrays
  //     Find number of modules within the time window
  // 5.- Set the relevant trigger

  // 1.- Get loaders and pointers
  AliRunLoader* runLoader = AliRunLoader::Instance();
  AliACORDELoader* loader = 
    (AliACORDELoader* )runLoader->GetLoader( "ACORDELoader" );
  loader->LoadDigits("READ");
  TTree* acordeDigitsTree = loader->TreeD();
  if (!acordeDigitsTree) return;
  TClonesArray* acordeDigits = new TClonesArray("AliACORDEdigit",1000);
  TBranch* digitBranch = acordeDigitsTree->GetBranch("ACORDEdigit");
  digitBranch->SetAddress(&acordeDigits);

  // 2.- Loop over all entries
  //     set temporal variables to default values

  Int_t MultiMin = AliACORDEConstants::Instance()->MultiMuonThreshold();
  Float_t time_window = AliACORDEConstants::Instance()->MultiMuonWindow();
  Int_t MinTimeModule = -1;
  Float_t MinTime = 1e10;
  Float_t ModuleTimes[60];
  for (Int_t i=0; i<60; i++) ModuleTimes[i] = -1.0;

  //     start the loop
  Int_t nEntries = (Int_t)acordeDigitsTree->GetEntries();
  cout << " ===AliACORDETrigger=== nEntries  " <<nEntries << endl; 
  for (Int_t e=0; e<nEntries; e++) {
    acordeDigitsTree->GetEvent(e);
    // 3.- Loop over all digits in an entrie
    //     Fill temporal arrays
    //     Find module with lowest time
    Int_t nDigits = acordeDigits->GetEntriesFast();
    cout << " ===AliACORDETrigger=== nDigits  " <<nDigits << endl; 
    for (Int_t d=0; d<nDigits; d++) {
      AliACORDEdigit* digit = (AliACORDEdigit*)acordeDigits->At(d);
      Int_t module = digit->GetModule();
      Float_t mod_time = digit->GetTime();
      ModuleTimes[module-1]=mod_time;
      if (mod_time < MinTime) {
	MinTime = mod_time;
	MinTimeModule = module;
      }
    } // end of loop over digits
  } // end of loop over events in digits tree

  // 4.- Loop over temporal arrays
  //     Find number of modules within the time window
  if (MinTimeModule == -1) return;
  for (Int_t i=0; i<60; i++) {
    if (ModuleTimes[i]<0) continue;
    Float_t diff = ModuleTimes[i]-MinTime;
    if (diff<time_window) {
      fMultiMuon++;
      fModuleFired[i]=kTRUE;
    }
  }
  cout << " fSingleMuon " << fSingleMuon
       << " MinTime " << MinTime
       << " fMultiMuon " << fMultiMuon << endl;
  // 5.- Set the relevant trigger
  fSingleMuon = MinTimeModule;
  SetInput( "0ASL" );
  if (fMultiMuon>=MultiMin) SetInput( "0AMU" );
  return;
}
