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

////////////////////////////////////
//
// MUON event reconstructor in ALICE
//
// This class contains as data:
// * the parameters for the event reconstruction
// * a pointer to the array of hits to be reconstructed (the event)
// * a pointer to the array of segments made with these hits inside each station
// * a pointer to the array of reconstructed tracks
//
// It contains as methods, among others:
// * MakeEventToBeReconstructed to build the array of hits to be reconstructed
// * MakeSegments to build the segments
// * MakeTracks to build the tracks
//
////////////////////////////////////

#include "AliMUONClusterReconstructor.h"
#include "AliMUON.h"
#include "AliMUONDigit.h"
#include "AliMUONConstants.h"
#include "AliMUONData.h"
#include "AliMUONClusterFinderVS.h"
#include "AliMUONClusterInput.h"
#include "AliMUONRawCluster.h"
#include "AliRun.h" // for gAlice
#include "AliRunLoader.h"
#include "AliLoader.h"

const Int_t AliMUONClusterReconstructor::fgkDefaultPrintLevel = 0;

ClassImp(AliMUONClusterReconstructor) // Class implementation in ROOT context

//__________________________________________________________________________
AliMUONClusterReconstructor::AliMUONClusterReconstructor(AliLoader* loader)
  : TObject()
{
  // Standard Constructor
 
  fDebug           = 0;
  fNCh             = 0;
  fNTrackingCh     = 0;
  fChambers        = 0;
  fMUONData        = 0;
  fChambers = new TObjArray(AliMUONConstants::NCh());

  fPrintLevel = fgkDefaultPrintLevel;

  // initialize loader's
  fLoader = loader;

  // initialize container
  fMUONData  = new AliMUONData(fLoader,"MUON","MUON");

  // Loading AliRun master
  AliRunLoader* runloader = fLoader->GetRunLoader();
  if (runloader->GetAliRun() == 0x0) runloader->LoadgAlice();
  gAlice = runloader->GetAliRun();

  // getting MUON
  fMUON = (AliMUON*) gAlice->GetDetector("MUON");
}

//__________________________________________________________________________
AliMUONClusterReconstructor::AliMUONClusterReconstructor()
  : TObject(),
    fNCh(0),
    fNTrackingCh(0),
    fMUONData(0),
    fMUON(0),
    fChambers(0),
    fPrintLevel(fgkDefaultPrintLevel),
    fDebug(0),
    fLoader(0)
{
  // Default Constructor
}

//____________________________________________________________________
void AliMUONClusterReconstructor::SetReconstructionModel(Int_t id, AliMUONClusterFinderVS *reconst)
{
  // take infos chambers from AliMUON
  AliMUONChamber* pCh = 0;
  pCh = &(fMUON->Chamber(id));

  fChambers->AddAt(pCh, id);

  // Set ClusterFinder for chamber id
  ((AliMUONChamber*) fChambers->At(id))->SetReconstructionModel(reconst);
}
//_______________________________________________________________________
AliMUONClusterReconstructor::AliMUONClusterReconstructor (const AliMUONClusterReconstructor& rhs)
  : TObject(rhs)
{
// Protected copy constructor

  Fatal("AliMUONClusterReconstructor", "Not implemented.");
}

//_______________________________________________________________________
AliMUONClusterReconstructor & 
AliMUONClusterReconstructor::operator=(const AliMUONClusterReconstructor& rhs)
{
// Protected assignement operator

  if (this == &rhs) return *this;

  Fatal("operator=", "Not implemented.");
    
  return *this;  
}

//__________________________________________________________________________
AliMUONClusterReconstructor::~AliMUONClusterReconstructor(void)
{
  if (fChambers){
    fChambers->Clear(); // Sets pointers to 0 sinche it is not the owner
    delete fChambers;
  } 
  if (fMUONData)
    delete fMUONData;

  return;
}
//____________________________________________________________________
void AliMUONClusterReconstructor::Digits2Clusters()
{
//
//  Perform cluster finding
//
    TClonesArray *dig1, *dig2;
    Int_t ndig, k;
    dig1 = new TClonesArray("AliMUONDigit",1000);
    dig2 = new TClonesArray("AliMUONDigit",1000);
    AliMUONDigit *digit;
// Loop on chambers and on cathode planes
//
//    fMUONData->ResetRawClusters();        
    TClonesArray * muonDigits;

    for (Int_t ich = 0; ich < 10; ich++) {
 	AliMUONChamber* iChamber = (AliMUONChamber*) fChambers->At(ich);
	AliMUONClusterFinderVS* rec = iChamber->ReconstructionModel();
	//AliMUONClusterFinderAZ* rec = (AliMUONClusterFinderAZ*)iChamber->ReconstructionModel();

	fMUONData->ResetDigits();
	fMUONData->GetCathode(0);
	//TClonesArray *
	muonDigits = fMUONData->Digits(ich); 
	ndig=muonDigits->GetEntriesFast();
	if(fDebug)
	  printf("1 Found %d digits in %p chamber %d\n", ndig, (void*)muonDigits,ich);
	TClonesArray &lhits1 = *dig1;
	Int_t n = 0;
	for (k = 0; k < ndig; k++) {
	    digit = (AliMUONDigit*) muonDigits->UncheckedAt(k);
	    if (rec->TestTrack(digit->Track(0)))
	      new(lhits1[n++]) AliMUONDigit(*digit);
	}
	fMUONData->ResetDigits();
	fMUONData->GetCathode(1);
	muonDigits =  fMUONData->Digits(ich);  
	ndig=muonDigits->GetEntriesFast();
	if(fDebug)
	  printf("\n 2 Found %d digits in %p %d", ndig, (void*)muonDigits, ich);
	TClonesArray &lhits2 = *dig2;
	n=0;
	
	for (k=0; k<ndig; k++) {
	    digit= (AliMUONDigit*) muonDigits->UncheckedAt(k);
	    if (rec->TestTrack(digit->Track(0)))
	      new(lhits2[n++]) AliMUONDigit(*digit);
	}

	if (rec) {	 
	    AliMUONClusterInput::Instance()->SetDigits(ich, dig1, dig2);
	    rec->FindRawClusters();
	}
	// copy into the container
	TClonesArray* tmp = rec->GetRawClusters();
	for (Int_t id = 0; id < tmp->GetEntriesFast(); id++) {
	  AliMUONRawCluster* pClus = (AliMUONRawCluster*) tmp->At(id);
	  fMUONData->AddRawCluster(ich, *pClus);
	}
	dig1->Delete();
	dig2->Delete();
    } // for ich
    delete dig1;
    delete dig2;
}
