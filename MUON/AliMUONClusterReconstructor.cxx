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
// MUON cluster reconstructor for MUON
//
// Should implement a virtual class ClusterFinder to chose between VS and AZ method
////////////////////////////////////

#include "AliMUONClusterReconstructor.h"
#include "AliRun.h" // for gAlice
#include "AliRunLoader.h"
#include "AliLoader.h"

#include "AliMUONDigit.h"
#include "AliMUONConstants.h"
#include "AliMUONData.h"
#include "AliMUONClusterFinderVS.h"
#include "AliMUONClusterInput.h"
#include "AliMUONRawCluster.h"
#include "AliRawReader.h" // for raw data
#include "AliLog.h"


const Int_t AliMUONClusterReconstructor::fgkDefaultPrintLevel = 0;

ClassImp(AliMUONClusterReconstructor) // Class implementation in ROOT context

//__________________________________________________________________________
AliMUONClusterReconstructor::AliMUONClusterReconstructor(AliLoader* loader)
  : TObject(),
    fMUONData(0),
    fPrintLevel(fgkDefaultPrintLevel),
    fDebug(0)
{
  // Standard Constructor

  // initialize loader's
  fLoader = loader;

  // initialize container
  fMUONData  = new AliMUONData(fLoader,"MUON","MUON");

  // reconstruction model
  fRecModel = new AliMUONClusterFinderVS();
  //fRecModel = new AliMUONClusterFinderAZ();

}

//__________________________________________________________________________
AliMUONClusterReconstructor::AliMUONClusterReconstructor()
  : TObject(),
    fMUONData(0),
    fPrintLevel(fgkDefaultPrintLevel),
    fDebug(0),
    fLoader(0)
{
  // Default Constructor
}

//_______________________________________________________________________
AliMUONClusterReconstructor::AliMUONClusterReconstructor (const AliMUONClusterReconstructor& rhs)
  : TObject(rhs)
{
// Protected copy constructor

  AliFatal("Not implemented.");
}

//_______________________________________________________________________
AliMUONClusterReconstructor & 
AliMUONClusterReconstructor::operator=(const AliMUONClusterReconstructor& rhs)
{
// Protected assignement operator

  if (this == &rhs) return *this;

  AliFatal("Not implemented.");
    
  return *this;  
}

//__________________________________________________________________________
AliMUONClusterReconstructor::~AliMUONClusterReconstructor(void)
{

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
    TClonesArray * muonDigits;

    for (Int_t ich = 0; ich < 10; ich++) {

	fMUONData->ResetDigits();
	fMUONData->GetCathode(0);
	//TClonesArray *
	muonDigits = fMUONData->Digits(ich); 
	ndig=muonDigits->GetEntriesFast();
	AliDebug(1,Form("1 Found %d digits in %p chamber %d", ndig, (void*)muonDigits,ich));
	TClonesArray &lhits1 = *dig1;
	Int_t n = 0;
	for (k = 0; k < ndig; k++) {
	    digit = (AliMUONDigit*) muonDigits->UncheckedAt(k);
	    if (fRecModel->TestTrack(digit->Track(0)))
	      new(lhits1[n++]) AliMUONDigit(*digit);
	}
	fMUONData->ResetDigits();
	fMUONData->GetCathode(1);
	muonDigits =  fMUONData->Digits(ich);  
	ndig=muonDigits->GetEntriesFast();
	AliDebug(1,Form("2 Found %d digits in %p %d", ndig, (void*)muonDigits, ich));
	TClonesArray &lhits2 = *dig2;
	n=0;
	
	for (k=0; k<ndig; k++) {
	    digit= (AliMUONDigit*) muonDigits->UncheckedAt(k);
	    if (fRecModel->TestTrack(digit->Track(0)))
	      new(lhits2[n++]) AliMUONDigit(*digit);
	}

	if (fRecModel) {	 
	    AliMUONClusterInput::Instance()->SetDigits(ich, dig1, dig2);
	    fRecModel->FindRawClusters();
	}
	// copy into the container
	TClonesArray* tmp = fRecModel->GetRawClusters();
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

//____________________________________________________________________
void AliMUONClusterReconstructor::Digits2Clusters(AliRawReader* /*rawReader*/)
{

//  Perform cluster finding form raw data

   AliFatal("clusterization not implemented for raw data input");
}
//_______________________________________________________________________
void AliMUONClusterReconstructor::Trigger2Trigger() 
{
// copy trigger from TreeD to TreeR

  fMUONData->SetTreeAddress("GLT");
  fMUONData->GetTriggerD();
}
//_______________________________________________________________________
void AliMUONClusterReconstructor::Trigger2Trigger(AliRawReader* /*rawReader*/) 
{
// call the Trigger Algorithm from raw data and fill TreeR 

   AliFatal("Trigger not implemented for raw data input");

}
