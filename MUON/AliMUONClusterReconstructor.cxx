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

#include "AliMUON.h"
#include "AliMUONDigit.h"
#include "AliMUONConstants.h"
#include "AliMUONData.h"
#include "AliMUONClusterFinderVS.h"
#include "AliMUONClusterInput.h"
#include "AliMUONRawCluster.h"
#include "AliLog.h"

ClassImp(AliMUONClusterReconstructor) // Class implementation in ROOT context

//__________________________________________________________________________
  AliMUONClusterReconstructor::AliMUONClusterReconstructor(AliLoader* loader, AliMUONData* data)
  : TObject()
{
  // Standard Constructor

  // initialize loader's
  fLoader = loader;

  // initialize container
  if (data == 0x0)
    fMUONData = new AliMUONData(fLoader,"MUON","MUON");
  else
    fMUONData = data;
  
  // reconstruction model
  fRecModel = new AliMUONClusterFinderVS();
  //fRecModel = new AliMUONClusterFinderAZ();

}

//__________________________________________________________________________
AliMUONClusterReconstructor::AliMUONClusterReconstructor()
  : TObject(),
    fMUONData(0),
    fRecModel(0),
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
  if (fRecModel)
    delete fRecModel;

  return;
}
//____________________________________________________________________
void AliMUONClusterReconstructor::Digits2Clusters()
{

    TClonesArray *dig1, *dig2, *digAll;
    Int_t ndig, k, idDE, idDEprev;
    dig1 = new TClonesArray("AliMUONDigit",1000);
    dig2 = new TClonesArray("AliMUONDigit",1000);
    digAll = new TClonesArray("AliMUONDigit",2000);

    AliMUONDigit* digit;

    TArrayI id(200); // contains the different IdDE
   
  
// Loop on chambers and on cathode planes     
    TClonesArray* muonDigits;
    Int_t n2;
    Int_t n1;
  
    for (Int_t ich = 0; ich < AliMUONConstants::NTrackingCh(); ich++) {
 
      id.Reset();
      n1 = 0;
      n2 = 0;
      //cathode 0 & 1
      fMUONData->ResetDigits();
      fMUONData->GetDigits();
      muonDigits = fMUONData->Digits(ich); 
      ndig = muonDigits->GetEntriesFast();
      TClonesArray &lDigit = *digAll;

      idDEprev = 0;
      muonDigits->Sort();
      for (k = 0; k < ndig; k++) {

	digit = (AliMUONDigit*) muonDigits->UncheckedAt(k);
	new(lDigit[n1++]) AliMUONDigit(*digit);
	idDE = digit->DetElemId();
	if (idDE != idDEprev) {
	  id.AddAt(idDE,n2++);
	}
	idDEprev = idDE;
      }

      Int_t idSize = n2;

      // loop over id DE
      for (idDE = 0; idDE < idSize; idDE++) {
	TClonesArray &lhits1 = *dig1;
	TClonesArray &lhits2 = *dig2;
	dig1->Clear();
	dig2->Clear();
	n1 = n2 = 0;

	for (k = 0; k < digAll->GetEntriesFast(); k++) {
	  digit = (AliMUONDigit*) digAll->UncheckedAt(k);
	  //	  printf("digit idDE %d\n", digit->DetElemId());
	  if (id[idDE] == digit->DetElemId()) {
	    if (digit->Cathode() == 0)
	      new(lhits1[n1++]) AliMUONDigit(*digit);
	    else 
	      new(lhits2[n2++]) AliMUONDigit(*digit);
	  }
	}

	// cluster finder
	if (fRecModel) {
	  AliMUONClusterInput::Instance()->SetDigits(ich, id[idDE], dig1, dig2);
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
      } // idDE
      digAll->Delete();
    } // for ich
    delete dig1;
    delete dig2;
    delete digAll;
}

//_______________________________________________________________________
void AliMUONClusterReconstructor::Trigger2Trigger() 
{
// copy trigger from TreeD to TreeR
  fMUONData->GetTriggerD();
}
