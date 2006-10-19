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

// -----------------------------------
// Class AliMUONClusterReconstructor
// ----------------------------------
// MUON cluster reconstructor for MUON
// Should implement a virtual class ClusterFinder to choose between VS and AZ method

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

#include "AliMpDEManager.h"

#include "AliLog.h"

ClassImp(AliMUONClusterReconstructor) // Class implementation in ROOT context
 
//__________________________________________________________________________
AliMUONClusterReconstructor::AliMUONClusterReconstructor(AliMUONData* data)
: TObject(),
  fMUONData(data),
  fRecModel(new AliMUONClusterFinderVS()),
  fDigitsCath0(new TClonesArray("AliMUONDigit",1000)),
  fDigitsCath1(new TClonesArray("AliMUONDigit",1000))
{
/// Standard Constructor

  fDigitsCath0->SetOwner(kTRUE); 
  fDigitsCath1->SetOwner(kTRUE);
}

//__________________________________________________________________________
AliMUONClusterReconstructor::~AliMUONClusterReconstructor(void)
{
/// Destructor

  delete fRecModel;
  delete fDigitsCath0;
  delete fDigitsCath1;
}

//______________________________________________________________________________
void
AliMUONClusterReconstructor::ClusterizeOneDE(Int_t detElemId)
{
/// Clusterize one detection element, and let fMUONData know about
/// the results.
  
  if ( fDigitsCath0->GetEntriesFast() || fDigitsCath1->GetEntriesFast() )
  {
    Int_t iChamber = AliMpDEManager::GetChamberId(detElemId);
    AliMUONClusterInput::Instance()->SetDigits(iChamber, detElemId,
                                               fDigitsCath0,fDigitsCath1);
    AliDebug(3,Form("ClusterizeOneDE iChamber=%d DE=%d",iChamber,detElemId));
    StdoutToAliDebug(3,cout << "DigitsCath0=" << endl;
                     fDigitsCath0->Print();
                     cout << "DigitsCath1=" << endl;
                     fDigitsCath1->Print(););
    fRecModel->FindRawClusters();
    
    // copy results into the output container
    TClonesArray* tmp = fRecModel->GetRawClusters();
    for (Int_t id = 0; id < tmp->GetEntriesFast(); ++id) 
    {
      AliMUONRawCluster* pClus = (AliMUONRawCluster*) tmp->At(id);
      fMUONData->AddRawCluster(iChamber, *pClus);
    }        
    
    // Reset the arrays
    fDigitsCath0->Clear("C");
    fDigitsCath1->Clear("C");
  }
}

//____________________________________________________________________
void AliMUONClusterReconstructor::Digits2Clusters(Int_t chBeg)
{
/// Clusterize all the tracking chamber digits.
///
/// For each chamber, we loop *once* on that chamber digits, and store them
/// in 2 temporary arrays (one pair of arrays per detection element, 
/// one array per cathode). Once a pair of arrays is full (i.e. all the digits
/// of that detection element have been stored), we clusterize this DE, and
/// move to the next one.
  
  if (!fRecModel)
  {
    AliWarning("No reco model defined. Nothing to do...");
    return;
  }
  
  Int_t iChamber(-1);
  Int_t currentDE(-1);
  
  // Loop on chambers 
  for ( iChamber = chBeg; iChamber < AliMUONConstants::NTrackingCh(); ++iChamber ) 
  {
    TClonesArray* muonDigits = fMUONData->Digits(iChamber); 
    
    Int_t ndig = muonDigits->GetEntriesFast();
    if (!ndig) continue;
    
    muonDigits->Sort(); // the sort *must* be per DE (at least), otherwise
                        // the following logic with currentDE will fail.
    
    currentDE = -1; // initialize the DE counter (that is used to track 
                    // when we change of DE in the following loop over
                    // all digits) to an invalid value.

    for ( Int_t k = 0; k < ndig; ++k ) 
    {
      AliMUONDigit* digit = (AliMUONDigit*) muonDigits->UncheckedAt(k);
      if ( ! digit->Signal() > 0 ) continue; // skip void digits.
      
      if ( digit->DetElemId() != currentDE )
      {
        AliDebug(3,Form("Switching DE from %d to %d",currentDE,digit->DetElemId()));
        // we get to a new DE, so clusterize the previous one before
        // moving on.
        ClusterizeOneDE(currentDE);
        currentDE = digit->DetElemId();
      }
      
      // Add the digit to the array with the right cathode number.
      if (digit->Cathode() == 0)
      {
        new((*fDigitsCath0)[fDigitsCath0->GetLast()+1]) AliMUONDigit(*digit);
      }
      else 
      {
        new((*fDigitsCath1)[fDigitsCath1->GetLast()+1]) AliMUONDigit(*digit);
      }
    } // end of loop on chamber digits
    
    // As the above logic is based on detecting a change in DE number,
    // the last DE of each chamber has not been clusterized, so we do 
    // it here.
    ClusterizeOneDE(currentDE);
  } // end of loop over chambers
}

//_______________________________________________________________________
void 
AliMUONClusterReconstructor::SetRecoModel(AliMUONClusterFinderVS* rec)
{ 
/// Set reconstruction model

  delete fRecModel; 
  fRecModel = rec;
} 

//_______________________________________________________________________
void AliMUONClusterReconstructor::Trigger2Trigger() 
{
/// Copy trigger from TreeD to TreeR

  fMUONData->SetTreeAddress("GLT");
  fMUONData->GetTriggerD();
}
