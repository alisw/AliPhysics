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
//_________________________________________________________________________
// Version of AliPHOSv0_hits which allows for keeping all hits in TreeH
//                  
//*-- Author: Gines MARTINEZ (SUBATECH)


// --- ROOT system ---

#include "TBRIK.h"
#include "TNode.h"
#include "TRandom.h"

// --- Standard library ---

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <strstream.h>

// --- AliRoot header files ---

#include "AliPHOSv0_hits.h"
#include "AliPHOSHit.h"
#include "AliPHOSDigit.h"
#include "AliPHOSReconstructioner.h"
#include "AliRun.h"
#include "AliConst.h"

ClassImp(AliPHOSv0_hits)

//____________________________________________________________________________
AliPHOSv0_hits::AliPHOSv0_hits()
{
  // ctor
  fNTmpHits = 0 ; 
  fTmpHits  = 0 ; 
}

//____________________________________________________________________________
AliPHOSv0_hits::AliPHOSv0_hits(const char *name, const char *title):
  AliPHOSv0(name,title)
{
   fHits= new TClonesArray("AliPHOSHit",1000) ;
}

//____________________________________________________________________________
AliPHOSv0_hits::~AliPHOSv0_hits()
{
  // dtor
  fTmpHits->Delete() ; 
  delete fTmpHits ;
  fTmpHits = 0 ; 

  fEmcClusters->Delete() ; 
  delete fEmcClusters ; 
  fEmcClusters = 0 ; 

  fPpsdClusters->Delete() ;
  delete fPpsdClusters ;
  fPpsdClusters = 0 ; 

  fTrackSegments->Delete() ; 
  delete fTrackSegments ;
  fTrackSegments = 0 ; 
}

//____________________________________________________________________________
void AliPHOSv0_hits::AddHit(Int_t primary, Int_t Id, Float_t * hits)
{
  // Add a hit to the hit list.
  // In this version of AliPHOSv0, a PHOS hit is real geant 
  // hits in a single crystal or in a single PPSD gas cell

  //  cout << "Primary particle is " << primary << endl;
  //cout << "Vol Id is " << Id << endl;
  //cout << "hits is " << hits[0] << "  " << hits[1] << "  " << hits[2] << "   " << hits[3] <<endl;

  cout << " Adding a hit number " << fNhits << endl ;

  TClonesArray &ltmphits = *fHits ;
  AliPHOSHit *newHit ;

  //  fHits->Print("");

  newHit = new AliPHOSHit(primary, Id, hits) ;

  // We DO want to save in TreeH the raw hits 
  //  TClonesArray &lhits = *fHits;
  cout << " Adding a hit to fHits TCloneArray number " << fNhits << endl ;     
  new(ltmphits[fNhits]) AliPHOSHit(*newHit) ;
  fNhits++ ;

  cout << " Added a hit to fHits TCloneArray number " << fNhits << endl ; 
  // We do not want to save in TreeH the raw hits 
  //   new(lhits[fNhits]) AliPHOSHit(*newHit) ;    
  //   fNhits++ ;

  // Please note that the fTmpHits array must survive up to the
  // end of the events, so it does not appear e.g. in ResetHits() (
  // which is called at the end of each primary).  

  delete newHit;

}




//___________________________________________________________________________
void AliPHOSv0_hits::FinishEvent()
{
  // Makes the digits from the sum of summed hit in a single crystal or PPSD gas cell
  // Adds to the energy the electronic noise
  // Keeps digits with energy above fDigitThreshold

  // Save the cumulated hits instead of raw hits (need to create the branch myself)
  // It is put in the Digit Tree because the TreeH is filled after each primary
  // and the TreeD at the end of the event.
  
  Int_t i ;
  Int_t relid[4];
  Int_t j ; 
  TClonesArray &lDigits = *fDigits ;
  AliPHOSHit  * hit ;
  AliPHOSDigit * newdigit ;
  AliPHOSDigit * curdigit ;
  Bool_t deja = kFALSE ; 
  
  for ( i = 0 ; i < fNhits ; i++ ) {
    hit = (AliPHOSHit*)fHits->At(i) ;
    newdigit = new AliPHOSDigit( hit->GetPrimary(), hit->GetId(), Digitize( hit->GetEnergy() ) ) ;
    deja =kFALSE ;
    for ( j = 0 ; j < fNdigits ;  j++) { 
      curdigit = (AliPHOSDigit*) lDigits[j] ;
      if ( *curdigit == *newdigit) {
	*curdigit = *curdigit + *newdigit ; 
	deja = kTRUE ; 
      }
    }
    if ( !deja ) {
      new(lDigits[fNdigits]) AliPHOSDigit(* newdigit) ;
      fNdigits++ ;  
    }
 
    delete newdigit ;    
  } 
  
  // Noise induced by the PIN diode of the PbWO crystals

  Float_t energyandnoise ;
  for ( i = 0 ; i < fNdigits ; i++ ) {
    newdigit =  (AliPHOSDigit * ) fDigits->At(i) ;
    fGeom->AbsToRelNumbering(newdigit->GetId(), relid) ;

    if (relid[1]==0){   // Digits belong to EMC (PbW0_4 crystals)
      energyandnoise = newdigit->GetAmp() + Digitize(gRandom->Gaus(0., fPinElectronicNoise)) ;

      if (energyandnoise < 0 ) 
	energyandnoise = 0 ;

      if ( newdigit->GetAmp() < fDigitThreshold ) // if threshold not surpassed, remove digit from list
	fDigits->RemoveAt(i) ; 
    }
  }
  
  fDigits->Compress() ;  

  fNdigits =  fDigits->GetEntries() ; 
  for (i = 0 ; i < fNdigits ; i++) { 
    newdigit = (AliPHOSDigit *) fDigits->At(i) ; 
    newdigit->SetIndexInList(i) ; 
  }

}

