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
// Version of AliPHOSv1 which keeps all hits in TreeH
// AddHit, StepManager,and FinishEvent are redefined 
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

#include "AliPHOSv2.h"
#include "AliPHOSHit.h"
#include "AliPHOSDigit.h"
#include "AliPHOSReconstructioner.h"
#include "AliRun.h"
#include "AliConst.h"

ClassImp(AliPHOSv2)

//____________________________________________________________________________
AliPHOSv2::AliPHOSv2()
{
  // default ctor
  fNTmpHits = 0 ; 
  fTmpHits  = 0 ; 
}

//____________________________________________________________________________
AliPHOSv2::AliPHOSv2(const char *name, const char *title):
AliPHOSv1(name,title)
{
  // ctor
  
   fHits= new TClonesArray("AliPHOSHit",1000) ;
}

//____________________________________________________________________________
AliPHOSv2::~AliPHOSv2()
{
  // dtor
  if ( fTmpHits) {
    fTmpHits->Delete() ; 
    delete fTmpHits ;
    fTmpHits = 0 ; 
  }
  
  if ( fEmcRecPoints ) { 
    fEmcRecPoints->Delete() ; 
    delete fEmcRecPoints ; 
    fEmcRecPoints = 0 ; 
  }
  
  if ( fPpsdRecPoints ) { 
    fPpsdRecPoints->Delete() ;
    delete fPpsdRecPoints ;
    fPpsdRecPoints = 0 ; 
  }
 
  if ( fTrackSegments ) {
    fTrackSegments->Delete() ; 
    delete fTrackSegments ;
    fTrackSegments = 0 ; 
  }
}

//____________________________________________________________________________
void AliPHOSv2::AddHit(Int_t shunt, Int_t primary, Int_t tracknumber, Int_t Id, Float_t * hits, Int_t pid)
{
  // Add a hit to the hit list.
  // A PHOS hit is the sum of all hits in a single crystal
  //   or in a single PPSD gas cell

  Int_t hitCounter ;
  TClonesArray &ltmphits = *fTmpHits ;
  AliPHOSHit *newHit ;
  AliPHOSHit *curHit ;
  Bool_t deja = kFALSE ;

  // In any case, fills the fTmpHit TClonesArray (with "accumulated hits")

  newHit = new AliPHOSHit(shunt, primary, tracknumber, Id, hits, pid) ;

  // We do want to save in TreeH the raw hits 
  TClonesArray &lhits = *fHits;

  for ( hitCounter = 0 ; hitCounter < fNTmpHits && !deja ; hitCounter++ ) {
    curHit = (AliPHOSHit*) ltmphits[hitCounter] ;
  if( *curHit == *newHit ) {
    *curHit = *curHit + *newHit ;
    deja = kTRUE ;
    }
  }
         
  if ( !deja ) {
    new(ltmphits[fNTmpHits]) AliPHOSHit(*newHit) ;
    fNTmpHits++ ;
  }

  // We do want to save in TreeH the raw hits 
  new(lhits[fNhits]) AliPHOSHit(*newHit) ;    
  fNhits++ ;

  // Please note that the fTmpHits array must survive up to the
  // end of the events, so it does not appear e.g. in ResetHits() (
  // which is called at the end of each primary).  

  delete newHit;

}


