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
// Implementation version v1 of EMCAL Manager class 
// An object of this class does not produce digits
// It is the one to use if you do want to produce outputs in TREEH 
//                  
//*-- Author: Sahal Yacoob (LBL /UCT)
//*--       : Jennifer Klay (LBL)


// --- ROOT system ---
#include "TPGON.h"
#include "TTUBS.h"
#include "TNode.h"
#include "TRandom.h"
#include "TTree.h"
#include "TGeometry.h"


// --- Standard library ---

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <strstream.h>
#include <iostream.h>
#include <math.h>
// --- AliRoot header files ---

#include "AliEMCALv1.h"
#include "AliEMCALHit.h"
#include "AliEMCALGeometry.h"
#include "AliConst.h"
#include "AliRun.h"
#include "AliMC.h"

ClassImp(AliEMCALv1)


//______________________________________________________________________
AliEMCALv1::AliEMCALv1():AliEMCALv0(){
  // ctor
}
//______________________________________________________________________
AliEMCALv1::AliEMCALv1(const char *name, const char *title):
    AliEMCALv0(name,title){
    // Standard Creator.

    fHits= new TClonesArray("AliEMCALHit",1000);
    gAlice->AddHitList(fHits);

    fNhits = 0;

    fIshunt     =  1; // All hits are associated with primary particles
}
//______________________________________________________________________
AliEMCALv1::~AliEMCALv1(){
    // dtor

    if ( fHits) {
	fHits->Delete();
	delete fHits;
	fHits = 0;
    }
}
//______________________________________________________________________
void AliEMCALv1::AddHit(Int_t shunt, Int_t primary, Int_t tracknumber, 
			Int_t id, Float_t * hits,TLorentzVector *p){
    // Add a hit to the hit list.
    // A PHOS hit is the sum of all hits in a single crystal
    //   or in a single PPSD gas cell
    Int_t hitCounter;
    
    AliEMCALHit *newHit;
    AliEMCALHit *curHit;
    Bool_t deja = kFALSE;

    newHit = new AliEMCALHit(shunt, primary, tracknumber, id, hits, p);

    for ( hitCounter = fNhits-1; hitCounter >= 0 && !deja; hitCounter-- ) {
	curHit = (AliEMCALHit*) (*fHits)[hitCounter];
	// We add hits with the same primary, while GEANT treats
	// primaries succesively
	if(curHit->GetPrimary() != primary) break;
	if( *curHit == *newHit ) {
	    *curHit = *curHit + *newHit;
	    deja = kTRUE;
	} // end if
    } // end for hitCounter

    if ( !deja ) {
	new((*fHits)[fNhits]) AliEMCALHit(*newHit);
	fNhits++;
    } // end if

    delete newHit;
}
//______________________________________________________________________
void AliEMCALv1::StepManager(void){
    // Accumulates hits as long as the track stays in a single
    // crystal or PPSD gas Cell
    Int_t          id[2];           // (layer, phi, Eta) indices
    Int_t          absid;
    // position wrt MRS and energy deposited
    Float_t        xyze[4]={0.,0.,0.,0.};
    TLorentzVector pos; // Lorentz vector of the track current position.
    TLorentzVector mom; // Lorentz vector of the track current momentum.
    Int_t tracknumber =  gAlice->CurrentTrack();
    static Int_t primary;

    gMC->CurrentVolOffID(1, id[0]); // get the POLY copy number;
    gMC->CurrentVolID(id[1]); // get the phi number inside the layer
    primary = gAlice->GetPrimary(tracknumber);
    if(gMC->IsTrackEntering()&&id[0]==1) primary = tracknumber;
//    tracknumber = primary;
    gMC->TrackPosition(pos);
    gMC->TrackMomentum(mom);
    xyze[0] = pos[0];
    xyze[1] = pos[1];
    xyze[2] = pos[2];
    xyze[3] = gMC->Edep();

    if(xyze[3] > 0.){// Track is inside the crystal and deposits some energy
	absid = (id[0]-1)*(fGeom->GetNPhi()) + id[1];
	AddHit(fIshunt, primary,tracknumber, absid, xyze, &mom);
    } // there is deposited energy

}
