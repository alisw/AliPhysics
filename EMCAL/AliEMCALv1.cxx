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
    fSamplingFraction = 12.9 ; 
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
void AliEMCALv1::AddHit(Int_t shunt, Int_t primary, Int_t tracknumber, Int_t iparent, Float_t ienergy, 
			Int_t id, Float_t * hits,Float_t * p){
    // Add a hit to the hit list.
    // An EMCAL hit is the sum of all hits in a single segment 
    //   originating from the same enterring particle 
    Int_t hitCounter;
    
    AliEMCALHit *newHit;
    AliEMCALHit *curHit;
    Bool_t deja = kFALSE;

    newHit = new AliEMCALHit(shunt, primary, tracknumber, iparent, ienergy, id, hits, p);
    for ( hitCounter = fNhits-1; hitCounter >= 0 && !deja; hitCounter-- ) {
	curHit = (AliEMCALHit*) (*fHits)[hitCounter];
	// We add hits with the same tracknumber, while GEANT treats
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
    Float_t        pmom[4]={0.,0.,0.,0.};
    TLorentzVector pos; // Lorentz vector of the track current position.
    TLorentzVector mom; // Lorentz vector of the track current momentum.
    Int_t tracknumber =  gAlice->CurrentTrack();
    Int_t primary = 0;
    static Int_t iparent = 0;
    static Float_t ienergy = 0;
    Int_t copy = 0;


 if(gMC->IsTrackEntering() && (strcmp(gMC->CurrentVolName(),"XALU") == 0)) // This Particle in enterring the Calorimeter 
 {
    iparent = tracknumber;
    gMC->TrackMomentum(mom);
    ienergy = mom[3]; 
}
 if(gMC->CurrentVolID(copy) == gMC->VolId("XPHI") ) { // We are in a Scintillator Layer 

    gMC->CurrentVolOffID(1, id[0]); // get the POLY copy number;
    gMC->CurrentVolID(id[1]); // get the phi number inside the layer
    primary = gAlice->GetPrimary(tracknumber);
    gMC->TrackPosition(pos);
    gMC->TrackMomentum(mom);
    xyze[0] = pos[0];
    xyze[1] = pos[1];
    xyze[2] = pos[2];
    xyze[3] = fSamplingFraction*(gMC->Edep()); // Correct for sampling calorimeter
    pmom[0] = mom[0];
    pmom[1] = mom[1];
    pmom[2] = mom[2];
    pmom[3] = mom[3];
    
    if(xyze[3] > 0.){// Track is inside the crystal and deposits some energy
	absid = (id[0]-1)*(fGeom->GetNPhi()) + id[1];
        if((absid/fGeom->GetNPhi()) < (2*fGeom->GetNZ()))
        {xyze[3] = 5*xyze[3]/6 ;}  //                                             Preshower readout must be scaled
        AddHit(fIshunt, primary,tracknumber, iparent, ienergy, absid, xyze, pmom);
    } // there is deposited energy
}
}
