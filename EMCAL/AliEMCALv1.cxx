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
//*-- Implementation version v1 of EMCAL Manager class 
//*-- An object of this class does not produce digits
//*-- It is the one to use if you do want to produce outputs in TREEH 
//*--                  
//*-- Author: Sahal Yacoob (LBL /UCT)
//*--       : Jennifer Klay (LBL)
// This Class not stores information on all particles prior to EMCAL entry - in order to facilitate analysis.
// This is done by setting fIShunt =2, and flagging all parents of particles entering the EMCAL.

// 15/02/2002 .... Yves Schutz
//  1. fSamplingFraction and fLayerToPreshowerRatio have been removed
//  2. Timing signal is collected and added to hit

// --- ROOT system ---
#include <TParticle.h>
#include <TVirtualMC.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliEMCALv1.h"
#include "AliEMCALHit.h"
#include "AliEMCALGeometry.h"
#include "AliRun.h"
#include "AliMC.h"

ClassImp(AliEMCALv1)


//______________________________________________________________________
AliEMCALv1::AliEMCALv1()
  : AliEMCALv0(), 
    fCurPrimary(-1), 
    fCurParent(-1), 
    fCurTrack(-1),
    fTimeCut(30e-09)
{
  // default ctor
}

//______________________________________________________________________
AliEMCALv1::AliEMCALv1(const char *name, const char *title)
  : AliEMCALv0(name,title), 
    fCurPrimary(-1), 
    fCurParent(-1), 
    fCurTrack(-1), 
    fTimeCut(30e-09)
{
    // Standard Creator.

    fHits= new TClonesArray("AliEMCALHit",1000);
    //    gAlice->GetMCApp()->AddHitList(fHits); // 20-dec-04 - advice of Andreas

    fNhits = 0;
    fIshunt     =  2; // All hits are associated with particles entering the calorimeter
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
    // An EMCAL hit is the sum of all hits in a tower section
    //   originating from the same entering particle 
    Int_t hitCounter;
    
    AliEMCALHit *newHit;
    AliEMCALHit *curHit;
    Bool_t deja = kFALSE;

    newHit = new AliEMCALHit(shunt, primary, tracknumber, iparent, ienergy, id, hits, p);
     for ( hitCounter = fNhits-1; hitCounter >= 0 && !deja; hitCounter-- ) {
	curHit = (AliEMCALHit*) (*fHits)[hitCounter];
	// We add hits with the same tracknumber, while GEANT treats
	// primaries succesively
	if(curHit->GetPrimary() != primary) 
	  break;
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
  // Accumulates hits as long as the track stays in a tower

  Int_t          id[2];           // (phi, Eta) indices
  // position wrt MRS and energy deposited
  Float_t        xyzte[5]={0.,0.,0.,0.,0.};// position wrt MRS, time and energy deposited
  Float_t        pmom[4]={0.,0.,0.,0.};
  TLorentzVector pos; // Lorentz vector of the track current position.
  TLorentzVector mom; // Lorentz vector of the track current momentum.
  Int_t tracknumber =  gAlice->GetMCApp()->GetCurrentTrackNumber();
  static Float_t ienergy = 0;
  Int_t copy = 0;
  
  AliEMCALGeometry * geom = GetGeometry() ; 

  static Int_t idXPHI = gMC->VolId("XPHI");
  if(gMC->CurrentVolID(copy) == idXPHI ) { // We are in a Scintillator Layer 
    Float_t depositedEnergy ; 
    
    if( ((depositedEnergy = gMC->Edep()) > 0.)  && (gMC->TrackTime() < fTimeCut)){// Track is inside a scintillator and deposits some energy
       if (fCurPrimary==-1) 
	fCurPrimary=gAlice->GetMCApp()->GetPrimary(tracknumber);

      if (fCurParent==-1 || tracknumber != fCurTrack) {
	// Check parentage
	Int_t parent=tracknumber;
	if (fCurParent != -1) {
	  while (parent != fCurParent && parent != -1) {
	    TParticle *part=gAlice->GetMCApp()->Particle(parent);
	    parent=part->GetFirstMother();
	  }
	}
	if (fCurParent==-1 || parent==-1) {
	  Int_t parent=tracknumber;
	  TParticle *part=gAlice->GetMCApp()->Particle(parent);
	  while (parent != -1 && geom->IsInEMCAL(part->Vx(),part->Vy(),part->Vz())) {
	    parent=part->GetFirstMother();
	    if (parent!=-1) 
	      part=gAlice->GetMCApp()->Particle(parent);
	  } 
	  fCurParent=parent;
	  if (fCurParent==-1)
	    Error("StepManager","Cannot find parent");
	  else {
	    TParticle *part=gAlice->GetMCApp()->Particle(fCurParent);
	    ienergy = part->Energy(); 
	  }
	  while (parent != -1) {
	    part=gAlice->GetMCApp()->Particle(parent);
	    part->SetBit(kKeepBit);
	    parent=part->GetFirstMother();
	  }
	}
	fCurTrack=tracknumber;
      }    
      gMC->TrackPosition(pos);
      xyzte[0] = pos[0];
      xyzte[1] = pos[1];
      xyzte[2] = pos[2];
      xyzte[3] = gMC->TrackTime() ;       
      
      gMC->TrackMomentum(mom);
      pmom[0] = mom[0];
      pmom[1] = mom[1];
      pmom[2] = mom[2];
      pmom[3] = mom[3];
      
      gMC->CurrentVolOffID(1, id[0]); // get the POLY copy number;
      gMC->CurrentVolID(id[1]); // get the phi number inside the layer
      
      Int_t tower = (id[0]-1) % geom->GetNZ() + 1 + (id[1] - 1) * geom->GetNZ() ;  
      Int_t layer = static_cast<Int_t>((id[0]-1)/(geom->GetNZ())) + 1 ; 
      Int_t absid = tower; 
      Int_t nlayers = geom->GetNECLayers();
      if ((layer > nlayers)||(layer<1)) 
        Fatal("StepManager", "Wrong calculation of layer number: layer = %d > %d\n", layer, nlayers) ;

      Float_t lightYield =  depositedEnergy ;
      // Apply Birk's law (copied from G3BIRK)

      if (gMC->TrackCharge()!=0) { // Check
	  Float_t birkC1Mod = 0;
	if (fBirkC0==1){ // Apply correction for higher charge states
	  if (TMath::Abs(gMC->TrackCharge())>=2)
	    birkC1Mod=fBirkC1*7.2/12.6;
	  else
	    birkC1Mod=fBirkC1;
	}
	Float_t dedxcm;
	if (gMC->TrackStep()>0) 
	  dedxcm=1000.*gMC->Edep()/gMC->TrackStep();
	else
	  dedxcm=0;
	lightYield=lightYield/(1.+birkC1Mod*dedxcm+fBirkC2*dedxcm*dedxcm);
      } 

      // use sampling fraction to get original energy --HG
      xyzte[4] = lightYield * geom->GetSampling();
        
      if (gDebug == 2) 
	printf("StepManager: id0 = %d, id1 = %d, absid = %d tower = %d layer = %d energy = %f\n", id[0], id[1], absid, tower, layer, xyzte[4]) ;

      AddHit(fIshunt, fCurPrimary,tracknumber, fCurParent, ienergy, absid,  xyzte, pmom);
    } // there is deposited energy
  }
}

void AliEMCALv1::RemapTrackHitIDs(Int_t *map) {
  // remap track index numbers for primary and parent indices
  // (Called by AliStack::PurifyKine)
  if (Hits()==0)
    return;
  TIter hitIter(Hits());
  Int_t iHit=0;
  while (AliEMCALHit *hit=dynamic_cast<AliEMCALHit*>(hitIter()) ) {
    if (map[hit->GetIparent()]==-99)
      cout << "Remapping, found -99 for parent id " << hit->GetIparent() << ", " << map[hit->GetIparent()] << ", iHit " << iHit << endl;
    hit->SetIparent(map[hit->GetIparent()]);
    if (map[hit->GetPrimary()]==-99)
      cout << "Remapping, found -99 for primary id " << hit->GetPrimary() << ", " << map[hit->GetPrimary()] << ", iHit " << iHit << endl;
    hit->SetPrimary(map[hit->GetPrimary()]);
    iHit++;
  }
}

void AliEMCALv1::FinishPrimary() {
  // finish primary
  fCurPrimary=-1;
  fCurParent=-1;
  fCurTrack=-1;
}
