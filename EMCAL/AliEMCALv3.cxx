/**************************************************************************
 * Copyright(c) 1998-2004, ALICE Experiment at CERN, All rights reserved. *
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
//*-- Implementation version v2 of EMCAL Manager class; SHASHLYK version
//*-- An object of this class does not produce digits
//*-- It is the one to use if you do want to produce outputs in TREEH 
//*--                  
//*-- Author : Aleksei Pavlinov (WSU)

// This Class not stores information on all particles prior to EMCAL entry - in order to facilitate analysis.
// This is done by setting fIShunt =2, and flagging all parents of particles entering the EMCAL.

// --- ROOT system ---
#include "AliEMCALv3.h"
#include "AliEMCALHitv1.h"

#include "TParticle.h"
#include "TVirtualMC.h"
#include "TBrowser.h"
#include "TH2.h"
#include <cassert>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliEMCALGeometry.h"
#include "AliRun.h"
#include "AliHeader.h"
#include "AliMC.h"
#include "AliStack.h"
// for TRD1,2 case

ClassImp(AliEMCALv3)

  //extern "C" void gdaxis_(float &x0, float &y0, float &z0, float &axsiz);


//______________________________________________________________________
AliEMCALv3::AliEMCALv3():AliEMCALv1(){
  // ctor

}

//______________________________________________________________________
AliEMCALv3::AliEMCALv3(const char *name, const char *title): AliEMCALv1(name,title) {
    // Standard Creator.

    fHits= new TClonesArray("AliEMCALHitv1",2000);
    gAlice->GetMCApp()->AddHitList(fHits);

    fNhits    = 0;
    fIshunt   = 2; // All hits are associated with particles entering the calorimeter
    fTimeCut  = 30e-09;

    fGeometry = GetGeometry(); 
}

//______________________________________________________________________
AliEMCALv3::AliEMCALv3(const AliEMCALv3 & emcal):AliEMCALv1(emcal)
{
  fGeometry = emcal.fGeometry;
}

//______________________________________________________________________
AliEMCALv3::~AliEMCALv3(){
    // dtor

    if ( fHits) {
	fHits->Delete();
	delete fHits;
	fHits = 0;
    }
}

//______________________________________________________________________
void AliEMCALv3::AddHit(Int_t shunt, Int_t primary, Int_t tracknumber, Int_t iparent, Float_t ienergy, 
			Int_t id, Float_t * hits,Float_t * p){
    // Add a hit to the hit list.
    // An EMCAL hit is the sum of all hits in a tower section
    //   originating from the same entering particle 
    static Int_t hitCounter;
    static AliEMCALHitv1 *newHit, *curHit;
    static Bool_t deja;

    deja = kFALSE;

    newHit = new AliEMCALHitv1(shunt, primary, tracknumber, iparent, ienergy, id, hits, p);
    for ( hitCounter = fNhits-1; hitCounter >= 0 && !deja; hitCounter-- ) {
	curHit = (AliEMCALHitv1*) (*fHits)[hitCounter];
	// We add hits with the same tracknumber, while GEANT treats
	// primaries succesively
	if(curHit->GetPrimary() != primary) 
	  break;
	if( *curHit == *newHit ) {
	    *curHit = *curHit + *newHit;
	    deja = kTRUE;
	    //            break; // 30-aug-04 by PAI 
	} // end if
    } // end for hitCounter
    
    if ( !deja ) {
	new((*fHits)[fNhits]) AliEMCALHitv1(*newHit);
	fNhits++;
    }
    //    printf(" fNhits %i \n", fNhits); 
    delete newHit;
}
//______________________________________________________________________
void AliEMCALv3::StepManager(void){
  // Accumulates hits as long as the track stays in a tower

  // position wrt MRS and energy deposited
  static Float_t        xyzte[5]={0.,0.,0.,0.,0.};// position wrt MRS, time and energy deposited
  static Float_t        pmom[4]={0.,0.,0.,0.};
  static TLorentzVector pos; // Lorentz vector of the track current position.
  static TLorentzVector mom; // Lorentz vector of the track current momentum.
  static Float_t ienergy = 0;
  static TString curVolName;
  static int supModuleNumber, moduleNumber, yNumber, xNumber, absid;
  static int keyGeom=0;
  static char *vn = "SX"; // 15-mar-05
  static int nSMOP[7]={1,3,5,7,9,11}; // 30-mar-05
  static int nSMON[7]={2,4,6,8,10,12};
  static Float_t depositedEnergy=0.0; 

  if(keyGeom == 0) {
    keyGeom = 2;
    if(gMC->VolId("PBMO")==0 || gMC->VolId("WSUC")==1) {
      vn      = "SCMX";   // old TRD2(TRD1) or WSUC
      keyGeom = 1;
    }    
    printf("AliEMCALv3::StepManager():  keyGeom %i : Sensetive volume %s \n", 
    keyGeom, vn); 
    if(gMC->VolId("WSUC")==1) printf(" WSUC - cosmic ray stand geometry \n");
  }
  Int_t tracknumber =  gAlice->GetMCApp()->GetCurrentTrackNumber();

  curVolName = gMC->CurrentVolName();
  if(curVolName.Contains(vn)) { // We are in a scintillator layer
    //    printf(" keyGeom %i : Sensetive volume %s (%s) \n", keyGeom, curVolName.Data(), vn); 
    
    if( ((depositedEnergy = gMC->Edep()) > 0.)  && (gMC->TrackTime() < fTimeCut)){// Track is inside a scintillator and deposits some energy
      //       Info("StepManager "," entry %i DE %f",++ientry, depositedEnergy); // for testing
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
	  while (parent != -1 && fGeometry->IsInEMCAL(part->Vx(),part->Vy(),part->Vz())) {
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
      
      //      if(ientry%200 > 0) return; // testing
      supModuleNumber = moduleNumber = yNumber = xNumber = absid = 0;
      if(keyGeom >= 1) { // old style
        gMC->CurrentVolOffID(4, supModuleNumber);
        gMC->CurrentVolOffID(3, moduleNumber);
        gMC->CurrentVolOffID(1, yNumber);
        gMC->CurrentVolOffID(0, xNumber); // really x number now
        if(strcmp(gMC->CurrentVolOffName(4),"SM10")==0) supModuleNumber += 10; // 13-oct-05
      } else {
        gMC->CurrentVolOffID(5, supModuleNumber);
        gMC->CurrentVolOffID(4, moduleNumber);
        gMC->CurrentVolOffID(1, yNumber);
        gMC->CurrentVolOffID(0, xNumber);
        if     (strcmp(gMC->CurrentVolOffName(5),"SMOP")==0) supModuleNumber = nSMOP[supModuleNumber-1];
        else if(strcmp(gMC->CurrentVolOffName(5),"SMON")==0) supModuleNumber = nSMON[supModuleNumber-1];
        else   assert(0); // something wrong
      }
      absid = fGeometry->GetAbsCellId(supModuleNumber, moduleNumber, yNumber, xNumber);
    
      if (absid == -1) Fatal("StepManager()", "Wrong id ") ;

      Float_t lightYield =  depositedEnergy ;
      // Apply Birk's law (copied from G3BIRK)

      if (gMC->TrackCharge()!=0) { // Check
	  Float_t BirkC1_mod = 0;
	if (fBirkC0==1){ // Apply correction for higher charge states
	  if (TMath::Abs(gMC->TrackCharge())>=2) BirkC1_mod=fBirkC1*7.2/12.6;
	  else                                    BirkC1_mod=fBirkC1;
	}

	Float_t dedxcm;
	if (gMC->TrackStep()>0)  dedxcm=1000.*gMC->Edep()/gMC->TrackStep();
	else                     dedxcm=0;
	lightYield=lightYield/(1.+BirkC1_mod*dedxcm+fBirkC2*dedxcm*dedxcm);
      } 

      xyzte[4] = lightYield; // 15-dec-04
      Float_t xd[5]; // see Gmtod
      gMC->Gmtod(xyzte, xd, 1);
      xd[3] = xyzte[3];
      xd[4] = xyzte[4];
      //      printf(" x %7.2f y %7.2f z %7.2f de %f\n", xd[0], xd[1], xd[2], xd[4]);
        
      if (gDebug == -2) 
      printf("#sm %2i #m %3i #x %1i #z %1i -> absid %i : xyzte[4] = %f\n",
      supModuleNumber,moduleNumber,yNumber,xNumber,absid, xd[4]);

      AddHit(fIshunt, fCurPrimary,tracknumber, fCurParent, ienergy, absid,  xd, pmom);
    } // there is deposited energy
  }
}

//_____________________________________________
void AliEMCALv3::Browse(TBrowser* b)
{
  TObject::Browse(b);
}
