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
//*-- Author : Alexei Pavlinov (WSU)

// This Class not stores information on all particles prior to EMCAL entry - in order to facilitate analysis.
// This is done by setting fIShunt =2, and flagging all parents of particles entering the EMCAL.

#include <cassert>
// --- ROOT system ---
#include <TBrowser.h>
#include <TClonesArray.h>
#include <TH2.h>
#include <TParticle.h>
#include <TROOT.h>
#include <TVirtualMC.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliEMCALv2.h"
#include "AliEMCALHit.h"
#include "AliEMCALGeometry.h"
#include "AliRun.h"
#include "AliHeader.h"
#include "AliMC.h"
#include "AliStack.h"
#include "AliTrackReference.h"
// for TRD1 case only; May 31,2006

ClassImp(AliEMCALv2)

//______________________________________________________________________
AliEMCALv2::AliEMCALv2()
  : AliEMCALv1()
{
  // ctor
}

//______________________________________________________________________
AliEMCALv2::AliEMCALv2(const char *name, const char *title)
  : AliEMCALv1(name,title)
{
    // Standard Creator.

    //fHits= new TClonesArray("AliEMCALHit",1000); //Already done in ctor of v1
    gAlice->GetMCApp()->AddHitList(fHits);

    fNhits    = 0;
    fIshunt   = 2; // All hits are associated with particles entering the calorimeter
    fTimeCut  = 30e-09;

    fGeometry = GetGeometry(); 
}

//______________________________________________________________________
AliEMCALv2::~AliEMCALv2(){
    // dtor

  //Already done in dtor of v1
//  if ( fHits ) {
//    fHits->Clear();
//    delete fHits;
//    fHits = 0;
//  }
}

//______________________________________________________________________
void AliEMCALv2::AddHit(Int_t shunt, Int_t primary, Int_t tracknumber, Int_t iparent, Float_t ienergy, 
			Int_t id, Float_t * hits,Float_t * p){
    // Add a hit to the hit list.
    // An EMCAL hit is the sum of all hits in a tower section
    //   originating from the same entering particle 
    static Int_t hitCounter;
    static AliEMCALHit *newHit, *curHit;
    static Bool_t deja ;
  
    deja = kFALSE;

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
	    //            break; // 30-aug-04 by PAI 
	} // end if
    } // end for hitCounter
    
    if ( !deja ) {
	new((*fHits)[fNhits]) AliEMCALHit(*newHit);
	fNhits++;
    }
    //    printf(" fNhits %i \n", fNhits); 
    delete newHit;
}

//______________________________________________________________________
void AliEMCALv2::StepManager(void){
  // Accumulates hits as long as the track stays in a tower

  // position wrt MRS and energy deposited
  static Float_t        xyzte[5]={0.,0.,0.,0.,0.};// position wrt MRS, time and energy deposited
  static Float_t        pmom[4]={0.,0.,0.,0.};
  static TLorentzVector pos;  // Lorentz vector of the track current position.
  static TLorentzVector mom;  // Lorentz vector of the track current momentum.
  static Float_t ienergy = 0; // part->Energy();
  static TString curVolName="";
  static int supModuleNumber=-1, moduleNumber=-1, yNumber=-1, xNumber=-1, absid=-1;
  static int keyGeom=1;  //real TRD1 geometry
  static const char *vn = "SCMX"; // Apr 13, 2006 - only TRD1 case now
  static int nSMOP[7]={1,3,5,7,9,11}; // 30-mar-05
  static int nSMON[7]={2,4,6,8,10,12};
  static Float_t depositedEnergy=0.0; 

  if(keyGeom == 0) {
    keyGeom = 2;
    if(gMC->VolId("PBMO")==0 || gMC->VolId("WSUC")==1) {
      vn      = "SCMX";   // old TRD2(TRD1) or WSUC
      keyGeom = 1;
    }    
    printf("AliEMCALv2::StepManager():  keyGeom %i : Sensetive volume %s \n", 
    keyGeom, vn); 
    if(gMC->VolId("WSUC")==1) printf(" WSUC - cosmic ray stand geometry \n");
  }
  Int_t tracknumber =  gAlice->GetMCApp()->GetCurrentTrackNumber();
  Int_t parent=0;
  TParticle* part=0;

  curVolName = gMC->CurrentVolName();
  if(curVolName.Contains(vn) || curVolName.Contains("SCX")) { // We are in a scintillator layer; SCX for 3X3
    
    if( ((depositedEnergy = gMC->Edep()) > 0.)  && (gMC->TrackTime() < fTimeCut)){// Track is inside a scintillator and deposits some energy
      //       Info("StepManager "," entry %i DE %f",++ientry, depositedEnergy); // for testing
       if (fCurPrimary==-1) 
	fCurPrimary=gAlice->GetMCApp()->GetPrimary(tracknumber);

      if (fCurParent==-1 || tracknumber != fCurTrack) {
	// Check parentage
	parent=tracknumber;

	if (fCurParent != -1) {
	  while (parent != fCurParent && parent != -1) {
	    //TParticle *part=gAlice->GetMCApp()->Particle(parent);
	    part=gAlice->GetMCApp()->Particle(parent);
	    parent=part->GetFirstMother();
	  }
	}
	if (fCurParent==-1 || parent==-1) {
	  //Int_t parent=tracknumber;
	  //TParticle *part=gAlice->GetMCApp()->Particle(parent);
	  parent=tracknumber;
	  part=gAlice->GetMCApp()->Particle(parent);
	  while (parent != -1 && fGeometry->IsInEMCAL(part->Vx(),part->Vy(),part->Vz())) {
	    parent=part->GetFirstMother();
	    if (parent!=-1) 
	      part=gAlice->GetMCApp()->Particle(parent);
	  } 
	  fCurParent=parent;
	  if (fCurParent==-1)
	    Error("StepManager","Cannot find parent");
	  else {
	    //TParticle *part=gAlice->GetMCApp()->Particle(fCurParent);
	    part=gAlice->GetMCApp()->Particle(fCurParent);
	    ienergy = part->Energy(); 

	    //Add reference to parent in TR tree. 	
	    AddTrackReference(tracknumber, AliTrackReference::kEMCAL);

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
      if(keyGeom >= 1) { // TRD1 case now
        gMC->CurrentVolOffID(4, supModuleNumber);
        gMC->CurrentVolOffID(3, moduleNumber);
        gMC->CurrentVolOffID(1, yNumber);
        gMC->CurrentVolOffID(0, xNumber); // really x number now
        if(strcmp(gMC->CurrentVolOffName(4),"SM10")==0) supModuleNumber += 10; // 13-oct-05
	// Nov 10,2006
        if(strcmp(gMC->CurrentVolOffName(0),vn) != 0) { // 3X3 case
          if     (strcmp(gMC->CurrentVolOffName(0),"SCX1")==0) xNumber=1;
          else if(strcmp(gMC->CurrentVolOffName(0),"SCX2")==0) xNumber=2;
          else if(strcmp(gMC->CurrentVolOffName(0),"SCX3")==0) xNumber=3;
          else Fatal("StepManager()", "Wrong name of sensitive volume in 3X3 case : %s ", gMC->CurrentVolOffName(0));
	}
      } else {
        gMC->CurrentVolOffID(5, supModuleNumber);
        gMC->CurrentVolOffID(4, moduleNumber);
        gMC->CurrentVolOffID(1, yNumber);
        gMC->CurrentVolOffID(0, xNumber);
        if     (strcmp(gMC->CurrentVolOffName(5),"SMOP")==0) supModuleNumber = nSMOP[supModuleNumber-1];
        else if(strcmp(gMC->CurrentVolOffName(5),"SMON")==0) supModuleNumber = nSMON[supModuleNumber-1];
        else   assert(0); // something wrong
      }
		
      // Due to problem with index ordering conventions the calcultation of absid is no more like this:	
      //absid = fGeometry->GetAbsCellId(smNumber, moduleNumber-1, yNumber-1, xNumber-1);
      
      //Swap A side in Phi and C side in Eta due to wrong indexing.
      Int_t iphi = -1;
      Int_t ieta = -1;
      Int_t smNumber = supModuleNumber-1;
      Int_t smType   = 1;
      fGeometry->GetCellPhiEtaIndexInSModule(smNumber,moduleNumber-1,yNumber-1,xNumber-1, iphi, ieta);
      if (smNumber%2 == 0) {
	ieta = ((fGeometry->GetCentersOfCellsEtaDir()).GetSize()-1)-ieta;// 47-ieta, revert the ordering on A side in order to keep convention.
      }
      else {  
	if(smNumber >= 10) smType = 2 ; //half supermodule
	iphi= ((fGeometry->GetCentersOfCellsPhiDir()).GetSize()/smType-1)-iphi;//23-iphi, revert the ordering on C side in order to keep convention.
      }
      
      //Once we know the indexes, calculate the absolute ID
      absid = fGeometry->GetAbsCellIdFromCellIndexes(smNumber, iphi, ieta);
      
      if (absid < 0) {
        printf(" supModuleNumber %i : moduleNumber %i : yNumber %i : xNumber %i \n",
	       supModuleNumber, moduleNumber, yNumber, xNumber); 
	Fatal("StepManager()", "Wrong id : %i ", absid) ; 
      }

      Float_t lightYield =  depositedEnergy ;
      // Apply Birk's law (copied from G3BIRK)

      if (gMC->TrackCharge()!=0) { // Check
	  Float_t birkC1Mod = 0;
	if (fBirkC0==1){ // Apply correction for higher charge states
	  if (TMath::Abs(gMC->TrackCharge())>=2) birkC1Mod = fBirkC1*7.2/12.6;
	  else                                   birkC1Mod = fBirkC1;
	}

	Float_t dedxcm=0.;
	if (gMC->TrackStep()>0)  dedxcm=1000.*gMC->Edep()/gMC->TrackStep();
	else                     dedxcm=0;
	lightYield=lightYield/(1.+birkC1Mod*dedxcm+fBirkC2*dedxcm*dedxcm);
      } 

      // use sampling fraction to get original energy --HG
      //      xyzte[4] = lightYield * fGeometry->GetSampling();
      xyzte[4] = lightYield; // 15-dec-04
        
      if (gDebug == -2) 
      printf("#sm %2i #m %3i #x %1i #z %1i -> absid %i : xyzte[4] = %f\n",
      supModuleNumber,moduleNumber,yNumber,xNumber,absid, xyzte[4]);

      AddHit(fIshunt, fCurPrimary,tracknumber, fCurParent, ienergy, absid,  xyzte, pmom);
    } // there is deposited energy
  }
}
