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
// Implementation version v1 of PHOS Manager class 
//---
// Layout EMC + PPSD has name GPS2:
// Produces cumulated hits (no hits) and digits                  
//---
// Layout EMC + CPV  has name IHEP:
// Produces hits for CPV, cumulated hits and digits                  
//*-- Author: Yves Schutz (SUBATECH)


// --- ROOT system ---

#include "TBRIK.h"
#include "TNode.h"
#include "TRandom.h"
#include "TTree.h"


// --- Standard library ---

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <strstream.h>

// --- AliRoot header files ---

#include "AliPHOSv1.h"
#include "AliPHOSHit.h"
#include "AliPHOSDigit.h"
#include "AliPHOSReconstructioner.h"
#include "AliRun.h"
#include "AliConst.h"
#include "AliMC.h"

ClassImp(AliPHOSv1)

//____________________________________________________________________________
AliPHOSv1::AliPHOSv1()
{
  // ctor
  fNTmpHits = 0 ; 
  fTmpHits  = 0 ;

  // Create an empty array of AliPHOSCPVModule to satisfy
  // AliPHOSv1::Streamer when reading root file

  if ( NULL==(fCPVModules=new TClonesArray("AliPHOSCPVModule",0)) ) {
    Error("AliPHOSv1","Can not create array of CPV modules");
    exit(1);
  }

}

//____________________________________________________________________________
AliPHOSv1::AliPHOSv1(const char *name, const char *title):
AliPHOSv0(name,title) 
{
  // ctor : title is used to identify the layout
  //        GPS2 = 5 modules (EMC + PPSD)
  //        IHEP = 5 modules (EMC + CPV )
  // We use 2 arrays of hits :
  //
  //   - fHits (the "normal" one), which retains the hits associated with
  //     the current primary particle being tracked
  //     (this array is reset after each primary has been tracked).
  //
  //   - fTmpHits, which retains all the hits of the current event. It 
  //     is used for the digitization part.

  fPinElectronicNoise = 0.010 ;
  fDigitThreshold      = 0.1 ;   // 1 GeV 

  // We do not want to save in TreeH the raw hits
  // But save the cumulated hits instead (need to create the branch myself)
  // It is put in the Digit Tree because the TreeH is filled after each primary
  // and the TreeD at the end of the event (branch is set in FinishEvent() ).
  
  fTmpHits= new TClonesArray("AliPHOSHit",1000) ;

  fNTmpHits = fNhits = 0 ;

  fDigits = new TClonesArray("AliPHOSDigit",1000) ;

  fIshunt     =  1 ; // All hits are associated with primary particles
  
  // Create array of CPV modules for the IHEP's version of CPV

  if ( strcmp(fGeom->GetName(),"IHEP") == 0 ) {
    // Create array of AliPHOSCPVmodule of the size of PHOS modules number

    if ( NULL==(fCPVModules=new TClonesArray("AliPHOSCPVModule",fGeom->GetNModules())) ) {
      Error("AliPHOSv1","Can not create array of CPV modules");
      exit(1);
    }
    TClonesArray &lcpvmodule = *fCPVModules;
    for (Int_t i=0; i<fGeom->GetNModules(); i++) new(lcpvmodule[i]) AliPHOSCPVModule();
  }
  else {
    // Create an empty array of AliPHOSCPVModule to satisfy
    // AliPHOSv1::Streamer when writing root file

    fCPVModules=new TClonesArray("AliPHOSCPVModule",0);

  }
}

//____________________________________________________________________________
AliPHOSv1::AliPHOSv1(AliPHOSReconstructioner * Reconstructioner, const char *name, const char *title):
  AliPHOSv0(name,title)
{
  // ctor : title is used to identify the layout
  //        GPS2 = 5 modules (EMC + PPSD)   
  // We use 2 arrays of hits :
  //
  //   - fHits (the "normal" one), which retains the hits associated with
  //     the current primary particle being tracked
  //     (this array is reset after each primary has been tracked).
  //
  //   - fTmpHits, which retains all the hits of the current event. It 
  //     is used for the digitization part.

  fPinElectronicNoise = 0.010 ;

  // We do not want to save in TreeH the raw hits
  //fHits   = new TClonesArray("AliPHOSHit",100) ;

  fDigits = new TClonesArray("AliPHOSDigit",1000) ;
  fTmpHits= new TClonesArray("AliPHOSHit",1000) ;

  fNTmpHits = fNhits = 0 ;

  fIshunt     =  1 ; // All hits are associated with primary particles
 
  // gets an instance of the geometry parameters class  
  fGeom =  AliPHOSGeometry::GetInstance(title, "") ; 

  if (fGeom->IsInitialized() ) 
    cout << "AliPHOS" << Version() << " : PHOS geometry intialized for " << fGeom->GetName() << endl ;
  else
    cout << "AliPHOS" << Version() << " : PHOS geometry initialization failed !" << endl ;   

  // Defining the PHOS Reconstructioner
 
 fReconstructioner = Reconstructioner ;

}

//____________________________________________________________________________
AliPHOSv1::~AliPHOSv1()
{
  // dtor

  if ( fTmpHits) {
    fTmpHits->Delete() ; 
    delete fTmpHits ;
    fTmpHits = 0 ; 
  }

//   if ( fEmcRecPoints ) {
//     fEmcRecPoints->Delete() ; 
//     delete fEmcRecPoints ; 
//     fEmcRecPoints = 0 ; 
//   }

//   if ( fPpsdRecPoints ) { 
//     fPpsdRecPoints->Delete() ;
//     delete fPpsdRecPoints ;
//     fPpsdRecPoints = 0 ; 
//   }
  
//   if ( fTrackSegments ) {
//     fTrackSegments->Delete() ; 
//     delete fTrackSegments ;
//     fTrackSegments = 0 ; 
//   }

}

//____________________________________________________________________________
void AliPHOSv1::AddHit(Int_t shunt, Int_t primary, Int_t tracknumber, Int_t Id, Float_t * hits)
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

  newHit = new AliPHOSHit(shunt, primary, tracknumber, Id, hits) ;

  // We do not want to save in TreeH the raw hits 
  //  TClonesArray &lhits = *fHits;

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

  // We do not want to save in TreeH the raw hits 
  //   new(lhits[fNhits]) AliPHOSHit(*newHit) ;    
  //   fNhits++ ;

  // Please note that the fTmpHits array must survive up to the
  // end of the events, so it does not appear e.g. in ResetHits() (
  // which is called at the end of each primary).  

  delete newHit;

}

//___________________________________________________________________________
Int_t AliPHOSv1::Digitize(Float_t Energy)
{
  // Applies the energy calibration
  
  Float_t fB = 100000000. ;
  Float_t fA = 0. ;
  Int_t chan = Int_t(fA + Energy*fB ) ;
  return chan ;
}

//___________________________________________________________________________
void AliPHOSv1::FinishEvent()
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
  AliPHOSHit   * hit ;
  AliPHOSDigit * newdigit ;
  AliPHOSDigit * curdigit ;
  Bool_t deja = kFALSE ; 
  
  for ( i = 0 ; i < fNTmpHits ; i++ ) {
    hit = (AliPHOSHit*)fTmpHits->At(i) ;

    // Assign primary number only if contribution is significant
    if( hit->GetEnergy() > fDigitThreshold)
      newdigit = new AliPHOSDigit( hit->GetPrimary(), hit->GetId(), Digitize( hit->GetEnergy() ) ) ;
    else
      newdigit = new AliPHOSDigit( -1               , hit->GetId(), Digitize( hit->GetEnergy() ) ) ;
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

//      fGeom->AbsToRelNumbering(newdigit->GetId(), relid) ;
//      printf("FinishEvent(): relid=(%d,%d,%d,%d) Amp=%d\n",
//  	   relid[0],relid[1],relid[2],relid[3], newdigit->GetAmp());
  }

}

//___________________________________________________________________________
void AliPHOSv1::MakeBranch(Option_t* opt)
{  
  // Create new branche in the current Root Tree in the digit Tree
  AliDetector::MakeBranch(opt) ;
  
  char branchname[10];
  sprintf(branchname,"%s",GetName());
  char *cdD = strstr(opt,"D");
  if (fDigits && gAlice->TreeD() && cdD) {
    gAlice->TreeD()->Branch(branchname, &fDigits, fBufferSize);
  }

  // Create new branche PHOSCH in the current Root Tree in the digit Tree for accumulated Hits
  if ( ! (gAlice->IsLegoRun()) ) { // only when not in lego plot mode 
    if ( fTmpHits && gAlice->TreeD()  && cdD) {
      char branchname[10] ;
      sprintf(branchname, "%sCH", GetName()) ;
      gAlice->TreeD()->Branch(branchname, &fTmpHits, fBufferSize) ;
    }   
  }

  // Create new branches CPV<i> for hits in CPV modules for IHEP geometry
  // Yuri Kharlov, 28 September 2000.

  if ( strcmp(fGeom->GetName(),"IHEP") == 0 ) {
    for( Int_t i=0; i<fGeom->GetNModules(); i++ ) GetCPVModule(i).MakeBranch(i+1);
  }

}

//_____________________________________________________________________________
void AliPHOSv1::Reconstruction(AliPHOSReconstructioner * Reconstructioner)
{ 
  // 1. Reinitializes the existing RecPoint, TrackSegment, and RecParticles Lists and 
  // 2. Creates TreeR with a branch for each list
  // 3. Steers the reconstruction processes
  // 4. Saves the 3 lists in TreeR
  // 5. Write the Tree to File
  
  fReconstructioner = Reconstructioner ;
  
  char branchname[10] ;
  
  // 1.

  //  gAlice->MakeTree("R") ; 
  Int_t splitlevel = 0 ; 
  
  fEmcRecPoints->Delete() ; 

  if ( fEmcRecPoints && gAlice->TreeR() ) {
    sprintf(branchname,"%sEmcRP",GetName()) ;
    gAlice->TreeR()->Branch(branchname, "TObjArray", &fEmcRecPoints, fBufferSize, splitlevel) ; 
  }

  fPpsdRecPoints->Delete() ; 

  if ( fPpsdRecPoints && gAlice->TreeR() ) {
    sprintf(branchname,"%sPpsdRP",GetName()) ;
    gAlice->TreeR()->Branch(branchname, "TObjArray", &fPpsdRecPoints, fBufferSize, splitlevel) ;
  }

  fTrackSegments->Delete() ; 

  if ( fTrackSegments && gAlice->TreeR() ) { 
    sprintf(branchname,"%sTS",GetName()) ;
    gAlice->TreeR()->Branch(branchname, &fTrackSegments, fBufferSize) ;
  }

  fRecParticles->Delete() ; 

  if      (strcmp(fGeom->GetName(),"GPS2") == 0) {
    if ( fRecParticles && gAlice->TreeR() ) { 
      sprintf(branchname,"%sRP",GetName()) ;
      gAlice->TreeR()->Branch(branchname, &fRecParticles, fBufferSize) ;
    }
  }
  
  // 3.
  if      (strcmp(fGeom->GetName(),"GPS2") == 0)
    fReconstructioner->Make(fDigits, fEmcRecPoints, fPpsdRecPoints, fTrackSegments, fRecParticles);
  else if (strcmp(fGeom->GetName(),"IHEP") == 0)
    fReconstructioner->Make(fDigits, fEmcRecPoints, fPpsdRecPoints);

  // 4. Expand or Shrink the arrays to the proper size
  
  Int_t size ;
  
  size = fEmcRecPoints->GetEntries() ;
  fEmcRecPoints->Expand(size) ;
 
  size = fPpsdRecPoints->GetEntries() ;
  fPpsdRecPoints->Expand(size) ;

  size = fTrackSegments->GetEntries() ;
  fTrackSegments->Expand(size) ;

  size = fRecParticles->GetEntries() ;
  fRecParticles->Expand(size) ;

  gAlice->TreeR()->Fill() ;
  // 5.

  gAlice->TreeR()->Write(0,TObject::kOverwrite) ;
 
  // Deleting reconstructed objects
  ResetReconstruction();
  
}

//____________________________________________________________________________
void AliPHOSv1::ResetHits() 
{ 
  // Reset hit tree for CPV in IHEP geometry
  // Yuri Kharlov, 28 September 2000

  AliDetector::ResetHits();
  if ( strcmp(fGeom->GetName(),"IHEP") == 0 ) {
    for (Int_t i=0; i<fGeom->GetNModules(); i++) ((AliPHOSCPVModule*)(*fCPVModules)[i]) -> Clear();
  }
}  
//____________________________________________________________________________
void AliPHOSv1::ResetDigits() 
{ 
  // May sound strange, but cumulative hits are store in digits Tree
  AliDetector::ResetDigits();
  if(  fTmpHits ) {
    fTmpHits->Delete();
    fNTmpHits = 0 ;
  }
}  
//____________________________________________________________________________
void AliPHOSv1::ResetReconstruction() 
{ 
  // Deleting reconstructed objects

  if ( fEmcRecPoints )   fEmcRecPoints->Delete();
  if ( fPpsdRecPoints )  fPpsdRecPoints->Delete();
  if ( fTrackSegments )  fTrackSegments->Delete();
  if ( fRecParticles )   fRecParticles->Delete();
  
}

//____________________________________________________________________________
void AliPHOSv1::SetTreeAddress()
{ 
  //  TBranch *branch;
  AliPHOS::SetTreeAddress();

 //  //Branch address for TreeR: RecPpsdRecPoint
//   TTree *treeR = gAlice->TreeR();
//   if ( treeR && fPpsdRecPoints ) {
//     branch = treeR->GetBranch("PHOSPpsdRP");
//     if (branch) branch->SetAddress(&fPpsdRecPoints) ;
//  }

  // Set branch address for the Hits Tree for hits in CPV modules for IHEP geometry
  // Yuri Kharlov, 28 September 2000.

  if ( strcmp(fGeom->GetName(),"IHEP") == 0 ) {
    for( Int_t i=0; i<fGeom->GetNModules(); i++ ) GetCPVModule(i).SetTreeAddress(i+1);
  }

}

//____________________________________________________________________________

void AliPHOSv1::StepManager(void)
{
  // Accumulates hits as long as the track stays in a single crystal or PPSD gas Cell

  Int_t          relid[4] ;      // (box, layer, row, column) indices
  Int_t          absid    ;      // absolute cell ID number
  Float_t        xyze[4]  ;      // position wrt MRS and energy deposited
  TLorentzVector pos      ;      // Lorentz vector of the track current position
  Int_t          copy     ;

  Int_t tracknumber =  gAlice->CurrentTrack() ; 
  Int_t primary     =  gAlice->GetPrimary( gAlice->CurrentTrack() ); 
  TString name      =  fGeom->GetName() ; 

  if ( name == "GPS2" ) {                                       // ======> CPV is a GPS' PPSD

    if( gMC->CurrentVolID(copy) == gMC->VolId("GCEL") ) // We are inside a gas cell 
    {
      gMC->TrackPosition(pos) ;
      xyze[0] = pos[0] ;
      xyze[1] = pos[1] ;
      xyze[2] = pos[2] ;
      xyze[3] = gMC->Edep() ; 

      if ( xyze[3] != 0 ) { // there is deposited energy 
       	gMC->CurrentVolOffID(5, relid[0]) ;  // get the PHOS Module number
       	gMC->CurrentVolOffID(3, relid[1]) ;  // get the Micromegas Module number 
      // 1-> Geom->GetNumberOfModulesPhi() *  fGeom->GetNumberOfModulesZ() upper                         
      //  >  fGeom->GetNumberOfModulesPhi()  *  fGeom->GetNumberOfModulesZ() lower
       	gMC->CurrentVolOffID(1, relid[2]) ;  // get the row number of the cell
        gMC->CurrentVolID(relid[3]) ;        // get the column number 

	// get the absolute Id number

       	fGeom->RelToAbsNumbering(relid, absid) ; 

	// add current hit to the hit list      
	AddHit(fIshunt, primary, tracknumber, absid, xyze);

      } // there is deposited energy 
    } // We are inside the gas of the CPV  
  } // GPS2 configuration

  else if ( name == "IHEP" ) {                                  // ======> CPV is a IHEP's one

    // Yuri Kharlov, 28 September 2000

    if( gMC->CurrentVolID(copy) == gMC->VolId("CPVQ") &&
	gMC->IsTrackEntering() &&
	gMC->TrackCharge() != 0) {

      // Charged track has just entered to the CPV sensitive plane
      
      AliPHOSv1 &phos = *(AliPHOSv1*)gAlice->GetModule("PHOS");
      
      Int_t moduleNumber;
      gMC->CurrentVolOffID(3,moduleNumber);
      moduleNumber--;
      
      // Current position of the hit in the CPV module ref. system

      gMC -> TrackPosition(pos);
      Float_t xyzm[3], xyzd[3], xyd[2];
      Int_t i;
      for (i=0; i<3; i++) xyzm[i] = pos[i];
      gMC -> Gmtod (xyzm, xyzd, 1);    // transform coordinate from master to daughter system
      xyd[0]  = xyzd[0];
      xyd[1]  =-xyzd[2];
      
      // Current momentum of the hit's track in the CPV module ref. system
      
      TLorentzVector  pmom;
      gMC -> TrackMomentum(pmom);
      Float_t pm[3], pd[3];
      for (i=0; i<3; i++) pm[i]   = pmom[i];
      gMC -> Gmtod (pm, pd, 2);        // transform 3-momentum from master to daughter system
      pmom[0] = pd[0];
      pmom[1] =-pd[1];
      pmom[2] =-pd[2];

      // Current particle type of the hit's track

      Int_t ipart = gMC->TrackPid();

      // Add the current particle in the list of the CPV hits.

      phos.GetCPVModule(moduleNumber).AddHit(pmom,xyd,ipart);

      if (fDebugLevel == 1) {
	printf("CPV hit added to module #%2d: p = (% .4f, % .4f, % .4f, % .4f) GeV,\n",
	       moduleNumber+1,pmom.Px(),pmom.Py(),pmom.Pz(),pmom.E());
	printf( "                            xy = (%8.4f, %8.4f) cm, ipart = %d\n",
		xyd[0],xyd[1],ipart);
      }

      // Digitize the current CPV hit:

      // 1. find pad response and
      
      TClonesArray *cpvDigits = new TClonesArray("AliPHOSCPVDigit",0);   // array of digits for current hit
      CPVDigitize(pmom,xyd,moduleNumber,cpvDigits);
      
      Float_t xmean = 0;
      Float_t zmean = 0;
      Float_t qsum  = 0;
      Int_t   idigit,ndigits;

      // 2. go through the current digit list and sum digits in pads

      ndigits = cpvDigits->GetEntriesFast();
      for (idigit=0; idigit<ndigits-1; idigit++) {
	AliPHOSCPVDigit  *cpvDigit1 = (AliPHOSCPVDigit*) cpvDigits->UncheckedAt(idigit);
	Float_t x1 = cpvDigit1->GetXpad() ;
	Float_t z1 = cpvDigit1->GetYpad() ;
	for (Int_t jdigit=idigit+1; jdigit<ndigits; jdigit++) {
	  AliPHOSCPVDigit  *cpvDigit2 = (AliPHOSCPVDigit*) cpvDigits->UncheckedAt(jdigit);
	  Float_t x2 = cpvDigit2->GetXpad() ;
	  Float_t z2 = cpvDigit2->GetYpad() ;
	  if (x1==x2 && z1==z2) {
	    Float_t qsum = cpvDigit1->GetQpad() + cpvDigit2->GetQpad() ;
	    cpvDigit2->SetQpad(qsum) ;
	    cpvDigits->RemoveAt(idigit) ;
	  }
	}
      }
      cpvDigits->Compress() ;

      // 3. add digits to temporary hit list fTmpHits

      ndigits = cpvDigits->GetEntriesFast();
      for (idigit=0; idigit<ndigits; idigit++) {
	AliPHOSCPVDigit  *cpvDigit = (AliPHOSCPVDigit*) cpvDigits->UncheckedAt(idigit);
	relid[0] = moduleNumber + 1 ;                             // CPV (or PHOS) module number
	relid[1] =-1 ;                                            // means CPV
	relid[2] = cpvDigit->GetXpad() ;                          // column number of a pad
	relid[3] = cpvDigit->GetYpad() ;                          // row    number of a pad
	
	// get the absolute Id number
	fGeom->RelToAbsNumbering(relid, absid) ; 

	// add current digit to the temporary hit list
	xyze[0] = 0. ;
	xyze[1] = 0. ;
	xyze[2] = 0. ;
	xyze[3] = cpvDigit->GetQpad() ;                           // amplitude in a pad
	primary = -1;                                             // No need in primary for CPV
	AddHit(fIshunt, primary, tracknumber, absid, xyze);

	if (cpvDigit->GetQpad() > 0.02) {
	  xmean += cpvDigit->GetQpad() * (cpvDigit->GetXpad() + 0.5);
	  zmean += cpvDigit->GetQpad() * (cpvDigit->GetYpad() + 0.5);
	  qsum  += cpvDigit->GetQpad();
	}
      }
      delete cpvDigits;
    }
  } // end of IHEP configuration
  
  if(gMC->CurrentVolID(copy) == gMC->VolId("PXTL") ) { //  We are inside a PBWO crystal
    gMC->TrackPosition(pos) ;
    xyze[0] = pos[0] ;
    xyze[1] = pos[1] ;
    xyze[2] = pos[2] ;
    xyze[3] = gMC->Edep() ;
    
    if ( xyze[3] != 0 ) {
      gMC->CurrentVolOffID(10, relid[0]) ; // get the PHOS module number ;
      relid[1] = 0   ;                    // means PBW04
      gMC->CurrentVolOffID(4, relid[2]) ; // get the row number inside the module
      gMC->CurrentVolOffID(3, relid[3]) ; // get the cell number inside the module
      
      // get the absolute Id number
      
      fGeom->RelToAbsNumbering(relid, absid) ; 
      
      // add current hit to the hit list
      
      AddHit(fIshunt, primary,tracknumber, absid, xyze);
      
    } // there is deposited energy
  } // we are inside a PHOS Xtal
}

//____________________________________________________________________________
void AliPHOSv1::CPVDigitize (TLorentzVector p, Float_t *zxhit, Int_t moduleNumber, TClonesArray *cpvDigits)
{
  // ------------------------------------------------------------------------
  // Digitize one CPV hit:
  // On input take exact 4-momentum p and position zxhit of the hit,
  // find the pad response around this hit and
  // put the amplitudes in the pads into array digits
  //
  // Author: Yuri Kharlov (after Serguei Sadovsky)
  // 2 October 2000
  // ------------------------------------------------------------------------

  const Float_t kCelWr  = fGeom->GetPadSizePhi()/2;  // Distance between wires (2 wires above 1 pad)
  const Float_t kDetR   = 0.1;     // Relative energy fluctuation in track for 100 e-
  const Float_t kdEdx   = 4.0;     // Average energy loss in CPV;
  const Int_t   kNgamz  = 5;       // Ionization size in Z
  const Int_t   kNgamx  = 9;       // Ionization size in Phi
  const Float_t kNoise = 0.03;    // charge noise in one pad

  Float_t rnor1,rnor2;

  // Just a reminder on axes notation in the CPV module:
  // axis Z goes along the beam
  // axis X goes across the beam in the module plane
  // axis Y is a normal to the module plane showing from the IP

  Float_t hitX  = zxhit[0];
  Float_t hitZ  =-zxhit[1];
  Float_t pX    = p.Px();
  Float_t pZ    =-p.Pz();
  Float_t pNorm = p.Py();
  Float_t eloss = kdEdx;

  Float_t dZY   = pZ/pNorm * fGeom->GetCPVGasThickness();
  Float_t dXY   = pX/pNorm * fGeom->GetCPVGasThickness();
  gRandom->Rannor(rnor1,rnor2);
  eloss *= (1 + kDetR*rnor1) *
           TMath::Sqrt((1 + ( pow(dZY,2) + pow(dXY,2) ) / pow(fGeom->GetCPVGasThickness(),2)));
  Float_t zhit1 = hitZ + fGeom->GetCPVActiveSize(1)/2 - dZY/2;
  Float_t xhit1 = hitX + fGeom->GetCPVActiveSize(0)/2 - dXY/2;
  Float_t zhit2 = zhit1 + dZY;
  Float_t xhit2 = xhit1 + dXY;

  Int_t   iwht1 = (Int_t) (xhit1 / kCelWr);           // wire (x) coordinate "in"
  Int_t   iwht2 = (Int_t) (xhit2 / kCelWr);           // wire (x) coordinate "out"

  Int_t   nIter;
  Float_t zxe[3][5];
  if (iwht1==iwht2) {                      // incline 1-wire hit
    nIter = 2;
    zxe[0][0] = (zhit1 + zhit2 - dZY*0.57735) / 2;
    zxe[1][0] = (iwht1 + 0.5) * kCelWr;
    zxe[2][0] =  eloss/2;
    zxe[0][1] = (zhit1 + zhit2 + dZY*0.57735) / 2;
    zxe[1][1] = (iwht1 + 0.5) * kCelWr;
    zxe[2][1] =  eloss/2;
  }
  else if (TMath::Abs(iwht1-iwht2) != 1) { // incline 3-wire hit
    nIter = 3;
    Int_t iwht3 = (iwht1 + iwht2) / 2;
    Float_t xwht1 = (iwht1 + 0.5) * kCelWr; // wire 1
    Float_t xwht2 = (iwht2 + 0.5) * kCelWr; // wire 2
    Float_t xwht3 = (iwht3 + 0.5) * kCelWr; // wire 3
    Float_t xwr13 = (xwht1 + xwht3) / 2;   // center 13
    Float_t xwr23 = (xwht2 + xwht3) / 2;   // center 23
    Float_t dxw1  = xhit1 - xwr13;
    Float_t dxw2  = xhit2 - xwr23;
    Float_t egm1  = TMath::Abs(dxw1) / ( TMath::Abs(dxw1) + TMath::Abs(dxw2) + kCelWr );
    Float_t egm2  = TMath::Abs(dxw2) / ( TMath::Abs(dxw1) + TMath::Abs(dxw2) + kCelWr );
    Float_t egm3  =           kCelWr / ( TMath::Abs(dxw1) + TMath::Abs(dxw2) + kCelWr );
    zxe[0][0] = (dXY*(xwr13-xwht1)/dXY + zhit1 + zhit1) / 2;
    zxe[1][0] =  xwht1;
    zxe[2][0] =  eloss * egm1;
    zxe[0][1] = (dXY*(xwr23-xwht1)/dXY + zhit1 + zhit2) / 2;
    zxe[1][1] =  xwht2;
    zxe[2][1] =  eloss * egm2;
    zxe[0][2] =  dXY*(xwht3-xwht1)/dXY + zhit1;
    zxe[1][2] =  xwht3;
    zxe[2][2] =  eloss * egm3;
  }
  else {                                   // incline 2-wire hit
    nIter = 2;
    Float_t xwht1 = (iwht1 + 0.5) * kCelWr;
    Float_t xwht2 = (iwht2 + 0.5) * kCelWr;
    Float_t xwr12 = (xwht1 + xwht2) / 2;
    Float_t dxw1  = xhit1 - xwr12;
    Float_t dxw2  = xhit2 - xwr12;
    Float_t egm1  = TMath::Abs(dxw1) / ( TMath::Abs(dxw1) + TMath::Abs(dxw2) );
    Float_t egm2  = TMath::Abs(dxw2) / ( TMath::Abs(dxw1) + TMath::Abs(dxw2) );
    zxe[0][0] = (zhit1 + zhit2 - dZY*egm1) / 2;
    zxe[1][0] =  xwht1;
    zxe[2][0] =  eloss * egm1;
    zxe[0][1] = (zhit1 + zhit2 + dZY*egm2) / 2;
    zxe[1][1] =  xwht2;
    zxe[2][1] =  eloss * egm2;
  }

  // Finite size of ionization region

  Int_t nCellZ  = fGeom->GetNumberOfPadsZ();
  Int_t nCellX  = fGeom->GetNumberOfPadsPhi();
  Int_t nz3     = (kNgamz+1)/2;
  Int_t nx3     = (kNgamx+1)/2;
  cpvDigits->Expand(nIter*kNgamx*kNgamz);
  TClonesArray &ldigits = *(TClonesArray *)cpvDigits;

  for (Int_t iter=0; iter<nIter; iter++) {

    Float_t zhit = zxe[0][iter];
    Float_t xhit = zxe[1][iter];
    Float_t qhit = zxe[2][iter];
    Float_t zcell = zhit / fGeom->GetPadSizeZ();
    Float_t xcell = xhit / fGeom->GetPadSizePhi();
    if ( zcell<=0      || xcell<=0 ||
	 zcell>=nCellZ || xcell>=nCellX) return;
    Int_t izcell = (Int_t) zcell;
    Int_t ixcell = (Int_t) xcell;
    Float_t zc = zcell - izcell - 0.5;
    Float_t xc = xcell - ixcell - 0.5;
    for (Int_t iz=1; iz<=kNgamz; iz++) {
      Int_t kzg = izcell + iz - nz3;
      if (kzg<=0 || kzg>nCellZ) continue;
      Float_t zg = (Float_t)(iz-nz3) - zc;
      for (Int_t ix=1; ix<=kNgamx; ix++) {
	Int_t kxg = ixcell + ix - nx3;
	if (kxg<=0 || kxg>nCellX) continue;
	Float_t xg = (Float_t)(ix-nx3) - xc;
	
	// Now calculate pad response
	Float_t qpad = CPVPadResponseFunction(qhit,zg,xg);
	qpad += kNoise*rnor2;
	if (qpad<0) continue;
	
	// Fill the array with pad response ID and amplitude
	new(ldigits[cpvDigits->GetEntriesFast()]) AliPHOSCPVDigit(kxg,kzg,qpad);
      }
    }
  }
}

//____________________________________________________________________________
Float_t AliPHOSv1::CPVPadResponseFunction(Float_t qhit, Float_t zhit, Float_t xhit) {
  // ------------------------------------------------------------------------
  // Calculate the amplitude in one CPV pad using the
  // cumulative pad response function
  // Author: Yuri Kharlov (after Serguei Sadovski)
  // 3 October 2000
  // ------------------------------------------------------------------------

  Double_t dz = fGeom->GetPadSizeZ()   / 2;
  Double_t dx = fGeom->GetPadSizePhi() / 2;
  Double_t z  = zhit * fGeom->GetPadSizeZ();
  Double_t x  = xhit * fGeom->GetPadSizePhi();
  Double_t amplitude = qhit *
    (CPVCumulPadResponse(z+dz,x+dx) - CPVCumulPadResponse(z+dz,x-dx) -
     CPVCumulPadResponse(z-dz,x+dx) + CPVCumulPadResponse(z-dz,x-dx));
  return (Float_t)amplitude;
}

//____________________________________________________________________________
Double_t AliPHOSv1::CPVCumulPadResponse(Double_t x, Double_t y) {
  // ------------------------------------------------------------------------
  // Cumulative pad response function
  // It includes several terms from the CF decomposition in electrostatics
  // Note: this cumulative function is wrong since omits some terms
  //       but the cell amplitude obtained with it is correct because
  //       these omitting terms cancel
  // Author: Yuri Kharlov (after Serguei Sadovski)
  // 3 October 2000
  // ------------------------------------------------------------------------

  const Double_t kA=1.0;
  const Double_t kB=0.7;

  Double_t r2       = x*x + y*y;
  Double_t xy       = x*y;
  Double_t cumulPRF = 0;
  for (Int_t i=0; i<=4; i++) {
    Double_t b1 = (2*i + 1) * kB;
    cumulPRF += TMath::Power(-1,i) * TMath::ATan( xy / (b1*TMath::Sqrt(b1*b1 + r2)) );
  }
  cumulPRF *= kA/(2*TMath::Pi());
  return cumulPRF;
}
