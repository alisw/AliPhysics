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
// Layout EMC + PPSD has name GPS2  
// Produces cumulated hits (no hits) and digits                  
//*-- Author: Yves Schutz (SUBATECH)


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

#include "AliPHOSv1.h"
#include "AliPHOSHit.h"
#include "AliPHOSDigit.h"
#include "AliPHOSReconstructioner.h"
#include "AliRun.h"
#include "AliConst.h"

ClassImp(AliPHOSv1)

//____________________________________________________________________________
AliPHOSv1::AliPHOSv1()
{
  // ctor
  fNTmpHits = 0 ; 
  fTmpHits  = 0 ; 


}

//____________________________________________________________________________
AliPHOSv1::AliPHOSv1(const char *name, const char *title):
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
  fDigitThreshold      = 0.1 ;   // 1 GeV 

  // We do not want to save in TreeH the raw hits
  // But save the cumulated hits instead (need to create the branch myself)
  // It is put in the Digit Tree because the TreeH is filled after each primary
 // and the TreeD at the end of the event (branch is set in FinishEvent() ).
  
  fTmpHits= new TClonesArray("AliPHOSHit",1000) ;

  fNTmpHits = fNhits = 0 ;

  fDigits = new TClonesArray("AliPHOSDigit",1000) ;


  fIshunt     =  1 ; // All hits are associated with primary particles
 
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
  AliPHOSHit  * hit ;
  AliPHOSDigit * newdigit ;
  AliPHOSDigit * curdigit ;
  Bool_t deja = kFALSE ; 
  
  for ( i = 0 ; i < fNTmpHits ; i++ ) {
    hit = (AliPHOSHit*)fTmpHits->At(i) ;

    // Assign primary number only if contribution is significant
    if( hit->GetEnergy() > fDigitThreshold)
      newdigit = new AliPHOSDigit( hit->GetPrimary(), hit->GetId(), Digitize( hit->GetEnergy() ) ) ;
    else
      newdigit = new AliPHOSDigit( -1 , hit->GetId(), Digitize( hit->GetEnergy() ) ) ;
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

  if ( fRecParticles && gAlice->TreeR() ) { 
     sprintf(branchname,"%sRP",GetName()) ;
     gAlice->TreeR()->Branch(branchname, &fRecParticles, fBufferSize) ;
  }

  
  // 3.
  fReconstructioner->Make(fDigits, fEmcRecPoints, fPpsdRecPoints, fTrackSegments, fRecParticles);

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
}

//____________________________________________________________________________

void AliPHOSv1::StepManager(void)
{
  // Accumulates hits as long as the track stays in a single crystal or PPSD gas Cell

  Int_t          relid[4] ;      // (box, layer, row, column) indices
  Float_t        xyze[4] ;       // position wrt MRS and energy deposited
  TLorentzVector pos ;
  Int_t copy ;

  Int_t tracknumber =  gAlice->CurrentTrack() ; 
  Int_t primary =  gAlice->GetPrimary( gAlice->CurrentTrack() ); 
  TString name = fGeom->GetName() ; 
  if ( name == "GPS2" ) { // the CPV is a PPSD
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

	Int_t absid ; 
       	fGeom->RelToAbsNumbering(relid, absid) ; 

	// add current hit to the hit list      
	AddHit(fIshunt, primary, tracknumber, absid, xyze);

      } // there is deposited energy 
     } // We are inside the gas of the CPV  
   } // GPS2 configuration
  
   if(gMC->CurrentVolID(copy) == gMC->VolId("PXTL") )  //  We are inside a PBWO crystal
     {
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
	  
          Int_t absid ; 
          fGeom->RelToAbsNumbering(relid, absid) ; 
 
	  // add current hit to the hit list

	  AddHit(fIshunt, primary,tracknumber, absid, xyze);
	  
       } // there is deposited energy
     } // we are inside a PHOS Xtal
}

