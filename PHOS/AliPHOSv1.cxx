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
// Produces cumulated hits
//---
// Layout EMC + CPV  has name IHEP:
// Produces hits for CPV, cumulated hits
//---
// Layout EMC + CPV + PPSD has name GPS:
// Produces hits for CPV, cumulated hits
//---
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

  // Create an empty array of AliPHOSCPVModule to satisfy
  // AliPHOSv1::Streamer when reading root file
 
  fReconstructioner  = 0;
  fTrackSegmentMaker = 0;

  fHits = new TClonesArray("AliPHOSHit",1000) ;
 
  //  if ( 0==(fEMCModules=new TClonesArray("AliPHOSCPVModule",0)) ) {
  //    Error("AliPHOSv1","Can not create array of EMC modules");
  //    exit(1);
  //  }

  //  if ( 0==(fCPVModules=new TClonesArray("AliPHOSCPVModule",0)) ) {
  //    Error("AliPHOSv1","Can not create array of CPV modules");
  //    exit(1);
  //  }

}

//____________________________________________________________________________
AliPHOSv1::AliPHOSv1(const char *name, const char *title):
AliPHOSv0(name,title) 
{
  // ctor : title is used to identify the layout
  //        GPS2 = 5 modules (EMC + PPSD)
  //        IHEP = 5 modules (EMC + CPV )
  //        MIXT = 4 modules (EMC + CPV ) and 1 module (EMC + PPSD)
  //
  // We store hits :
  //   - fHits (the "normal" one), which retains the hits associated with
  //     the current primary particle being tracked
  //     (this array is reset after each primary has been tracked).
  //

  fPinElectronicNoise = 0.010 ;
  fDigitThreshold      = 0.01 ;   // 1 GeV 
  fDigitizeA= 0. ;             
  fDigitizeB = 10000000. ;    


  // We do not want to save in TreeH the raw hits
  // But save the cumulated hits instead (need to create the branch myself)
  // It is put in the Digit Tree because the TreeH is filled after each primary
  // and the TreeD at the end of the event (branch is set in FinishEvent() ).
  
  fHits= new TClonesArray("AliPHOSHit",1000) ;

  fNhits = 0 ;

  fReconstructioner  = 0;
  fTrackSegmentMaker = 0;

  fIshunt     =  1 ; // All hits are associated with primary particles
  
  // Create array of EMC modules of the size of PHOS modules number
  
//  if ( 0==(fEMCModules=new TClonesArray("AliPHOSCPVModule",fGeom->GetNModules())) ) {
//    Error("AliPHOSv1","Can not create array of EMC modules");
//    exit(1);
//  }
//  TClonesArray &lemcmodule = *fEMCModules;
//  for (Int_t i=0; i<fGeom->GetNModules(); i++) new(lemcmodule[i]) AliPHOSCPVModule();

  // Create array of CPV modules for the IHEP's version of CPV

//  if ( strcmp(fGeom->GetName(),"IHEP") == 0 || strcmp(fGeom->GetName(),"MIXT") == 0 ) {
//    // Create array of CPV modules of the size of PHOS modules number

//    if ( 0==(fCPVModules=new TClonesArray("AliPHOSCPVModule",fGeom->GetNCPVModules())) ) {
//      Error("AliPHOSv1","Can not create array of CPV modules");
//      exit(1);
//    }
//    TClonesArray &lcpvmodule = *fCPVModules;
//    for (Int_t i=0; i<fGeom->GetNCPVModules(); i++) new(lcpvmodule[i]) AliPHOSCPVModule();
//  }
//  else {
//    // Create an empty array of AliPHOSCPVModule to satisfy
//    // AliPHOSv1::Streamer when writing root file

//    fCPVModules=new TClonesArray("AliPHOSCPVModule",0);

//  }
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

  fDigits = 0 ;
  fHits= new TClonesArray("AliPHOSHit",1000) ;

  fNhits = 0 ;

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

  if ( fHits) {
    fHits->Delete() ; 
    delete fHits ;
    fHits = 0 ; 
  }

  if ( fSDigits) {
    fSDigits->Delete() ; 
    delete fSDigits ;
    fSDigits = 0 ; 
  }

  if ( fDigits) {
    fDigits->Delete() ; 
    delete fDigits ;
    fDigits = 0 ; 
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
void AliPHOSv1::AddHit(Int_t shunt, Int_t primary, Int_t tracknumber, Int_t Id, Float_t * hits, Int_t trackpid, TLorentzVector p, Float_t * lpos)
{
  // Add a hit to the hit list.
  // A PHOS hit is the sum of all hits in a single crystal
  //   or in a single PPSD gas cell

  Int_t hitCounter ;
  AliPHOSHit *newHit ;
  AliPHOSHit *curHit ;
  Bool_t deja = kFALSE ;

  newHit = new AliPHOSHit(shunt, primary, tracknumber, Id, hits, trackpid, p, lpos) ;

  for ( hitCounter = fNhits-1 ; hitCounter >= 0 && !deja ; hitCounter-- ) {
    curHit = (AliPHOSHit*) (*fHits)[hitCounter] ;
    if( *curHit == *newHit ) {
      *curHit = *curHit + *newHit ;
      deja = kTRUE ;
    }
  }
         
  if ( !deja ) {
    new((*fHits)[fNhits]) AliPHOSHit(*newHit) ;
    fNhits++ ;
  }

  delete newHit;
}

//____________________________________________________________________________
void AliPHOSv1::Hits2SDigits(){
  //Collects all hits in the same active volume into digit

  Int_t i ;
  Int_t j ; 
  AliPHOSHit   * hit ;
  AliPHOSDigit * newdigit ;
  AliPHOSDigit * curdigit ;
  Bool_t deja = kFALSE ; 
  

  Int_t itrack ;
  for (itrack=0; itrack<gAlice->GetNtrack(); itrack++){
        
    //=========== Get the Hits Tree for the Primary track itrack
    gAlice->ResetHits();    
    gAlice->TreeH()->GetEvent(itrack);
      

    for ( i = 0 ; i < fHits->GetEntries() ; i++ ) {
      hit = (AliPHOSHit*)fHits->At(i) ;
    
      // Assign primary number only if contribution is significant
      if( hit->GetEnergy() > fDigitThreshold)
	newdigit = new AliPHOSDigit( hit->GetPrimary(), hit->GetId(), Digitize( hit->GetEnergy() ) ) ;
      else
	newdigit = new AliPHOSDigit( -1               , hit->GetId(), Digitize( hit->GetEnergy() ) ) ;
      deja =kFALSE ;

      for ( j = 0 ; j < fnSdigits ;  j++) { 
	curdigit = (AliPHOSDigit*) fSDigits->At(j) ;
	if ( *curdigit == *newdigit) {
	  *curdigit = *curdigit + *newdigit ; 
	  deja = kTRUE ; 
	}
      }

      if ( !deja ) {
	new((*fSDigits)[fnSdigits]) AliPHOSDigit(* newdigit) ;
	fnSdigits++ ;  
      }
      
      delete newdigit ;    
    } 

  } // loop over tracks

  fSDigits->Sort() ;

  fnSdigits = fSDigits->GetEntries() ;
  fSDigits->Expand(fnSdigits) ;

  for (i = 0 ; i < fnSdigits ; i++) { 
    AliPHOSDigit * digit = (AliPHOSDigit *) fSDigits->At(i) ; 
    digit->SetIndexInList(i) ;     
  }

  gAlice->TreeS()->Fill() ;
  
  TTree *ts = gAlice->TreeS();
  TDirectory *dir = ts->GetDirectory() ; 
  cout << "dir name is " << dir->GetName() << endl ;  
  //  ts->Write(0,TObject::kOverwrite) ;

gAlice->TreeS()->GetBranch("PHOS")->Fill();
gAlice->TreeS()->GetBranch("PHOS")->Write();


}
//____________________________________________________________________________
void AliPHOSv1::SDigits2Digits(){
  //Adds noise to the summable digits and removes everething below thresholds
  //Note, that sDigits should be SORTED in accordance with abs ID.


  gAlice->TreeS()->GetEvent(0) ;

  cout << "fSdigits " <<  fSDigits << " " << fSDigits->GetEntries() << endl ;
  

  // Noise induced by the PIN diode of the PbWO crystals
  Int_t iCurSDigit = 0 ;
  //we assume, that there is al least one EMC digit...
  Int_t idCurSDigit = ((AliPHOSDigit *)fSDigits->At(0))->GetId() ;

  cout << "fDigits " << fDigits << " " <<idCurSDigit<<  endl ;

  Int_t absID ;
  for(absID = 1; absID < fGeom->GetNModules()*fGeom->GetNPhi()*fGeom->GetNZ(); absID++){
    
    cout << "absID " << absID << " " << idCurSDigit << endl ; 

    Float_t noise = gRandom->Gaus(0., fPinElectronicNoise) ; 
    if(absID < idCurSDigit ){ 
      cout << "In < idC  " << noise << " " << fDigitThreshold << endl ;
      if(noise >fDigitThreshold ){
	cout << "noise " << absID << " " << noise << endl;
	new((*fDigits)[fNdigits]) AliPHOSDigit( -1,absID,Digitize(noise) ) ;
	cout << "FN " << fNdigits << endl ;
	fNdigits++ ;  
      }
    }
    else{ //add noise and may be remove the true hit
      cout << "correcting digit " << iCurSDigit << endl ;
      Float_t signal = noise + Calibrate(((AliPHOSDigit *)fSDigits->At(iCurSDigit))->GetAmp()) ;
      cout << "signal " << signal << endl ;
      if( signal >fDigitThreshold ){
	cout << "signal " << signal << endl ;
	AliPHOSDigit * digit = (AliPHOSDigit*) fSDigits->At(iCurSDigit) ;
	new((*fDigits)[fNdigits]) AliPHOSDigit( *digit ) ;
	((AliPHOSDigit *)fDigits->At(fNdigits))->SetAmp(Digitize(signal));
	cout << "fNdigits " << fNdigits << endl ;
	fNdigits++ ;  
      }

      if(iCurSDigit < fSDigits->GetEntries()-1){
	iCurSDigit++ ;
	idCurSDigit = ((AliPHOSDigit*)fSDigits->At(iCurSDigit))->GetId() ;
      }
      else
	idCurSDigit = 10000000; //no real hits left    
    }
    
  }  

  //remove PPSD/CPV digits below thresholds
  Int_t idigit ;
  for(idigit = iCurSDigit; idigit < fSDigits->GetEntries() ; idigit++){  //loop over CPV/PPSD digits
    
    AliPHOSDigit * digit = (AliPHOSDigit *) fSDigits->At(idigit) ; 
    Float_t ene = Calibrate(digit->GetAmp()) ;
    
    Int_t relid[4] ; 
    fGeom->AbsToRelNumbering(digit->GetId(), relid) ; 
    if ( relid[0] > fGeom->GetNCPVModules() ){ //ppsd
      if ( ( (relid[1] > 0) && (ene > fPpsdEnergyThreshold)) ||    //PPSD digit
	   ( (relid[1] < 0) && (ene > fCpvEnergyThreshold ) ) )	   //CPV digit 
	new((*fDigits)[fNdigits]) AliPHOSDigit( *digit ) ;
	fNdigits++ ;
    }
  }    
    
  fDigits->Compress() ;  
  
  fNdigits = fDigits->GetEntries() ;
  fDigits->Expand(fNdigits) ;

  Int_t i ;
  for (i = 0 ; i < fNdigits ; i++) { 
    AliPHOSDigit * digit = (AliPHOSDigit *) fDigits->At(i) ; 
    digit->SetIndexInList(i) ;     
  }


  gAlice->TreeD()->Fill() ;

  gAlice->TreeD()->Write(0,TObject::kOverwrite) ;
 
}

//___________________________________________________________________________
void AliPHOSv1::MakeBranch(Option_t* opt, char *file)
{ 


  char *cH ; 
  // Create new branche in the current Root Tree in the digit Tree
  AliDetector::MakeBranch(opt) ;


  cH = strstr(opt,"S");
  //Create a branch for SDigits
  if( cH ){
    char branchname[20];
    sprintf(branchname,"%s",GetName());  
    if(fSDigits)
      fSDigits->Clear();
    else
      fSDigits = new TClonesArray("AliPHOSDigit",1000);
    fnSdigits = 0 ;
    cout << " AliPHOSv1::MakeBranch : " << file << endl ; 
    gAlice->MakeBranchInTree(gAlice->TreeS(),branchname,&fSDigits,fBufferSize,file);  
  }
  
  cH = strstr(opt,"D");
  //Create a branch for Digits
  if( cH ){
    char branchname[20];
    sprintf(branchname,"%s",GetName());  
    if(fSDigits)
      fDigits->Clear();
    else
      fDigits = new TClonesArray("AliPHOSDigit",1000);
    fNdigits = 0 ;
    
    gAlice->MakeBranchInTree(gAlice->TreeD(),branchname,&fSDigits,fBufferSize,file);  
  }

  cH = strstr(opt,"R");
  //Create a branch for Reconstruction
  if( cH ){
    char branchname[20];

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
    
  // 1.

  //  gAlice->MakeTree("R") ; 
  
  MakeBranch("R") ;
  
  // 3.

  fReconstructioner->Make(fDigits, fEmcRecPoints, fPpsdRecPoints, fTrackSegments, fRecParticles);

  printf("Reconstruction: %d %d %d %d\n",
	 fEmcRecPoints->GetEntries(),fPpsdRecPoints->GetEntries(),
	 fTrackSegments->GetEntries(),fRecParticles->GetEntries());

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
  // Reset hit tree for EMC and CPV
  // Yuri Kharlov, 28 September 2000

  AliDetector::ResetHits();
  //  for (Int_t i=0; i<fGeom->GetNModules(); i++) ((AliPHOSCPVModule*)(*fEMCModules)[i]) -> Clear();
  //  if ( strcmp(fGeom->GetName(),"IHEP") == 0 || strcmp(fGeom->GetName(),"MIXT") == 0 ) {
  //    for (Int_t i=0; i<fGeom->GetNCPVModules(); i++) ((AliPHOSCPVModule*)(*fCPVModules)[i]) -> Clear();
  //  }
  
}  
//____________________________________________________________________________
void AliPHOSv1::ResetReconstruction() 
{ 
  // Deleting reconstructed objects

  if ( fEmcRecPoints  )  fEmcRecPoints ->Delete();
  if ( fPpsdRecPoints )  fPpsdRecPoints->Delete();
  if ( fTrackSegments )  fTrackSegments->Delete();
  if ( fRecParticles  )  fRecParticles ->Delete();
  
}

//____________________________________________________________________________
//void AliPHOSv1::SDigits2Digits() {
//  // Adds the noise to the summable digits and keeps digits above a threshold 
//  // To make a digit
//}

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

  // Set branch address for the Hits Tree for hits in EMC modules
  // Yuri Kharlov, 23 November 2000.

  //  for( Int_t i=0; i<fGeom->GetNModules(); i++ ) GetEMCModule(i).SetTreeAddress("EMC",i+1);

  // Set branch address for the Hits Tree for hits in CPV modules for IHEP geometry
  // Yuri Kharlov, 28 September 2000.

  //  if ( strcmp(fGeom->GetName(),"IHEP") == 0 || strcmp(fGeom->GetName(),"MIXT") == 0 ) {
  //    for( Int_t i=0; i<fGeom->GetNCPVModules(); i++ ) GetCPVModule(i).SetTreeAddress("CPV",i+1);
  //  }

}

//____________________________________________________________________________

void AliPHOSv1::StepManager(void)
{
  // Accumulates hits as long as the track stays in a single crystal or PPSD gas Cell

//    if (gMC->IsTrackEntering())
//      cout << "Track enters the volume " << gMC->CurrentVolName() << endl;
//    if (gMC->IsTrackExiting())
//      cout << "Track leaves the volume " << gMC->CurrentVolName() << endl;

  Int_t          relid[4] ;      // (box, layer, row, column) indices
  Int_t          absid    ;      // absolute cell ID number
  Float_t        xyze[4]={0,0,0,0}  ; // position wrt MRS and energy deposited
  TLorentzVector pos      ;      // Lorentz vector of the track current position
  TLorentzVector pmom     ;      //momentum of the particle initiated hit
  Float_t        xyd[3]   ;      //local posiiton of the entering
  Bool_t         entered = kFALSE    ;  
  Int_t          copy     ;

  Int_t tracknumber =  gAlice->CurrentTrack() ; 
  Int_t primary     =  gAlice->GetPrimary( gAlice->CurrentTrack() ); 
  TString name      =  fGeom->GetName() ; 
  Int_t trackpid    =  0  ; 

  if( gMC->IsTrackEntering() ){ // create hit with position and momentum of new particle, 
                                // but may be without energy deposition

    // Current position of the hit in the local ref. system
      gMC -> TrackPosition(pos);
      Float_t xyzm[3], xyzd[3] ;
      Int_t i;
      for (i=0; i<3; i++) xyzm[i] = pos[i];
      gMC -> Gmtod (xyzm, xyzd, 1);    // transform coordinate from master to daughter system
      xyd[0]  = xyzd[0];
      xyd[1]  =-xyzd[1];
      xyd[2]  =-xyzd[2];

      
      // Current momentum of the hit's track in the local ref. system
      gMC -> TrackMomentum(pmom);
      Float_t pm[3], pd[3];
      for (i=0; i<3; i++) pm[i]   = pmom[i];
      gMC -> Gmtod (pm, pd, 2);        // transform 3-momentum from master to daughter system
      pmom[0] = pd[0];
      pmom[1] =-pd[1];
      pmom[2] =-pd[2];

      trackpid = gMC->TrackPid();
      entered = kTRUE ;      // Mark to create hit even withou energy deposition

  }


  if ( name == "GPS2" || name == "MIXT" ) {            // ======> CPV is a GPS' PPSD

    if( gMC->CurrentVolID(copy) == gMC->VolId("GCEL") ) // We are inside a gas cell 
    {
      gMC->TrackPosition(pos) ;
      xyze[0] = pos[0] ;
      xyze[1] = pos[1] ;
      xyze[2] = pos[2] ;
      xyze[3] = gMC->Edep() ; 

      if ( (xyze[3] != 0) || entered ) { // there is deposited energy or new particle entering  PPSD
       	gMC->CurrentVolOffID(5, relid[0]) ;  // get the PHOS Module number
	if ( name == "MIXT" && strcmp(gMC->CurrentVolOffName(5),"PHO1") == 0 ){
	  relid[0] += fGeom->GetNModules() - fGeom->GetNPPSDModules();
	}
       	gMC->CurrentVolOffID(3, relid[1]) ;  // get the Micromegas Module number 
      // 1-> fGeom->GetNumberOfModulesPhi() * fGeom->GetNumberOfModulesZ() upper
      //   > fGeom->GetNumberOfModulesPhi() * fGeom->GetNumberOfModulesZ() lower
       	gMC->CurrentVolOffID(1, relid[2]) ;  // get the row number of the cell
        gMC->CurrentVolID(relid[3]) ;        // get the column number 

	// get the absolute Id number

       	fGeom->RelToAbsNumbering(relid, absid) ; 

	// add current hit to the hit list      
	  AddHit(fIshunt, primary, tracknumber, absid, xyze, trackpid, pmom, xyd);


      } // there is deposited energy 
    } // We are inside the gas of the CPV  
  } // GPS2 configuration

  if ( name == "IHEP" || name == "MIXT" ) {       // ======> CPV is a IHEP's one

    // Yuri Kharlov, 28 September 2000

    if( gMC->CurrentVolID(copy) == gMC->VolId("CPVQ") &&
	entered &&
	gMC->TrackCharge() != 0) {      
      
      // Digitize the current CPV hit:

      // 1. find pad response and
      
      Int_t moduleNumber;
      gMC->CurrentVolOffID(3,moduleNumber);
      moduleNumber--;


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
	AddHit(fIshunt, primary, tracknumber, absid, xyze, trackpid, pmom, xyd);

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

  
    if ( (xyze[3] != 0) || entered ) {  // Track is inside the crystal and deposits some energy or just entered 

      gMC->CurrentVolOffID(10, relid[0]) ; // get the PHOS module number ;

      if ( name == "MIXT" && strcmp(gMC->CurrentVolOffName(10),"PHO1") == 0 )
	relid[0] += fGeom->GetNModules() - fGeom->GetNPPSDModules();      

      relid[1] = 0   ;                    // means PBW04
      gMC->CurrentVolOffID(4, relid[2]) ; // get the row number inside the module
      gMC->CurrentVolOffID(3, relid[3]) ; // get the cell number inside the module
      
      // get the absolute Id number
      fGeom->RelToAbsNumbering(relid, absid) ; 

      // add current hit to the hit list
	AddHit(fIshunt, primary,tracknumber, absid, xyze, trackpid,pmom, xyd);


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

//    cout << "CPVDigitize: YVK : "<<hitX<<" "<<hitZ<<" | "<<pX<<" "<<pZ<<" "<<pNorm<<endl;

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

  Int_t nCellZ  = fGeom->GetNumberOfCPVPadsZ();
  Int_t nCellX  = fGeom->GetNumberOfCPVPadsPhi();
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

