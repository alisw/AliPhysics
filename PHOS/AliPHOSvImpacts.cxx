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
// Implementation version vImpacts of PHOS Manager class.
// This class inherits from v1 and adds impacts storing.
// Impacts stands for exact values of track coming to the detectors
// EMC, CPV or PPSD.
// Impacts are written to the same tree as hits are
// but in separate branches.
//---
//*-- Author: Yuri Kharlov (IHEP, Protvino/SUBATECH, Nantes)


// --- ROOT system ---

#include "TTree.h"

// --- Standard library ---

// --- AliRoot header files ---

#include "AliRun.h"
#include "AliMC.h"
#include "AliPHOSvImpacts.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSImpact.h"

ClassImp(AliPHOSvImpacts)

//____________________________________________________________________________
AliPHOSvImpacts::AliPHOSvImpacts()
{
  // ctor
}

//____________________________________________________________________________
AliPHOSvImpacts::AliPHOSvImpacts(const char *name, const char *title):
AliPHOSv1(name,title) 
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
  //     This part inherits from AliPHOSv1
  //
  // We store impacts :
  //  - fEMCImpacts, fCPVImpacts, fPPSDImpacts which are
  //    TList of EMC, CPV and PPSD modules respectively, each
  //    modules contains TClonesArray of AliPHOSImpacts
  
  fEMCImpacts  = new TList();
  fCPVImpacts  = new TList();
  fPPSDImpacts = new TList();

  Int_t nPHOSModules = fGeom->GetNModules();
  Int_t nCPVModules  = fGeom->GetNCPVModules();
  Int_t nPPSDModules = fGeom->GetNPPSDModules();

  Int_t iPHOSModule;
  TClonesArray * impacts;
  for (iPHOSModule=0; iPHOSModule<nPHOSModules; iPHOSModule++) {
    fEMCImpacts->Add(new TClonesArray("AliPHOSImpact",200)) ;
    fNEMCImpacts[iPHOSModule] = 0;
    impacts = (TClonesArray *)fEMCImpacts->At(iPHOSModule);
  }
  for (iPHOSModule=0; iPHOSModule<nCPVModules; iPHOSModule++) {
    fCPVImpacts->Add(new TClonesArray("AliPHOSImpact",200)) ;
    fNCPVImpacts[iPHOSModule] = 0;
    impacts = (TClonesArray *)fCPVImpacts->At(iPHOSModule);
  }
  for (iPHOSModule=0; iPHOSModule<nPPSDModules; iPHOSModule++) {
    fPPSDImpacts->Add(new TClonesArray("AliPHOSImpact",200)) ;
    fNPPSDImpacts[iPHOSModule] = 0;
    impacts = (TClonesArray *)fPPSDImpacts->At(iPHOSModule);
  }

}

//____________________________________________________________________________
AliPHOSvImpacts::~AliPHOSvImpacts()
{
  // dtor

  // Delete hits
  if ( fHits ) {
    fHits->Delete() ; 
    delete fHits ;
    fHits = 0 ; 
  }

  // Delete impacts in EMC, CPV and PPSD
  if ( fEMCImpacts ) {
    fEMCImpacts->Delete() ; 
    delete fEMCImpacts ;
    fEMCImpacts = 0 ; 
  }
  if ( fCPVImpacts ) {
    fCPVImpacts->Delete() ; 
    delete fCPVImpacts ;
    fCPVImpacts = 0 ; 
  }
  if ( fPPSDImpacts ) {
    fPPSDImpacts->Delete() ; 
    delete fPPSDImpacts ;
    fPPSDImpacts = 0 ; 
  }
}

//____________________________________________________________________________
void AliPHOSvImpacts::AddImpact( char* det, Int_t shunt, Int_t primary, Int_t track, Int_t module,
			   Int_t pid, TLorentzVector p, Float_t *xyz)
{
  // Add an impact to the impact list.

  TClonesArray * impacts = 0;
  Int_t         nImpacts = 0;

  if (strcmp(det,"EMC ")==0) {
    impacts = (TClonesArray *)fEMCImpacts->At(module);
    nImpacts= fNEMCImpacts[module];
    fNEMCImpacts[module]++ ;
  }
  else if (strcmp(det,"CPV ")==0) {
    impacts = (TClonesArray *)fCPVImpacts->At(module);
    nImpacts= fNCPVImpacts[module];
    fNCPVImpacts[module]++ ;
  }
  else if (strcmp(det,"PPSD")==0) {
    impacts = (TClonesArray *)fPPSDImpacts->At(module);
    nImpacts= fNPPSDImpacts[module];
    fNPPSDImpacts[module]++ ;
  }

  new((*impacts)[nImpacts]) AliPHOSImpact(shunt,primary,track,pid,p,xyz) ;

  if (fDebug==1) {
    printf("Module %d %s: ",module,det);
    ((AliPHOSImpact*)(impacts->At(nImpacts)))->Print();
  }
}

//____________________________________________________________________________
void AliPHOSvImpacts::MakeBranch(Option_t *opt, const char *file)
{  
  // Create new branch in the current Hits Root Tree containing
  // a list of PHOS impacts (exact values of track coming to detector)

  AliDetector::MakeBranch(opt,file);
  
  Int_t bufferSize = 32000 ;
  Int_t splitlevel = 0 ;
  gAlice->TreeH()->Branch("PHOSEmcImpacts" , "TList", &fEMCImpacts , bufferSize, splitlevel);
  gAlice->TreeH()->Branch("PHOSCpvImpacts" , "TList", &fCPVImpacts , bufferSize, splitlevel);
  gAlice->TreeH()->Branch("PHOSPpsdImpacts", "TList", &fPPSDImpacts, bufferSize, splitlevel);
  
}

//____________________________________________________________________________
void AliPHOSvImpacts::ResetHits() 
{              
  // Reset impact branches for EMC, CPV and PPSD

  AliDetector::ResetHits();

  Int_t i;
  for (i=0; i<fGeom->GetNModules(); i++) {
    ((TClonesArray*)fEMCImpacts->At(i)) -> Clear();
    fNEMCImpacts[i] = 0 ;
  }

  if ( strcmp(fGeom->GetName(),"IHEP") == 0 || strcmp(fGeom->GetName(),"MIXT") == 0 ) {
    for (i=0; i<fGeom->GetNCPVModules(); i++) {
      ((TClonesArray*)fCPVImpacts->At(i)) -> Clear();
      fNCPVImpacts[i] = 0 ;
    }
  }

  if ( strcmp(fGeom->GetName(),"GPS2") == 0 || strcmp(fGeom->GetName(),"MIXT") == 0 ) {
    for (i=0; i<fGeom->GetNPPSDModules(); i++) {
      ((TClonesArray*)fPPSDImpacts->At(i)) -> Clear();
      fNPPSDImpacts[i] = 0 ;
    }
  }
  
}

//_____________________________________________________________________________
void AliPHOSvImpacts::StepManager(void)
{
  // Find impacts (tracks which enter the EMC, CPV or PPSD)
  // and add them to to the impact lists

  AliPHOSv1::StepManager();

  Float_t xyzm[3], xyzd[3], pm[3], pd[3];
  TLorentzVector pmom     ;           // Lorentz momentum of the particle initiated hit
  TLorentzVector pos      ;           // Lorentz vector of the track current position
  Int_t          copy     ;

  Int_t tracknumber =  gAlice->CurrentTrack() ; 
  Int_t primary     =  gAlice->GetPrimary( gAlice->CurrentTrack() ); 
  TString name      =  fGeom->GetName() ; 

  // Add impact to EMC

  if( gMC->CurrentVolID(copy) == gMC->VolId("PXTL") &&
      gMC->IsTrackEntering() ) {
    gMC->TrackMomentum(pmom);
    gMC->TrackPosition(pos) ;

    Int_t i;
    for (i=0; i<3; i++) xyzm[i] = pos[i];

    for (i=0; i<3; i++) {
      xyzm[i] = pos[i] ;
      pm[i]   = pmom[i];
    }
    gMC -> Gmtod (xyzm, xyzd, 1);    // transform coordinate from master to daughter system
    gMC -> Gmtod (pm,   pd,   2);    // transform 3-momentum from master to daughter system

    // Select tracks coming to the crystal from up or down sides
    if (pd[1]<0 && xyzd[1] >  fGeom->GetCrystalSize(1)/2-0.001 ||
	pd[1]>0 && xyzd[1] < -fGeom->GetCrystalSize(1)/2+0.001) {
      Int_t pid = gMC->TrackPid();
      Int_t module;
      gMC->CurrentVolOffID(10,module);
      if ( name == "MIXT" && strcmp(gMC->CurrentVolOffName(10),"PHO1") == 0 )
	module += fGeom->GetNModules() - fGeom->GetNPPSDModules();
      module--;
      AddImpact("EMC ",fIshunt, primary,tracknumber, module, pid, pmom, xyzm);
    }
  }

  // Add impact to CPV

  if( (name == "IHEP" || name == "MIXT") &&
      gMC->CurrentVolID(copy) == gMC->VolId("PCPQ") &&
      gMC->IsTrackEntering() ) {
    gMC->TrackMomentum(pmom);
    gMC->TrackPosition(pos) ;

    Int_t i;
    for (i=0; i<3; i++) xyzm[i] = pos[i];

    for (i=0; i<3; i++) {
      xyzm[i] = pos[i] ;
      pm[i]   = pmom[i];
    }
    Int_t pid = gMC->TrackPid();
    Int_t module;
    gMC->CurrentVolOffID(3,module);
    module--;
    AddImpact("CPV ",fIshunt, primary,tracknumber, module, pid, pmom, xyzm);
  }

  // Add impact to PPSD

  if( (name == "GPS2" || name == "MIXT") &&
      gMC->CurrentVolID(copy) == gMC->VolId("PPCE") &&
      gMC->IsTrackEntering() ) {
    gMC->TrackMomentum(pmom);
    gMC->TrackPosition(pos) ;

    Int_t i;
    for (i=0; i<3; i++) xyzm[i] = pos[i];

    for (i=0; i<3; i++) {
      xyzm[i] = pos[i] ;
      pm[i]   = pmom[i];
    }
    gMC -> Gmtod (xyzm, xyzd, 1);    // transform coordinate from master to daughter system
    gMC -> Gmtod (pm,   pd,   2);    // transform 3-momentum from master to daughter system

    // Select tracks coming to the crystal from up or down sides
    if (pd[1]<0 && xyzd[1] >  (fGeom->GetConversionGap() +  fGeom->GetAvalancheGap())/2-0.001 ||
	pd[1]>0 && xyzd[1] < -(fGeom->GetConversionGap() +  fGeom->GetAvalancheGap())/2+0.001) {
      Int_t pid = gMC->TrackPid();
      Int_t module;
      gMC->CurrentVolOffID(5,module);
      module--;
      AddImpact("PPSD",fIshunt, primary,tracknumber, module, pid, pmom, xyzm);
    }
  }
}
