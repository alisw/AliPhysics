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
//*--
//*-- Yves Schutz (SUBATECH) 
// Reconstruction class. Redesigned from the old AliReconstructionner class and 
// derived from STEER/AliReconstructor. 
// 
// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---
#include "AliESD.h"
#include "AliPHOSReconstructor.h"
#include "AliPHOSClusterizerv1.h"
#include "AliPHOSTrackSegmentMakerv1.h"
#include "AliPHOSPIDv1.h"
#include "AliPHOSGetter.h"


ClassImp(AliPHOSReconstructor)

//____________________________________________________________________________
  AliPHOSReconstructor::AliPHOSReconstructor() : fDebug(kFALSE) 
{
  // ctor
} 


//____________________________________________________________________________
void AliPHOSReconstructor::Reconstruct(AliRunLoader* runLoader) const 
{
  // method called by AliReconstruction; 
  // Only the clusterization is performed,; the rest of the reconstruction is done in FillESD because the track
  // segment maker needs access to the AliESD object to retrieve the tracks reconstructed by 
  // the global tracking.
 
  TString headerFile(runLoader->GetFileName()) ; 
  TString branchName("Default") ;  
  
  AliPHOSClusterizerv1 clu(headerFile, branchName);
  clu.SetEventRange(0, -1) ; // do all the events
  if ( Debug() ) 
    clu.ExecuteTask("deb all") ; 
  else 
    clu.ExecuteTask("") ;  
}

//____________________________________________________________________________
void AliPHOSReconstructor::FillESD(AliRunLoader* runLoader, AliESD* esd) const
{
  // Called by AliReconstruct after Reconstruct() and global tracking and vertxing 
  //Creates the tracksegments and Recparticles
  
  TString headerFile(runLoader->GetFileName()) ; 
  TString branchName("Default") ;  

  AliPHOSTrackSegmentMakerv1 tsMaker(headerFile, branchName);
  AliPHOSPIDv1 pID(headerFile, branchName);

  AliPHOSGetter *gime = AliPHOSGetter::Instance( (runLoader->GetFileName()).Data() ) ;
  Int_t eventNumber = gime->EventNumber();
  // do current event; the loop over events is done by AliReconstruction::Run()
  tsMaker.SetEventRange(eventNumber, eventNumber) ; 
  pID.SetEventRange(eventNumber, eventNumber) ; 
  if ( Debug() ) {
   tsMaker.ExecuteTask("deb all") ;
   pID.ExecuteTask("deb all") ;
  }
  else {
    tsMaker.ExecuteTask("") ;
    pID.ExecuteTask("") ;
  }
  
  // Creates AliESDtrack from AliPHOSRecParticles 
  gime->Event(eventNumber, "P") ; 
  TClonesArray *recParticles = gime->RecParticles();
  Int_t nOfRecParticles = recParticles->GetEntries();
  for (Int_t recpart = 0 ; recpart < nOfRecParticles ; recpart++) {
    AliPHOSRecParticle * rp = dynamic_cast<AliPHOSRecParticle*>(recParticles->At(recpart));
    if (Debug()) 
      rp->Print();
    AliESDtrack * et = new AliESDtrack() ; 
    // fills the ESDtrack
    Double_t xyz[3];
    for (Int_t ixyz=0; ixyz<3; ixyz++) xyz[ixyz] = rp->GetPos()[ixyz];
    et->SetPHOSposition(xyz) ; 
    et->SetPHOSsignal  (rp->Energy()) ; 
    et->SetPHOSpid     (rp->GetPID()) ;
    // add the track to the esd object
    esd->AddTrack(et);
    delete et;
  }
}
