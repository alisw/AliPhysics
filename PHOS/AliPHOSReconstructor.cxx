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
#include "AliPHOSTracker.h"
#include "AliRawReader.h"

 
ClassImp(AliPHOSReconstructor)

Bool_t AliPHOSReconstructor::fgDebug = kFALSE ; 

//____________________________________________________________________________
  AliPHOSReconstructor::AliPHOSReconstructor() 
{
  // ctor

} 

//____________________________________________________________________________
  AliPHOSReconstructor::~AliPHOSReconstructor()
{
  // dtor

} 

//____________________________________________________________________________
void AliPHOSReconstructor::Reconstruct(AliRunLoader* runLoader) const
{
  // method called by AliReconstruction; 
  // Only the clusterization is performed,; the rest of the reconstruction is done in FillESD because the track
  // segment maker needs access to the AliESD object to retrieve the tracks reconstructed by 
  // the global tracking.
 
  TString headerFile(runLoader->GetFileName()) ; 
  TString branchName(runLoader->GetEventFolder()->GetName()) ;  
  
  AliPHOSClusterizerv1 clu(headerFile, branchName);
  clu.SetEventRange(0, -1) ; // do all the events
  if ( Debug() ) 
    clu.ExecuteTask("deb all") ; 
  else 
    clu.ExecuteTask("") ;  

}

//____________________________________________________________________________
void AliPHOSReconstructor::Reconstruct(AliRunLoader* runLoader, AliRawReader* rawreader) const
{
  // method called by AliReconstruction; 
  // Only the clusterization is performed,; the rest of the reconstruction is done in FillESD because the track
  // segment maker needs access to the AliESD object to retrieve the tracks reconstructed by 
  // the global tracking.
  // Here we reconstruct from Raw Data

  rawreader->Reset() ; 
  TString headerFile(runLoader->GetFileName()) ; 
  TString branchName(runLoader->GetEventFolder()->GetName()) ;  
  
  AliPHOSClusterizerv1 clu(headerFile, branchName);
  clu.SetEventRange(0, -1) ; // do all the events
  clu.SetRawReader(rawreader);
  if ( Debug() ) 
    clu.ExecuteTask("deb all") ; 
  else 
    clu.ExecuteTask("") ;

}

//____________________________________________________________________________
void AliPHOSReconstructor::FillESD(AliRunLoader* runLoader, AliESD* esd) const
{
  // This function creates AliESDtracks from AliPHOSRecParticles
  //         and
  // writes them to the ESD

  Int_t eventNumber = runLoader->GetEventNumber() ;

  AliPHOSGetter::Instance()->Event(eventNumber, "P") ; 
  TClonesArray *recParticles = AliPHOSGetter::Instance()->RecParticles();
  Int_t nOfRecParticles = recParticles->GetEntries();
  esd->SetNumberOfPHOSClusters(nOfRecParticles) ; 
  esd->SetFirstPHOSCluster(esd->GetNumberOfTracks()) ; 

  for (Int_t recpart = 0 ; recpart < nOfRecParticles ; recpart++) {
    AliPHOSRecParticle * rp = dynamic_cast<AliPHOSRecParticle*>(recParticles->At(recpart));
    if (Debug()) 
      rp->Print();
    AliESDCaloCluster * ec = new AliESDCaloCluster() ; 
//     AliESDtrack * et = new AliESDtrack() ; 
    // fills the ESDCaloCluster
    Float_t xyz[3];
    for (Int_t ixyz=0; ixyz<3; ixyz++) 
      xyz[ixyz] = rp->GetPos()[ixyz];
    ec->SetGlobalPosition(xyz);
    ec->SetClusterEnergy(rp->Energy());
    ec->SetPid          (rp->GetPID()) ;
//     ec->SetPrimaryIndex (rp->GetPrimaryIndex()); // for dry events only
    // add the track to the esd object
    esd->AddCaloCluster(ec);
    delete ec;
  }
}

AliTracker* AliPHOSReconstructor::CreateTracker(AliRunLoader* runLoader) const
{
// creates the PHOS tracker
  if (!runLoader) return NULL; 
  return new AliPHOSTracker(runLoader);
}

