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
#include "AliEMCALReconstructor.h"

#include "AliESD.h"
#include "AliRunLoader.h"
#include "AliEMCALLoader.h"
#include "AliEMCALClusterizerv1.h"
#include "AliEMCALRecPoint.h"
#include "AliRawReaderFile.h"

ClassImp(AliEMCALReconstructor)

//____________________________________________________________________________
  AliEMCALReconstructor::AliEMCALReconstructor() : fDebug(kFALSE) 
{
  // ctor

} 


//____________________________________________________________________________
AliEMCALReconstructor::~AliEMCALReconstructor()
{
  // dtor
} 

//____________________________________________________________________________
void AliEMCALReconstructor::Reconstruct(AliRunLoader* runLoader) const 
{
  // method called by AliReconstruction; 
  // Only the clusterization is performed,; the rest of the reconstruction is done in FillESD because the track
  // segment maker needs access to the AliESD object to retrieve the tracks reconstructed by 
  // the global tracking.
 
  TString headerFile(runLoader->GetFileName()) ; 
  TString branchName(runLoader->GetEventFolder()->GetName() ) ;  
  
  AliEMCALClusterizerv1 clu(headerFile, branchName);
  clu.SetEventRange(0, -1) ; // do all the events
  if ( Debug() ) 
    clu.ExecuteTask("deb all") ; 
  else 
    clu.ExecuteTask("") ;  

}

//____________________________________________________________________________
void AliEMCALReconstructor::Reconstruct(AliRunLoader* runLoader, AliRawReader* rawreader) const 
{
  // method called by AliReconstruction; 
  // Only the clusterization is performed,; the rest of the reconstruction is done in FillESD because the track
  // segment maker needs access to the AliESD object to retrieve the tracks reconstructed by 
  // the global tracking.
  // Here we reconstruct from Raw Data
  
  rawreader->Reset() ; 
  TString headerFile(runLoader->GetFileName()) ; 
  TString branchName(runLoader->GetEventFolder()->GetName()) ;  
  
  AliEMCALClusterizerv1 clu(headerFile, branchName);
  clu.SetEventRange(0, -1) ; // do all the events
  if ( Debug() ) 
    clu.ExecuteTask("deb all") ; 
  else 
    clu.ExecuteTask("") ;  

}

//____________________________________________________________________________
void AliEMCALReconstructor::FillESD(AliRunLoader* runLoader, AliESD* esd) const
{
  // Called by AliReconstruct after Reconstruct() and global tracking and vertxing 
  
  Int_t eventNumber = runLoader->GetEventNumber() ;

  TString headerFile(runLoader->GetFileName()) ; 
  TString branchName(runLoader->GetEventFolder()->GetName()) ;  

  // PID is not implemented. Skipping for now
  //AliEMCALPIDv1 pid(headerFile, branchName);

  // do current event; the loop over events is done by AliReconstruction::Run()
  /*
  pid.SetEventRange(eventNumber, eventNumber) ; 
  if ( Debug() ) 
    pid.ExecuteTask("deb all") ;
  else 
    pid.ExecuteTask("") ;
  */

  // Creates AliESDtrack from AliEMCALRecPoints 
  AliRunLoader *rl = AliRunLoader::GetRunLoader();
  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(rl->GetDetectorLoader("EMCAL"));
  rl->LoadRecPoints();
  rl->GetEvent(eventNumber);
  TObjArray *clusters = emcalLoader->RecPoints();
  Int_t nClusters = clusters->GetEntries();
  esd->SetNumberOfEMCALParticles(nClusters) ; 
  esd->SetFirstEMCALParticle(esd->GetNumberOfTracks()) ; 
  for (Int_t iClust = 0 ; iClust < nClusters ; iClust++) {
    const AliEMCALRecPoint * clust = emcalLoader->RecPoint(iClust);
    if (Debug()) 
      clust->Print();
    AliESDtrack * et = new AliESDtrack() ; 
    // fills the ESDtrack
    Double_t xyz[3];
    TVector3 gpos;
    clust->GetGlobalPosition(gpos);
    for (Int_t ixyz=0; ixyz<3; ixyz++) 
      xyz[ixyz] = gpos[ixyz];
    et->SetEMCALposition(xyz) ; 
    et->SetEMCALsignal  (clust->GetEnergy()) ; 
    // add the track to the esd object
    esd->AddTrack(et);
    delete et;
  }
}
