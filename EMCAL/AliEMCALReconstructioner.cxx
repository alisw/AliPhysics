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
//*-- Author: Gines Martinez & Yves Schutz (SUBATECH) 
//*-- Compleetely redesigned by Dmitri Peressounko (SUBATECH & RRC KI) March 2001
/////////////////////////////////////////////////////////////////////////////////////
//  Wrapping class for reconstruction. Allows to produce reconstruction from 
//  different steps: from previously produced hits,sdigits, etc. Each new reconstruction
//  flow (e.g. digits, made from them RecPoints,
//  subsequently made RecParticles) are distinguished by the title of created branches. One can 
//  use this title as a comment, see use case below. 
//  Thanks to getters, one can set 
//  parameters to reconstruction briks. The full set of parameters is saved in the 
//  corresponding branch: e.g. parameters of clusterizer are stored in branch 
//  TreeR::AliEMCALClusterizer with the same title as the branch containing the RecPoints. //  TTree does not support overwriting, therefore one can not produce several 
//  branches with the same names and titles - use different titles.
//
//  Use case: 
//
//  root [0] AliEMCALReconstructioner * r = new AliEMCALReconstructioner("galice.root")
//              //  Set the header file
//  root [1] r->ExecuteTask() 
//              //  Make full chain of reconstruction
//
//              // One can specify the title for each branch 
//  root [2] r->SetBranchFileName("RecPoints","RecPoints1") ;
//      
//             // One can change parameters of reconstruction algorithms
//  root [3] r->GetClusterizer()->SetEmcLocalMaxCut(0.02)
//
//             // One can specify the starting point of the reconstruction and title of all 
//             // branches produced in this pass
//  root [4] r->StartFrom("AliEMCALClusterizer","Local max cut 0.02") 
//             // means that will use already generated Digits and produce only RecPoints, 
//             // and RecParticles 
//
//             // And finally one can call ExecuteTask() with the following options
//  root [5] r->ExecuteTask("debug all timing")
//             // deb     - prints the numbers of RecPoints and RecParticles
//             // deb all - prints in addition list of RecPoints and RecParticles
//             // timing  - prints benchmarking results
///////////////////////////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---
#include "AliESD.h"
#include "AliESDCaloTrack.h"
#include "AliEMCALReconstructioner.h"
#include "AliEMCALClusterizerv1.h"
#include "AliEMCALPIDv1.h"
#include "AliEMCALGetter.h"

ClassImp(AliEMCALReconstructioner)

//____________________________________________________________________________
  AliEMCALReconstructioner::AliEMCALReconstructioner():TTask("AliEMCALReconstructioner","")
{
  // ctor
  fClusterizer = 0 ;
  fPID         = 0 ; 

  fIsInitialized = kFALSE ;

} 

//____________________________________________________________________________
AliEMCALReconstructioner::AliEMCALReconstructioner(const char* evFoldName,const char * branchName):
TTask("AliEMCALReconstructioner",evFoldName)
{
  // ctor
  
  AliEMCALGetter::Instance(evFoldName) ; 

  fRecPointBranch=branchName ; 
  fClusterizer = new AliEMCALClusterizerv1(evFoldName, GetTitle());
  Add(fClusterizer);
  
  fRecPartBranch=branchName ; 
  fPID         = new AliEMCALPIDv1(evFoldName, GetTitle());
  Add(fPID);
  
  fIsInitialized = kTRUE ;
} 
//____________________________________________________________________________
void AliEMCALReconstructioner::Exec(Option_t *opt)
{
  //check, if the names of branches, which should be made coincide with already
  //existing
  if (!opt) 
    return ; 
  if(!fIsInitialized)
    Init() ;
}
//____________________________________________________________________________
void AliEMCALReconstructioner:: Clusters2Tracks(Int_t ievent, AliESD *event)
{
  // Convert EMCAL reconstructed particles into ESD object for event# ievent.
  // ESD object is returned as an argument event

  if(!fIsInitialized) Init() ;

  fClusterizer->SetEventRange(ievent,ievent);
  fClusterizer->ExecuteTask();
  
  fPID        ->SetEventRange(ievent,ievent);
  fPID        ->ExecuteTask();

  AliEMCALGetter *gime = AliEMCALGetter::Instance();
  TClonesArray *recParticles = gime->RecParticles();
  Int_t nOfRecParticles = recParticles->GetEntries();
  for (Int_t recpart=0; recpart<nOfRecParticles; recpart++) {
    AliESDCaloTrack *ct = new AliESDCaloTrack((AliEMCALRecParticle*)recParticles->At(recpart));
    event->AddCaloTrack(ct);
  }
  
}
//____________________________________________________________________________
 void AliEMCALReconstructioner::Init()
{
  // initiliaze Reconstructioner if necessary: we can not do this in default constructor

  if(!fIsInitialized){

    fRecPointBranch="Default" ; 
    fClusterizer = new AliEMCALClusterizerv1(GetTitle(),fRecPointBranch.Data());
    Add(fClusterizer) ;

    fRecPartBranch="Default"; 
    fPID         = new AliEMCALPIDv1(GetTitle(),fRecPartBranch.Data()) ;
    Add(fPID) ;
    
    fIsInitialized = kTRUE ;
    
  }
} 
//____________________________________________________________________________
AliEMCALReconstructioner::~AliEMCALReconstructioner()
{
  // Delete data members if any
} 

//____________________________________________________________________________
void AliEMCALReconstructioner::Print()const {
  // Print reconstructioner data  

  TString message ; 
  message  = "-----------------AliEMCALReconstructioner---------------\n" ;
  message += " Reconstruction of the header file %s\n" ;
  message += " with the following modules:\n" ;
  
  if(fClusterizer->IsActive()){
    message += "   (+)   %s to branch %s\n" ;
  }

  if(fPID->IsActive()){
    message += "   (+)   %s to branch %s\n" ;  
  }
  Info("Print", message.Data(), 
       GetTitle(), 
       fClusterizer->GetName(), fRecPointBranch.Data(), 
       fPID->GetName(), fRecPartBranch.Data() ) ; 
}

//____________________________________________________________________________
void AliEMCALReconstructioner::SetEventRange(Int_t first, Int_t last)
{
  // Set the event range to process
  fFirstEvent=first; 
  fLastEvent=last; 
  fClusterizer->SetEventRange(fFirstEvent, fLastEvent) ; 
  fPID->SetEventRange(fFirstEvent, fLastEvent) ;
}
