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
//  flow (e.g. digits, made from them RecPoints, subsequently made TrackSegments, 
//  subsequently made RecParticles) are distinguished by the title of created branches. One can 
//  use this title as a comment, see use case below. 
//  Thanks to getters, one can set 
//  parameters to reconstruction briks. The full set of parameters is saved in the 
//  corresponding branch: e.g. parameters of clusterizer are stored in branch 
//  TreeR::AliPHOSClusterizer with the same title as the branch containing the RecPoints. //  TTree does not support overwriting, therefore one can not produce several 
//  branches with the same names and titles - use different titles.
//
//  Use case: 
//
//  root [0] AliPHOSReconstructor * r = new AliPHOSReconstructor("galice.root")
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
//  root [4] r->StartFrom("AliPHOSClusterizer","Local max cut 0.02") 
//             // means that will use already generated Digits and produce only RecPoints, 
//             // TS and RecParticles 
//
//             // And finally one can call ExecuteTask() with the following options
//  root [5] r->ExecuteTask("debug all timing")
//            // deb     - prints the numbers of RecPoints, TrackSegments, RecParticles 
//            // deb all - prints in addition list of RecPoints, TrackSegments, RecParticles  
//            // timing  - prints benchmarking results
///////////////////////////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---
#include "AliESD.h"
#include "AliESDCaloTrack.h"
#include "AliPHOSReconstructor.h"
#include "AliPHOSClusterizerv1.h"
#include "AliPHOSTrackSegmentMakerv1.h"
#include "AliPHOSPIDv1.h"
#include "AliPHOSGetter.h"


ClassImp(AliPHOSReconstructor)

//____________________________________________________________________________
  AliPHOSReconstructor::AliPHOSReconstructor():TTask("AliPHOSReconstructor","")
{
  // ctor
  fClusterizer = 0 ;
  fTSMaker     = 0 ;
  fPID         = 0 ; 
  fFirstEvent  = 0 ; 
  fLastEvent   = -1 ; 
  fIsInitialized = kFALSE ;

} 

//____________________________________________________________________________
AliPHOSReconstructor::AliPHOSReconstructor(const char* evFoldName,const char * branchName,const TString taskName):
TTask("AliPHOSReconstructor",evFoldName)
{
  // Create a PHOS reconstructioner for the tasks defined by taskName
  // "C" - clusterization
  // "T" - track segment making
  // "P" - PID

  AliPHOSGetter::Instance(evFoldName) ; 

  if (taskName.Contains("C")) {
    fRecPointBranch=branchName ; 
    fClusterizer = new AliPHOSClusterizerv1(evFoldName, GetTitle());
    Add(fClusterizer);
  }
  
  if (taskName.Contains("T")) {
    fTSBranch=branchName ; 
    fTSMaker     = new AliPHOSTrackSegmentMakerv1(evFoldName, GetTitle());
    Add(fTSMaker) ;
  }
  
  if (taskName.Contains("P")) {
    fRecPartBranch=branchName ; 
    fPID         = new AliPHOSPIDv1(evFoldName, GetTitle());
    Add(fPID);
  }
  
  fIsInitialized = kTRUE ;
} 
//____________________________________________________________________________
void AliPHOSReconstructor::Exec(Option_t *opt)
{
  //check, if the names of branches, which should be made coincide with already
  //existing
  if (!opt) 
    return ; 
  if(!fIsInitialized)
    Init() ;
}
//____________________________________________________________________________
void AliPHOSReconstructor:: Clusters2Tracks(Int_t ievent, AliESD *event)
{
  // Convert PHOS reconstructed particles into ESD object for event# ievent.
  // ESD object is returned as an argument event

  if(!fIsInitialized) Init() ;

  fClusterizer->SetEventRange(ievent,ievent);
  fClusterizer->ExecuteTask();

  fTSMaker    ->SetEventRange(ievent,ievent);
  fTSMaker    ->ExecuteTask();
  
  fPID        ->SetEventRange(ievent,ievent);
  fPID        ->ExecuteTask();

  AliPHOSGetter *gime = AliPHOSGetter::Instance();
  TClonesArray *recParticles = gime->RecParticles();
  Int_t nOfRecParticles = recParticles->GetEntries();
  for (Int_t recpart=0; recpart<nOfRecParticles; recpart++) {
    AliESDCaloTrack *ct = new AliESDCaloTrack((AliPHOSRecParticle*)recParticles->At(recpart));
    event->AddCaloTrack(ct);
  }
  
}
//____________________________________________________________________________
 void AliPHOSReconstructor::Init()
{
  // initiliaze Reconstructioner if necessary: we can not do this in default constructor

  if(!fIsInitialized){
    
    fRecPointBranch="Default" ; 
    fClusterizer = new AliPHOSClusterizerv1(GetTitle(),fRecPointBranch.Data());
    Add(fClusterizer) ;

    fTSBranch="Default" ; 
    fTSMaker     = new AliPHOSTrackSegmentMakerv1(GetTitle(),fTSBranch.Data());
    Add(fTSMaker) ;


    fRecPartBranch="Default"; 
    fPID         = new AliPHOSPIDv1(GetTitle(),fRecPartBranch.Data()) ;
    Add(fPID) ;
    
    fIsInitialized = kTRUE ;
    
  }
} 
//____________________________________________________________________________
AliPHOSReconstructor::~AliPHOSReconstructor()
{
  // Delete data members if any
} 

void AliPHOSReconstructor::Print()const {
  // Print reconstructioner data  

  TString message ; 
  message  = "-----------------AliPHOSReconstructor---------------\n" ;
  message += " Reconstruction of the header file %s\n" ;
  message += " with the following modules:\n" ;

  if(fClusterizer->IsActive()){
    message += "   (+)   %s to branch %s\n" ;
  }

  if(fTSMaker->IsActive()){
    message += "   (+)   %s to branch %s\n" ; 
  }

  if(fPID->IsActive()){
    message += "   (+)   %s to branch %s\n" ;  
  }
  Info("Print", message.Data(), 
       GetTitle(), 
       fClusterizer->GetName(), fRecPointBranch.Data(), 
       fTSMaker->GetName(), fTSBranch.Data() , 
       fPID->GetName(), fRecPartBranch.Data() ) ; 
}

//____________________________________________________________________________
void AliPHOSReconstructor::SetEventRange(Int_t first, Int_t last)
{
  // Set the event range to process
  fFirstEvent=first; 
  fLastEvent=last; 
  fClusterizer->SetEventRange(fFirstEvent, fLastEvent) ; 
  fTSMaker->SetEventRange(fFirstEvent, fLastEvent) ;
  fPID->SetEventRange(fFirstEvent, fLastEvent) ;
}
