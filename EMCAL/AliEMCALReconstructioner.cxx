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
//  TreeR::AliEMCALClusterizer with the same title as the branch containing the RecPoints. 
//  TTree does not support overwriting, therefore one can not produce several 
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
//             // TS and RecParticles 
//
//             // And finally one can call ExecuteTask() with the following options
//  root [5] r->ExecuteTask("debug all timing")
//             // deb     - prints the numbers of produced SDigits, Digits etc.
//             // deb all - prints in addition list of made SDigits, digits etc.
//             // timing  - prints benchmarking results
///////////////////////////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

#include "TClonesArray.h"
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

// --- Standard library ---

// --- AliRoot header files ---
#include "AliRun.h"
#include "AliEMCALReconstructioner.h"
#include "AliEMCALClusterizerv1.h"
#include "AliEMCALDigitizer.h"
#include "AliEMCALSDigitizer.h"
//#include "AliEMCALTrackSegmentMakerv1.h"
//#include "AliEMCALPIDv1.h"
//#include "AliEMCALFastRecParticle.h"
//#include "AliEMCALCpvRecPoint.h"

ClassImp(AliEMCALReconstructioner)

//____________________________________________________________________________
  AliEMCALReconstructioner::AliEMCALReconstructioner():TTask("AliEMCALReconstructioner","")
{
  // ctor
  fToSplit = kFALSE ;
  fDigitizer   = 0 ;
  fClusterizer = 0 ;
  //fTSMaker     = 0 ;
  //fPID         = 0 ; 
  fSDigitizer  = 0 ;
  fHeaderFileName = "galice.root" ;

  fIsInitialized = kFALSE ;

} 

//____________________________________________________________________________
AliEMCALReconstructioner::AliEMCALReconstructioner(const char* headerFile,const char * branchName,Bool_t toSplit):
TTask("AliEMCALReconstructioner","")
{
  // ctor
  
  fHeaderFileName = headerFile ;
  fToSplit = toSplit ;
  fSDigitsBranch= branchName; 
  fSDigitizer  = new AliEMCALSDigitizer(fHeaderFileName.Data(),fSDigitsBranch.Data(),toSplit) ; 
  Add(fSDigitizer) ;

  fDigitsBranch=branchName ; 
  fDigitizer   = new AliEMCALDigitizer(fHeaderFileName.Data(),fDigitsBranch.Data(),toSplit) ; 
  Add(fDigitizer) ;


  fRecPointBranch=branchName ; 
  fClusterizer = new AliEMCALClusterizerv1(fHeaderFileName.Data(),fRecPointBranch.Data(),toSplit) ; 
  Add(fClusterizer) ;
  

//   fTSBranch=branchName ; 
//   fTSMaker     = new AliEMCALTrackSegmentMakerv1(fHeaderFileName.Data(),fTSBranch.Data(),toSplit) ;
//   Add(fTSMaker) ;
  
  
//   fRecPartBranch=branchName ; 
//   fPID         = new AliEMCALPIDv1(fHeaderFileName.Data(),fRecPartBranch.Data(),toSplit) ;
//   Add(fPID) ;
  
  fIsInitialized = kTRUE ;
  
} 
//____________________________________________________________________________
void AliEMCALReconstructioner::Exec(Option_t *option)
{
  //chesk, if the names of branches, which should be made conicide with already
  //existing
  if(!fIsInitialized)
    Init() ;
}

//____________________________________________________________________________
 void AliEMCALReconstructioner::Init()
{
  // initiliaze Reconstructioner if necessary: we can not do this in default constructor

  if(!fIsInitialized){
    // Initialisation

    fSDigitsBranch="Default" ; 
    fSDigitizer  = new AliEMCALSDigitizer(fHeaderFileName.Data(),fSDigitsBranch.Data(),fToSplit) ; 
    Add(fSDigitizer) ;

    fDigitsBranch="Default" ; 
    fDigitizer   = new AliEMCALDigitizer(fHeaderFileName.Data(),fDigitsBranch.Data(),fToSplit) ; 
    Add(fDigitizer) ;

    fRecPointBranch="Default" ; 
    fClusterizer = new AliEMCALClusterizerv1(fHeaderFileName.Data(),fRecPointBranch.Data(),fToSplit) ; 
    Add(fClusterizer) ;

//     fTSBranch="Default" ; 
//     fTSMaker     = new AliEMCALTrackSegmentMakerv1(fHeaderFileName.Data(),fTSBranch.Data(),fToSplit) ;
//     Add(fTSMaker) ;


//     fRecPartBranch="Default" ; 
//     fPID         = new AliEMCALPIDv1(fHeaderFileName.Data(),fRecPartBranch.Data(),fToSplit) ;
//     Add(fPID) ;
    
    fIsInitialized = kTRUE ;
  }
} 
//____________________________________________________________________________
AliEMCALReconstructioner::~AliEMCALReconstructioner()
{
  // Delete data members if any

} 

//____________________________________________________________________________

void AliEMCALReconstructioner::Print(Option_t * option)const {
  // Print reconstructioner data  

  printf("\n") ; 

  printf(" Reconstruction of the header file "); 
  printf(fHeaderFileName.Data());
  printf("\n with the following modules:\n");

  if(fSDigitizer->IsActive()){
     printf("   (+)   "); 
     printf(fSDigitizer->GetName()); 
     printf(" to branch : "); 
     printf(fSDigitsBranch.Data()); 
  }
  if(fDigitizer->IsActive()){
    printf("\n   (+)   "); 
    printf(fDigitizer->GetName()); 
    printf(" to branch : "); 
    printf(fDigitsBranch.Data()); 
  }
  
  if(fClusterizer->IsActive()){
    printf("\n   (+)   "); 
    printf(fClusterizer->GetName()); 
    printf(" to branch : "); 
    printf(fRecPointBranch.Data()); 
  } 
}
