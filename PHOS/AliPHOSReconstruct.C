// YS Subatech Mai 2002
// YK Subatech 6 Aug 2002

// PHOS Reconstruction chain:
// Hits -> SDigits -> Digits -> RecPoints -> TrackSegments -> RecParticles

//Root
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TString.h"

//AliRoot
#include "STEER/AliRun.h"
#include "PHOS/AliPHOSSDigitizer.h"
#include "PHOS/AliPHOSDigitizer.h"
#include "PHOS/AliPHOSClusterizerv1.h"
#include "PHOS/AliPHOSTrackSegmentMakerv1.h"
#include "PHOS/AliPHOSPIDv1.h"
#include "EMCAL/AliEMCALSDigitizer.h"
#include "EMCAL/AliEMCALDigitizer.h"
#include "EMCAL/AliEMCALClusterizerv1.h"
#endif

void PHOSHits2SDigits( Bool_t split=kFALSE, TString fileName = "galice.root") {

  // usage : 
  // 1. write SDigits in the same file as Hits --------------- (OK)
  //root [0] .L Reconstruct.C++
  //root [1] SDigits2Digits()
  // 2. write SDigits in a separate file, one per detector, from Hits --------------- (OK)
  //root [0] .L Reconstruct.C++
  //root [1] SDigits2Digits(kTRUE) // SDigits saved in [DET}.SDigits.root (DET=PHOS, EMCAL)

  delete gAlice ; 
  gAlice = 0 ; 
  
  AliPHOSSDigitizer * sdp = new AliPHOSSDigitizer(fileName) ; 
  if (split) 
    sdp->SetSplitFile() ;
  sdp->ExecuteTask("deb") ; 

  delete sdp ;
}

//________________________________________________________________________
void PHOSSDigits2Digits( Bool_t split=kFALSE, TString fileName = "galice.root") {
  
 // usage : 
  // 1. write SDigits in the same file as SDigits --------------- (OK)
  //root [0] .L Reconstruct.C++
  //root [1] SDigits2Digits()
  // 2. write SDigits in a separate file, one per detector, from SDigits --------------- (OK)
  //root [0] .L Reconstruct.C++
  //root [1] SDigitsDigits(kTRUE) // Digits saved in [DET}.Digits.root (DET=PHOS, EMCAL)

  delete gAlice ; 
  gAlice = 0 ; 
  
  // PHOS
  AliPHOSDigitizer * dp = 0 ; 
 
  if (split) {
    dp = new AliPHOSDigitizer("PHOS.SDigits.root") ; 
    dp->SetSplitFile() ; } 
  else 
    dp = new AliPHOSDigitizer(fileName) ; 
  
  dp->ExecuteTask("deb") ; 
  
  delete dp ;
}

//________________________________________________________________________
void PHOSDigits2RecPoints( Bool_t split=kFALSE, TString fileName = "galice.root") {
  
 // usage : 
  // 1. write RecPoints in the same file as Digits --------------- OK 
  //root [0] .L Reconstruct.C++
  //root [1] Digits2RecPoints()
  // 2. write RecPoints in a separate file, one per detector, from Digits --------------- OK 
  //root [0] .L Reconstruct.C++
  //root [1] Digits2RecPoints(kTRUE) // RecPoints saved in [DET}.RecPoints.root (DET=PHOS, EMCAL)

  delete gAlice ; 
  gAlice = 0 ; 
 
  AliPHOSClusterizer * cp = 0 ; 
 
  if (split) {
    cp = new AliPHOSClusterizerv1("PHOS.Digits.root") ; 
    cp->SetSplitFile() ; } 
  else 
    cp = new AliPHOSClusterizerv1(fileName) ; 
  
  cp->ExecuteTask("deb") ; 
  
  delete cp ;
}

//________________________________________________________________________
void PHOSRecPoints2TrackSegments( Bool_t split=kFALSE, TString fileName = "galice.root") {
  
 // usage : 
  // 1. write TrackSegments in the same file as RecPoints --------------- (OK) 
  //root [0] .L Reconstruct.C++
  //root [1] RecPoints2TrackSegments()
  // 2. write TrackSegments in a separate file, one per detector, from RecPoints --------------- (Not needed) 
  //root [0] .L Reconstruct.C++
  //root [1] RecPoints2TrackSegments(kTRUE) // TrackSegments saved in [DET}.RecData.root (DET=PHOS, EMCAL)

  delete gAlice ; 
  gAlice = 0 ; 
  
  AliPHOSTrackSegmentMaker * tmp = 0 ; 
 
  if (split)
    tmp = new AliPHOSTrackSegmentMakerv1("PHOS.RecData.root") ; 
  else 
    tmp = new AliPHOSTrackSegmentMakerv1(fileName) ; 
  
  tmp->ExecuteTask("deb") ; 
  
  delete tmp ;
}

//________________________________________________________________________
void PHOSTrackSegments2RecParticles( Bool_t split=kFALSE, TString fileName = "galice.root") {
  
 // usage : 
  // 1. write RecParticles in the same file as TrackSegments ---------------  (OK)
  //root [0] .L Reconstruct.C++
  //root [1] TrackSegments2RecParticles()
  // 2. write RecParticles in a separate file, one per detector, from TrackSegments --------------- (Not needed) 
  //root [0] .L Reconstruct.C++
  //root [1] TrackSegments2RecParticles(kTRUE) // RecParticles saved in [DET}.RecData.root (DET=PHOS, EMCAL)

  delete gAlice ; 
  gAlice = 0 ; 
  
  AliPHOSPID * pp = 0 ; 
 
  if (split) 
    pp = new AliPHOSPIDv1("PHOS.RecData.root") ; 
  else 
    pp = new AliPHOSPIDv1(fileName) ; 
  
  pp->ExecuteTask("deb") ; 
  
  delete pp ;
}

//________________________________________________________________________
void PHOSDigits2RecParticles( Bool_t split=kFALSE, TString fileName = "galice.root") {
  
 // usage : 
  // 1. write RecPoints, TrackSegments and RecParticles in the same file as Digits --------------- (OK)
  //root [0] .L Reconstruct.C++
  //root [1] Digits2RecParticles()
  // 2. write RecPoints , TrackSegments and RecParticles in a separate file, one per detector, from Digits --------------- (OK)
  //root [0] .L Reconstruct.C++
  //root [1] Digits2RecParticles(kTRUE) // TrackSegments saved in [DET}.RecData.root (DET=PHOS, EMCAL)

 
  delete gAlice ; 
  gAlice = 0 ; 
 
  // PHOS
  AliPHOSClusterizer * cp = 0 ; 
 
  if (split) {
    cp = new AliPHOSClusterizerv1("PHOS.Digits.root") ; 
    cp->SetSplitFile() ; } 
  else 
    cp = new AliPHOSClusterizerv1(fileName) ; 
  
  cp->ExecuteTask("deb") ; 

  if (split) 
    delete cp ;
  
  AliPHOSTrackSegmentMaker * tmp = 0 ; 
  
  if (split) 
    tmp = new AliPHOSTrackSegmentMakerv1("PHOS.RecData.root") ; 
  else 
    tmp = new AliPHOSTrackSegmentMakerv1(fileName) ; 
  
  tmp->ExecuteTask("deb") ; 
  
  AliPHOSPID * pp = 0 ; 
 
  if (split) 
    pp = new AliPHOSPIDv1("PHOS.RecData.root") ; 
  else 
    pp = new AliPHOSPIDv1(fileName) ; 
  
  pp->ExecuteTask("deb") ; 
  
  delete tmp; 
  delete pp ; 
}

//________________________________________________________________________
void PHOSHits2Digits (Bool_t split=kFALSE, TString fileName = "galice.root") {
  // usage : 
  // 1. write (S)Digits in the same file as Hits --------------- (OK)
  //root [0] .L Reconstruct.C++
  //root [1] Hits2Digits()
  // 2. write (S)Digits in a separate file, one per detector, from Hits --------------- (OK)
  //root [0] .L Reconstruct.C++
  //root [1] Hits2Digits(kTRUE) // SDigits saved in [DET}.SDigits.root (DET=PHOS, EMCAL)
                                // Digits  saved in [DET}.Digits.root  (DET=PHOS, EMCAL)

  delete gAlice ; 
  gAlice = 0 ; 
  
  //PHOS
  AliPHOSSDigitizer * sdp = new AliPHOSSDigitizer(fileName) ; 
  if (split) 
    sdp->SetSplitFile() ;
  sdp->ExecuteTask("deb") ; 

  if (split) 
    delete sdp ; 

  AliPHOSDigitizer * dp = 0 ; 
 
  if (split) {
    dp = new AliPHOSDigitizer("PHOS.SDigits.root") ; 
    dp->SetSplitFile() ; } 
  else 
    dp = new AliPHOSDigitizer(fileName) ; 
  
  dp->ExecuteTask("deb") ; 

  if (split) 
    delete dp ; 

  if (!split) { 
    delete sdp ; 
    delete dp ; 
  }
}


 
